#include "core/localization/solid_relocalizer.h"

#include "utils/compute_profiling.h"
#include "utils/thread_scheduling.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <future>
#include <iomanip>
#include <limits>
#include <stdexcept>

#include <glog/logging.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <pcl/registration/icp.h>
#include <yaml-cpp/yaml.h>

#include "core/maps/map_frame.h"
#include "utils/compute_budget.h"

namespace lightning::loc {
namespace {

template <typename T>
T ReadOr(const YAML::Node& node, const char* key, const T& fallback) {
    return node && node[key] ? node[key].as<T>() : fallback;
}

SolidDescriptorOptions ReadDescriptorOptions(
    const YAML::Node& node, SolidDescriptorOptions options = {}) {
    options.range_bins = ReadOr<int>(node, "range_bins", options.range_bins);
    options.angle_bins = ReadOr<int>(node, "angle_bins", options.angle_bins);
    options.elevation_bins =
        ReadOr<int>(node, "elevation_bins", options.elevation_bins);
    options.min_distance = ReadOr<double>(node, "min_distance", options.min_distance);
    options.max_distance = ReadOr<double>(node, "max_distance", options.max_distance);
    options.elevation_min_deg =
        ReadOr<double>(node, "elevation_min_deg", options.elevation_min_deg);
    options.elevation_max_deg =
        ReadOr<double>(node, "elevation_max_deg", options.elevation_max_deg);
    return options;
}

YAML::Node WriteDescriptorOptions(const SolidDescriptorOptions& options) {
    YAML::Node node;
    node["range_bins"] = options.range_bins;
    node["angle_bins"] = options.angle_bins;
    node["elevation_bins"] = options.elevation_bins;
    node["min_distance"] = options.min_distance;
    node["max_distance"] = options.max_distance;
    node["elevation_min_deg"] = options.elevation_min_deg;
    node["elevation_max_deg"] = options.elevation_max_deg;
    return node;
}

std::vector<double> ToVector(const VecXd& values) {
    return std::vector<double>(values.data(), values.data() + values.size());
}

VecXd ReadVector(const YAML::Node& node, int expected_size, const char* label) {
    if (!node || !node.IsSequence()) {
        throw std::runtime_error(std::string("missing SOLiD ") + label);
    }
    const auto values = node.as<std::vector<double>>();
    if (static_cast<int>(values.size()) != expected_size ||
        std::any_of(values.begin(), values.end(),
                    [](double value) { return !std::isfinite(value); })) {
        throw std::runtime_error(std::string("invalid SOLiD ") + label);
    }
    return Eigen::Map<const VecXd>(values.data(), values.size());
}

SE3 ReadPose(const YAML::Node& entry) {
    if (!entry || !entry["pose_xyzw"]) {
        throw std::runtime_error("SOLiD database entry is missing pose_xyzw");
    }
    const auto values = entry["pose_xyzw"].as<std::vector<double>>();
    if (values.size() != 7 ||
        std::any_of(values.begin(), values.end(),
                    [](double value) { return !std::isfinite(value); })) {
        throw std::runtime_error("invalid SOLiD database pose_xyzw");
    }
    Quatd quaternion(values[6], values[3], values[4], values[5]);
    if (quaternion.norm() < 1e-9) {
        throw std::runtime_error("SOLiD database pose has zero quaternion");
    }
    return SE3(quaternion.normalized(), Vec3d(values[0], values[1], values[2]));
}

double Yaw(const SE3& pose) {
    const Mat3d rotation = pose.rotationMatrix();
    return std::atan2(rotation(1, 0), rotation(0, 0));
}

double AngleDistance(double left, double right) {
    return std::abs(std::atan2(std::sin(left - right), std::cos(left - right)));
}

PointCloudType ConvertCloud(const pcl::PointCloud<pcl::PointXYZI>& source) {
    PointCloudType output;
    output.reserve(source.size());
    for (const auto& point : source.points) {
        PointType converted;
        converted.x = point.x;
        converted.y = point.y;
        converted.z = point.z;
        converted.intensity = point.intensity;
        output.push_back(converted);
    }
    return output;
}

}  // namespace

SolidRelocalizer::SolidRelocalizer() : descriptor_engine_(descriptor_options_) {}

bool SolidRelocalizer::BuildDatabase(const std::string& config_path,
                                     const std::string& map_path) {
    namespace fs = std::filesystem;
    try {
        const YAML::Node root = YAML::LoadFile(config_path);
        const YAML::Node config = root["relocalization"];
        const YAML::Node solid = config["solid"];
        const std::string source_subdirectory = ReadOr<std::string>(
            solid, "source_database_subdirectory", "btc_relocalization");
        const std::string output_subdirectory = ReadOr<std::string>(
            solid, "database_subdirectory", "solid_relocalization");
        const SolidDescriptorOptions descriptor_options =
            ReadDescriptorOptions(solid["descriptor"]);
        const SolidDescriptorEngine engine(descriptor_options);
        const fs::path source_directory = fs::path(map_path) / source_subdirectory;
        const fs::path source_manifest_path = source_directory / "database.yaml";
        const YAML::Node source_manifest = YAML::LoadFile(source_manifest_path.string());
        const YAML::Node source_entries = source_manifest["entries"];
        if (!source_entries || !source_entries.IsSequence() || source_entries.size() == 0) {
            LOG(ERROR) << "BTC source database has no entries: " << source_manifest_path;
            return false;
        }

        YAML::Node manifest;
        manifest["schema_version"] = 1;
        manifest["source"] = "SOLiD descriptors generated from optimized BTC submap clouds";
        manifest["source_database_subdirectory"] = source_subdirectory;
        if (source_manifest["map_frame"]) manifest["map_frame"] = source_manifest["map_frame"];
        manifest["descriptor"] = WriteDescriptorOptions(descriptor_options);
        YAML::Node entries(YAML::NodeType::Sequence);
        for (std::size_t index = 0; index < source_entries.size(); ++index) {
            const YAML::Node source_entry = source_entries[index];
            if (!source_entry["cloud"] || !source_entry["pose_xyzw"]) {
                LOG(ERROR) << "invalid BTC source entry " << index;
                return false;
            }
            const fs::path cloud_path =
                source_directory / source_entry["cloud"].as<std::string>();
            pcl::PointCloud<pcl::PointXYZI> cloud;
            if (pcl::io::loadPCDFile(cloud_path.string(), cloud) != 0 || cloud.empty()) {
                LOG(ERROR) << "failed to load SOLiD source cloud: " << cloud_path;
                return false;
            }
            const auto descriptor = engine.Compute(ConvertCloud(cloud));
            if (!descriptor) {
                LOG(ERROR) << "failed to compute SOLiD descriptor for " << cloud_path;
                return false;
            }
            YAML::Node entry;
            entry["descriptor_id"] = static_cast<int>(index);
            entry["source_descriptor_id"] =
                ReadOr<int>(source_entry, "descriptor_id", static_cast<int>(index));
            if (source_entry["timestamp"]) entry["timestamp"] = source_entry["timestamp"];
            entry["pose_xyzw"] = source_entry["pose_xyzw"];
            entry["cloud"] = source_entry["cloud"];
            entry["accepted_points"] = static_cast<unsigned long long>(
                descriptor->accepted_points);
            entry["range_descriptor"] = ToVector(descriptor->range);
            entry["angle_descriptor"] = ToVector(descriptor->angle);
            entries.push_back(entry);
        }
        manifest["entries"] = entries;

        const fs::path output_directory = fs::path(map_path) / output_subdirectory;
        std::error_code error;
        fs::create_directories(output_directory, error);
        if (error) {
            LOG(ERROR) << "failed to create SOLiD database directory: " << error.message();
            return false;
        }
        const fs::path manifest_path = output_directory / "database.yaml";
        const fs::path temporary_path = output_directory / "database.yaml.tmp";
        std::ofstream stream(temporary_path);
        if (!stream) {
            LOG(ERROR) << "failed to open SOLiD manifest: " << temporary_path;
            return false;
        }
        stream << std::setprecision(17) << manifest;
        stream.close();
        if (!stream) {
            LOG(ERROR) << "failed to write SOLiD manifest: " << temporary_path;
            return false;
        }
        fs::rename(temporary_path, manifest_path, error);
        if (error) {
            LOG(ERROR) << "failed to publish SOLiD manifest: " << error.message();
            return false;
        }
        LOG(INFO) << "saved SOLiD relocalization database: entries=" << entries.size()
                  << ", path=" << manifest_path;
        return true;
    } catch (const std::exception& error) {
        LOG(ERROR) << "failed to build SOLiD database: " << error.what();
        return false;
    }
}

bool SolidRelocalizer::Init(const std::string& config_path,
                            const std::string& map_path,
                            const SE3& T_imu_lidar) {
    ready_ = false;
    entries_.clear();
    ResetQuery();
    T_imu_lidar_ = T_imu_lidar;
    try {
        const YAML::Node root = YAML::LoadFile(config_path);
        const YAML::Node config = root["relocalization"];
        const YAML::Node solid = config["solid"];
        options_.database_subdirectory = ReadOr<std::string>(
            solid, "database_subdirectory", options_.database_subdirectory);
        options_.source_database_subdirectory = ReadOr<std::string>(
            solid, "source_database_subdirectory", options_.source_database_subdirectory);
        options_.compute_backend =
            ReadOr<std::string>(solid, "compute_backend", options_.compute_backend);
        if (options_.compute_backend != "cpu" && options_.compute_backend != "auto") {
            LOG(ERROR) << "SOLiD compute_backend=" << options_.compute_backend
                       << " is unavailable in this build";
            return false;
        }
        if (solid && solid["query_submap_sizes"]) {
            options_.query_submap_sizes =
                solid["query_submap_sizes"].as<std::vector<int>>();
        }
        std::sort(options_.query_submap_sizes.begin(), options_.query_submap_sizes.end());
        options_.query_submap_sizes.erase(
            std::unique(options_.query_submap_sizes.begin(), options_.query_submap_sizes.end()),
            options_.query_submap_sizes.end());
        options_.query_stride = ReadOr<int>(solid, "query_stride", options_.query_stride);
        options_.top_k = ReadOr<int>(solid, "top_k", options_.top_k);
        options_.retrieval_pool_size = ReadOr<int>(
            solid, "retrieval_pool_size", options_.retrieval_pool_size);
        options_.min_similarity =
            ReadOr<double>(solid, "min_similarity", options_.min_similarity);
        options_.candidate_dedup_radius = ReadOr<double>(
            solid, "candidate_dedup_radius", options_.candidate_dedup_radius);
        options_.candidate_dedup_yaw_deg = ReadOr<double>(
            solid, "candidate_dedup_yaw_deg", options_.candidate_dedup_yaw_deg);
        options_.max_points_per_submap = ReadOr<int>(
            solid, "max_points_per_submap", options_.max_points_per_submap);
        options_.min_points_per_submap = ReadOr<int>(
            solid, "min_points_per_submap", options_.min_points_per_submap);
        options_.downsample_leaf_size = ReadOr<double>(
            solid, "downsample_leaf_size", options_.downsample_leaf_size);
        options_.refine_with_icp =
            ReadOr<bool>(solid, "refine_with_icp", options_.refine_with_icp);
        options_.icp_max_iterations = ReadOr<int>(
            solid, "icp_max_iterations", options_.icp_max_iterations);
        options_.icp_max_correspondence_distance = ReadOr<double>(
            solid, "icp_max_correspondence_distance",
            options_.icp_max_correspondence_distance);
        options_.icp_max_fitness_score = ReadOr<double>(
            solid, "icp_max_fitness_score", options_.icp_max_fitness_score);
        options_.icp_max_translation_correction = ReadOr<double>(
            solid, "icp_max_translation_correction",
            options_.icp_max_translation_correction);
        options_.icp_batch_size =
            ReadOr<int>(solid, "icp_batch_size", options_.icp_batch_size);
        options_.icp_workers =
            ReadOr<int>(solid, "icp_workers", options_.icp_workers);
        compute::ComputeBudget compute_budget;
        compute_budget.solid_icp_workers = options_.icp_workers;
        std::string compute_budget_error;
        if (!compute::LoadComputeBudget(root, compute_budget, &compute_budget_error)) {
            LOG(ERROR) << "invalid compute budget: " << compute_budget_error;
            return false;
        }
        options_.icp_workers = compute_budget.solid_icp_workers;
        options_.icp_worker_nice = compute_budget.solid_worker_nice;
        std::string affinity_error;
        if (!threading::ParseCpuList(compute_budget.solid_cpu_affinity,
                                     options_.icp_cpu_affinity,
                                     &affinity_error)) {
            LOG(ERROR) << "invalid SOLiD worker CPU affinity: " << affinity_error;
            return false;
        }
        if (solid && solid["icp_yaw_hypothesis_offsets_deg"]) {
            options_.icp_yaw_hypothesis_offsets_deg =
                solid["icp_yaw_hypothesis_offsets_deg"].as<std::vector<double>>();
        }
        if (solid && solid["debug_candidate_ids"]) {
            options_.debug_candidate_ids =
                solid["debug_candidate_ids"].as<std::vector<int>>();
        }
        if (options_.query_submap_sizes.empty() || options_.query_stride <= 0 ||
            options_.top_k <= 0 || options_.retrieval_pool_size < options_.top_k ||
            options_.min_similarity < -1.0 ||
            options_.min_similarity > 1.0 || options_.candidate_dedup_radius < 0.0 ||
            options_.candidate_dedup_yaw_deg < 0.0 ||
            options_.max_points_per_submap <= 0 ||
            options_.min_points_per_submap <= 0 ||
            options_.downsample_leaf_size < 0.0 ||
            options_.icp_max_iterations <= 0 ||
            options_.icp_max_correspondence_distance <= 0.0 ||
            options_.icp_max_fitness_score < 0.0 ||
            options_.icp_max_translation_correction < 0.0 ||
            options_.icp_batch_size <= 0 || options_.icp_workers <= 0 ||
            options_.icp_yaw_hypothesis_offsets_deg.empty() ||
            std::any_of(options_.icp_yaw_hypothesis_offsets_deg.begin(),
                        options_.icp_yaw_hypothesis_offsets_deg.end(),
                        [](double value) { return !std::isfinite(value); }) ||
            std::any_of(options_.query_submap_sizes.begin(),
                        options_.query_submap_sizes.end(),
                        [](int size) { return size <= 0; })) {
            LOG(ERROR) << "invalid SOLiD relocalization configuration";
            return false;
        }

        const std::filesystem::path database_directory =
            std::filesystem::path(map_path) / options_.database_subdirectory;
        const std::filesystem::path manifest_path = database_directory / "database.yaml";
        const YAML::Node manifest = YAML::LoadFile(manifest_path.string());
        if (ReadOr<int>(manifest, "schema_version", 0) != 1) {
            LOG(ERROR) << "unsupported SOLiD database schema: " << manifest_path;
            return false;
        }
        map_frame::ExportOptions export_options;
        std::string map_frame_error;
        if (!map_frame::ReadExportOptions(root, export_options, map_frame_error)) {
            LOG(ERROR) << map_frame_error;
            return false;
        }
        if (export_options.normalize_start_ground_z) {
            map_frame::Metadata package_metadata;
            map_frame::Metadata database_metadata;
            if (!map_frame::LoadMetadata(map_path, package_metadata, map_frame_error) ||
                !map_frame::ReadTransformReference(
                    manifest["map_frame"], database_metadata, map_frame_error) ||
                !map_frame::SameTransform(package_metadata, database_metadata)) {
                LOG(ERROR) << "SOLiD database map frame mismatch: " << map_frame_error;
                return false;
            }
        }
        descriptor_options_ = ReadDescriptorOptions(manifest["descriptor"]);
        descriptor_engine_ = SolidDescriptorEngine(descriptor_options_);
        const YAML::Node yaml_entries = manifest["entries"];
        if (!yaml_entries || !yaml_entries.IsSequence() || yaml_entries.size() == 0) {
            LOG(ERROR) << "SOLiD database has no entries: " << manifest_path;
            return false;
        }
        entries_.reserve(yaml_entries.size());
        for (std::size_t index = 0; index < yaml_entries.size(); ++index) {
            const YAML::Node yaml_entry = yaml_entries[index];
            if (ReadOr<int>(yaml_entry, "descriptor_id", -1) !=
                static_cast<int>(index)) {
                LOG(ERROR) << "SOLiD descriptor ids must be dense";
                return false;
            }
            SolidDescriptor descriptor;
            descriptor.range = ReadVector(
                yaml_entry["range_descriptor"], descriptor_options_.range_bins,
                "range descriptor");
            descriptor.angle = ReadVector(
                yaml_entry["angle_descriptor"], descriptor_options_.angle_bins,
                "angle descriptor");
            descriptor.accepted_points = ReadOr<std::size_t>(
                yaml_entry, "accepted_points", 0);
            const std::string cloud_name = ReadOr<std::string>(
                yaml_entry, "cloud", "submap_" +
                    std::string(6 - std::min<std::size_t>(6, std::to_string(index).size()), '0') +
                    std::to_string(index) + ".pcd");
            entries_.push_back(DatabaseEntry{
                static_cast<int>(index), ReadPose(yaml_entry), std::move(descriptor),
                (std::filesystem::path(map_path) /
                 options_.source_database_subdirectory / cloud_name).string(),
                nullptr});
        }
        ready_ = !entries_.empty();
        LOG(INFO) << "loaded SOLiD relocalization database: entries=" << entries_.size()
                  << ", query_sizes=" << options_.query_submap_sizes.front() << ".."
                  << options_.query_submap_sizes.back()
                  << ", retrieval_pool_size=" << options_.retrieval_pool_size
                  << ", validation_top_k=" << options_.top_k
                  << ", min_similarity=" << options_.min_similarity
                  << ", refine_with_icp=" << options_.refine_with_icp
                  << ", icp_batch_size=" << options_.icp_batch_size
                  << ", icp_workers=" << options_.icp_workers
                  << ", icp_worker_nice=" << options_.icp_worker_nice
                  << ", icp_cpu_affinity="
                  << threading::FormatCpuList(options_.icp_cpu_affinity)
                  << ", icp_yaw_hypotheses="
                  << options_.icp_yaw_hypothesis_offsets_deg.size()
                  << ", compute_backend=cpu, path=" << manifest_path;
        return ready_;
    } catch (const std::exception& error) {
        LOG(ERROR) << "failed to initialize SOLiD relocalization: " << error.what();
        entries_.clear();
        ready_ = false;
        return false;
    }
}

void SolidRelocalizer::ResetQuery() {
    query_frames_.clear();
    frames_since_query_ = 0;
    icp_batch_cursor_ = 0;
}

pcl::PointCloud<pcl::PointXYZI>::Ptr SolidRelocalizer::BuildQuerySubmap(
    std::size_t frame_count) const {
    pcl::PointCloud<pcl::PointXYZI>::Ptr combined(
        new pcl::PointCloud<pcl::PointXYZI>);
    if (query_frames_.empty() || frame_count == 0) return combined;
    frame_count = std::min(frame_count, query_frames_.size());
    const auto first = query_frames_.end() -
                       static_cast<std::ptrdiff_t>(frame_count);
    std::size_t available_points = 0;
    for (auto frame = first; frame != query_frames_.end(); ++frame) {
        if (frame->cloud) available_points += frame->cloud->size();
    }
    const std::size_t maximum =
        static_cast<std::size_t>(options_.max_points_per_submap);
    const std::size_t stride =
        std::max<std::size_t>(1, (available_points + maximum - 1) / maximum);
    combined->reserve(std::min(available_points, maximum));
    const SE3 T_current_lidar_odom = query_frames_.back().T_odom_lidar.inverse();
    std::size_t source_index = 0;
    for (auto frame = first; frame != query_frames_.end(); ++frame) {
        if (!frame->cloud) continue;
        const SE3 T_current_lidar_frame_lidar =
            T_current_lidar_odom * frame->T_odom_lidar;
        for (const auto& source : frame->cloud->points) {
            if (source_index++ % stride != 0 || !std::isfinite(source.x) ||
                !std::isfinite(source.y) || !std::isfinite(source.z)) {
                continue;
            }
            const Vec3d position = T_current_lidar_frame_lidar *
                                   Vec3d(source.x, source.y, source.z);
            pcl::PointXYZI point;
            point.x = static_cast<float>(position.x());
            point.y = static_cast<float>(position.y());
            point.z = static_cast<float>(position.z());
            point.intensity = source.intensity;
            combined->push_back(point);
        }
    }
    if (options_.downsample_leaf_size <= 0.0 || combined->empty()) return combined;
    pcl::PointCloud<pcl::PointXYZI>::Ptr filtered(
        new pcl::PointCloud<pcl::PointXYZI>);
    pcl::VoxelGrid<pcl::PointXYZI> voxel;
    const float leaf = static_cast<float>(options_.downsample_leaf_size);
    voxel.setLeafSize(leaf, leaf, leaf);
    voxel.setInputCloud(combined);
    voxel.filter(*filtered);
    return filtered;
}

bool SolidRelocalizer::RefineCandidateWithIcp(
    const pcl::PointCloud<pcl::PointXYZI>::ConstPtr& query,
    DatabaseEntry& entry, RelocalizationCandidate& candidate,
    double& fitness_score) {
    fitness_score = std::numeric_limits<double>::infinity();
    if (!query || query->empty()) return false;
    if (!entry.cloud) {
        entry.cloud.reset(new pcl::PointCloud<pcl::PointXYZI>);
        if (pcl::io::loadPCDFile(entry.cloud_path, *entry.cloud) != 0 ||
            entry.cloud->empty()) {
            LOG(ERROR) << "failed to load SOLiD ICP cloud: " << entry.cloud_path;
            entry.cloud.reset();
            return false;
        }
    }

    const SE3 T_world_current_lidar = candidate.T_world_imu * T_imu_lidar_;
    const SE3 initial = entry.T_world_lidar.inverse() * T_world_current_lidar;
    pcl::IterativeClosestPoint<pcl::PointXYZI, pcl::PointXYZI> icp;
    icp.setInputSource(query);
    icp.setInputTarget(entry.cloud);
    icp.setMaximumIterations(options_.icp_max_iterations);
    icp.setMaxCorrespondenceDistance(options_.icp_max_correspondence_distance);
    icp.setTransformationEpsilon(1e-4);
    icp.setEuclideanFitnessEpsilon(1e-4);
    pcl::PointCloud<pcl::PointXYZI> aligned;
    icp.align(aligned, initial.matrix().cast<float>());
    if (!icp.hasConverged()) return false;
    fitness_score = icp.getFitnessScore(options_.icp_max_correspondence_distance);
    if (!std::isfinite(fitness_score) ||
        fitness_score > options_.icp_max_fitness_score) {
        return false;
    }

    const Eigen::Matrix4d matrix = icp.getFinalTransformation().cast<double>();
    const SE3 refined(Quatd(matrix.block<3, 3>(0, 0)).normalized(),
                      matrix.block<3, 1>(0, 3));
    const SE3 correction = initial.inverse() * refined;
    if (!correction.translation().allFinite() ||
        correction.translation().norm() >
            options_.icp_max_translation_correction) {
        return false;
    }

    const SE3 T_map_current_lidar = entry.T_world_lidar * refined;
    const SE3 T_map_odom_raw =
        T_map_current_lidar * query_frames_.back().T_odom_lidar.inverse();
    const SE3 T_map_odom_planar(
        Quatd(AngAxisd(Yaw(T_map_odom_raw), Vec3d::UnitZ())),
        T_map_odom_raw.translation());
    candidate.T_world_imu =
        T_map_odom_planar * query_frames_.back().T_odom_lidar *
        T_imu_lidar_.inverse();
    return true;
}

std::optional<RelocalizationResult> SolidRelocalizer::AddFrame(
    const CloudPtr& cloud, const SE3& T_odom_imu, double timestamp) {
    if (!ready_ || !cloud || cloud->empty() || !std::isfinite(timestamp)) {
        return std::nullopt;
    }
    QueryFrame frame;
    frame.cloud.reset(new PointCloudType(*cloud));
    frame.T_odom_lidar = T_odom_imu * T_imu_lidar_;
    frame.timestamp = timestamp;
    query_frames_.push_back(std::move(frame));
    ++frames_since_query_;
    const std::size_t maximum_window =
        static_cast<std::size_t>(options_.query_submap_sizes.back());
    while (query_frames_.size() > maximum_window) query_frames_.pop_front();
    const std::size_t minimum_window =
        static_cast<std::size_t>(options_.query_submap_sizes.front());
    if (query_frames_.size() < minimum_window) return std::nullopt;
    const bool first_query = query_frames_.size() == minimum_window &&
                             frames_since_query_ >= minimum_window;
    if (!first_query &&
        frames_since_query_ < static_cast<std::size_t>(options_.query_stride)) {
        return std::nullopt;
    }
    frames_since_query_ = 0;

    RelocalizationResult result;
    result.attempted = true;
    result.timestamp = query_frames_.back().timestamp;
    const bool profiling_enabled = profiling::ComputeProfilingEnabled();
    profiling::Stopwatch search_profile_timer(profiling_enabled);
    const auto begin = std::chrono::steady_clock::now();
    std::vector<RelocalizationCandidate> hypotheses;
    std::vector<std::pair<int, pcl::PointCloud<pcl::PointXYZI>::Ptr>> query_clouds;
    bool descriptor_generated = false;
    bool score_candidate_found = false;
    for (const int configured_size : options_.query_submap_sizes) {
        if (configured_size > static_cast<int>(query_frames_.size())) continue;
        const auto query_cloud = BuildQuerySubmap(configured_size);
        query_clouds.emplace_back(configured_size, query_cloud);
        result.point_count = std::max(result.point_count, query_cloud->size());
        if (query_cloud->size() <
            static_cast<std::size_t>(options_.min_points_per_submap)) {
            continue;
        }
        const auto query_descriptor = descriptor_engine_.Compute(ConvertCloud(*query_cloud));
        if (!query_descriptor) continue;
        descriptor_generated = true;
        result.descriptor_count = std::max<std::size_t>(
            result.descriptor_count,
            query_descriptor->range.size() + query_descriptor->angle.size());

        struct ScoredEntry {
            int index = -1;
            double score = -1.0;
        };
        std::vector<ScoredEntry> scored;
        scored.reserve(entries_.size());
        for (std::size_t index = 0; index < entries_.size(); ++index) {
            const double score = SolidDescriptorEngine::Similarity(
                *query_descriptor, entries_[index].descriptor);
            if (std::isfinite(score)) {
                scored.push_back(ScoredEntry{static_cast<int>(index), score});
            }
        }
        std::stable_sort(scored.begin(), scored.end(), [](const auto& left, const auto& right) {
            if (left.score != right.score) return left.score > right.score;
            return left.index < right.index;
        });
        if (!scored.empty() && scored.front().score >= options_.min_similarity) {
            score_candidate_found = true;
        }
        for (const int debug_id : options_.debug_candidate_ids) {
            const auto found = std::find_if(
                scored.begin(), scored.end(), [&](const ScoredEntry& value) {
                    return entries_[value.index].descriptor_id == debug_id;
                });
            if (found == scored.end()) continue;
            const std::size_t rank =
                static_cast<std::size_t>(std::distance(scored.begin(), found));
            const auto debug_yaw =
                SolidDescriptorEngine::EstimateCandidateFromQueryYaw(
                    *query_descriptor, entries_[found->index].descriptor);
            LOG(INFO) << "SOLID_RETRIEVAL_DIAGNOSTIC query_frames="
                      << configured_size << ", candidate=" << debug_id
                      << ", rank=" << rank << ", score=" << found->score
                      << ", yaw_deg="
                      << (debug_yaw ? *debug_yaw * 180.0 / M_PI
                                    : std::numeric_limits<double>::quiet_NaN());
        }
        const std::size_t limit = std::min<std::size_t>(
            scored.size(), static_cast<std::size_t>(options_.retrieval_pool_size));
        for (std::size_t rank = 0; rank < limit; ++rank) {
            if (scored[rank].score < options_.min_similarity) break;
            const auto& entry = entries_[scored[rank].index];
            const auto yaw = SolidDescriptorEngine::EstimateCandidateFromQueryYaw(
                *query_descriptor, entry.descriptor);
            if (!yaw) continue;
            const SE3 T_source_lidar_current_lidar(
                Quatd(AngAxisd(*yaw, Vec3d::UnitZ())), Vec3d::Zero());
            // SOLiD estimates heading in the sensor frame.  Directly applying
            // that local-Z rotation to a gravity-aligned but tilted Livox pose
            // can introduce artificial roll/pitch.  Project the resulting
            // map<-odom rotation onto world Z, preserving the LIO gravity
            // estimate while retaining the descriptor's heading hypothesis.
            const SE3 T_map_odom_raw =
                entry.T_world_lidar * T_source_lidar_current_lidar *
                query_frames_.back().T_odom_lidar.inverse();
            const SE3 T_map_odom_planar(
                Quatd(AngAxisd(Yaw(T_map_odom_raw), Vec3d::UnitZ())),
                T_map_odom_raw.translation());
            RelocalizationCandidate candidate;
            candidate.candidate_id = entry.descriptor_id;
            candidate.score = scored[rank].score;
            candidate.query_submap_size = configured_size;
            candidate.point_count = query_cloud->size();
            candidate.descriptor_count =
                query_descriptor->range.size() + query_descriptor->angle.size();
            candidate.T_world_imu =
                T_map_odom_planar * query_frames_.back().T_odom_lidar *
                T_imu_lidar_.inverse();
            hypotheses.push_back(std::move(candidate));
        }
    }
    const double maximum_yaw_distance =
        options_.candidate_dedup_yaw_deg * 3.14159265358979323846 / 180.0;
    std::stable_sort(hypotheses.begin(), hypotheses.end(), [](const auto& left, const auto& right) {
        if (left.score != right.score) return left.score > right.score;
        return left.candidate_id < right.candidate_id;
    });
    for (const auto& hypothesis : hypotheses) {
        const bool duplicate = std::any_of(
            result.candidates.begin(), result.candidates.end(),
            [&](const auto& selected) {
                return (selected.T_world_imu.translation() -
                        hypothesis.T_world_imu.translation())
                               .norm() <= options_.candidate_dedup_radius &&
                       AngleDistance(Yaw(selected.T_world_imu),
                                     Yaw(hypothesis.T_world_imu)) <=
                           maximum_yaw_distance;
            });
        if (!duplicate) result.candidates.push_back(hypothesis);
        if (result.candidates.size() >= static_cast<std::size_t>(options_.top_k)) break;
    }
    if (options_.refine_with_icp) {
        const std::size_t available = result.candidates.size();
        const std::size_t batch_size = std::min<std::size_t>(
            available, static_cast<std::size_t>(options_.icp_batch_size));
        std::vector<RelocalizationCandidate> batch;
        batch.reserve(batch_size);
        if (available > 0) {
            const std::size_t start = icp_batch_cursor_ % available;
            for (std::size_t offset = 0; offset < batch_size; ++offset) {
                batch.push_back(result.candidates[(start + offset) % available]);
            }
            icp_batch_cursor_ = (start + batch_size) % available;
        }

        // Load target clouds on this thread so parallel ICP workers only read
        // shared PCD data and never race while publishing a lazy cache entry.
        for (const auto& candidate : batch) {
            auto& entry = entries_[static_cast<std::size_t>(candidate.candidate_id)];
            if (entry.cloud) continue;
            entry.cloud.reset(new pcl::PointCloud<pcl::PointXYZI>);
            if (pcl::io::loadPCDFile(entry.cloud_path, *entry.cloud) != 0 ||
                entry.cloud->empty()) {
                LOG(ERROR) << "failed to load SOLiD ICP cloud: " << entry.cloud_path;
                entry.cloud.reset();
            }
        }

        struct IcpResult {
            bool accepted = false;
            RelocalizationCandidate candidate;
            double fitness_score = std::numeric_limits<double>::infinity();
            std::size_t batch_index = 0;
            std::size_t yaw_index = 0;
        };
        struct IcpTask {
            RelocalizationCandidate candidate;
            std::size_t batch_index = 0;
            std::size_t yaw_index = 0;
        };
        std::vector<IcpTask> icp_tasks;
        icp_tasks.reserve(batch.size() *
                          options_.icp_yaw_hypothesis_offsets_deg.size());
        for (std::size_t batch_index = 0; batch_index < batch.size();
             ++batch_index) {
            for (std::size_t yaw_index = 0;
                 yaw_index < options_.icp_yaw_hypothesis_offsets_deg.size();
                 ++yaw_index) {
                auto candidate = batch[batch_index];
                const double radians =
                    options_.icp_yaw_hypothesis_offsets_deg[yaw_index] * M_PI / 180.0;
                if (std::abs(radians) > 1e-12) {
                    const Mat3d rotation =
                        AngAxisd(radians, Vec3d::UnitZ()).toRotationMatrix() *
                        candidate.T_world_imu.rotationMatrix();
                    candidate.T_world_imu = SE3(
                        Quatd(rotation).normalized(),
                        candidate.T_world_imu.translation());
                }
                icp_tasks.push_back(
                    IcpTask{std::move(candidate), batch_index, yaw_index});
            }
        }
        std::vector<IcpResult> icp_results(icp_tasks.size());
        const std::size_t worker_count = std::min<std::size_t>(
            icp_tasks.size(), static_cast<std::size_t>(options_.icp_workers));
        std::vector<std::future<void>> workers;
        workers.reserve(worker_count);
        std::atomic<bool> scheduling_warning_logged{false};
        for (std::size_t worker = 0; worker < worker_count; ++worker) {
            workers.emplace_back(std::async(std::launch::async, [&, worker]() {
                std::string scheduling_error;
                if (!threading::ConfigureCurrentThread(
                        "solid_icp_" + std::to_string(worker),
                        options_.icp_worker_nice, options_.icp_cpu_affinity,
                        &scheduling_error) &&
                    !scheduling_warning_logged.exchange(true)) {
                    LOG(WARNING) << "failed to apply SOLiD worker scheduling policy: "
                                 << scheduling_error;
                }
                for (std::size_t index = worker; index < icp_tasks.size();
                     index += worker_count) {
                    auto candidate = icp_tasks[index].candidate;
                    const auto query = std::find_if(
                        query_clouds.begin(), query_clouds.end(),
                        [&](const auto& value) {
                            return value.first == candidate.query_submap_size;
                        });
                    if (query == query_clouds.end()) continue;
                    auto& entry = entries_[static_cast<std::size_t>(
                        candidate.candidate_id)];
                    if (!entry.cloud) continue;
                    double fitness_score = std::numeric_limits<double>::infinity();
                    if (!RefineCandidateWithIcp(
                            query->second, entry, candidate, fitness_score)) {
                        continue;
                    }
                    icp_results[index] =
                        IcpResult{true, std::move(candidate), fitness_score,
                                  icp_tasks[index].batch_index,
                                  icp_tasks[index].yaw_index};
                }
            }));
        }
        for (auto& worker : workers) worker.get();

        std::vector<RelocalizationCandidate> refined_candidates;
        refined_candidates.reserve(batch.size());
        for (std::size_t batch_index = 0; batch_index < batch.size();
             ++batch_index) {
            IcpResult* best = nullptr;
            for (auto& icp_result : icp_results) {
                if (!icp_result.accepted ||
                    icp_result.batch_index != batch_index) {
                    continue;
                }
                if (!best || icp_result.fitness_score < best->fitness_score) {
                    best = &icp_result;
                }
            }
            if (!best) continue;
            auto& candidate = best->candidate;
            LOG_IF(INFO, std::find(options_.debug_candidate_ids.begin(),
                                   options_.debug_candidate_ids.end(),
                                   candidate.candidate_id) !=
                             options_.debug_candidate_ids.end())
                << "SOLID_ICP_DIAGNOSTIC candidate=" << candidate.candidate_id
                << ", query_frames=" << candidate.query_submap_size
                << ", fitness=" << best->fitness_score
                << ", yaw_offset_deg="
                << options_.icp_yaw_hypothesis_offsets_deg[best->yaw_index]
                << ", pose=" << candidate.T_world_imu.translation().transpose();
            refined_candidates.push_back(std::move(candidate));
        }
        result.candidates = std::move(refined_candidates);
    }
    // Report the complete SOLiD search cost, including candidate
    // de-duplication, PCD loading and ICP refinement. Previously this timer
    // stopped before ICP and hid the dominant relocalization latency.
    result.search_time_ms = std::chrono::duration<double, std::milli>(
                                std::chrono::steady_clock::now() - begin)
                                .count();
    const profiling::TimingSample search_timing = search_profile_timer.Stop();
    const auto emit_profile = [&]() {
        if (!profiling_enabled) return;
        LOG(INFO) << "COMPUTE_BENCH_EVENT module=global_relocalization backend=solid"
                  << " timestamp_s=" << result.timestamp
                  << " query_frames=" << query_frames_.size()
                  << " query_points=" << result.point_count
                  << " candidates=" << result.candidates.size()
                  << " accepted=" << result.accepted
                  << ' ' << profiling::FormatTimingSample("search", search_timing);
    };
    result.candidate_found = score_candidate_found;
    if (result.candidates.empty()) {
        result.reason = descriptor_generated ? "score_below_threshold"
                                             : "no_query_descriptor";
        emit_profile();
        return result;
    }
    const auto& best = result.candidates.front();
    result.accepted = true;
    result.candidate_found = true;
    result.candidate_id = best.candidate_id;
    result.score = best.score;
    result.point_count = best.point_count;
    result.descriptor_count = best.descriptor_count;
    result.T_world_imu = best.T_world_imu;
    result.reason = "accepted";
    emit_profile();
    return result;
}

}  // namespace lightning::loc
