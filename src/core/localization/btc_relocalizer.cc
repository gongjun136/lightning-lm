#include "core/localization/btc_relocalizer.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <stdexcept>

#include <glog/logging.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <yaml-cpp/yaml.h>

#include "core/maps/map_frame.h"

namespace lightning::loc {
namespace {

template <typename T>
T ReadOr(const YAML::Node& node, const char* key, const T& fallback) {
    return node && node[key] ? node[key].as<T>() : fallback;
}

void ReadDescriptorConfig(const YAML::Node& node, ConfigSetting& config) {
    config.useful_corner_num_ = ReadOr<int>(node, "useful_corner_num", config.useful_corner_num_);
    config.plane_merge_normal_thre_ =
        ReadOr<float>(node, "plane_merge_normal_threshold", config.plane_merge_normal_thre_);
    config.plane_merge_dis_thre_ =
        ReadOr<float>(node, "plane_merge_distance_threshold", config.plane_merge_dis_thre_);
    config.plane_detection_thre_ =
        ReadOr<float>(node, "plane_detection_threshold", config.plane_detection_thre_);
    config.voxel_size_ = ReadOr<float>(node, "voxel_size", config.voxel_size_);
    config.voxel_init_num_ = ReadOr<int>(node, "voxel_init_points", config.voxel_init_num_);
    config.proj_plane_num_ = ReadOr<int>(node, "projection_plane_count", config.proj_plane_num_);
    config.proj_image_resolution_ =
        ReadOr<float>(node, "projection_resolution", config.proj_image_resolution_);
    config.proj_image_high_inc_ =
        ReadOr<float>(node, "projection_height_increment", config.proj_image_high_inc_);
    config.proj_dis_min_ = ReadOr<float>(node, "projection_min_distance", config.proj_dis_min_);
    config.proj_dis_max_ = ReadOr<float>(node, "projection_max_distance", config.proj_dis_max_);
    config.summary_min_thre_ = ReadOr<float>(node, "summary_min_threshold", config.summary_min_thre_);
    config.line_filter_enable_ = ReadOr<int>(node, "line_filter", config.line_filter_enable_);
    config.touch_filter_enable_ = ReadOr<int>(node, "touch_filter", config.touch_filter_enable_);
    config.descriptor_near_num_ =
        ReadOr<float>(node, "descriptor_near_count", config.descriptor_near_num_);
    config.descriptor_min_len_ =
        ReadOr<float>(node, "descriptor_min_length", config.descriptor_min_len_);
    config.descriptor_max_len_ =
        ReadOr<float>(node, "descriptor_max_length", config.descriptor_max_len_);
    config.non_max_suppression_radius_ =
        ReadOr<float>(node, "non_max_suppression_radius", config.non_max_suppression_radius_);
    config.std_side_resolution_ =
        ReadOr<float>(node, "triangle_side_resolution", config.std_side_resolution_);
    config.skip_near_num_ = ReadOr<int>(node, "skip_near_descriptors", config.skip_near_num_);
    config.candidate_num_ = ReadOr<int>(node, "candidate_count", config.candidate_num_);
    config.candidate_min_votes_ =
        ReadOr<int>(node, "candidate_min_votes", config.candidate_min_votes_);
    config.verification_threads_ =
        ReadOr<int>(node, "verification_threads", config.verification_threads_);
    config.rough_dis_threshold_ =
        ReadOr<float>(node, "rough_distance_threshold", config.rough_dis_threshold_);
    config.similarity_threshold_ =
        ReadOr<float>(node, "similarity_threshold", config.similarity_threshold_);
    config.icp_threshold_ =
        ReadOr<float>(node, "internal_icp_threshold", config.icp_threshold_);
    config.normal_threshold_ = ReadOr<float>(node, "normal_threshold", config.normal_threshold_);
    config.dis_threshold_ =
        ReadOr<float>(node, "plane_distance_threshold", config.dis_threshold_);
}

SE3 ReadPose(const YAML::Node& entry) {
    if (!entry || !entry["pose_xyzw"]) {
        throw std::runtime_error("BTC database entry is missing pose_xyzw");
    }
    const auto values = entry["pose_xyzw"].as<std::vector<double>>();
    if (values.size() != 7) {
        throw std::runtime_error("BTC database pose_xyzw must have 7 values");
    }
    for (const double value : values) {
        if (!std::isfinite(value)) throw std::runtime_error("BTC database pose contains non-finite value");
    }
    Quatd quaternion(values[6], values[3], values[4], values[5]);
    if (quaternion.norm() < 1e-9) throw std::runtime_error("BTC database pose has zero quaternion");
    quaternion.normalize();
    return SE3(quaternion, Vec3d(values[0], values[1], values[2]));
}

double Yaw(const SE3& pose) {
    const Mat3d rotation = pose.rotationMatrix();
    return std::atan2(rotation(1, 0), rotation(0, 0));
}

double AngleDistance(double left, double right) {
    return std::abs(std::atan2(std::sin(left - right), std::cos(left - right)));
}

}  // namespace

BtcRelocalizer::BtcRelocalizer() : manager_(descriptor_config_) {}

bool BtcRelocalizer::Init(const std::string& config_path, const std::string& map_path,
                          const SE3& T_imu_lidar) {
    ready_ = false;
    entries_.clear();
    query_frames_.clear();
    frames_since_query_ = 0;
    T_imu_lidar_ = T_imu_lidar;

    try {
        const YAML::Node root = YAML::LoadFile(config_path);
        const YAML::Node config = root["relocalization"];
        map_frame::ExportOptions export_options;
        std::string map_frame_error;
        if (!map_frame::ReadExportOptions(root, export_options, map_frame_error)) {
            LOG(ERROR) << map_frame_error;
            return false;
        }
        map_frame::Metadata package_metadata;
        if (export_options.normalize_start_ground_z &&
            !map_frame::LoadMetadata(map_path, package_metadata, map_frame_error)) {
            LOG(ERROR) << map_frame_error;
            return false;
        }
        options_.enabled = ReadOr<bool>(config, "enabled", false);
        if (!options_.enabled) {
            LOG(INFO) << "BTC relocalization is disabled";
            return true;
        }
        options_.database_subdirectory =
            ReadOr<std::string>(config, "database_subdirectory", options_.database_subdirectory);
        options_.min_btc_score = ReadOr<double>(config, "min_btc_score", options_.min_btc_score);

        const std::filesystem::path database_directory =
            std::filesystem::path(map_path) / options_.database_subdirectory;
        const std::filesystem::path manifest_path = database_directory / "database.yaml";
        if (!std::filesystem::exists(manifest_path)) {
            LOG(ERROR) << "BTC relocalization database not found: " << manifest_path;
            return false;
        }

        const YAML::Node manifest = YAML::LoadFile(manifest_path.string());
        const int schema_version = ReadOr<int>(manifest, "schema_version", 0);
        if (schema_version != 1 && schema_version != 2 && schema_version != 3) {
            LOG(ERROR) << "unsupported BTC relocalization database schema: " << manifest_path;
            return false;
        }
        if (export_options.normalize_start_ground_z) {
            if (schema_version < 2) {
                LOG(ERROR) << "normalized map requires a map-frame-aware BTC database: "
                           << manifest_path;
                return false;
            }
            map_frame::Metadata database_metadata;
            if (!map_frame::ReadTransformReference(
                    manifest["map_frame"], database_metadata, map_frame_error) ||
                !map_frame::SameTransform(package_metadata, database_metadata)) {
                LOG(ERROR) << "BTC database map frame does not match the global map: "
                           << (map_frame_error.empty() ? "transform mismatch" : map_frame_error);
                return false;
            }
        }
        if (config && config["query_submap_sizes"]) {
            options_.query_submap_sizes = config["query_submap_sizes"].as<std::vector<int>>();
        } else {
            options_.query_submap_sizes = {ReadOr<int>(
                config, "query_submap_size", ReadOr<int>(manifest, "descriptor_submap_size", 10))};
        }
        std::sort(options_.query_submap_sizes.begin(), options_.query_submap_sizes.end());
        options_.query_submap_sizes.erase(
            std::unique(options_.query_submap_sizes.begin(), options_.query_submap_sizes.end()),
            options_.query_submap_sizes.end());
        options_.query_submap_size = options_.query_submap_sizes.empty()
                                             ? 0
                                             : options_.query_submap_sizes.back();
        options_.query_stride = ReadOr<int>(config, "query_stride", options_.query_stride);
        options_.top_k = ReadOr<int>(config, "top_k", options_.top_k);
        options_.candidate_dedup_radius =
            ReadOr<double>(config, "candidate_dedup_radius", options_.candidate_dedup_radius);
        options_.candidate_dedup_yaw_deg =
            ReadOr<double>(config, "candidate_dedup_yaw_deg", options_.candidate_dedup_yaw_deg);
        options_.max_points_per_submap = ReadOr<int>(
            config, "max_points_per_submap", ReadOr<int>(manifest, "max_points_per_submap", 50000));
        options_.min_points_per_submap =
            ReadOr<int>(config, "min_points_per_submap", options_.min_points_per_submap);
        options_.downsample_leaf_size = ReadOr<double>(
            config, "downsample_leaf_size", ReadOr<double>(manifest, "downsample_leaf_size", 0.20));
        if (options_.query_submap_size <= 0 || options_.query_stride <= 0 || options_.top_k <= 0 ||
            options_.candidate_dedup_radius < 0.0 || options_.candidate_dedup_yaw_deg < 0.0 ||
            options_.max_points_per_submap <= 0 || options_.min_points_per_submap <= 0 ||
            options_.min_btc_score < 0.0 ||
            std::any_of(options_.query_submap_sizes.begin(), options_.query_submap_sizes.end(),
                        [](int size) { return size <= 0; })) {
            LOG(ERROR) << "invalid BTC relocalization configuration";
            return false;
        }

        descriptor_config_ = ConfigSetting{};
        ReadDescriptorConfig(manifest["descriptor"], descriptor_config_);
        descriptor_config_.candidate_num_ =
            ReadOr<int>(config, "search_candidate_count", descriptor_config_.candidate_num_);
        descriptor_config_.candidate_min_votes_ = ReadOr<int>(
            config, "search_candidate_min_votes", descriptor_config_.candidate_min_votes_);
        descriptor_config_.rough_dis_threshold_ = ReadOr<float>(
            config, "search_rough_distance_threshold", descriptor_config_.rough_dis_threshold_);
        descriptor_config_.similarity_threshold_ = ReadOr<float>(
            config, "search_similarity_threshold", descriptor_config_.similarity_threshold_);
        descriptor_config_.icp_threshold_ = ReadOr<float>(
            config, "search_internal_icp_threshold", descriptor_config_.icp_threshold_);
        if (descriptor_config_.candidate_num_ <= 0 ||
            descriptor_config_.candidate_min_votes_ <= 0 ||
            descriptor_config_.rough_dis_threshold_ <= 0.0F ||
            descriptor_config_.similarity_threshold_ < 0.0F ||
            descriptor_config_.similarity_threshold_ > 1.0F ||
            descriptor_config_.icp_threshold_ < 0.0F) {
            LOG(ERROR) << "invalid BTC relocalization search configuration";
            return false;
        }
        manager_ = STDescManager(descriptor_config_);

        const YAML::Node yaml_entries = manifest["entries"];
        if (!yaml_entries || !yaml_entries.IsSequence() || yaml_entries.size() == 0) {
            LOG(ERROR) << "BTC relocalization database has no entries: " << manifest_path;
            return false;
        }
        entries_.reserve(yaml_entries.size());
        for (std::size_t index = 0; index < yaml_entries.size(); ++index) {
            const YAML::Node yaml_entry = yaml_entries[index];
            const int descriptor_id = ReadOr<int>(yaml_entry, "descriptor_id", -1);
            if (descriptor_id != static_cast<int>(index)) {
                LOG(ERROR) << "BTC descriptor ids must be dense; entry=" << index
                           << ", id=" << descriptor_id;
                return false;
            }
            if (!yaml_entry["cloud"]) {
                LOG(ERROR) << "BTC database entry " << index << " is missing its cloud";
                return false;
            }
            const std::filesystem::path cloud_path =
                database_directory / yaml_entry["cloud"].as<std::string>();
            pcl::PointCloud<pcl::PointXYZI>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZI>);
            if (pcl::io::loadPCDFile(cloud_path.string(), *cloud) != 0 || cloud->empty()) {
                LOG(ERROR) << "failed to load BTC relocalization submap: " << cloud_path;
                return false;
            }
            std::vector<STD> descriptors;
            manager_.GenerateSTDescs(cloud, descriptors, descriptor_id);
            manager_.AddSTDescs(descriptors);
            entries_.push_back(DatabaseEntry{descriptor_id, ReadPose(yaml_entry)});
        }
        // This is a previous-session database, so no entry is temporally
        // adjacent to the live query.  Keep the descriptor frame id equal to
        // STDescManager's next frame id and disable the online near-frame gate,
        // matching Voxel-SLAM's multi-session relocalization convention.
        manager_.config_setting_.skip_near_num_ = -1;

        ready_ = !entries_.empty();
        LOG(INFO) << "loaded BTC relocalization database: entries=" << entries_.size()
                  << ", max_query_submap_size=" << options_.query_submap_size
                  << ", query_stride=" << options_.query_stride
                  << ", top_k=" << options_.top_k
                  << ", min_score=" << options_.min_btc_score
                  << ", search_candidates=" << descriptor_config_.candidate_num_
                  << ", search_min_votes=" << descriptor_config_.candidate_min_votes_
                  << ", search_rough_distance=" << descriptor_config_.rough_dis_threshold_
                  << ", search_similarity=" << descriptor_config_.similarity_threshold_
                  << ", search_internal_icp=" << descriptor_config_.icp_threshold_
                  << ", path=" << database_directory;
        return ready_;
    } catch (const std::exception& error) {
        LOG(ERROR) << "failed to initialize BTC relocalization: " << error.what();
        ready_ = false;
        entries_.clear();
        return false;
    }
}

void BtcRelocalizer::ResetQuery() {
    query_frames_.clear();
    frames_since_query_ = 0;
}

pcl::PointCloud<pcl::PointXYZI>::Ptr BtcRelocalizer::BuildQuerySubmap(
    std::size_t frame_count) const {
    pcl::PointCloud<pcl::PointXYZI>::Ptr combined(new pcl::PointCloud<pcl::PointXYZI>);
    if (query_frames_.empty() || frame_count == 0) return combined;
    frame_count = std::min(frame_count, query_frames_.size());
    const auto first = query_frames_.end() - static_cast<std::ptrdiff_t>(frame_count);

    std::size_t available_points = 0;
    for (auto frame = first; frame != query_frames_.end(); ++frame) {
        if (frame->cloud) available_points += frame->cloud->size();
    }
    const std::size_t maximum = static_cast<std::size_t>(options_.max_points_per_submap);
    const std::size_t stride = std::max<std::size_t>(1, (available_points + maximum - 1) / maximum);
    combined->reserve(std::min(available_points, maximum));

    const SE3 T_current_lidar_odom = query_frames_.back().T_odom_lidar.inverse();
    std::size_t source_index = 0;
    for (auto frame = first; frame != query_frames_.end(); ++frame) {
        if (!frame->cloud) continue;
        const SE3 T_current_lidar_frame_lidar = T_current_lidar_odom * frame->T_odom_lidar;
        for (const auto& source : frame->cloud->points) {
            if (source_index++ % stride != 0) continue;
            if (!std::isfinite(source.x) || !std::isfinite(source.y) || !std::isfinite(source.z)) continue;
            const Vec3d transformed = T_current_lidar_frame_lidar * Vec3d(source.x, source.y, source.z);
            pcl::PointXYZI point;
            point.x = static_cast<float>(transformed.x());
            point.y = static_cast<float>(transformed.y());
            point.z = static_cast<float>(transformed.z());
            point.intensity = source.intensity;
            combined->push_back(point);
        }
    }

    if (options_.downsample_leaf_size <= 0.0 || combined->empty()) return combined;
    pcl::PointCloud<pcl::PointXYZI>::Ptr filtered(new pcl::PointCloud<pcl::PointXYZI>);
    pcl::VoxelGrid<pcl::PointXYZI> voxel;
    const float leaf = static_cast<float>(options_.downsample_leaf_size);
    voxel.setLeafSize(leaf, leaf, leaf);
    voxel.setInputCloud(combined);
    voxel.filter(*filtered);
    return filtered;
}

std::optional<BtcRelocalizationResult> BtcRelocalizer::AddFrame(
    const CloudPtr& cloud, const SE3& T_odom_imu, double timestamp) {
    if (!ready_ || !cloud || cloud->empty() || !std::isfinite(timestamp)) return std::nullopt;

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
    if (query_frames_.size() < minimum_window) {
        return std::nullopt;
    }
    const bool first_query =
        query_frames_.size() == minimum_window && frames_since_query_ >= minimum_window;
    if (!first_query &&
        frames_since_query_ < static_cast<std::size_t>(options_.query_stride)) {
        return std::nullopt;
    }
    frames_since_query_ = 0;

    BtcRelocalizationResult result;
    result.attempted = true;
    result.timestamp = query_frames_.back().timestamp;
    const auto search_begin = std::chrono::steady_clock::now();
    bool generated_descriptors = false;
    bool raw_candidate_found = false;
    std::vector<BtcRelocalizationCandidate> hypotheses;
    for (const int configured_size : options_.query_submap_sizes) {
        const std::size_t window = static_cast<std::size_t>(configured_size);
        if (window > query_frames_.size()) continue;
        auto query_cloud = BuildQuerySubmap(window);
        result.point_count = std::max(result.point_count, query_cloud->size());
        if (query_cloud->size() < static_cast<std::size_t>(options_.min_points_per_submap)) continue;

        STDescManager query_manager(descriptor_config_);
        std::vector<STD> descriptors;
        const int query_id = static_cast<int>(manager_.current_frame_id_);
        query_manager.GenerateSTDescs(query_cloud, descriptors, query_id);
        result.descriptor_count = std::max(result.descriptor_count, descriptors.size());
        generated_descriptors = generated_descriptors || !descriptors.empty();
        if (descriptors.empty() || query_manager.plane_cloud_vec_.empty()) continue;

        const auto candidates = manager_.SearchLoopTopK(
            descriptors, query_manager.plane_cloud_vec_.back(),
            static_cast<std::size_t>(options_.top_k));
        raw_candidate_found = raw_candidate_found || !candidates.empty();
        for (const auto& candidate : candidates) {
            if (candidate.candidate_id < 0 ||
                candidate.candidate_id >= static_cast<int>(entries_.size()) ||
                !std::isfinite(candidate.score) || candidate.score < options_.min_btc_score) {
                continue;
            }
            const SE3 T_candidate_lidar_current_lidar(
                Quatd(candidate.transform.second).normalized(), candidate.transform.first);
            const SE3 T_world_current_lidar =
                entries_[candidate.candidate_id].T_world_lidar *
                T_candidate_lidar_current_lidar;
            BtcRelocalizationCandidate hypothesis;
            hypothesis.candidate_id = candidate.candidate_id;
            hypothesis.score = candidate.score;
            hypothesis.query_submap_size = configured_size;
            hypothesis.point_count = query_cloud->size();
            hypothesis.descriptor_count = descriptors.size();
            hypothesis.rough_match_count = candidate.rough_match_count;
            hypothesis.spatial_coverage = candidate.spatial_coverage;
            hypothesis.T_world_imu = T_world_current_lidar * T_imu_lidar_.inverse();
            hypotheses.push_back(std::move(hypothesis));
        }
    }
    result.search_time_ms = std::chrono::duration<double, std::milli>(
                                std::chrono::steady_clock::now() - search_begin)
                                .count();
    result.candidate_found = raw_candidate_found;

    std::stable_sort(hypotheses.begin(), hypotheses.end(), [](const auto& left, const auto& right) {
        if (left.score != right.score) return left.score > right.score;
        if (left.rough_match_count != right.rough_match_count) {
            return left.rough_match_count > right.rough_match_count;
        }
        return left.candidate_id < right.candidate_id;
    });
    const double maximum_yaw_distance =
        options_.candidate_dedup_yaw_deg * M_PI / 180.0;
    for (const auto& hypothesis : hypotheses) {
        const bool duplicate = std::any_of(
            result.candidates.begin(), result.candidates.end(), [&](const auto& selected) {
                return (selected.T_world_imu.translation() -
                        hypothesis.T_world_imu.translation())
                               .norm() <= options_.candidate_dedup_radius &&
                       AngleDistance(Yaw(selected.T_world_imu), Yaw(hypothesis.T_world_imu)) <=
                           maximum_yaw_distance;
            });
        if (!duplicate) result.candidates.push_back(hypothesis);
        if (result.candidates.size() >= static_cast<std::size_t>(options_.top_k)) break;
    }

    if (result.candidates.empty()) {
        result.reason = raw_candidate_found ? "score_below_threshold"
                                            : (generated_descriptors ? "no_candidate"
                                                                     : "no_query_descriptors");
        return result;
    }

    const auto& best = result.candidates.front();
    result.candidate_found = true;
    result.candidate_id = best.candidate_id;
    result.score = best.score;
    result.point_count = best.point_count;
    result.descriptor_count = best.descriptor_count;
    result.T_world_imu = best.T_world_imu;
    result.accepted = true;
    result.reason = "accepted";
    return result;
}

}  // namespace lightning::loc
