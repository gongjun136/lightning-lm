#include "core/localization/btc_relocalizer.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <stdexcept>

#include <glog/logging.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <yaml-cpp/yaml.h>

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

}  // namespace

BtcRelocalizer::BtcRelocalizer() : manager_(descriptor_config_) {}

bool BtcRelocalizer::Init(const std::string& config_path, const std::string& map_path,
                          const SE3& T_imu_lidar) {
    ready_ = false;
    entries_.clear();
    query_frames_.clear();
    T_imu_lidar_ = T_imu_lidar;

    try {
        const YAML::Node root = YAML::LoadFile(config_path);
        const YAML::Node config = root["relocalization"];
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
        if (ReadOr<int>(manifest, "schema_version", 0) != 1) {
            LOG(ERROR) << "unsupported BTC relocalization database schema: " << manifest_path;
            return false;
        }
        options_.query_submap_size = ReadOr<int>(
            config, "query_submap_size", ReadOr<int>(manifest, "descriptor_submap_size", 10));
        options_.max_points_per_submap = ReadOr<int>(
            config, "max_points_per_submap", ReadOr<int>(manifest, "max_points_per_submap", 50000));
        options_.min_points_per_submap =
            ReadOr<int>(config, "min_points_per_submap", options_.min_points_per_submap);
        options_.downsample_leaf_size = ReadOr<double>(
            config, "downsample_leaf_size", ReadOr<double>(manifest, "downsample_leaf_size", 0.20));
        if (options_.query_submap_size <= 0 || options_.max_points_per_submap <= 0 ||
            options_.min_points_per_submap <= 0 || options_.min_btc_score < 0.0) {
            LOG(ERROR) << "invalid BTC relocalization configuration";
            return false;
        }

        descriptor_config_ = ConfigSetting{};
        ReadDescriptorConfig(manifest["descriptor"], descriptor_config_);
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
                  << ", query_submap_size=" << options_.query_submap_size
                  << ", min_score=" << options_.min_btc_score
                  << ", path=" << database_directory;
        return ready_;
    } catch (const std::exception& error) {
        LOG(ERROR) << "failed to initialize BTC relocalization: " << error.what();
        ready_ = false;
        entries_.clear();
        return false;
    }
}

void BtcRelocalizer::ResetQuery() { query_frames_.clear(); }

pcl::PointCloud<pcl::PointXYZI>::Ptr BtcRelocalizer::BuildQuerySubmap() const {
    pcl::PointCloud<pcl::PointXYZI>::Ptr combined(new pcl::PointCloud<pcl::PointXYZI>);
    if (query_frames_.empty()) return combined;

    std::size_t available_points = 0;
    for (const auto& frame : query_frames_) {
        if (frame.cloud) available_points += frame.cloud->size();
    }
    const std::size_t maximum = static_cast<std::size_t>(options_.max_points_per_submap);
    const std::size_t stride = std::max<std::size_t>(1, (available_points + maximum - 1) / maximum);
    combined->reserve(std::min(available_points, maximum));

    const SE3 T_current_lidar_odom = query_frames_.back().T_odom_lidar.inverse();
    std::size_t source_index = 0;
    for (const auto& frame : query_frames_) {
        if (!frame.cloud) continue;
        const SE3 T_current_lidar_frame_lidar = T_current_lidar_odom * frame.T_odom_lidar;
        for (const auto& source : frame.cloud->points) {
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
    if (query_frames_.size() < static_cast<std::size_t>(options_.query_submap_size)) {
        return std::nullopt;
    }

    BtcRelocalizationResult result;
    result.attempted = true;
    result.timestamp = query_frames_.back().timestamp;
    auto query_cloud = BuildQuerySubmap();
    result.point_count = query_cloud->size();
    query_frames_.clear();
    if (query_cloud->size() < static_cast<std::size_t>(options_.min_points_per_submap)) {
        result.reason = "too_few_query_points";
        return result;
    }

    STDescManager query_manager(descriptor_config_);
    std::vector<STD> descriptors;
    const int query_id = static_cast<int>(manager_.current_frame_id_);
    query_manager.GenerateSTDescs(query_cloud, descriptors, query_id);
    result.descriptor_count = descriptors.size();
    if (descriptors.empty() || query_manager.plane_cloud_vec_.empty()) {
        result.reason = "no_query_descriptors";
        return result;
    }

    std::pair<int, double> search_result(-1, 0.0);
    std::pair<Eigen::Vector3d, Eigen::Matrix3d> transform;
    std::vector<std::pair<STD, STD>> matches;
    manager_.SearchLoop(descriptors, search_result, transform, matches,
                        query_manager.plane_cloud_vec_.back());
    if (search_result.first < 0 || search_result.first >= static_cast<int>(entries_.size())) {
        result.reason = "no_candidate";
        return result;
    }

    result.candidate_found = true;
    result.candidate_id = search_result.first;
    result.score = search_result.second;
    if (!std::isfinite(result.score) || result.score < options_.min_btc_score) {
        result.reason = "score_below_threshold";
        return result;
    }

    const SE3 T_candidate_lidar_current_lidar(
        Quatd(transform.second).normalized(), transform.first);
    const SE3 T_world_current_lidar =
        entries_[result.candidate_id].T_world_lidar * T_candidate_lidar_current_lidar;
    result.T_world_imu = T_world_current_lidar * T_imu_lidar_.inverse();
    result.accepted = true;
    result.reason = "accepted";
    return result;
}

}  // namespace lightning::loc
