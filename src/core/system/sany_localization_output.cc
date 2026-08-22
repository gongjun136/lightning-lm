#include "core/system/sany_localization_output.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

#include <sensor_msgs/point_cloud2_iterator.hpp>
#include <yaml-cpp/yaml.h>

namespace lightning::sany_output {
namespace {

constexpr double kDegreesPerRadian = 180.0 / M_PI;
constexpr double kRadiansPerDegree = M_PI / 180.0;

builtin_interfaces::msg::Time ToRosStamp(double seconds) {
    const std::int64_t nanoseconds = static_cast<std::int64_t>(std::llround(seconds * 1e9));
    builtin_interfaces::msg::Time stamp;
    stamp.sec = static_cast<std::int32_t>(nanoseconds / 1000000000LL);
    stamp.nanosec = static_cast<std::uint32_t>(nanoseconds % 1000000000LL);
    return stamp;
}

double ToSec(const builtin_interfaces::msg::Time& stamp) {
    return static_cast<double>(stamp.sec) + static_cast<double>(stamp.nanosec) * 1e-9;
}

}  // namespace

bool LoadFixedMapTransform(const YAML::Node& root, FixedMapTransform& transform,
                           std::string& error) {
    transform = FixedMapTransform{};
    error.clear();
    try {
        const YAML::Node node = root["output"] ? root["output"]["fixed_map_transform"]
                                                 : YAML::Node();
        if (!node) return true;
        transform.enabled = node["enabled"] ? node["enabled"].as<bool>() : false;
        if (!transform.enabled) return true;

        const std::string convention =
            node["convention"] ? node["convention"].as<std::string>() : "";
        if (convention != "target_from_localization") {
            error = "output.fixed_map_transform.convention must be target_from_localization";
            return false;
        }
        transform.source_frame = node["source_frame"]
                                     ? node["source_frame"].as<std::string>()
                                     : "localization_map";
        transform.target_frame = node["target_frame"]
                                     ? node["target_frame"].as<std::string>()
                                     : "map";
        const auto translation = node["translation_xyz"].as<std::vector<double>>();
        const auto quaternion = node["quaternion_xyzw"].as<std::vector<double>>();
        if (translation.size() != 3 || quaternion.size() != 4) {
            error = "output.fixed_map_transform requires 3D translation_xyz and XYZW quaternion_xyzw";
            return false;
        }
        for (const double value : translation) {
            if (!std::isfinite(value)) {
                error = "output.fixed_map_transform translation contains a non-finite value";
                return false;
            }
        }
        for (const double value : quaternion) {
            if (!std::isfinite(value)) {
                error = "output.fixed_map_transform quaternion contains a non-finite value";
                return false;
            }
        }
        Quatd q(quaternion[3], quaternion[0], quaternion[1], quaternion[2]);
        if (q.norm() < 1e-9) {
            error = "output.fixed_map_transform quaternion has zero norm";
            return false;
        }
        q.normalize();
        transform.T_target_localization =
            SE3(q, Vec3d(translation[0], translation[1], translation[2]));
        return true;
    } catch (const std::exception& exception) {
        error = std::string("invalid output.fixed_map_transform: ") + exception.what();
        return false;
    }
}

SE3 TransformPoseForOutput(const SE3& localization_pose,
                           const FixedMapTransform& transform) {
    return transform.enabled ? transform.T_target_localization * localization_pose
                             : localization_pose;
}

geometry_msgs::msg::TransformStamped TransformTfForOutput(
    const geometry_msgs::msg::TransformStamped& localization_tf,
    const FixedMapTransform& transform) {
    if (!transform.enabled) return localization_tf;
    Quatd q(localization_tf.transform.rotation.w, localization_tf.transform.rotation.x,
            localization_tf.transform.rotation.y, localization_tf.transform.rotation.z);
    if (q.norm() < 1e-9) q = Quatd::Identity();
    q.normalize();
    const SE3 localization_pose(
        q, Vec3d(localization_tf.transform.translation.x,
                 localization_tf.transform.translation.y,
                 localization_tf.transform.translation.z));
    const SE3 output_pose = TransformPoseForOutput(localization_pose, transform);
    geometry_msgs::msg::TransformStamped output = localization_tf;
    output.header.frame_id = transform.target_frame;
    output.transform.translation.x = output_pose.translation().x();
    output.transform.translation.y = output_pose.translation().y();
    output.transform.translation.z = output_pose.translation().z();
    const Quatd output_q = output_pose.unit_quaternion();
    output.transform.rotation.x = output_q.x();
    output.transform.rotation.y = output_q.y();
    output.transform.rotation.z = output_q.z();
    output.transform.rotation.w = output_q.w();
    return output;
}

geosun_msgs::msg::PosRes MakePosResMessage(const SE3& map_rear_axle_pose, double vehicle_speed, double stamp,
                                           const std::string& frame_id) {
    geosun_msgs::msg::PosRes message;
    message.header.stamp = ToRosStamp(stamp);
    message.header.frame_id = frame_id;
    message.f8enh = {map_rear_axle_pose.translation().x(), map_rear_axle_pose.translation().y(),
                     map_rear_axle_pose.translation().z()};
    message.f8vehiclespeed = vehicle_speed;

    const Mat3d rotation = map_rear_axle_pose.so3().matrix();
    const double roll = std::atan2(rotation(2, 1), rotation(2, 2));
    const double pitch = std::asin(std::clamp(-rotation(2, 0), -1.0, 1.0));
    double yaw = std::atan2(rotation(1, 0), rotation(0, 0));
    if (yaw < 0.0) yaw += 2.0 * M_PI;
    message.f8pry = {roll * kDegreesPerRadian, pitch * kDegreesPerRadian, yaw * kDegreesPerRadian};
    return message;
}

geometry_msgs::msg::PoseStamped MakePoseMessage(const geosun_msgs::msg::PosRes& position) {
    geometry_msgs::msg::PoseStamped message;
    message.header = position.header;
    message.pose.position.x = position.f8enh[0];
    message.pose.position.y = position.f8enh[1];
    message.pose.position.z = position.f8enh[2];

    const double roll = position.f8pry[0] * kRadiansPerDegree;
    const double pitch = position.f8pry[1] * kRadiansPerDegree;
    const double yaw = position.f8pry[2] * kRadiansPerDegree;
    const Quatd quaternion =
        Eigen::AngleAxisd(yaw, Vec3d::UnitZ()) * Eigen::AngleAxisd(pitch, Vec3d::UnitY()) *
        Eigen::AngleAxisd(roll, Vec3d::UnitX());
    message.pose.orientation.x = quaternion.x();
    message.pose.orientation.y = quaternion.y();
    message.pose.orientation.z = quaternion.z();
    message.pose.orientation.w = quaternion.w();
    return message;
}

lightning::msg::VehiclePose MakeVehiclePoseMessage(
    const geosun_msgs::msg::PosRes& position) {
    lightning::msg::VehiclePose message;
    message.header = position.header;
    message.x = position.f8enh[0];
    message.y = position.f8enh[1];
    message.z = position.f8enh[2];
    message.roll = position.f8pry[0];
    message.pitch = position.f8pry[1];
    message.yaw = position.f8pry[2];
    message.speed = position.f8vehiclespeed;
    return message;
}

sensor_msgs::msg::PointCloud2 MakeCloudMessage(const CloudPtr& cloud, double begin_time, double end_time,
                                               const SE3& T_output_lidar, const std::string& frame_id) {
    sensor_msgs::msg::PointCloud2 message;
    message.header.stamp = ToRosStamp(end_time);
    message.header.frame_id = frame_id;
    sensor_msgs::PointCloud2Modifier modifier(message);
    modifier.setPointCloud2Fields(7, "x", 1, sensor_msgs::msg::PointField::FLOAT32, "y", 1,
                                  sensor_msgs::msg::PointField::FLOAT32, "z", 1,
                                  sensor_msgs::msg::PointField::FLOAT32, "intensity", 1,
                                  sensor_msgs::msg::PointField::FLOAT32, "tag", 1,
                                  sensor_msgs::msg::PointField::UINT8, "line", 1,
                                  sensor_msgs::msg::PointField::UINT8, "timestamp", 1,
                                  sensor_msgs::msg::PointField::FLOAT64);
    const std::size_t size = cloud ? cloud->size() : 0;
    modifier.resize(size);
    message.is_dense = cloud ? cloud->is_dense : true;

    sensor_msgs::PointCloud2Iterator<float> x(message, "x"), y(message, "y"), z(message, "z"),
        intensity(message, "intensity");
    sensor_msgs::PointCloud2Iterator<std::uint8_t> tag(message, "tag"), line(message, "line");
    sensor_msgs::PointCloud2Iterator<double> timestamp(message, "timestamp");
    if (cloud) {
        for (const auto& point : cloud->points) {
            const Vec3d transformed = T_output_lidar * point.getVector3fMap().cast<double>();
            *x = static_cast<float>(transformed.x());
            *y = static_cast<float>(transformed.y());
            *z = static_cast<float>(transformed.z());
            *intensity = point.intensity;
            *tag = 0;
            *line = point.lidar_id;
            *timestamp = begin_time + point.time * 1e-3;
            ++x;
            ++y;
            ++z;
            ++intensity;
            ++tag;
            ++line;
            ++timestamp;
        }
    }
    return message;
}

SE3 MakeLivoxLidarTransform(const SO3& initial_lidar_rotation) {
    return SE3(initial_lidar_rotation, Vec3d::Zero());
}

SE3 MakeMapLivoxPose(const SE3& map_lidar_pose, const SO3& initial_lidar_rotation) {
    return map_lidar_pose * MakeLivoxLidarTransform(initial_lidar_rotation).inverse();
}

SE3 MakeRearAxleLidarTransform(const SO3& initial_lidar_rotation,
                               const Vec3d& primary_lidar_position_in_body) {
    const SE3 T_rear_livox(SO3(), primary_lidar_position_in_body);
    return T_rear_livox * MakeLivoxLidarTransform(initial_lidar_rotation);
}

SE3 MakeMapRearAxlePose(const SE3& map_lidar_pose, const SO3& initial_lidar_rotation,
                        const Vec3d& primary_lidar_position_in_body) {
    return map_lidar_pose *
           MakeRearAxleLidarTransform(initial_lidar_rotation, primary_lidar_position_in_body).inverse();
}

LocalizationPublicationGate::LocalizationPublicationGate(std::size_t lost_frame_threshold)
    : lost_frame_threshold_(std::max<std::size_t>(1, lost_frame_threshold)) {}

void LocalizationPublicationGate::SetLostFrameThreshold(std::size_t lost_frame_threshold) {
    std::lock_guard<std::mutex> lock(mutex_);
    lost_frame_threshold_ = std::max<std::size_t>(1, lost_frame_threshold);
    consecutive_lost_frames_ = 0;
    has_valid_match_ = false;
    last_lidar_match_stamp_ = 0.0;
    last_valid_lidar_match_stamp_ = 0.0;
    stale_latched_ = false;
}

void LocalizationPublicationGate::SetMaxLidarMatchAge(double max_age_sec) {
    std::lock_guard<std::mutex> lock(mutex_);
    max_lidar_match_age_sec_ = std::max(0.0, max_age_sec);
}

void LocalizationPublicationGate::ObserveLidarMatch(bool valid, double sensor_stamp) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (sensor_stamp > 0.0) {
        last_lidar_match_stamp_ = std::max(last_lidar_match_stamp_, sensor_stamp);
    }
    if (valid) {
        has_valid_match_ = true;
        consecutive_lost_frames_ = 0;
        stale_latched_ = false;
        if (sensor_stamp > 0.0) {
            last_valid_lidar_match_stamp_ =
                std::max(last_valid_lidar_match_stamp_, sensor_stamp);
        }
    } else if (has_valid_match_) {
        ++consecutive_lost_frames_;
    }
}

bool LocalizationPublicationGate::PoseOutputsEnabled() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return has_valid_match_ && consecutive_lost_frames_ < lost_frame_threshold_;
}

bool LocalizationPublicationGate::MapOutputsEnabled(double current_sensor_stamp) const {
    std::lock_guard<std::mutex> lock(mutex_);
    const bool stale = max_lidar_match_age_sec_ > 0.0 && current_sensor_stamp > 0.0 &&
                       last_lidar_match_stamp_ > 0.0 &&
                       current_sensor_stamp - last_lidar_match_stamp_ > max_lidar_match_age_sec_;
    if (stale) stale_latched_ = true;
    return has_valid_match_ && consecutive_lost_frames_ < lost_frame_threshold_ &&
           !stale_latched_;
}

bool LocalizationPublicationGate::LidarMatchStale(double current_sensor_stamp) const {
    std::lock_guard<std::mutex> lock(mutex_);
    const bool stale = has_valid_match_ && max_lidar_match_age_sec_ > 0.0 &&
                       current_sensor_stamp > 0.0 && last_lidar_match_stamp_ > 0.0 &&
                       current_sensor_stamp - last_lidar_match_stamp_ > max_lidar_match_age_sec_;
    if (stale) stale_latched_ = true;
    return stale;
}

double LocalizationPublicationGate::LidarMatchAgeSec(double current_sensor_stamp) const {
    std::lock_guard<std::mutex> lock(mutex_);
    if (current_sensor_stamp <= 0.0 || last_lidar_match_stamp_ <= 0.0) return -1.0;
    return std::max(0.0, current_sensor_stamp - last_lidar_match_stamp_);
}

double LocalizationPublicationGate::LastLidarMatchStamp() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return last_lidar_match_stamp_;
}

double LocalizationPublicationGate::LastValidLidarMatchStamp() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return last_valid_lidar_match_stamp_;
}

std::size_t LocalizationPublicationGate::ConsecutiveLostFrames() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return consecutive_lost_frames_;
}

LocalizationTelemetryState::LocalizationTelemetryState(std::size_t lost_frame_threshold,
                                                       std::size_t path_capacity,
                                                       double path_sample_period)
    : lost_frame_threshold_(std::max<std::size_t>(1, lost_frame_threshold)),
      path_capacity_(std::max<std::size_t>(1, path_capacity)),
      path_sample_period_(std::max(0.0, path_sample_period)) {}

void LocalizationTelemetryState::Start() {
    std::lock_guard<std::mutex> lock(mutex_);
    current_status_ = lightning::msg::LocalizationStatus::STATUS_INITIALIZING;
}

void LocalizationTelemetryState::ObserveLocalization(loc::LocalizationStatus status,
                                                      std::size_t consecutive_lost_frames) {
    std::lock_guard<std::mutex> lock(mutex_);
    current_status_ = static_cast<std::uint8_t>(status);
    if (status == loc::LocalizationStatus::GOOD) {
        has_good_localization_ = true;
        localization_lost_latched_ = false;
        return;
    }
    if (status == loc::LocalizationStatus::FAIL ||
        (has_good_localization_ && consecutive_lost_frames >= lost_frame_threshold_)) {
        localization_lost_latched_ = true;
    }
}

void LocalizationTelemetryState::ObserveLocalizationStale(bool stale) {
    if (!stale) return;
    std::lock_guard<std::mutex> lock(mutex_);
    if (!has_good_localization_ || localization_lost_latched_) return;
    current_status_ = lightning::msg::LocalizationStatus::STATUS_FOLLOWING_DR;
}

void LocalizationTelemetryState::ObservePose(const geometry_msgs::msg::PoseStamped& pose) {
    const double timestamp = ToSec(pose.header.stamp);
    if (timestamp <= 0.0) return;

    std::lock_guard<std::mutex> lock(mutex_);
    if (current_status_ != lightning::msg::LocalizationStatus::STATUS_GOOD ||
        localization_lost_latched_) {
        return;
    }
    constexpr double kTimestampTolerance = 1e-6;
    if (last_path_sample_time_ > 0.0 &&
        timestamp - last_path_sample_time_ < path_sample_period_ - kTimestampTolerance) {
        return;
    }
    path_poses_.push_back(pose);
    last_path_sample_time_ = timestamp;
    while (path_poses_.size() > path_capacity_) path_poses_.pop_front();
}

lightning::msg::LocalizationStatus LocalizationTelemetryState::MakeLocalizationStatus(
    const builtin_interfaces::msg::Time& stamp) const {
    std::lock_guard<std::mutex> lock(mutex_);
    lightning::msg::LocalizationStatus message;
    message.header.stamp = stamp;
    message.status = current_status_;
    return message;
}

lightning::msg::FaultStatus LocalizationTelemetryState::MakeFaultStatus(
    const builtin_interfaces::msg::Time& stamp) const {
    std::lock_guard<std::mutex> lock(mutex_);
    lightning::msg::FaultStatus message;
    message.header.stamp = stamp;
    if (localization_lost_latched_) {
        message.level = lightning::msg::FaultStatus::LEVEL_P0;
        message.fault_type = static_cast<std::int32_t>(LocalizationFaultType::LOCALIZATION_LOST);
        message.description = "Localization lost; global relocalization in progress";
    } else if (current_status_ == lightning::msg::LocalizationStatus::STATUS_FOLLOWING_DR) {
        message.level = lightning::msg::FaultStatus::LEVEL_P1;
        message.fault_type = static_cast<std::int32_t>(LocalizationFaultType::LOCALIZATION_DEGRADED);
        message.description = "Localization degraded; following dead reckoning";
    } else {
        message.level = lightning::msg::FaultStatus::LEVEL_NO_FAULT;
        message.fault_type = static_cast<std::int32_t>(LocalizationFaultType::NONE);
    }
    return message;
}

nav_msgs::msg::Path LocalizationTelemetryState::MakePath(
    const builtin_interfaces::msg::Time& stamp) const {
    std::lock_guard<std::mutex> lock(mutex_);
    nav_msgs::msg::Path message;
    message.header.stamp = stamp;
    if (!path_poses_.empty()) message.header.frame_id = path_poses_.back().header.frame_id;
    message.poses.assign(path_poses_.begin(), path_poses_.end());
    return message;
}

bool LocalizationTelemetryState::OfflineHealthPublishDue(double sensor_time) {
    if (sensor_time <= 0.0) return false;
    std::lock_guard<std::mutex> lock(mutex_);
    constexpr double kPeriod = 0.1;
    constexpr double kTolerance = 1e-6;
    if (last_offline_health_publish_time_ > 0.0 &&
        sensor_time - last_offline_health_publish_time_ < kPeriod - kTolerance) {
        return false;
    }
    last_offline_health_publish_time_ = sensor_time;
    return true;
}

bool LocalizationTelemetryState::OfflinePathPublishDue(double sensor_time, double interval) {
    if (sensor_time <= 0.0 || interval <= 0.0) return false;
    std::lock_guard<std::mutex> lock(mutex_);
    if (path_poses_.empty()) return false;
    if (last_offline_path_publish_time_ <= 0.0) {
        last_offline_path_publish_time_ = sensor_time;
        return false;
    }
    constexpr double kTolerance = 1e-6;
    if (sensor_time - last_offline_path_publish_time_ < interval - kTolerance) return false;
    last_offline_path_publish_time_ = sensor_time;
    return true;
}

std::size_t LocalizationTelemetryState::PathSize() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return path_poses_.size();
}

}  // namespace lightning::sany_output
