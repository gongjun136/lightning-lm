#pragma once

#include <cstddef>
#include <cstdint>
#include <deque>
#include <iosfwd>
#include <mutex>
#include <string>

#include <builtin_interfaces/msg/time.hpp>
#include <geometry_msgs/msg/pose_stamped.hpp>
#include <geometry_msgs/msg/transform_stamped.hpp>
#include <nav_msgs/msg/path.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <yaml-cpp/yaml.h>

#include "common/eigen_types.h"
#include "common/point_def.h"
#include "core/localization/localization_result.h"
#include "geosun_msgs/msg/pos_res.hpp"
#include "lightning/msg/fault_status.hpp"
#include "lightning/msg/localization_status.hpp"
#include "lightning/msg/vehicle_pose.hpp"

namespace lightning::sany_output {

struct FixedMapTransform {
    bool enabled = false;
    std::string source_frame = "localization_map";
    std::string target_frame = "map";
    SE3 T_target_localization;
};

/// Read output.fixed_map_transform.  The stored transform convention is
/// T_target_localization and is left-multiplied only at the online ROS output boundary.
bool LoadFixedMapTransform(const YAML::Node& root, FixedMapTransform& transform,
                           std::string& error);

SE3 TransformPoseForOutput(const SE3& localization_pose,
                           const FixedMapTransform& transform);

geometry_msgs::msg::TransformStamped TransformTfForOutput(
    const geometry_msgs::msg::TransformStamped& localization_tf,
    const FixedMapTransform& transform);

geosun_msgs::msg::PosRes MakePosResMessage(const SE3& map_rear_axle_pose, double vehicle_speed, double stamp,
                                           const std::string& frame_id);

geometry_msgs::msg::PoseStamped MakePoseMessage(const geosun_msgs::msg::PosRes& position);
lightning::msg::VehiclePose MakeVehiclePoseMessage(
    const geosun_msgs::msg::PosRes& position);

/// Append one strictly monotonic PoseStamped sample in TUM format. The caller
/// owns stream lifetime and flushing.
bool WriteTumPoseLine(std::ostream& stream,
                      const geometry_msgs::msg::PoseStamped& pose,
                      double& last_timestamp);

sensor_msgs::msg::PointCloud2 MakeCloudMessage(const CloudPtr& cloud, double begin_time, double end_time,
                                               const SE3& T_output_lidar, const std::string& frame_id);

SE3 MakeLivoxLidarTransform(const SO3& initial_lidar_rotation);
SE3 MakeMapLivoxPose(const SE3& map_lidar_pose, const SO3& initial_lidar_rotation);
SE3 MakeRearAxleLidarTransform(const SO3& initial_lidar_rotation,
                               const Vec3d& primary_lidar_position_in_body);
SE3 MakeMapRearAxlePose(const SE3& map_lidar_pose, const SO3& initial_lidar_rotation,
                        const Vec3d& primary_lidar_position_in_body);

class LocalizationPublicationGate {
   public:
    explicit LocalizationPublicationGate(std::size_t lost_frame_threshold = 5);

    void SetLostFrameThreshold(std::size_t lost_frame_threshold);
    void SetMaxLidarMatchAge(double max_age_sec);
    void ObserveLidarMatch(bool valid, double sensor_stamp = 0.0);
    bool PoseOutputsEnabled() const;
    bool MapOutputsEnabled(double current_sensor_stamp = 0.0) const;
    bool LidarMatchStale(double current_sensor_stamp) const;
    double LidarMatchAgeSec(double current_sensor_stamp) const;
    double LastLidarMatchStamp() const;
    double LastValidLidarMatchStamp() const;
    std::size_t ConsecutiveLostFrames() const;

   private:
    mutable std::mutex mutex_;
    std::size_t lost_frame_threshold_ = 5;
    std::size_t consecutive_lost_frames_ = 0;
    bool has_valid_match_ = false;
    double max_lidar_match_age_sec_ = 0.0;
    double last_lidar_match_stamp_ = 0.0;
    double last_valid_lidar_match_stamp_ = 0.0;
    mutable bool stale_latched_ = false;
};

class FrameDecimator {
   public:
    explicit FrameDecimator(std::size_t interval) : interval_(interval) {}

    bool Tick() { return interval_ > 0 && ++frame_count_ % interval_ == 0; }

   private:
    std::size_t interval_;
    std::size_t frame_count_ = 0;
};

enum class LocalizationFaultType : std::int32_t {
    NONE = 0,
    LOCALIZATION_DEGRADED = 1,
    LOCALIZATION_LOST = 2,
};

class LocalizationTelemetryState {
   public:
    explicit LocalizationTelemetryState(std::size_t lost_frame_threshold = 5,
                                        std::size_t path_capacity = 500,
                                        double path_sample_period = 0.1);

    void Start();
    void ObserveLocalization(loc::LocalizationStatus status,
                             std::size_t consecutive_lost_frames);
    void ObserveLocalizationStale(bool stale);
    void ObservePose(const geometry_msgs::msg::PoseStamped& pose);

    lightning::msg::LocalizationStatus MakeLocalizationStatus(
        const builtin_interfaces::msg::Time& stamp) const;
    lightning::msg::FaultStatus MakeFaultStatus(
        const builtin_interfaces::msg::Time& stamp) const;
    nav_msgs::msg::Path MakePath(const builtin_interfaces::msg::Time& stamp) const;

    bool OfflineHealthPublishDue(double sensor_time);
    bool OfflinePathPublishDue(double sensor_time, double interval = 2.0);
    std::size_t PathSize() const;

   private:
    mutable std::mutex mutex_;
    std::size_t lost_frame_threshold_ = 5;
    std::size_t path_capacity_ = 500;
    double path_sample_period_ = 0.1;
    std::uint8_t current_status_ = lightning::msg::LocalizationStatus::STATUS_IDLE;
    bool has_good_localization_ = false;
    bool localization_lost_latched_ = false;
    std::deque<geometry_msgs::msg::PoseStamped> path_poses_;
    double last_path_sample_time_ = 0.0;
    double last_offline_health_publish_time_ = 0.0;
    double last_offline_path_publish_time_ = 0.0;
};

}  // namespace lightning::sany_output
