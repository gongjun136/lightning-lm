#pragma once

#include <cstddef>
#include <chrono>
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
#include "lightning/msg/fault_status.hpp"
#include "lightning/msg/localization_status.hpp"
#include "lightning/msg/vehicle_pose.hpp"

namespace lightning::sany_output {

struct FixedMapTransform {
    bool enabled = false;  ///< 是否启用业务输出地图的固定坐标变换。
    std::string source_frame = "localization_map";  ///< 固定变换的源定位地图坐标系名称。
    std::string target_frame = "map";  ///< 固定变换的目标业务地图坐标系名称。
    SE3 T_target_localization;  ///< 将定位地图坐标变换到目标业务地图坐标的固定变换。
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

lightning::msg::VehiclePose MakeVehiclePoseMessage(const SE3& map_rear_axle_pose, double vehicle_speed, double stamp,
                                           const std::string& frame_id);

geometry_msgs::msg::PoseStamped MakePoseMessage(const lightning::msg::VehiclePose& position);

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

/** @see @ref output_contracts */
class LocalizationPublicationGate {
   public:
    using Clock = std::chrono::steady_clock;
    explicit LocalizationPublicationGate(std::size_t lost_frame_threshold = 5);

    void SetLostFrameThreshold(std::size_t lost_frame_threshold);
    void SetMaxLidarMatchAge(double max_age_sec);
    void ObserveLidarMatch(bool valid, double sensor_stamp = 0.0,
                           Clock::time_point now = Clock::now());
    bool PoseOutputsEnabled(double current_sensor_stamp = 0.0,
                            Clock::time_point now = Clock::now()) const;
    bool MapOutputsEnabled(double current_sensor_stamp = 0.0,
                           Clock::time_point now = Clock::now()) const;
    bool LidarMatchStale(double current_sensor_stamp,
                        Clock::time_point now = Clock::now()) const;
    double LidarMatchAgeSec(double current_sensor_stamp) const;
    double LastLidarMatchStamp() const;
    double LastValidLidarMatchStamp() const;
    std::size_t ConsecutiveLostFrames() const;

   private:
    bool MatchStaleLocked(double current_sensor_stamp, Clock::time_point now) const;
    mutable std::mutex mutex_;  ///< 保护匹配时效、连续丢帧与输出门控状态。
    std::size_t lost_frame_threshold_ = 5;  ///< 判定定位失锁所需的连续无效匹配帧数。
    std::size_t consecutive_lost_frames_ = 0;  ///< 当前连续无效地图匹配的帧数。
    bool has_valid_match_ = false;  ///< 是否已观察到有效地图匹配。
    double max_lidar_match_age_sec_ = 0.0;  ///< 有效地图匹配允许的最大年龄（s）。
    double last_lidar_match_stamp_ = 0.0;  ///< 最近一次地图匹配输入的时间戳（s）。
    double last_valid_lidar_match_stamp_ = 0.0;  ///< 最近有效地图匹配的传感器时间戳（s）。
    mutable double latest_sensor_stamp_ = 0.0;  ///< 门控已观察到的最新传感器时间戳（s）。
    Clock::time_point last_valid_match_arrival_{};  ///< 最近有效匹配到达的单调时钟时刻。
    mutable bool stale_latched_ = false;  ///< 地图匹配过期后保持的失效状态。
};

class FrameDecimator {
   public:
    explicit FrameDecimator(std::size_t interval) : interval_(interval) {}

    bool Tick() { return interval_ > 0 && ++frame_count_ % interval_ == 0; }

   private:
    std::size_t interval_;  ///< 地图点云降频发布的帧间隔。
    std::size_t frame_count_ = 0;  ///< 已经调用抽帧器的帧数。
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
    mutable std::mutex mutex_;  ///< 保护健康状态与可视化轨迹缓存。
    std::size_t lost_frame_threshold_ = 5;  ///< 判定定位失锁所需的连续无效匹配帧数。
    std::size_t path_capacity_ = 500;  ///< 可视化轨迹最多保留的位姿数。
    double path_sample_period_ = 0.1;  ///< 可视化轨迹的采样间隔（s）。
    std::uint8_t current_status_ = lightning::msg::LocalizationStatus::STATUS_IDLE;  ///< 最近一次定位模块报告的运行状态。
    bool has_good_localization_ = false;  ///< 定位是否曾达到有效状态。
    bool localization_lost_latched_ = false;  ///< 失锁后保持的故障状态，供健康消息使用。
    std::deque<geometry_msgs::msg::PoseStamped> path_poses_;  ///< 按时间保存的后轴可视化位姿路径。
    double last_path_sample_time_ = 0.0;  ///< 最近一次路径采样的传感器时间戳（s）。
    double last_offline_health_publish_time_ = 0.0;  ///< 离线诊断最近一次发布的传感器时间戳（s）。
    double last_offline_path_publish_time_ = 0.0;  ///< 离线路径最近一次发布的传感器时间戳（s）。
};

}  // namespace lightning::sany_output
