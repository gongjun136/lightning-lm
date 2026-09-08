#pragma once

#include <atomic>
#include <chrono>
#include <cstdint>
#include <deque>
#include <mutex>
#include <map>
#include <shared_mutex>
#include <string>

#include "geometry_msgs/msg/transform_stamped.hpp"
#include "std_msgs/msg/int32.hpp"

#include "common/imu.h"
#include "common/timestamp_gate.h"
#include "core/lio/laser_mapping.h"
#include "core/localization/localization_result.h"
#include "core/system/async_message_process.h"
#include "utils/compute_profiling.h"

/// 预声明
namespace lightning {
namespace ui {
class PangolinWindow;
}

namespace loc {

class LidarLoc;
class PGO;

/**
 * 实时定位接口实现
 */
class Localization {
   public:
    struct RuntimeStats {
        std::uint64_t sensor_queue_pending = 0;
        std::uint64_t sensor_queue_dropped = 0;
        std::uint64_t sensor_queue_processed = 0;
        std::uint64_t localization_queue_pending = 0;
        std::uint64_t localization_queue_dropped = 0;
        std::uint64_t localization_queue_processed = 0;
        std::uint64_t high_frequency_queue_pending = 0;
        std::uint64_t high_frequency_queue_dropped = 0;
        std::uint64_t high_frequency_queue_processed = 0;
        double latest_enqueued_sensor_stamp = 0.0;
        double latest_processed_sensor_stamp = 0.0;
        double current_sensor_lag_sec = 0.0;
        double max_sensor_lag_sec = 0.0;
        std::uint64_t severe_timestamp_rollback_count = 0;
        double worst_timestamp_rollback_sec = 0.0;
        std::uint64_t live_output_non_monotonic_drop_count = 0;
        double worst_live_output_timestamp_rollback_sec = 0.0;
        std::uint64_t relocalization_attempt_count = 0;
        std::uint64_t relocalization_accept_count = 0;
        bool last_relocalization_candidate_found = false;
        bool last_relocalization_accepted = false;
        int last_relocalization_candidate_id = -1;
        double last_relocalization_score = 0.0;
        double last_relocalization_search_time_ms = 0.0;
        std::string last_relocalization_reason;
        bool adaptive_lidar_load_enabled = false;
        int adaptive_lidar_degradation_step = 0;
        int selected_lidar_point_stride = 1;
        std::vector<int> selected_lidar_ids;
        std::vector<int> current_frame_lidar_ids;
        std::vector<int> current_frame_missing_lidar_ids;
        double lidar_correction_age_sec = 0.0;
        double last_lio_processing_ms = 0.0;
        std::uint64_t adaptive_lidar_stale_drop_count = 0;
        std::uint32_t cloud_publish_min_lidars = 0;
        bool cloud_publish_eligible = true;
        bool wheel_speed_dr_enabled = false;
        WheelSpeedDrStats wheel_speed_dr_stats;
    };

    struct Options {
        Options() {}

        bool online_mode_ = false;  // 在线模式还是离线模式
        bool with_ui_ = false;      // 是否带ui

        /// 参数
        SE3 T_body_lidar_;

        bool enable_lidar_odom_skip_ = false;  // 是否允许激光里程计跳帧
        int lidar_odom_skip_num_ = 1;          // 如果允许跳帧，跳多少帧
        bool enable_lidar_loc_skip_ = true;    // 是否允许激光定位跳帧
        bool enable_lidar_loc_rviz_ = false;   // 是否允许调试用rviz
        int lidar_loc_skip_num_ = 4;           // 如果允许跳帧，跳多少帧
        bool loc_on_kf_ = false;
    };

    Localization(Options options = Options());
    ~Localization() = default;

    /**
     * 初始化，读配置参数
     * @param yaml_path
     * @param global_map_path
     * @param init_reloc_pose
     */
    bool Init(const std::string& yaml_path, const std::string& global_map_path);

    /// 处理lidar消息
    void ProcessLidarMsg(const sensor_msgs::msg::PointCloud2::SharedPtr laser_msg);
    void ProcessLidarMsg(const sensor_msgs::msg::PointCloud2::SharedPtr laser_msg, int lidar_id);
    void ProcessLivoxLidarMsg(const livox_ros_driver2::msg::CustomMsg::SharedPtr laser_msg);

    /// 处理IMU消息
    void ProcessIMUMsg(IMUPtr imu);

    /// Observe signed longitudinal vehicle speed converted from CAN motor rpm.
    void ProcessWheelSpeed(double timestamp, double longitudinal_speed_mps,
                           double motor_torque_nm = 0.0);

    // void ProcessOdomMsg(const nav_msgs::msg::Odometry::SharedPtr odom_msg) override;

    /// 由外部设置pose，适用于手动重定位
    void SetExternalPose(const Eigen::Quaterniond& q, const Eigen::Vector3d& t);

    /// TODO: 其他初始化逻辑

    /// TODO: 处理odom消息

    /// 结束，保存临时地图
    void Finish();

    /// 异步处理函数
    void LidarOdomProcCloud(CloudPtr, int lidar_id);
    void DrainLioOutputs();
    struct LidarLocInput {
        CloudPtr localization_cloud;
        CloudPtr publication_cloud;
        MultiLidarFrameStats frame_stats;
        bool publication_eligible = true;
        std::chrono::steady_clock::time_point primary_received_at{};
        std::chrono::steady_clock::time_point enqueued_at{};
    };
    void LidarLocProcCloud(const LidarLocInput& input);

    using TFCallback = std::function<void(const geometry_msgs::msg::TransformStamped& odom)>;
    using LocalizationResultCallback = std::function<void(const LocalizationResult& result)>;
    using ProcessedCloudCallback = std::function<void(const CloudPtr& cloud, const LocalizationResult& result,
                                                       const MultiLidarFrameStats& stats, bool eligible)>;
    using LocStateCallback = std::function<void(const std_msgs::msg::Int32& state)>;
    using PointcloudBodyCallback = std::function<void(const sensor_msgs::msg::PointCloud2& pointcloud)>;
    using PointcloudWorldCallback = std::function<void(const sensor_msgs::msg::PointCloud2& pointcloud)>;

    void SetTFCallback(TFCallback&& callback);
    void SetLocalizationResultCallback(LocalizationResultCallback&& callback);
    void SetGlobalLocalizationResultCallback(LocalizationResultCallback&& callback);
    void SetProcessedCloudCallback(ProcessedCloudCallback&& callback);
    void SetLocStateCallback(LocStateCallback&& callback);

    bool IsMultiLidarEnabled() const;
    const MultiLidarConfig& GetMultiLidarConfig() const;
    SO3 GetInitialLidarRotation() const;
    RuntimeStats GetRuntimeStats() const;

    // void SetPathCallback(std::function<void(const nav_msgs::msg::Path& path)>&& callback);
    // void SetPointcloudWorldCallback(std::function<void(const sensor_msgs::msg::PointCloud2& pointcloud)>&& callback);
    // void SetPointcloudBodyCallback(std::function<void(const sensor_msgs::msg::PointCloud2& pointcloud)>&& callback);
    // void SetHealthDiagNormalCallback(interface::health_diag_normal_callback&& callback);

   private:
    friend class LocalizationLockingTestPeer;

    /// 模块  ========================================================================================================
    // Reinitialization needs exclusive ownership, while live callbacks and
    // workers only need the module pointers to remain valid. Shared ownership
    // keeps ROS input responsive while backend processing is in progress.
    mutable std::shared_mutex lifecycle_mutex_;
    // Preserve callback ordering and protect the shared preprocessor/skip
    // counter without coupling ROS input to backend processing latency.
    std::mutex input_mutex_;
    // Serialize ordered LIO/PGO state updates without blocking ROS callbacks.
    std::mutex processing_mutex_;
    // Keep the timestamp check and all live callbacks in the same critical
    // section so two producer threads cannot publish in reverse order.
    std::mutex live_output_dispatch_mutex_;
    MonotonicTimestampGate live_output_timestamp_gate_;
    Options options_;

    /// 预处理
    std::shared_ptr<PointCloudPreprocess> preprocess_ = nullptr;  // point cloud preprocess

    /// 前端
    std::shared_ptr<LaserMapping> lio_ = nullptr;
    Keyframe::Ptr lio_kf_ = nullptr;

    // ui
    std::shared_ptr<ui::PangolinWindow> ui_ = nullptr;

    // pose graph
    std::shared_ptr<PGO> pgo_ = nullptr;

    // lidar localization
    std::shared_ptr<LidarLoc> lidar_loc_;

    /// TODO async 处理
    struct SensorInput {
        IMUPtr imu;
        CloudPtr cloud;
        int lidar_id = 0;
        bool is_imu = false;
        std::chrono::steady_clock::time_point enqueued_at{};
    };
    void ProcessSensorInput(const SensorInput& input);
    void ProcessIMUData(IMUPtr imu);
    bool ShouldThrottleLidarInput(double timestamp);
    void ObserveSensorEnqueued(double timestamp);
    void ObserveSensorProcessed(double timestamp);
    void ObserveLioForStaticDetector(const NavState& state);
    void ObserveLidarLocForStaticDetector(const LocalizationResult& result);
    bool UpdateImuStaticState(const IMUPtr& imu);
    sys::AsyncMessageProcess<SensorInput> sensor_proc_;
    sys::AsyncMessageProcess<LidarLocInput> lidar_loc_proc_cloud_;  // 定位点云和同帧完整发布点云
    // ROS/DDS publication must never hold the ordered sensor/PGO processing
    // path. Capacity one intentionally keeps only the newest live pose.
    sys::AsyncMessageProcess<LocalizationResult> high_frequency_output_proc_;
    profiling::MultiStageTimingWindow lidar_input_timing_window_{
        {"preprocess", "dispatch", "outer"}};
    profiling::MultiStageTimingWindow imu_dr_timing_window_{
        {"lio_imu", "drain_lio", "state_static", "lidar_loc_dr", "pgo_dr", "outer"}};
    // Single sensor consumer owns these windows; queue age uses monotonic time,
    // not the sensor Header clock. It excludes processing after callback entry.
    profiling::MultiStageTimingWindow imu_queue_timing_window_{{"queue_wait"}};
    profiling::MultiStageTimingWindow lidar_queue_timing_window_{{"queue_wait"}};
    // Bounded primary-frame arrival ledger bridges preprocessing, sensor queue,
    // LIO and the asynchronous map matcher. Its clock never uses ROS wall time.
    void RememberPrimaryArrival(int lidar_id, std::uint64_t stamp,
                                std::chrono::steady_clock::time_point received_at);
    std::mutex arrival_mutex_;
    std::map<std::uint64_t, std::chrono::steady_clock::time_point> primary_arrivals_;
    profiling::MultiStageTimingWindow lidar_end_to_end_window_{{"primary_to_pgo"}, 10.0};
    std::size_t lidar_deadline_misses_ = 0;
    int lidar_odom_skip_cnt_ = 0;
    double online_sensor_max_lag_sec_ = 0.0;
    double online_sensor_resume_lag_sec_ = 0.0;
    std::atomic_bool lidar_overload_throttled_{false};

    struct StaticImuSample {
        double timestamp = 0.0;
        double gyro_norm = 0.0;
        double accel_norm = 0.0;
    };
    mutable std::mutex static_detector_mutex_;
    std::deque<StaticImuSample> static_imu_window_;
    double static_gyro_sum_ = 0.0;
    double static_gyro_sq_sum_ = 0.0;
    double static_accel_sum_ = 0.0;
    double static_accel_sq_sum_ = 0.0;
    double last_static_lio_stamp_ = 0.0;
    double last_static_lio_speed_ = 0.0;
    bool last_static_lio_reliable_ = false;
    double last_valid_lidar_loc_stamp_ = 0.0;
    struct StaticWheelSpeedSample {
        double timestamp = 0.0;
        double speed_mps = 0.0;
    };
    std::deque<StaticWheelSpeedSample> static_wheel_speed_history_;
    double last_wheel_speed_stamp_ = 0.0;
    bool wheel_speed_observed_ = false;
    bool wheel_timestamp_mismatch_reported_ = false;
    int static_exit_count_ = 0;
    bool imu_static_hold_enabled_ = false;
    bool imu_static_hold_active_ = false;

    /// 结果数据 =====================================================================================================
    LocalizationResult loc_result_;

    /// 框架相关
    TFCallback tf_callback_;
    LocalizationResultCallback localization_result_callback_;
    LocalizationResultCallback global_localization_result_callback_;
    ProcessedCloudCallback processed_cloud_callback_;
    LocStateCallback loc_state_callback_;
    PointcloudBodyCallback pointcloud_body_callback_;
    PointcloudWorldCallback pointcloud_world_callback_;

    /// 输入检查
    double last_imu_time_ = 0;
    double last_odom_time_ = 0;
    double last_cloud_time_ = 0;

    mutable std::mutex runtime_stats_mutex_;
    RuntimeStats runtime_stats_;
};
}  // namespace loc

}  // namespace lightning
