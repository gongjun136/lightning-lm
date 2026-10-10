//
// Created by xiang on 25-9-8.
//

#ifndef LIGHTNING_LOC_SYSTEM_H
#define LIGHTNING_LOC_SYSTEM_H

#include <atomic>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <map>
#include <tf2_ros/transform_broadcaster.h>
#include <geometry_msgs/msg/pose_stamped.hpp>
#include <geosun_msgs/msg/pos_res.hpp>
#include <geosun_msgs/msg/spe_thr_can4.hpp>
#include <nav_msgs/msg/path.hpp>
#include <rclcpp/rclcpp.hpp>
#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <memory>
#include <mutex>
#include <vector>

#include "livox_ros_driver2/msg/custom_msg.hpp"

#include "common/eigen_types.h"
#include "common/imu.h"
#include "common/keyframe.h"
#include "common/timestamp_gate.h"
#include "core/lio/multi_lidar_fusion.h"
#include "core/localization/localization_result.h"
#include "core/system/sany_localization_output.h"
#include "core/system/speed_smoothing_shadow.h"
#include "lightning/msg/debug_message.hpp"
#include "lightning/msg/pipeline_diagnostics.hpp"
#include "lightning/msg/vehicle_pose.hpp"
#include "utils/functional_safety_heartbeat.h"
#include "utils/compute_profiling.h"

namespace lightning {

namespace loc {
class Localization;
}

class LocSystem {
   public:
    struct Options {
        bool pub_tf_ = false;  // YAML system.pub_tf may explicitly enable /tf.
    };

    explicit LocSystem(Options options);
    ~LocSystem();

    /**
     * @brief 根据配置与地图路径装配在线定位系统，连接 ROS 输入输出与底层定位算法。
     *
     * @details
     * 由 @ref run_loc_online.cc "run_loc_online" 入口调用；成功后，入口通过 SetInitPose()
     * 设置初始位姿，再由 Spin() 开始处理输入。
     *
     * 主要步骤：
     *
     * 1. <b>读取配置</b>：确定地图路径、输出坐标及输入和发布约束。
     * 2. <b>初始化算法</b>：创建 ROS 节点与心跳，调用
     *    @ref lightning::loc::Localization::Init() "Localization::Init()"
     *    准备 LIO、地图匹配、PGO 和处理线程。
     * 3. <b>连接输入</b>：按单雷达或多雷达配置订阅 IMU、点云，并按需接入 CAN 轮速。
     * 4. <b>连接输出</b>：创建业务发布器和定时器，通过结果回调连接算法输出、业务发布与轨迹记录。
     *
     * @htmlonly[block]
     * <details>
     * <summary>关键代码：在线算法对象与实时结果回调</summary>
     * @endhtmlonly
     *
     * <b>在线算法对象</b>：创建并设置在线模式的 Localization 对象。
     *
     * @snippet{trimleft} loc_system.cc loc-system-init-online-object
     *
     * <b>实时结果回调</b>：算法产生结果后，交回本层 PublishLocalizationResult() 发布。
     *
     * @snippet{trimleft} loc_system.cc loc-system-init-live-result
     *
     * @htmlonly[block]
     * </details>
     * @endhtmlonly
     *
     * @param yaml_path 定位配置文件路径。
     * @param map_path_override 非空时覆盖 YAML 中的 `system.map_path`。
     * @return 装配成功返回 true，显式校验或初始化失败返回 false；有效定位结果由后续输入产生。
     * @note YAML 读取或字段转换异常可向调用方传播。
     * 后续处理流程见 @ref online_localization_flow "在线定位数据流"。
     */
    bool Init(const std::string& yaml_path, const std::string& map_path_override = "");

    /// 设置初始化位姿
    void SetInitPose(const SE3& pose);

    /// 处理IMU
    void ProcessIMU(const lightning::IMUPtr& imu);

    /// 处理点云
    void ProcessLidar(const sensor_msgs::msg::PointCloud2::SharedPtr& cloud);
    void ProcessLidar(const sensor_msgs::msg::PointCloud2::SharedPtr& cloud, int lidar_id);
    void ProcessLidar(const livox_ros_driver2::msg::CustomMsg::SharedPtr& cloud);

    /// 实时模式下的spin
    void Spin();
    /// Stop workers after draining all accepted sensor and localization data.
    void Finish();
    /// PGO global corrections, matching the offline trajectory definition.
    bool SaveTrajectoryTum(const std::string& path) const;
    /// Live high-frequency extrapolated output published to ROS.
    bool SaveHighFrequencyTrajectoryTum(const std::string& path) const;
    /// Open a live TUM stream containing the final rear-axle pose sent on /PosRes.
    bool StartPublishedTrajectoryTum(const std::string& path);

   private:
    void PublishLocalizationResult(const loc::LocalizationResult& result);
    void CaptureGlobalLocalizationResult(const loc::LocalizationResult& result);
    bool WriteTrajectoryTum(const std::string& path, const std::vector<NavState>& states,
                            const char* description) const;
    void RecordPublishedPoseTum(const geometry_msgs::msg::PoseStamped& pose);
    void ClosePublishedTrajectoryTum();
    void PublishProcessedCloud(const CloudPtr& cloud, const loc::LocalizationResult& result,
                               const MultiLidarFrameStats& stats, bool eligible);
    void PublishHealthStatus();
    void PublishPath();
    void RegisterLidarInput(int lidar_id, const std::string& topic);
    bool ObserveLidarInput(int lidar_id, double sensor_stamp);
    void ObserveImuInput(double sensor_stamp);
    void ObserveWheelSpeedInput(double sensor_stamp, double motor_rpm,
                                double motor_torque);

    struct InputTopicStats {
        std::string topic;
        std::uint64_t message_count = 0;
        std::uint64_t stale_drop_count = 0;
        double last_sensor_stamp = 0.0;
        std::chrono::steady_clock::time_point last_arrival;
        bool has_arrival = false;
    };

    Options options_;

    std::shared_ptr<loc::Localization> loc_ = nullptr;  // 定位接口

    std::atomic_bool loc_started_ = false;  // 是否开启定位
    std::atomic_bool map_loaded_ = false;   // 地图是否已载入

    /// 实时模式下的ros2 node, subscribers
    rclcpp::Node::SharedPtr node_;
    std::shared_ptr<tf2_ros::TransformBroadcaster> tf_broadcaster_ = nullptr;

    std::string imu_topic_;
    std::string cloud_topic_;
    std::string livox_topic_;
    std::string map_frame_ = "map";
    std::string rear_axle_frame_ = "rear_axle";
    bool wheel_speed_observation_enabled_ = true;
    bool speed_smoothing_shadow_enabled_ = false;
    sany_output::SpeedSmoothingShadow speed_smoothing_shadow_;
    std::string wheel_speed_topic_ = "/SpeThrCAN4_topic";
    double wheel_speed_scale_mps_per_rpm_ = 0.00120639253574024;
    double online_lidar_input_max_timestamp_lag_sec_ = 0.0;
    Vec3d primary_lidar_position_in_body_ = Vec3d::Zero();
    sany_output::FixedMapTransform fixed_map_transform_;
    sany_output::LocalizationPublicationGate publication_gate_;
    sany_output::FrameDecimator map_cloud_decimator_{10};
    std::unique_ptr<sany_output::LocalizationTelemetryState> telemetry_;
    mutable std::mutex trajectory_mutex_;
    std::vector<NavState> localization_states_;
    std::vector<NavState> global_localization_states_;
    mutable std::mutex published_tum_mutex_;
    std::ofstream published_tum_stream_;
    std::string published_tum_path_;
    double last_published_tum_timestamp_ = 0.0;
    double last_production_pose_audit_stamp_ = 0.0;
    std::uint64_t published_tum_pose_count_ = 0;
    std::uint64_t published_tum_write_error_count_ = 0;
    bool finished_ = false;
    mutable std::mutex input_stats_mutex_;
    MaximumLagTimestampGate lidar_input_timestamp_gate_;
    std::map<int, InputTopicStats> lidar_input_stats_;
    InputTopicStats imu_input_stats_;
    InputTopicStats wheel_speed_input_stats_;
    double last_wheel_speed_mps_ = 0.0;
    double last_motor_rpm_ = 0.0;
    double last_motor_torque_ = 0.0;
    std::atomic<double> last_localization_stamp_{0.0};
    std::atomic<double> last_posres_stamp_{0.0};
    std::atomic<std::uint64_t> cloud_publish_suppressed_count_{0};
    std::mutex posres_publish_mutex_;
    MonotonicTimestampGate posres_timestamp_gate_;
    std::atomic_bool map_outputs_ever_enabled_{false};
    std::atomic_bool map_outputs_enabled_last_{false};
    std::atomic_bool localization_ever_good_{false};
    std::atomic<std::uint32_t> pose_vel_comm_seq_{0};
    std::atomic<std::uint32_t> fault_status_comm_seq_{0};
    profiling::MultiStageTimingWindow output_timing_window_{
        {"pose_transform", "message_build", "ros_publish", "diagnostic_io",
         "record_tum", "outer"}};

    rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr imu_sub_ = nullptr;
    std::vector<rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr> cloud_subs_;
    rclcpp::Subscription<livox_ros_driver2::msg::CustomMsg>::SharedPtr livox_sub_ = nullptr;
    rclcpp::GenericSubscription::SharedPtr wheel_speed_sub_ = nullptr;

    std::unique_ptr<functional_safety::HeartbeatPublisher> heartbeat_;

    rclcpp::Publisher<geosun_msgs::msg::PosRes>::SharedPtr pos_res_pub_ = nullptr;
    rclcpp::Publisher<geometry_msgs::msg::PoseStamped>::SharedPtr pose_pub_ = nullptr;
    rclcpp::Publisher<lightning::msg::VehiclePose>::SharedPtr vehicle_pose_pub_ = nullptr;
    rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr inv_cloud_pub_ = nullptr;
    rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr map_cloud_pub_ = nullptr;
    rclcpp::Publisher<lightning::msg::FaultStatus>::SharedPtr fault_status_pub_ = nullptr;
    rclcpp::Publisher<lightning::msg::LocalizationStatus>::SharedPtr loc_status_pub_ = nullptr;
    rclcpp::Publisher<lightning::msg::PipelineDiagnostics>::SharedPtr pipeline_diagnostics_pub_ = nullptr;
    rclcpp::Publisher<lightning::msg::DebugMessage>::SharedPtr debug_message_pub_ = nullptr;
    rclcpp::Publisher<nav_msgs::msg::Path>::SharedPtr path_pub_ = nullptr;
    rclcpp::TimerBase::SharedPtr health_timer_ = nullptr;
    rclcpp::TimerBase::SharedPtr path_timer_ = nullptr;
};

};  // namespace lightning

#endif  // LIGHTNING_LOC_SYSTEM_H
