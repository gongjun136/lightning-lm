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
#include <geosun_msgs/msg/safety.hpp>
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
#include "utils/radar_dtc_state.h"
#include "core/localization/localization_dtc_monitor.h"
#include "utils/compute_profiling.h"

namespace lightning {

namespace loc {
class Localization;
}

/**
 * @brief 在线 LiDAR 定位的 ROS 编排层，连接传感器输入、定位算法与业务输出。
 * @details 持有 loc::Localization 和 ROS 节点，负责订阅、发布、心跳与轨迹记录。
 * 由 run_loc_online.cc 创建；算法内部的数据处理流程见 loc::Localization。
 * @see Init(), SetInitPose(), Spin(), Finish(),
 *      @ref online_localization_flow "在线定位数据流"
 */
class LocSystem {
   public:
    struct Options {
        bool pub_tf_ = false;  ///< 是否发布 TF，可由 YAML 的 system.pub_tf 显式启用。
    };

    explicit LocSystem(Options options);
    ~LocSystem();

    /**
     * @brief 根据配置与地图路径准备在线定位算法、ROS 输入与结果输出。
     *
     * @details
     * 由 @ref run_loc_online.cc "run_loc_online" 的 LiDAR 定位分支调用。
     * 成功后，入口通过 SetInitPose() 设置初始位姿，再调用 Spin() 处理 ROS 输入。
     *
     * @par 主要流程
     *
     * 以下名称与函数体中的“阶段 1”至“阶段 7”对应，可按阶段定位源码。
     *
     * 1. <b>重置状态并创建算法对象</b>：重置时间戳门控及部分运行状态，创建在线 Localization 对象。
     * 2. <b>加载地图与输出配置</b>：确定地图路径、输出坐标及输入/发布约束，准备轮速参数。
     * 3. <b>创建节点并初始化算法</b>：创建 ROS 节点与心跳，调用
     *    @ref lightning::loc::Localization::Init() "Localization::Init()"
     *    初始化底层定位算法。
     * 4. <b>建立诊断与 Safety 通道</b>：读取诊断配置，接入必需雷达的诊断状态，创建 Safety 发布器和诊断定时器。
     * 5. <b>注册传感器与轮速输入</b>：配置输入 QoS，按单/多雷达模式订阅 IMU 和点云，并按开关订阅 CAN 轮速。
     * 6. <b>建立业务输出并注册算法回调</b>：创建位姿、点云和状态发布器，绑定结果发布及轨迹记录回调。
     * 7. <b>完成初始化并报告 Idle</b>：更新心跳状态，返回初始化成功。
     *
     * 订阅、定时器和结果回调在这里注册；其回调体在后续被触发时执行。
     *
     * @htmlonly[block]
     * <details>
     * <summary>关键代码（点击展开）</summary>
     * @endhtmlonly
     *
     * <b>创建在线算法对象</b>
     *
     * @snippet core/system/loc_system.cc loc-system-init-online-object
     *
     * <b>连接结果发布</b>：将实时结果交给 PublishLocalizationResult()。
     *
     * @snippet core/system/loc_system.cc loc-system-init-live-result
     *
     * @htmlonly[block]
     * </details>
     * @endhtmlonly
     *
     * @param yaml_path 定位配置文件路径。
     * @param map_path_override 非空时覆盖 YAML 中的 `system.map_path`。
     * @return 初始化成功返回 true，显式校验或算法初始化失败返回 false。
     * @note 心跳与算法后台线程在初始化期间启动；ROS 订阅和定时器回调由后续 Spin() 调度。
     * @note Init() 成功不代表已有有效定位结果；SetInitPose() 负责放行定位输入处理。
     * @note YAML 读取或字段转换异常可向调用方传播。
     * @see SetInitPose(), Spin(), PublishLocalizationResult(),
     *      @ref online_localization_flow "在线定位数据流"
     */
    bool Init(const std::string& yaml_path, const std::string& map_path_override = "");

    /**
     * @brief 将初始位姿交给定位算法，并开启传感器处理。
     * @details 调用 loc::Localization::SetExternalPose() 后打开输入处理开关。
     * @param pose 用于初始化或重新定位的地图位姿。
     * @see Init(), Spin(), loc::Localization::SetExternalPose()
     */
    void SetInitPose(const SE3& pose);

    /**
     * @brief 记录 IMU 输入状态，并在定位启动后转交算法处理。
     * @details 由 ROS 订阅回调调用；SetInitPose() 之前只记录输入状态。
     * @see loc::Localization::ProcessIMUMsg()
     */
    void ProcessIMU(const lightning::IMUPtr& imu);

    /**
     * @brief 检查标准点云的输入时效，并在定位启动后转交算法。
     * @details 多雷达模式默认使用主雷达 ID；单雷达模式使用 ID 0。
     * @see ObserveLidarInput(), loc::Localization::ProcessLidarMsg()
     */
    void ProcessLidar(const sensor_msgs::msg::PointCloud2::SharedPtr& cloud);
    /**
     * @brief 按指定雷达 ID 检查并转发标准点云。
     * @param cloud 非空标准点云消息，使用消息头时间戳检查输入时效。
     * @param lidar_id 与多雷达配置对应的传感器 ID。
     * @see ObserveLidarInput(), loc::Localization::ProcessLidarMsg()
     */
    void ProcessLidar(const sensor_msgs::msg::PointCloud2::SharedPtr& cloud, int lidar_id);
    /**
     * @brief 检查 Livox 自定义点云的输入时效，并在定位启动后转交算法。
     * @details 输入统计使用 ID 0；算法侧通过 ProcessLivoxLidarMsg() 处理该消息格式。
     * @see loc::Localization::ProcessLivoxLidarMsg()
     */
    void ProcessLidar(const livox_ros_driver2::msg::CustomMsg::SharedPtr& cloud);

    /**
     * @brief 进入 ROS 事件循环，执行传感器回调和定时器，直到 ROS 关闭。
     * @see Init(), SetInitPose(), Finish()
     */
    void Spin();
    /**
     * @brief 等待算法队列处理结束，停止心跳并关闭已发布轨迹文件。
     * @details 入口在 ROS 输入停止后调用；析构也会调用，重复调用直接返回。
     * @see loc::Localization::Finish(), StartPublishedTrajectoryTum()
     */
    void Finish();
    /**
     * @brief 将 PGO 全局校正结果保存为 TUM 轨迹，与离线轨迹定义一致。
     * @return 文件打开成功且至少写入一个位姿时返回 true。
     * @see SaveHighFrequencyTrajectoryTum(), CaptureGlobalLocalizationResult()
     */
    bool SaveTrajectoryTum(const std::string& path) const;
    /**
     * @brief 将通过位姿发布门控的高频算法状态保存为 TUM 轨迹。
     * @details 记录算法坐标下的状态；业务输出坐标下的后轴轨迹由 StartPublishedTrajectoryTum() 记录。
     * @return 文件打开成功且至少写入一个位姿时返回 true。
     * @see SaveTrajectoryTum(), StartPublishedTrajectoryTum()
     */
    bool SaveHighFrequencyTrajectoryTum(const std::string& path) const;
    /**
     * @brief 打开 TUM 文件，在运行中记录实际发布到 /localization/pose_vel 的后轴位姿。
     * @details 创建父目录并截断目标文件；记录频率受运行开销策略控制。
     * @return 文件与表头初始化成功返回 true；路径、目录或文件操作失败返回 false。
     * @see PublishLocalizationResult(), Finish()
     */
    bool StartPublishedTrajectoryTum(const std::string& path);

   private:
    /**
     * @brief 将高频定位结果转换为后轴业务消息，经过门控后发布并记录轨迹。
     * @details 由实时结果回调调用，是算法坐标与业务输出坐标的连接处。
     * @par 主要流程
     * 以下名称与函数体中的阶段注释对应；未通过校验时提前返回。
     * 1. <b>校验结果与发布条件</b>：检查结果有效性、时间戳递增和位姿发布门控。
     * 2. <b>缓存状态并转换后轴位姿</b>：保存算法坐标轨迹，转换输出位姿及诊断协方差。
     * 3. <b>构造业务消息</b>：计算后轴纵向速度，填充位姿与通信头。
     * 4. <b>发布位姿并更新工作状态</b>：按发布器是否存在发送消息，业务发布后推进工作计数。
     * 5. <b>记录审计与耗时</b>：按运行策略记录影子速度、发布日志、TUM 及耗时统计。
     * @note 调用 publish() 不代表下游已经接收；缓存算法轨迹也不代表业务消息已发布。
     * @see PublishProcessedCloud(), RecordPublishedPoseTum(),
     *      @ref output_gate_reading "发布门控的源码导读",
     *      @ref output_pose_reading "后轴位姿的源码导读",
     *      @ref output_speed_reading "后轴速度的源码导读"
     */
    void PublishLocalizationResult(const loc::LocalizationResult& result);
    /// @brief 收集有效且时间递增的 PGO 全局校正结果，供 SaveTrajectoryTum() 导出。
    void CaptureGlobalLocalizationResult(const loc::LocalizationResult& result);
    bool WriteTrajectoryTum(const std::string& path, const std::vector<NavState>& states,
                            const char* description) const;
    void RecordPublishedPoseTum(const geometry_msgs::msg::PoseStamped& pose);
    void ClosePublishedTrajectoryTum();
    /**
     * @brief 更新匹配有效性与发布门控，并按条件发布后轴和地图坐标下的点云。
     * @details 雷达数量不满足要求时停止点云发布；地图点云还受匹配门控和降频控制。
     * @par 主要流程
     * 1. <b>更新匹配门控</b>：记录本次匹配有效性，更新遥测和定位状态事件。
     * 2. <b>检查点云发布条件</b>：检查雷达数量、发布器、点云与时间戳，失败时提前返回。
     * 3. <b>发布后轴点云</b>：按固定外参转换点云，发布到后轴坐标系。
     * 4. <b>按条件发布地图点云并记录耗时</b>：经过降频与地图门控后发布，并按开关记录耗时。
     * @note 阶段名称对应函数体；即使点云随后被抑制，阶段 1 仍会更新匹配门控。
     * @see PublishLocalizationResult(), PublishHealthStatus()
     */
    void PublishProcessedCloud(const CloudPtr& cloud, const loc::LocalizationResult& result,
                               const MultiLidarFrameStats& stats, bool eligible);
    /**
     * @brief 汇总输入、定位与队列状态，周期发布健康消息和流水线诊断。
     * @par 主要流程
     * 1. <b>采集输入统计</b>：检查输出资源，在输入锁下采集话题、时间戳与断流时长。
     * 2. <b>发布健康状态并更新工作状态</b>：发布故障与定位状态，按工作计数设置 Idle/Running。
     * 3. <b>汇总队列和定位状态</b>：读取算法运行统计与重定位结果。
     * 4. <b>报告地图输出门控变化</b>：记录失锁、匹配过期及地图输出恢复/停止。
     * 5. <b>补充降级与轮速诊断并发布</b>：填充输入拒绝、负载降级及轮速观测统计后发送。
     * @note 阶段名称对应函数体；Running 仅表示工作计数非零，不等价于当前诊断 PASS。
     */
    void PublishHealthStatus();
    void UpdateDtcHeartbeat();
    /**
     * @brief 按定位和雷达诊断及当前策略发布 Safety，并记录影子策略比较。
     * @par 主要流程
     * 1. <b>评估定位与雷达状态</b>：取得两类诊断快照，计算严格和有界授时策略结果。
     * 2. <b>按配置发布 Safety</b>：根据上报、授时告警和故障屏蔽开关选择实际 runstate。
     * 3. <b>记录影子评估</b>：满足开关与变化/周期条件时记录对比值，保存本次状态。
     * @note 阶段名称对应函数体；影子评估日志不替代阶段 2 已发布的值。
     */
    void PublishSafety();
    bool DtcLocalizationLost() const;
    void PublishPath();
    void RegisterLidarInput(int lidar_id, const std::string& topic);
    /**
     * @brief 记录点云到达情况，并检查相对最新传感器时间戳的滞后程度。
     * @return 时间戳有效且满足配置的输入时效门槛时返回 true。
     * @see ProcessLidar(), PublishHealthStatus()
     */
    bool ObserveLidarInput(int lidar_id, double sensor_stamp);
    void ObserveImuInput(double sensor_stamp);
    void ObserveWheelSpeedInput(double sensor_stamp, double motor_rpm,
                                double motor_torque);

    struct InputTopicStats {
        std::string topic;  ///< 该路输入的 ROS 话题名称。
        std::uint64_t message_count = 0;  ///< 该路输入累计收到的消息数。
        std::uint64_t stale_drop_count = 0;  ///< 该路点云因时间戳门控被拒绝的次数。
        double last_sensor_stamp = 0.0;  ///< 最近收到消息的传感器时间戳（s）。
        std::chrono::steady_clock::time_point last_arrival;  ///< 最近一次接收消息的单调时钟时刻，用于判断断流。
        bool has_arrival = false;  ///< 是否已接收到该路输入，用于区分未启动与断流。
    };

    Options options_;  ///< ROS 编排层的启动选项。

    std::shared_ptr<loc::Localization> loc_ = nullptr;  ///< 定位接口

    std::atomic_bool loc_started_ = false;  ///< 是否开启定位
    std::atomic_bool map_loaded_ = false;   ///< 地图是否已载入

    /// 实时模式下的ros2 node, subscribers
    rclcpp::Node::SharedPtr node_;
    std::shared_ptr<tf2_ros::TransformBroadcaster> tf_broadcaster_ = nullptr;  ///< 发布算法位姿对应 TF 的广播器。

    std::string imu_topic_;    ///< 单雷达配置使用的 IMU 输入话题。
    std::string cloud_topic_;  ///< 单雷达配置使用的 PointCloud2 输入话题。
    std::string livox_topic_;  ///< 单雷达配置使用的 Livox CustomMsg 输入话题。
    std::string map_frame_ = "map";              ///< 业务输出使用的地图坐标系名称。
    std::string rear_axle_frame_ = "rear_axle";  ///< 后轴坐标系点云使用的 frame ID。
    bool wheel_speed_observation_enabled_ = true;  ///< 是否订阅 CAN 轮速并送入定位滤波器。
    bool speed_smoothing_shadow_enabled_ = false;  ///< 是否计算用于诊断日志对照的平滑速度。
    sany_output::SpeedSmoothingShadow speed_smoothing_shadow_;  ///< 速度平滑影子实验的滤波状态。
    std::string wheel_speed_topic_ = "/SpeThrCAN4_topic";  ///< 提供电机转速与扭矩的 CAN 输入话题。
    double wheel_speed_scale_mps_per_rpm_ = 0.00120639253574024;  ///< 电机转速到车速的换算系数，单位 (m/s)/rpm。
    double online_lidar_input_max_timestamp_lag_sec_ = 0.0;  ///< 相对最新 LiDAR 时间戳允许的最大滞后（s）；0 关闭滞后限制。
    Vec3d primary_lidar_position_in_body_ = Vec3d::Zero();  ///< 主雷达在后轴坐标系中的安装位置，单位 m。
    sany_output::FixedMapTransform fixed_map_transform_;  ///< 定位地图到业务输出地图的固定坐标变换。
    sany_output::LocalizationPublicationGate publication_gate_;  ///< 根据匹配有效性、连续丢帧与匹配时效控制发布。
    sany_output::FrameDecimator map_cloud_decimator_{10};  ///< 地图坐标系点云的抽帧器，每 10 帧允许一次发布。
    std::unique_ptr<sany_output::LocalizationTelemetryState> telemetry_;  ///< 汇总定位健康状态并缓存待发布的轨迹路径。
    mutable std::mutex trajectory_mutex_;  ///< 保护高频与全局校正两组轨迹状态。
    std::vector<NavState> localization_states_;  ///< 高频定位状态轨迹，供 SaveHighFrequencyTrajectoryTum() 导出。
    std::vector<NavState> global_localization_states_;  ///< 全局校正状态轨迹，供 SaveTrajectoryTum() 导出。
    mutable std::mutex published_tum_mutex_;  ///< 保护已发布 TUM 文件流、时间戳与写入统计。
    std::ofstream published_tum_stream_;  ///< 记录实际发布的后轴位姿的 TUM 文件流。
    std::string published_tum_path_;  ///< 已发布位姿 TUM 轨迹文件的输出路径。
    double last_published_tum_timestamp_ = 0.0;  ///< TUM 写入门控使用的上一条位姿时间戳（s）。
    double last_production_pose_audit_stamp_ = 0.0;  ///< 最近一次生产位姿审计日志的传感器时间戳（s），用于限制日志频率。
    std::uint64_t published_tum_pose_count_ = 0;  ///< 成功写入已发布 TUM 轨迹的位姿数。
    std::uint64_t published_tum_write_error_count_ = 0;  ///< TUM 记录位姿时的失败次数。
    bool finished_ = false;  ///< 标记 Finish() 已完成，避免重复收尾。
    mutable std::mutex input_stats_mutex_;  ///< 保护传感器接收统计及最近 CAN 观测值。
    MaximumLagTimestampGate lidar_input_timestamp_gate_;  ///< 跨雷达共享的时间戳门控，拒绝无效或滞后超限的输入。
    std::map<int, InputTopicStats> lidar_input_stats_;  ///< 按雷达 ID 保存输入话题及接收统计。
    InputTopicStats imu_input_stats_;  ///< IMU 话题的接收与到达统计。
    InputTopicStats wheel_speed_input_stats_;  ///< CAN 轮速话题的接收与到达统计。
    double last_wheel_speed_mps_ = 0.0;  ///< 最近 CAN 观测换算得到的车速，单位 m/s。
    double last_motor_rpm_ = 0.0;  ///< 最近 CAN 观测的电机转速，单位 rpm。
    double last_motor_torque_ = 0.0;  ///< 最近 CAN 观测的电机扭矩。
    std::atomic<double> last_localization_stamp_{0.0};  ///< 最近通过输出时间戳检查的定位结果时间（s）。
    std::atomic<double> last_published_pose_stamp_{0.0};  ///< 最近实际发布的 pose_vel 时间戳（s）。
    std::atomic<std::uint64_t> cloud_publish_suppressed_count_{0};  ///< 因输入完整性要求而抑制点云发布的帧数。
    std::mutex pose_publish_mutex_;  ///< 串行保护位姿时间戳检查和业务消息发布。
    MonotonicTimestampGate pose_timestamp_gate_;  ///< 保证 pose_vel 业务输出时间戳严格递增。
    std::atomic_bool map_outputs_ever_enabled_{false};  ///< 地图输出门控是否曾经开启，用于状态变化诊断。
    std::atomic_bool map_outputs_enabled_last_{false};  ///< 上一次观察到的地图输出门控状态。
    std::atomic_bool localization_ever_good_{false};  ///< 定位是否曾达到有效状态，用于故障恢复诊断。
    std::atomic<std::uint32_t> pose_vel_comm_seq_{0};  ///< 位姿与速度消息的通信序列号。
    std::atomic<std::uint32_t> fault_status_comm_seq_{0};  ///< 故障状态消息的通信序列号。
    profiling::MultiStageTimingWindow output_timing_window_{
        {"pose_transform", "message_build", "ros_publish", "diagnostic_io",
         "record_tum", "outer"}};  ///< 位姿转换、消息构建、发布和记录各阶段的耗时窗口。

    rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr imu_sub_ = nullptr;  ///< 接收 IMU 消息的 ROS 订阅器。
    std::vector<rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr> cloud_subs_;  ///< 按雷达配置保存的 PointCloud2 订阅器。
    rclcpp::Subscription<livox_ros_driver2::msg::CustomMsg>::SharedPtr livox_sub_ = nullptr;  ///< 接收 Livox CustomMsg 的 ROS 订阅器。
    rclcpp::GenericSubscription::SharedPtr wheel_speed_sub_ = nullptr;  ///< 接收 CAN 转速与扭矩的通用 ROS 订阅器。

    std::unique_ptr<functional_safety::HeartbeatPublisher> heartbeat_;  ///< 独立发布功能安全心跳并记录有效工作进度。
    std::unique_ptr<loc::dtc::Monitor> dtc_monitor_;
    std::unique_ptr<functional_safety::RadarDtcState> radar_dtc_state_;
    std::string dtc_config_version_ = "localization-dtc-v3-optional-can";
    bool safety_fault_reporting_enabled_ = true;  ///< 关闭时 Safety.runstate 固定为正常值 0，源 DTC 不受影响。
    bool radar_time_sync_shadow_enabled_ = true;  ///< 仅记录反事实判定，不改变真实 Safety。
    bool radar_bounded_sync_warning_enabled_ = false;
    bool radar_fault_mask_enabled_ = false;
    int last_safety_published_ = -128;
    uint8_t last_safety_bounded_mask_ = 0;
    int last_safety_strict_ = -128, last_safety_shadow_ = -128;
    uint8_t last_safety_sync_mask_ = 0;
    uint64_t safety_shadow_ticks_ = 0;
    std::mutex dtc_publish_mutex_;
    std::array<uint8_t, 2> last_dtc_conditions_{{0, 0}};
    rclcpp::Publisher<geosun_msgs::msg::Safety>::SharedPtr safety_pub_;
    std::array<rclcpp::Subscription<livox_ros_driver2::msg::LidarDiagnosticState>::SharedPtr, 4> radar_availability_subs_;
    rclcpp::TimerBase::SharedPtr dtc_timer_, safety_timer_;

    rclcpp::Publisher<geometry_msgs::msg::PoseStamped>::SharedPtr pose_pub_ = nullptr;  ///< 发布后轴位姿 PoseStamped。
    rclcpp::Publisher<lightning::msg::VehiclePose>::SharedPtr vehicle_pose_pub_ = nullptr;  ///< 发布车辆位姿与速度消息。
    rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr inv_cloud_pub_ = nullptr;  ///< 发布后轴坐标系下的处理后点云。
    rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr map_cloud_pub_ = nullptr;  ///< 发布业务地图坐标系下的处理后点云。
    rclcpp::Publisher<lightning::msg::FaultStatus>::SharedPtr fault_status_pub_ = nullptr;  ///< 发布定位故障状态。
    rclcpp::Publisher<lightning::msg::LocalizationStatus>::SharedPtr loc_status_pub_ = nullptr;  ///< 发布定位运行状态。
    rclcpp::Publisher<lightning::msg::PipelineDiagnostics>::SharedPtr pipeline_diagnostics_pub_ = nullptr;  ///< 发布输入、队列与算法处理诊断。
    rclcpp::Publisher<lightning::msg::DebugMessage>::SharedPtr debug_message_pub_ = nullptr;  ///< 发布状态变化和异常调试事件。
    rclcpp::Publisher<nav_msgs::msg::Path>::SharedPtr path_pub_ = nullptr;  ///< 发布用于可视化的后轴轨迹路径。
    rclcpp::TimerBase::SharedPtr health_timer_ = nullptr;  ///< 周期触发健康状态与流水线诊断发布。
    rclcpp::TimerBase::SharedPtr path_timer_ = nullptr;  ///< 周期触发可视化轨迹路径发布。
};

};  // namespace lightning

#endif  // LIGHTNING_LOC_SYSTEM_H
