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
 * @brief 组织 LIO、地图匹配与 PGO，将传感器输入转换为定位结果。
 * @details 由 lightning::LocSystem 或离线入口使用；在线模式通过消息队列处理输入和输出，
 * 离线模式直接执行同一条算法链。ROS 订阅与业务发布由系统层负责。
 * @see Init(), ProcessIMUMsg(), ProcessLidarMsg(), Finish(),
 *      @ref localization_module "定位与重定位模块"
 */
class Localization {
   public:
    /// @brief 队列积压、处理滞后、重定位及雷达负载状态，供诊断接口读取。
    struct RuntimeStats {
        std::uint64_t sensor_queue_pending = 0;  ///< 传感器队列待执行与正在执行的消息数。
        std::uint64_t sensor_queue_dropped = 0;  ///< 传感器队列累计丢弃数。
        std::uint64_t sensor_queue_processed = 0;  ///< 传感器队列已完成的消息数。
        std::uint64_t localization_queue_pending = 0;  ///< 地图匹配队列待执行与正在执行的消息数。
        std::uint64_t localization_queue_dropped = 0;  ///< 地图匹配队列累计丢弃数。
        std::uint64_t localization_queue_processed = 0;  ///< 地图匹配队列已完成的消息数。
        std::uint64_t high_frequency_queue_pending = 0;  ///< 高频输出队列待执行与正在执行的消息数。
        std::uint64_t high_frequency_queue_dropped = 0;  ///< 高频输出队列累计丢弃数。
        std::uint64_t high_frequency_queue_processed = 0;  ///< 高频输出队列已完成的消息数。
        double latest_enqueued_sensor_stamp = 0.0;  ///< 已入队传感器数据的最新时间戳（s）。
        double latest_processed_sensor_stamp = 0.0;  ///< 已处理传感器数据的最新时间戳（s）。
        double current_sensor_lag_sec = 0.0;  ///< 最新入队与完成时间戳之间的处理滞后（s）。
        double max_sensor_lag_sec = 0.0;  ///< 运行期间观察到的最大传感器处理滞后（s）。
        std::uint64_t severe_timestamp_rollback_count = 0;  ///< 严重传感器时间戳回退的累计次数。
        double worst_timestamp_rollback_sec = 0.0;  ///< 最大的传感器时间戳回退量（s）。
        std::uint64_t live_output_non_monotonic_drop_count = 0;  ///< 实时输出因时间戳不递增而被拒绝的次数。
        double worst_live_output_timestamp_rollback_sec = 0.0;  ///< 实时输出最大的时间戳回退量（s）。
        std::uint64_t relocalization_attempt_count = 0;  ///< 全局重定位尝试总数。
        std::uint64_t relocalization_accept_count = 0;  ///< 全局重定位成功接受的次数。
        bool last_relocalization_candidate_found = false;  ///< 最近一次重定位是否找到候选。
        bool last_relocalization_accepted = false;  ///< 最近一次重定位候选是否通过验证。
        int last_relocalization_candidate_id = -1;  ///< 最近一次重定位候选的数据库 ID。
        double last_relocalization_score = 0.0;  ///< 最近一次重定位候选的检索分值。
        double last_relocalization_search_time_ms = 0.0;  ///< 最近一次重定位搜索耗时（ms）。
        std::string last_relocalization_reason;  ///< 最近一次重定位的接受或拒绝原因。
        bool adaptive_lidar_load_enabled = false;  ///< 是否启用多雷达自适应算力预算。
        int adaptive_lidar_degradation_step = 0;  ///< 当前自适应降载级别。
        int selected_lidar_point_stride = 1;  ///< 当前定位点云的点采样步长。
        std::vector<int> selected_lidar_ids;  ///< 当前选入定位计算的雷达 ID。
        std::vector<int> current_frame_lidar_ids;  ///< 当前组帧实际包含的雷达 ID。
        std::vector<int> current_frame_missing_lidar_ids;  ///< 当前组帧缺失的雷达 ID。
        double lidar_correction_age_sec = 0.0;  ///< 当前状态距离最近地图校正的时间（s）。
        double last_lio_processing_ms = 0.0;  ///< 最近一帧 LIO 处理耗时（ms）。
        std::uint64_t adaptive_lidar_stale_drop_count = 0;  ///< 自适应降载因时效不足而丢弃的帧数。
        std::uint32_t cloud_publish_min_lidars = 0;  ///< 点云发布所需的最少雷达数量。
        bool cloud_publish_eligible = true;  ///< 当前帧是否满足点云发布完整性要求。
        bool wheel_speed_dr_enabled = false;  ///< 是否启用轮速 DR 观测。
        WheelSpeedDrStats wheel_speed_dr_stats;  ///< 轮速输入、门控与滤波更新统计。
    };

    struct Options {
        Options() {}

        bool online_mode_ = false;  ///< 在线模式还是离线模式
        bool with_ui_ = false;      ///< 是否带ui

        /// 参数
        SE3 T_body_lidar_;

        bool enable_lidar_odom_skip_ = false;  ///< 是否允许激光里程计跳帧
        int lidar_odom_skip_num_ = 1;          ///< 如果允许跳帧，跳多少帧
        bool enable_lidar_loc_skip_ = true;    ///< 是否允许激光定位跳帧
        bool enable_lidar_loc_rviz_ = false;   ///< 是否允许调试用rviz
        int lidar_loc_skip_num_ = 4;           ///< 如果允许跳帧，跳多少帧
        bool loc_on_kf_ = false;  ///< 是否仅使用关键帧执行地图匹配。
    };

    Localization(Options options = Options());
    ~Localization() = default;

    /**
     * @brief 根据配置和地图创建定位模块，连接处理队列与结果回调。
     * @details 系统层在接收传感器数据前调用；再次初始化时先结束已启动的队列。
     * 初始位姿由后续 SetExternalPose() 设置。
     * @par 主要流程
     * 1. <b>装配算法</b>：初始化 LaserMapping、LidarLoc 和 PGO，并按配置创建 UI。
     * 2. <b>连接输入处理</b>：设置队列容量、跳帧与过载策略，绑定处理函数，在线模式启动线程。
     * 3. <b>连接结果输出</b>：区分全局校正与高频结果，设置输出队列和时间戳门控。
     * 4. <b>完成预处理配置</b>：设置雷达类型、点过滤与高度 ROI，供后续点云输入使用。
     * @htmlonly[block]
     * <details><summary>关键代码：输入与地图匹配的处理入口</summary>
     * @endhtmlonly
     * @snippet core/localization/localization.cpp localization-init-dispatch
     * @htmlonly[block]
     * </details>
     * @endhtmlonly
     * @param yaml_path 定位配置文件路径。
     * @param global_map_path 定位使用的地图包路径。
     * @return 模块和配置准备成功返回 true；显式配置校验或模块初始化失败返回 false。
     * @note YAML 读取或字段转换异常可向调用方传播。
     * @see SetExternalPose(), ProcessSensorInput(), LidarLocProcCloud(), Finish()
     */
    bool Init(const std::string& yaml_path, const std::string& global_map_path);

    /**
     * @brief 使用主雷达 ID 或单雷达 ID 0 接收标准点云。
     * @see ProcessLidarMsg(const sensor_msgs::msg::PointCloud2::SharedPtr, int)
     */
    void ProcessLidarMsg(const sensor_msgs::msg::PointCloud2::SharedPtr laser_msg);
    /**
     * @brief 预处理指定雷达的标准点云，并送入 LIO 数据处理链。
     * @details 先更新过载监测并按策略跳帧，再做点过滤；在线模式进入传感器队列，离线模式直接处理。
     * @param laser_msg 非空点云消息，保留传感器时间戳用于后续同步。
     * @param lidar_id 与多雷达配置对应的传感器 ID。
     * @see ProcessSensorInput(), LidarOdomProcCloud()
     */
    void ProcessLidarMsg(const sensor_msgs::msg::PointCloud2::SharedPtr laser_msg, int lidar_id);
    /**
     * @brief 将 Livox 自定义消息预处理为算法点云，并送入同一条 LIO 处理链。
     * @details 在线模式入队，离线模式直接调用 LidarOdomProcCloud()。
     * @see ProcessLidarMsg(), LidarOdomProcCloud()
     */
    void ProcessLivoxLidarMsg(const livox_ros_driver2::msg::CustomMsg::SharedPtr laser_msg);

    /**
     * @brief 接收 IMU 消息；在线模式与点云共用传感器队列，离线模式直接处理。
     * @details 在线入队后由 ProcessSensorInput() 调用 ProcessIMUData()。
     * @see ProcessSensorInput(), ProcessIMUData()
     */
    void ProcessIMUMsg(IMUPtr imu);

    /**
     * @brief 将 CAN 换算的车体纵向速度交给前端，并更新静止检测使用的轮速历史。
     * @param timestamp 与 IMU 同一时间基准的传感器时间，单位秒。
     * @param longitudinal_speed_mps 有符号纵向速度，单位 m/s。
     * @param motor_torque_nm 电机扭矩，单位 N·m。
     * @see ProcessIMUData()
     */
    void ProcessWheelSpeed(double timestamp, double longitudinal_speed_mps,
                           double motor_torque_nm = 0.0);

    // void ProcessOdomMsg(const nav_msgs::msg::Odometry::SharedPtr odom_msg) override;

    /**
     * @brief 为地图匹配设置外部初始位姿，并重置已有 PGO 状态。
     * @details 可用于启动定位或手动重新定位；不会自行接收传感器数据。
     * @see lightning::LocSystem::SetInitPose(), LidarLocProcCloud()
     */
    void SetExternalPose(const Eigen::Quaterniond& q, const Eigen::Vector3d& t);

    /**
     * @brief 依次处理完传感器、地图匹配和高频输出队列，再结束定位与 UI。
     * @details 调用者应先停止输入，并在对象销毁前调用；队列退出会等待已接收的工作结束。
     * @see Init(), lightning::sys::AsyncMessageProcess::Quit()
     */
    void Finish();

    /**
     * @brief 将预处理点云送入 LIO，并处理已具备同步数据的前端输出。
     * @details 在线模式由传感器队列调用；通过处理锁与 IMU 状态更新保持串行。
     * @see ProcessSensorInput(), DrainLioOutputs()
     */
    void LidarOdomProcCloud(CloudPtr, int lidar_id);
    /**
     * @brief 消费 LIO 已完成的帧，将里程计、点云及同帧状态连接到定位与 PGO。
     * @par 主要流程
     * 1. <b>取前端输出</b>：循环调用 LaserMapping::RunDetailed()，直到无可处理数据。
     * 2. <b>更新运动状态</b>：向地图匹配与 PGO 提交 LIO 里程计，并更新运行统计。
     * 3. <b>提交匹配帧</b>：按关键帧策略选择点云；在线模式入匹配队列，离线模式直接处理。
     * @see LidarOdomProcCloud(), ProcessIMUData(), LidarLocProcCloud()
     */
    void DrainLioOutputs();
    /// @brief 一帧地图匹配输入，关联定位点云、发布点云、雷达状态及处理时刻。
    struct LidarLocInput {
        CloudPtr localization_cloud;  ///< 用于地图匹配的当前帧点云。
        CloudPtr publication_cloud;  ///< 同帧用于上层发布的完整点云。
        MultiLidarFrameStats frame_stats;  ///< 该帧的多雷达完整性与点数统计。
        bool publication_eligible = true;  ///< 该帧是否满足点云发布的雷达数量与主雷达要求。
        std::chrono::steady_clock::time_point primary_received_at{};  ///< 主雷达帧到达的单调时钟时刻，用于端到端耗时。
        std::chrono::steady_clock::time_point enqueued_at{};  ///< 消息进入当前队列的单调时钟时刻，用于排队耗时。
    };
    /**
     * @brief 执行一帧地图匹配，将匹配结果连接到点云发布和 PGO 融合。
     * @details 在线模式由独立匹配队列调用，离线模式由 DrainLioOutputs() 直接调用。
     * @par 主要流程
     * 1. <b>匹配地图</b>：调用 LidarLoc::ProcessCloud()，读取结果并更新重定位统计。
     * 2. <b>处理恢复状态</b>：重定位接受后重置 PGO；静止地图冲突时标记无效并请求恢复。
     * 3. <b>交付结果</b>：调用点云回调，正常结果送入 PGO，并更新 UI、状态回调与耗时记录。
     * @see DrainLioOutputs(), SetProcessedCloudCallback(), SetLocalizationResultCallback()
     */
    void LidarLocProcCloud(const LidarLocInput& input);

    using TFCallback = std::function<void(const geometry_msgs::msg::TransformStamped& odom)>;
    using LocalizationResultCallback = std::function<void(const LocalizationResult& result)>;
    using ProcessedCloudCallback = std::function<void(const CloudPtr& cloud, const LocalizationResult& result,
                                                       const MultiLidarFrameStats& stats, bool eligible)>;
    using LocStateCallback = std::function<void(const std_msgs::msg::Int32& state)>;
    using PointcloudBodyCallback = std::function<void(const sensor_msgs::msg::PointCloud2& pointcloud)>;
    using PointcloudWorldCallback = std::function<void(const sensor_msgs::msg::PointCloud2& pointcloud)>;

    /// @brief 注册有效实时位姿的 TF 回调；在线模式在高频输出线程中执行。
    void SetTFCallback(TFCallback&& callback);
    /// @brief 注册时间戳门控后的实时结果回调，供系统层执行业务发布。
    void SetLocalizationResultCallback(LocalizationResultCallback&& callback);
    /// @brief 注册 PGO 全局校正回调，用于评估与轨迹记录，与实时发布路径分开。
    void SetGlobalLocalizationResultCallback(LocalizationResultCallback&& callback);
    /// @brief 注册同帧点云与匹配结果回调；在线模式在地图匹配线程中执行。
    void SetProcessedCloudCallback(ProcessedCloudCallback&& callback);
    /// @brief 注册每帧地图匹配后的定位状态回调。
    void SetLocStateCallback(LocStateCallback&& callback);

    bool IsMultiLidarEnabled() const;
    const MultiLidarConfig& GetMultiLidarConfig() const;
    SO3 GetInitialLidarRotation() const;
    SO3 GetImuToBodyRotation() const;
    /**
     * @brief 读取队列计数、传感器滞后与当前定位运行状态。
     * @note 各来源分别采样，结果不是所有线程在同一时刻的原子快照。
     * @see RuntimeStats, lightning::LocSystem::PublishHealthStatus()
     */
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
    std::mutex input_mutex_;  ///< 串行保护输入点云预处理及相关临时状态。
    // Serialize ordered LIO/PGO state updates without blocking ROS callbacks.
    std::mutex processing_mutex_;  ///< 串行保护 LIO 与 PGO 的有序状态更新。
    // Keep the timestamp check and all live callbacks in the same critical
    // section so two producer threads cannot publish in reverse order.
    std::mutex live_output_dispatch_mutex_;  ///< 保护实时输出时间戳检查及对应回调派发。
    MonotonicTimestampGate live_output_timestamp_gate_;  ///< 拒绝重复或回退时间戳的实时定位输出。
    Options options_;  ///< 定位模式、跳帧与可视化配置。

    /// 预处理
    std::shared_ptr<PointCloudPreprocess> preprocess_ = nullptr;  // point cloud preprocess

    /// 前端
    std::shared_ptr<LaserMapping> lio_ = nullptr;
    Keyframe::Ptr lio_kf_ = nullptr;  ///< 最近已交给定位流程的 LIO 关键帧。

    // ui
    std::shared_ptr<ui::PangolinWindow> ui_ = nullptr;  ///< 可选的定位状态与点云显示窗口。

    // pose graph
    std::shared_ptr<PGO> pgo_ = nullptr;  ///< 融合 LIO、DR 与地图匹配约束的位姿图模块。

    // lidar localization
    std::shared_ptr<LidarLoc> lidar_loc_;  ///< 地图匹配、跟踪与全局重定位模块。

    /// TODO async 处理
    struct SensorInput {
        IMUPtr imu;  ///< IMU 输入载荷，在 is_imu 为 true 时使用。
        CloudPtr cloud;  ///< 点云输入载荷，在 is_imu 为 false 时使用。
        int lidar_id = 0;  ///< 该输入对应的雷达配置 ID。
        bool is_imu = false;  ///< 区分统一传感器队列中的 IMU 与点云消息。
        std::chrono::steady_clock::time_point enqueued_at{};  ///< 消息进入当前队列的单调时钟时刻，用于排队耗时。
    };
    /**
     * @brief 从传感器队列取一条输入，分派到 IMU 或 LIO 点云处理，并更新完成进度。
     * @note queue_wait 记录入队到进入处理回调的耗时，不等同于操作系统调度等待时间。
     * @see ProcessIMUData(), LidarOdomProcCloud(),
     *      @ref resource_diagnostics "资源与调度诊断"
     */
    void ProcessSensorInput(const SensorInput& input);
    /**
     * @brief 校验 IMU、更新前端预测，并将可用航迹推算状态送入定位与 PGO。
     * @par 主要流程
     * 1. <b>检查输入</b>：过滤无效数据、重复或回退时间戳，串行进入算法处理。
     * 2. <b>更新前端</b>：提交 IMU，并先处理由该 IMU 补齐同步条件的 LIO 帧。
     * 3. <b>更新航迹推算</b>：位姿可用时处理静止状态，再向 LidarLoc 和 PGO 提交 DR 状态。
     * @see ProcessIMUMsg(), DrainLioOutputs(), ProcessWheelSpeed()
     */
    void ProcessIMUData(IMUPtr imu);
    /**
     * @brief 根据传感器处理滞后更新过载与恢复状态，供日志和诊断使用。
     * @note 当前实现始终返回 false，保持接收 LiDAR；积压旧输入由有界队列淘汰。
     * @see ProcessLidarMsg(), lightning::sys::AsyncMessageProcess::AddMessage()
     */
    bool ShouldThrottleLidarInput(double timestamp);
    void ObserveSensorEnqueued(double timestamp);
    void ObserveSensorProcessed(double timestamp);
    void ObserveLioForStaticDetector(const NavState& state);
    void ObserveLidarLocForStaticDetector(const LocalizationResult& result);
    bool UpdateImuStaticState(const IMUPtr& imu);
    sys::AsyncMessageProcess<SensorInput> sensor_proc_;  ///< 统一串行处理 IMU 与点云输入的有界队列。
    sys::AsyncMessageProcess<LidarLocInput> lidar_loc_proc_cloud_;  ///< 定位点云和同帧完整发布点云
    // ROS/DDS publication must never hold the ordered sensor/PGO processing
    // path. Capacity one intentionally keeps only the newest live pose.
    sys::AsyncMessageProcess<LocalizationResult> high_frequency_output_proc_;  ///< 异步发布高频位姿的队列，容量一保留最新结果。
    profiling::MultiStageTimingWindow lidar_input_timing_window_{
        {"preprocess", "dispatch", "outer"}};  ///< 点云预处理和分派耗时统计。
    profiling::MultiStageTimingWindow imu_dr_timing_window_{
        {"lio_imu", "drain_lio", "state_static", "lidar_loc_dr", "pgo_dr", "outer"}};  ///< IMU 预测、LIO 排出及 DR 处理耗时统计。
    // Single sensor consumer owns these windows; queue age uses monotonic time,
    // not the sensor Header clock. It excludes processing after callback entry.
    profiling::MultiStageTimingWindow imu_queue_timing_window_{{"queue_wait"}};  ///< IMU 入队到进入处理回调的等待耗时，使用单调时钟。
    profiling::MultiStageTimingWindow lidar_queue_timing_window_{{"queue_wait"}};  ///< 点云入队到进入处理回调的等待耗时，使用单调时钟。
    // Bounded primary-frame arrival ledger bridges preprocessing, sensor queue,
    // LIO and the asynchronous map matcher. Its clock never uses ROS wall time.
    void RememberPrimaryArrival(int lidar_id, std::uint64_t stamp,
                                std::chrono::steady_clock::time_point received_at);
    std::mutex arrival_mutex_;  ///< 保护主雷达帧到达时刻索引。
    std::map<std::uint64_t, std::chrono::steady_clock::time_point> primary_arrivals_;  ///< 按主雷达帧时间戳索引单调时钟接收时刻。
    profiling::MultiStageTimingWindow lidar_end_to_end_window_{{"primary_to_pgo"}, 10.0};  ///< 主雷达接收到地图匹配结果进入 PGO 的端到端耗时统计。
    std::size_t lidar_deadline_misses_ = 0;  ///< 主雷达端到端处理超出帧周期预算的次数。
    int lidar_odom_skip_cnt_ = 0;  ///< 激光里程计输入跳帧计数。
    double online_sensor_max_lag_sec_ = 0.0;  ///< 进入过载状态的传感器处理滞后阈值（s）。
    double online_sensor_resume_lag_sec_ = 0.0;  ///< 从过载状态恢复的传感器处理滞后阈值（s）。
    std::atomic_bool lidar_overload_throttled_{false};  ///< 过载诊断状态；当前实现仍接收点云，由队列容量淘汰旧输入。

    struct StaticImuSample {
        double timestamp = 0.0;  ///< 该观测的传感器时间戳（s）。
        Vec3d angular_velocity = Vec3d::Zero();  ///< IMU 角速度向量（rad/s），用于静止窗口均值计算。
        double gyro_norm = 0.0;  ///< IMU 角速度模长（rad/s）。
        double accel_norm = 0.0;  ///< IMU 加速度模长（m/s²）。
    };
    mutable std::mutex static_detector_mutex_;  ///< 保护静止检测窗口、统计量与轮速历史。
    std::deque<StaticImuSample> static_imu_window_;  ///< 静止判定使用的最近 IMU 样本窗口。
    double static_gyro_sum_ = 0.0;  ///< 窗口内角速度模长之和。
    double static_gyro_sq_sum_ = 0.0;  ///< 窗口内角速度模长的平方和，用于计算方差。
    Vec3d static_gyro_vector_sum_ = Vec3d::Zero();  ///< 窗口内角速度向量之和，用于计算均值。
    double static_accel_sum_ = 0.0;  ///< 窗口内加速度模长之和。
    double static_accel_sq_sum_ = 0.0;  ///< 窗口内加速度模长的平方和，用于计算方差。
    double last_static_lio_stamp_ = 0.0;  ///< 静止检测最近观察到的 LIO 时间戳（s）。
    double last_static_lio_speed_ = 0.0;  ///< 静止检测最近观察到的 LIO 速度模长（m/s）。
    bool last_static_lio_reliable_ = false;  ///< 最近 LIO 速度是否可信，可否参与静止判定。
    double last_valid_lidar_loc_stamp_ = 0.0;  ///< 最近有效地图匹配的时间戳（s）。
    struct StaticWheelSpeedSample {
        double timestamp = 0.0;  ///< 该观测的传感器时间戳（s）。
        double speed_mps = 0.0;  ///< 静止判定使用的轮速（m/s）。
    };
    std::deque<StaticWheelSpeedSample> static_wheel_speed_history_;  ///< 静止判定使用的带时间戳轮速历史。
    double last_wheel_speed_stamp_ = 0.0;  ///< 最近接收的轮速时间戳（s）。
    bool wheel_speed_observed_ = false;  ///< 是否已有可供静止判定的轮速观测。
    bool wheel_timestamp_mismatch_reported_ = false;  ///< 轮速与算法时间域不一致的告警是否已报告。
    int static_exit_count_ = 0;  ///< 连续满足静止退出条件的次数。
    bool imu_static_hold_enabled_ = false;  ///< 是否启用 IMU 高频输出的静止保持。
    bool imu_static_hold_active_ = false;  ///< 当前是否处于 IMU 高频输出静止保持状态。

    /// 结果数据 =====================================================================================================
    LocalizationResult loc_result_;

    /// 框架相关
    TFCallback tf_callback_;
    LocalizationResultCallback localization_result_callback_;  ///< 高频定位结果的上层输出回调。
    LocalizationResultCallback global_localization_result_callback_;  ///< PGO 全局校正结果的上层输出回调。
    ProcessedCloudCallback processed_cloud_callback_;  ///< 完整处理后点云、定位结果与组帧统计的输出回调。
    LocStateCallback loc_state_callback_;  ///< 定位状态变化的输出回调。
    PointcloudBodyCallback pointcloud_body_callback_;  ///< 机体坐标系点云的可选输出回调。
    PointcloudWorldCallback pointcloud_world_callback_;  ///< 世界坐标系点云的可选输出回调。

    /// 输入检查
    double last_imu_time_ = 0;
    double last_odom_time_ = 0;  ///< 最近接受的里程计时间戳（s）。
    double last_cloud_time_ = 0;  ///< 最近处理的点云时间戳（s）。

    mutable std::mutex runtime_stats_mutex_;  ///< 保护运行诊断统计。
    RuntimeStats runtime_stats_;  ///< 输入进度、队列、重定位与降载诊断快照。
};
}  // namespace loc

}  // namespace lightning
