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

// 前置声明。
namespace lightning {
namespace ui {
class PangolinWindow;
}

namespace loc {

class LidarLoc;
class PGO;

/**
 * @brief 组织 LIO、地图匹配与 PGO，将传感器输入转换为定位结果。
 * @details 在线入口通过 lightning::LocSystem 创建本对象，由系统层负责 ROS 订阅与业务发布。
 * 本类在线模式通过传感器、地图匹配和高频输出队列执行处理；同步模式直接调用相应算法函数。
 * 当前 run_loc_offline.cc 独立组装底层 LIO、地图匹配和 PGO，并未直接使用本类。
 * @note 输入结束后、对象销毁前必须调用 Finish()。回调应在开始输入前注册，
 * 其执行线程由处理路径决定；本类不保证所有公开接口可与 Init()/Finish() 并发调用。
 * @see Init(), ProcessIMUMsg(), ProcessLidarMsg(), Finish(),
 *      @ref localization_module "定位与重定位模块"
 */
class Localization {
   public:
    /**
     * @brief 供诊断接口读取的队列、传感器进度、重定位及雷达负载统计。
     * @details 队列完成数表示处理回调已返回，不代表数据被算法接受或已执行业务发布。
     * 时间戳使用传感器时间基准；GetRuntimeStats() 分别采样各来源，不提供跨线程原子快照。
     */
    struct RuntimeStats {
        std::uint64_t sensor_queue_pending = 0;  ///< 传感器队列待执行与正在执行的消息数。
        std::uint64_t sensor_queue_dropped = 0;  ///< 传感器队列容量淘汰及显式记录的丢弃数，不含策略跳帧。
        std::uint64_t sensor_queue_processed = 0;  ///< 传感器队列已完成的消息数。
        std::uint64_t localization_queue_pending = 0;  ///< 地图匹配队列待执行与正在执行的消息数。
        std::uint64_t localization_queue_dropped = 0;  ///< 地图匹配队列容量淘汰及显式记录的丢弃数，不含策略跳帧。
        std::uint64_t localization_queue_processed = 0;  ///< 地图匹配队列已完成的消息数。
        std::uint64_t high_frequency_queue_pending = 0;  ///< 高频输出队列待执行与正在执行的消息数。
        std::uint64_t high_frequency_queue_dropped = 0;  ///< 高频输出队列容量淘汰及显式记录的丢弃数，不含时间戳门控拒绝。
        std::uint64_t high_frequency_queue_processed = 0;  ///< 高频输出队列已完成的消息数。
        double latest_enqueued_sensor_stamp = 0.0;  ///< 在线入队时记录的最大正传感器时间戳（s），未记录时为 0。
        double latest_processed_sensor_stamp = 0.0;  ///< 传感器处理回调返回后记录的最大正时间戳（s），不等同于有效定位时间。
        double current_sensor_lag_sec = 0.0;  ///< 入队与完成时间戳水位差的非负值（s）；尚无完成记录时为 0，不是墙钟排队耗时。
        double max_sensor_lag_sec = 0.0;  ///< 运行期间观察到的最大传感器处理滞后（s）。
        std::uint64_t severe_timestamp_rollback_count = 0;  ///< 已完成输入比完成时间戳水位落后超过 1 s 的累计次数。
        double worst_timestamp_rollback_sec = 0.0;  ///< 上述严重回退中记录的最大回退量（s）。
        std::uint64_t live_output_non_monotonic_drop_count = 0;  ///< 实时输出因时间戳无效、重复或回退而被门控拒绝的次数。
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
        double lidar_correction_age_sec = 0.0;  ///< LIO 最近处理的点云帧末时刻相对最新传感器输入的滞后（s），不是地图校正间隔。
        double last_lio_processing_ms = 0.0;  ///< 最近一帧 LIO 处理耗时（ms）。
        std::uint64_t adaptive_lidar_stale_drop_count = 0;  ///< 自适应降载因时效不足而丢弃的帧数。
        std::uint32_t cloud_publish_min_lidars = 0;  ///< 点云发布所需的最少雷达数量。
        bool cloud_publish_eligible = true;  ///< 当前帧是否满足点云发布完整性要求。
        bool wheel_speed_dr_enabled = false;  ///< 是否启用轮速 DR 观测。
        WheelSpeedDrStats wheel_speed_dr_stats;  ///< 轮速输入、门控与滤波更新统计。
    };

    /// @brief 执行模式及处理选项；Init() 会从 YAML 更新可视化、跳帧和关键帧配置。
    struct Options {
        Options() {}

        bool online_mode_ = false;  ///< true 启动异步处理队列，false 在调用线程同步处理。
        bool with_ui_ = false;      ///< 是否创建 Pangolin 窗口；Init() 从 system.with_ui 读取。

        SE3 T_body_lidar_;  ///< 保留的外参字段；当前 Localization 实现未读写此字段，不参与坐标变换。

        bool enable_lidar_odom_skip_ = false;  ///< 是否在点云预处理后、提交 LIO 前按计数筛选输入帧。
        int lidar_odom_skip_num_ = 1;          ///< 里程计输入采纳周期 N：每 N 帧保留一帧，非正值按 1 处理。
        bool enable_lidar_loc_skip_ = true;    ///< 是否启用在线地图匹配队列的周期跳帧；同步路径不使用此队列策略。
        bool enable_lidar_loc_rviz_ = false;   ///< 从配置读取的保留开关；当前实现未据此创建 RViz 发布通道。
        int lidar_loc_skip_num_ = 4;           ///< 在线地图匹配采纳周期 N：开启跳帧时须为正，每 N 帧保留一帧。
        bool loc_on_kf_ = false;  ///< 是否仅使用关键帧执行地图匹配。
    };

    /// @brief 保存执行选项；算法对象及队列在 Init() 中装配。
    Localization(Options options = Options());
    /// @brief 默认析构；所有者须提前调用 Finish() 排空并回收工作线程。
    ~Localization() = default;

    /**
     * @brief 根据配置和地图创建定位模块，连接处理队列与结果回调。
     * @details 系统层在接收传感器数据前调用；再次初始化时先结束已启动的队列。
     * 初始位姿由后续 SetExternalPose() 设置。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>结束旧队列并重置门控</b>：持有生命周期锁准备模块；已有模块时先释放锁再等待旧队列结束。
     * 2. <b>装配 LIO 与地图定位</b>：读取配置，初始化 LIO、地图定位及可选 UI。
     * 3. <b>配置 PGO</b>：创建 PGO，设置 DR 平滑、外推与协方差过程噪声。
     * 4. <b>配置并启动输入队列</b>：设置跳帧、静止检测和有界队列，绑定处理函数；在线模式启动输入及匹配线程。
     * 5. <b>连接全局与实时输出</b>：定义带时间戳门控的实时回调，分开全局结果；在线高频输出走容量一队列。
     * 6. <b>完成点云预处理配置</b>：设置雷达类型、过滤参数和高度 ROI；在线模式下输入队列此时已经启动。
     * @htmlonly[block]
     * <details><summary>关键代码：输入与地图匹配的处理入口</summary>
     * @endhtmlonly
     * @snippet core/localization/localization.cpp localization-init-dispatch
     * @htmlonly[block]
     * </details>
     * @endhtmlonly
     * @param yaml_path 定位配置文件路径。
     * @param global_map_path 定位使用的地图包路径。
     * @pre 调用者尚未开始输入，重新初始化前应停止外部输入及其他公开接口的调用。
     * @return 模块和配置准备成功返回 true；显式配置校验或模块初始化失败返回 false。
     * @note YAML 读取或字段转换异常可向调用方传播。
     * @see SetExternalPose(), ProcessSensorInput(), LidarLocProcCloud(), Finish()
     */
    bool Init(const std::string& yaml_path, const std::string& global_map_path);

    /**
     * @brief 使用主雷达 ID 或单雷达 ID 0 接收标准点云。
     * @param laser_msg 非空标准点云；后续预处理和分派见指定雷达 ID 的重载。
     * @see ProcessLidarMsg(const sensor_msgs::msg::PointCloud2::SharedPtr, int)
     */
    void ProcessLidarMsg(const sensor_msgs::msg::PointCloud2::SharedPtr laser_msg);
    /**
     * @brief 预处理指定雷达的标准点云，并送入 LIO 数据处理链。
     * @details 先更新过载监测并做点过滤，再应用里程计跳帧策略；在线模式入队，离线模式直接处理。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>检查模块并串行进入输入处理</b>：记录到达时间和过载状态，锁定共享预处理器并检查模块是否就绪。
     * 2. <b>预处理并保留帧时间</b>：登记主雷达到达时间，转换标准点云并保留消息时间戳。
     * 3. <b>应用里程计跳帧策略</b>：在预处理之后推进跳帧计数；被跳过的帧不再入队。
     * 4. <b>分派点云并记录耗时</b>：在线加入统一传感器队列，离线释放输入锁后直接处理，再汇总本次耗时。
     * @param laser_msg 非空点云消息，保留传感器时间戳用于后续同步。
     * @param lidar_id 与多雷达配置对应的传感器 ID。
     * @see ProcessSensorInput(), LidarOdomProcCloud()
     */
    void ProcessLidarMsg(const sensor_msgs::msg::PointCloud2::SharedPtr laser_msg, int lidar_id);
    /**
     * @brief 将 Livox 自定义消息预处理为算法点云，并送入同一条 LIO 处理链。
     * @details 在线模式入队，离线模式直接调用 LidarOdomProcCloud()。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>检查模块并串行进入输入处理</b>：记录到达时间和过载状态，检查模块并锁定共享预处理器。
     * 2. <b>预处理并确定雷达身份</b>：转换 Livox 点云，保留帧时间，并选择主雷达或单雷达 ID。
     * 3. <b>应用里程计跳帧策略</b>：按计数决定是否将已预处理的帧继续送入 LIO。
     * 4. <b>分派点云并记录耗时</b>：在线入传感器队列，离线直接处理，并记录预处理和分派耗时。
     * @param laser_msg 非空 Livox 点云消息，时间戳用于后续多传感器同步。
     * @see ProcessLidarMsg(), LidarOdomProcCloud()
     */
    void ProcessLivoxLidarMsg(const livox_ros_driver2::msg::CustomMsg::SharedPtr laser_msg);

    /**
     * @brief 接收 IMU 消息；在线模式与点云共用传感器队列，离线模式直接处理。
     * @details 在线入队后由 ProcessSensorInput() 调用 ProcessIMUData()。
     * @param imu 内部 IMU 量测；实际处理时拒绝空指针、无效数值、重复或回退时间戳。
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
     * @param q 地图定位初始位姿的旋转，与当前地图定位坐标约定一致。
     * @param t 地图定位初始位姿的平移，单位 m。
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
     * @note 在线输入应经 ProcessLidarMsg()/ProcessLivoxLidarMsg() 入队，保持与 IMU 的处理次序。
     * @see ProcessSensorInput(), DrainLioOutputs()
     */
    void LidarOdomProcCloud(CloudPtr, int lidar_id);
    /**
     * @brief 消费 LIO 已完成的帧，将里程计、点云及同帧状态连接到定位与 PGO。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>取出可用 LIO 输出</b>：循环运行前端；无数据时结束，无输出的处理状态继续下一轮。
     * 2. <b>更新里程计状态</b>：将同帧 LO 交给静止检测、地图定位与 PGO。
     * 3. <b>整理同帧数据与统计</b>：关联定位点云、发布点云、雷达完整性与主雷达到达时间。
     * 4. <b>选择并提交地图匹配帧</b>：按关键帧策略筛选；在线入匹配队列，离线直接匹配，随后继续循环。
     * @pre 由已持有生命周期共享锁和 processing_mutex_ 的 IMU 或点云处理路径调用。
     * @see LidarOdomProcCloud(), ProcessIMUData(), LidarLocProcCloud()
     */
    void DrainLioOutputs();
    /// @brief 一帧地图匹配输入，关联定位点云、发布点云、雷达状态及处理时刻。
    struct LidarLocInput {
        CloudPtr localization_cloud;  ///< 用于地图匹配的当前帧点云。
        CloudPtr publication_cloud;  ///< 同帧用于上层发布的完整点云。
        MultiLidarFrameStats frame_stats;  ///< 该帧的多雷达完整性与点数统计。
        bool publication_eligible = true;  ///< 该帧是否满足点云发布的雷达数量与主雷达要求。
        std::chrono::steady_clock::time_point primary_received_at{};  ///< 可关联主雷达帧的单调时钟到达时刻；默认值表示未找到，不能计算端到端耗时。
        std::chrono::steady_clock::time_point enqueued_at{};  ///< 启用计算剖析时记录的单调时钟入队时刻；默认值表示未测量排队耗时。
    };
    /**
     * @brief 执行一帧地图匹配，将匹配结果连接到点云发布和 PGO 融合。
     * @details 在线模式由独立匹配队列调用，离线模式由 DrainLioOutputs() 直接调用。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>检查模块并执行地图匹配</b>：记录队列等待，在生命周期保护下调用 LidarLoc::ProcessCloud()。
     * 2. <b>读取结果并检查恢复状态</b>：记录重定位统计，接受重定位时重置 PGO，静止地图冲突时将结果置为无效。
     * 3. <b>交付点云与当前缓存结果</b>：把点云、结果及发布资格交给回调；冲突状态也在这里通知上层。
     * 4. <b>更新 PGO 或请求恢复</b>：冲突时重置 PGO 并请求全局重定位，否则向 PGO 提交地图匹配结果。
     * 5. <b>记录时延并通知状态</b>：更新端到端统计、UI 和状态回调，释放生命周期锁后记录计算耗时。
     * @note 此处未检查 ProcessCloud() 的返回值；点云过少而早退时，缓存结果可能仍来自上一帧，不能保证与本帧点云对应。
     * @see DrainLioOutputs(), SetProcessedCloudCallback(), SetLocalizationResultCallback()
     */
    void LidarLocProcCloud(const LidarLocInput& input);

    /// @brief 接收有效实时位姿对应的变换，由系统层决定是否广播 TF。
    using TFCallback = std::function<void(const geometry_msgs::msg::TransformStamped& odom)>;
    /// @brief 接收定位结果；实时结果和 PGO 全局结果通过不同的注册接口分别交付。
    using LocalizationResultCallback = std::function<void(const LocalizationResult& result)>;
    /// @brief 交付发布点云、当前地图匹配结果、组帧统计及发布资格；结果可能来自上一帧缓存。
    using ProcessedCloudCallback = std::function<void(const CloudPtr& cloud, const LocalizationResult& result,
                                                       const MultiLidarFrameStats& stats, bool eligible)>;
    /// @brief 接收地图匹配后的定位状态枚举值，并非仅在状态变化时触发。
    using LocStateCallback = std::function<void(const std_msgs::msg::Int32& state)>;
    /// @brief 保留的机体坐标系点云回调类型；当前本类没有启用对应派发路径。
    using PointcloudBodyCallback = std::function<void(const sensor_msgs::msg::PointCloud2& pointcloud)>;
    /// @brief 保留的世界坐标系点云回调类型；当前本类没有启用对应派发路径。
    using PointcloudWorldCallback = std::function<void(const sensor_msgs::msg::PointCloud2& pointcloud)>;

    // 以下设置器仅保存回调，应在开始输入前调用；回调执行与注册是两个时刻。
    /// @brief 注册有效实时位姿的 TF 回调；在线模式在高频输出线程中执行。
    void SetTFCallback(TFCallback&& callback);
    /// @brief 注册时间戳门控后的实时结果回调；在线高频输出线程调用，同步模式由生产者直接调用。
    void SetLocalizationResultCallback(LocalizationResultCallback&& callback);
    /// @brief 注册 PGO 全局校正回调，用于评估与轨迹记录，与实时发布路径分开。
    void SetGlobalLocalizationResultCallback(LocalizationResultCallback&& callback);
    /// @brief 注册发布点云与当前匹配结果回调；在线模式在地图匹配线程中执行，结果时效边界见 LidarLocProcCloud()。
    void SetProcessedCloudCallback(ProcessedCloudCallback&& callback);
    /// @brief 注册每帧地图匹配后的定位状态回调。
    void SetLocStateCallback(LocStateCallback&& callback);

    /// @brief 查询已创建的 LIO 是否启用多雷达；尚未创建 LIO 时返回 false。
    bool IsMultiLidarEnabled() const;
    /// @brief 返回 LIO 持有的多雷达配置引用；要求 LIO 已创建，且使用期间不重新初始化。
    const MultiLidarConfig& GetMultiLidarConfig() const;
    /// @brief 返回雷达坐标系到初始化车体坐标系的旋转，用于输出坐标转换；要求 LIO 已创建。
    SO3 GetInitialLidarRotation() const;
    /// @brief 返回 IMU 坐标系到初始化车体坐标系的旋转；要求 LIO 已创建。
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

    // 生命周期、输入顺序与输出派发。
    /**
     * @brief 协调算法模块初始化、替换及停止与正在进行的输入/算法处理。
     * @details Init() 持独占锁创建或替换模块；输入及算法处理路径持共享锁保持模块可用。
     * Finish() 先在不持此锁时排空并回收三个工作线程，再持独占锁结束地图定位与 UI。
     * 共享锁不负责串行化算法状态更新；该职责由处理队列及 processing_mutex_ 等锁承担。
     * @see Init(), Finish(), input_mutex_, processing_mutex_
     */
    mutable std::shared_mutex lifecycle_mutex_;
    std::mutex input_mutex_;  ///< 串行保护两种点云入口共享的预处理器、里程计跳帧计数及入队操作。
    std::mutex processing_mutex_;  ///< 串行化 IMU、LIO 点云及外部位姿处理；地图匹配提交 PGO 时使用模块内部锁。
    std::mutex live_output_dispatch_mutex_;  ///< 将实时结果的时间戳门控、缓存更新及 TF/结果回调放在同一临界区。
    MonotonicTimestampGate live_output_timestamp_gate_;  ///< 实时结果的严格递增门控，拒绝非有限、非正、重复或回退时间戳。
    Options options_;  ///< 定位模式、跳帧与可视化配置。

    // 算法模块及显示窗口。
    std::shared_ptr<PointCloudPreprocess> preprocess_ = nullptr;  ///< 两种点云消息共用的预处理器，使用时由 input_mutex_ 串行保护。

    std::shared_ptr<LaserMapping> lio_ = nullptr;  ///< 定位模式的 LIO 前端，提供里程计、IMU 预测及去畸变/投影点云。
    Keyframe::Ptr lio_kf_ = nullptr;  ///< 关键帧筛选模式下最近提交地图匹配的 LIO 关键帧，用于去重。

    std::shared_ptr<ui::PangolinWindow> ui_ = nullptr;  ///< 可选的定位状态与点云显示窗口。

    std::shared_ptr<PGO> pgo_ = nullptr;  ///< 融合 LIO、DR 与地图匹配约束的位姿图模块。

    std::shared_ptr<LidarLoc> lidar_loc_;  ///< 地图匹配、跟踪与全局重定位模块。

    /// @brief 统一传感器队列的单条输入，使用 is_imu 区分 IMU 与预处理点云。
    struct SensorInput {
        IMUPtr imu;  ///< IMU 输入载荷，在 is_imu 为 true 时使用。
        CloudPtr cloud;  ///< 点云输入载荷，在 is_imu 为 false 时使用。
        int lidar_id = 0;  ///< 该输入对应的雷达配置 ID。
        bool is_imu = false;  ///< 区分统一传感器队列中的 IMU 与点云消息。
        std::chrono::steady_clock::time_point enqueued_at{};  ///< 启用计算剖析时记录的单调时钟入队时刻；默认值表示未测量排队耗时。
    };
    /**
     * @brief 处理传感器队列交付的一条输入，分派到 IMU 或 LIO 点云处理并更新完成进度。
     * @details Init() 将本函数注册给 sensor_proc_；在线消费线程出队后调用，同步入口绕过此分派。
     * @note queue_wait 记录入队到进入处理回调的耗时，不等同于操作系统调度等待时间。
     * 完成进度只说明被调处理函数已返回，不保证输入被接受或产生定位结果。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>记录队列等待</b>：使用单调时钟统计入队到进入处理回调的等待。
     * 2. <b>分派 IMU 或点云</b>：按载荷类型调用 IMU 或 LIO 点云处理函数。
     * 3. <b>推进完成进度</b>：被调处理函数返回后更新进度，并按需记录因果追踪。
     * @see ProcessIMUData(), LidarOdomProcCloud(),
     *      @ref resource_diagnostics "资源与调度诊断"
     */
    void ProcessSensorInput(const SensorInput& input);
    /**
     * @brief 校验 IMU、更新前端预测，并将可用航迹推算状态送入定位与 PGO。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>串行检查 IMU 输入</b>：保护模块与有序算法状态，拒绝无效数值、重复时间和回退时间。
     * 2. <b>更新预测并排出 LIO 帧</b>：提交 IMU 后先排出已满足同步条件的 LIO 帧，保持 PGO 的传感器时间顺序。
     * 3. <b>检查预测位姿并更新静止状态</b>：取得 IMU 预测；位姿不可用时结束本次处理，否则应用保守静止检测。
     * 4. <b>分别交付地图定位与 PGO</b>：地图定位收到的 DR 清除停车标记以继续真实匹配，PGO 保留静止状态，再记录耗时。
     * @see ProcessIMUMsg(), DrainLioOutputs(), ProcessWheelSpeed()
     */
    void ProcessIMUData(IMUPtr imu);
    /**
     * @brief 根据传感器处理滞后更新过载与恢复状态，供日志和诊断使用。
     * @note 当前实现始终返回 false，保持接收 LiDAR；积压旧输入由有界队列淘汰。
     * @see ProcessLidarMsg(), lightning::sys::AsyncMessageProcess::AddMessage()
     */
    bool ShouldThrottleLidarInput(double timestamp);
    /// @brief 以正传感器时间戳更新在线入队水位及滞后统计，水位只增不减。
    void ObserveSensorEnqueued(double timestamp);
    /// @brief 在输入处理返回后更新完成水位，并记录超过 1 s 的回退；不代表算法接受该输入。
    void ObserveSensorProcessed(double timestamp);
    /// @brief 静止保持启用时，缓存最近 LIO 的时间、速度模长及可靠标志。
    void ObserveLioForStaticDetector(const NavState& state);
    /// @brief 静止保持启用时，缓存最近有效地图匹配时间供新鲜度检查。
    void ObserveLidarLocForStaticDetector(const LocalizationResult& result);
    /// @brief 根据 IMU 窗口、LIO、地图匹配新鲜度和因果轮速更新并返回静止保持状态。
    bool UpdateImuStaticState(const IMUPtr& imu);

    // 仅在线模式启动的异步队列，以及分阶段时延统计。
    sys::AsyncMessageProcess<SensorInput> sensor_proc_;  ///< 按入队顺序串行交付 IMU/点云给 ProcessSensorInput()，不按传感器时间重排。
    sys::AsyncMessageProcess<LidarLocInput> lidar_loc_proc_cloud_;  ///< 将定位/发布点云及同帧统计交给 LidarLocProcCloud()；等待容量默认 1，可配置。
    sys::AsyncMessageProcess<LocalizationResult> high_frequency_output_proc_;  ///< 将 PGO 高频结果交给实时输出回调；等待容量 1，保留最新待处理结果。
    profiling::MultiStageTimingWindow lidar_input_timing_window_{
        {"preprocess", "dispatch", "outer"}};  ///< 点云预处理和分派耗时统计。
    profiling::MultiStageTimingWindow imu_dr_timing_window_{
        {"lio_imu", "drain_lio", "state_static", "lidar_loc_dr", "pgo_dr", "outer"}};  ///< IMU 预测、LIO 排出及 DR 处理耗时统计。
    // 两个等待窗口由同一传感器消费线程更新，不含回调开始后的算法耗时。
    profiling::MultiStageTimingWindow imu_queue_timing_window_{{"queue_wait"}};  ///< IMU 入队到进入处理回调的等待耗时，使用单调时钟。
    profiling::MultiStageTimingWindow lidar_queue_timing_window_{{"queue_wait"}};  ///< 点云入队到进入处理回调的等待耗时，使用单调时钟。
    /// @brief 登记主雷达帧的单调时钟到达时刻，最多保留 256 个传感器时间戳索引。
    void RememberPrimaryArrival(int lidar_id, std::uint64_t stamp,
                                std::chrono::steady_clock::time_point received_at);
    std::mutex arrival_mutex_;  ///< 保护主雷达帧到达时刻索引。
    std::map<std::uint64_t, std::chrono::steady_clock::time_point> primary_arrivals_;  ///< 按主雷达帧时间戳索引单调时钟接收时刻。
    profiling::MultiStageTimingWindow lidar_end_to_end_window_{{"primary_to_pgo"}, 10.0};  ///< 主雷达接收到地图匹配结果进入 PGO 的端到端耗时统计。
    std::size_t lidar_deadline_misses_ = 0;  ///< 可关联到达时刻的有效跟踪帧端到端耗时达到 100 ms 的累计次数，不计本次重定位成功帧。
    int lidar_odom_skip_cnt_ = 0;  ///< 激光里程计输入跳帧计数。
    double online_sensor_max_lag_sec_ = 0.0;  ///< 在线积压诊断的进入阈值（s），0 禁用；当前不据此拒收点云。
    double online_sensor_resume_lag_sec_ = 0.0;  ///< 在线积压诊断的恢复阈值（s），启用时须小于进入阈值。
    std::atomic_bool lidar_overload_throttled_{false};  ///< 过载诊断状态；当前实现仍接收点云，由队列容量淘汰旧输入。

    // 静止检测窗口与轮速历史。
    /// @brief 静止检测窗口中的一条 IMU 样本，保留角速度方向以识别持续转动。
    struct StaticImuSample {
        double timestamp = 0.0;  ///< 该观测的传感器时间戳（s）。
        Vec3d angular_velocity = Vec3d::Zero();  ///< IMU 传感器坐标系的角速度向量（rad/s），用于静止窗口均值计算。
        double gyro_norm = 0.0;  ///< IMU 角速度模长（rad/s）。
        double accel_norm = 0.0;  ///< IMU 加速度模长（m/s²）。
    };
    mutable std::mutex static_detector_mutex_;  ///< 保护静止检测窗口、统计量与轮速历史。
    std::deque<StaticImuSample> static_imu_window_;  ///< 静止判定使用的最近 IMU 样本窗口。
    double static_gyro_sum_ = 0.0;  ///< IMU 窗口内角速度模长之和（rad/s）。
    double static_gyro_sq_sum_ = 0.0;  ///< IMU 窗口内角速度模长的平方和（rad²/s²），用于计算方差。
    Vec3d static_gyro_vector_sum_ = Vec3d::Zero();  ///< IMU 窗口内角速度向量之和（rad/s），用于计算向量均值。
    double static_accel_sum_ = 0.0;  ///< IMU 窗口内加速度模长之和（m/s²）。
    double static_accel_sq_sum_ = 0.0;  ///< IMU 窗口内加速度模长的平方和（m²/s⁴），用于计算方差。
    double last_static_lio_stamp_ = 0.0;  ///< 静止检测最近观察到的 LIO 时间戳（s）。
    double last_static_lio_speed_ = 0.0;  ///< 静止检测最近观察到的 LIO 速度模长（m/s）。
    bool last_static_lio_reliable_ = false;  ///< 最近 LIO 速度是否可信，可否参与静止判定。
    double last_valid_lidar_loc_stamp_ = 0.0;  ///< 静止检测最近记录的有效地图匹配时间戳（s），不是高频 PGO 输出时间。
    /// @brief 静止检测使用的有符号车体纵向轮速样本，用于按 IMU 时刻选择因果观测。
    struct StaticWheelSpeedSample {
        double timestamp = 0.0;  ///< 该观测的传感器时间戳（s）。
        double speed_mps = 0.0;  ///< 静止判定使用的轮速（m/s）。
    };
    std::deque<StaticWheelSpeedSample> static_wheel_speed_history_;  ///< 静止判定轮速历史，最多保留 2 s 且不超过 512 条；选取不晚于 IMU 的观测。
    double last_wheel_speed_stamp_ = 0.0;  ///< 静止检测历史最近接受的轮速时间戳（s），用于拒绝重复及回退样本。
    bool wheel_speed_observed_ = false;  ///< 是否已有可供静止判定的轮速观测。
    bool wheel_timestamp_mismatch_reported_ = false;  ///< 已报告缺少新鲜因果轮速的告警；可用轮速恢复后清除，不仅表示时间域不一致。
    int static_exit_count_ = 0;  ///< 静止保持期间连续满足退出条件的 IMU 样本数。
    bool imu_static_hold_enabled_ = false;  ///< 是否启用 IMU 高频输出的静止保持。
    bool imu_static_hold_active_ = false;  ///< 当前是否处于 IMU 高频输出静止保持状态。

    // 实时结果缓存与上层回调。
    LocalizationResult loc_result_;  ///< 最近通过实时输出时间戳门控的结果，更新与派发由 live_output_dispatch_mutex_ 保护。

    TFCallback tf_callback_;  ///< 实时结果通过时间戳门控且 valid_ 为 true 时调用，由上层广播 TF。
    LocalizationResultCallback localization_result_callback_;  ///< 时间戳门控后的实时结果回调；是否发布及有效性检查由系统层决定。
    LocalizationResultCallback global_localization_result_callback_;  ///< PGO 全局校正结果的回调，在 PGO 产出该结果的线程直接执行。
    ProcessedCloudCallback processed_cloud_callback_;  ///< 地图匹配路径交付本帧发布点云、缓存结果与组帧状态的回调，见 LidarLocProcCloud()。
    LocStateCallback loc_state_callback_;  ///< 地图匹配路径每次走到通知阶段时调用，状态未变化也可能触发。
    PointcloudBodyCallback pointcloud_body_callback_;  ///< 保留的机体坐标系点云回调字段；当前实现未注册或调用。
    PointcloudWorldCallback pointcloud_world_callback_;  ///< 保留的世界坐标系点云回调字段；当前实现未注册或调用。

    // 输入时间检查及运行诊断。
    double last_imu_time_ = 0;  ///< 最近通过数值和顺序检查的 IMU 传感器时间戳（s）；由 processing_mutex_ 保护，Init() 重置。
    double last_odom_time_ = 0;  ///< 保留的里程计时间戳字段；对应输入路径已注释，当前运行流程不更新。
    double last_cloud_time_ = 0;  ///< 保留的点云时间戳字段；当前运行流程未使用，不代表最近处理时间。

    mutable std::mutex runtime_stats_mutex_;  ///< 保护 runtime_stats_ 的复制和更新，不覆盖队列及 LIO 内部的独立统计。
    RuntimeStats runtime_stats_;  ///< 本类维护的进度、重定位和负载统计；GetRuntimeStats() 复制后再补齐队列等实时计数。
};
}  // namespace loc

}  // namespace lightning
