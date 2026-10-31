/**
 * @file pgo.h
 * @brief 地图定位与相对运动融合接口。
 */
#pragma once

#include "common/eigen_types.h"
#include "common/nav_state.h"
#include "core/localization/localization_result.h"

#include "pgo_impl.h"
#include "pose_extrapolator.h"
#include "smoother.h"

namespace lightning::loc {

/**
 * @brief 将地图匹配与里程计信息融合，并输出地图系定位结果。
 * @details 地图匹配帧经 ProcessLidarLoc()、ProcessPGOFrame() 进入窗口优化；
 * ProcessDR() 缓存相对运动并驱动外推或静止保持，不自行求解。
 * 当前优化问题包含地图匹配、激光里程计和先验因子，不加入 DR 因子。
 * @par 阅读顺序
 * 从输入接口和 ProcessPGOFrame() 开始，再阅读 PubResult() 与 ExtrapolateLocResult()。
 * @note 本类没有后台线程。输入处理在调用线程持数据锁同步执行，回调也可能在锁内同步执行；
 * 设置接口及直接调用 ProcessPGOFrame()、PubResult() 等不自行加锁，调用方须协调并发和重入。
 * @see PGOImpl::AddPGOFrame(), PGOImpl::RunOptimization(), @ref pose_graph_theory
 */
class PGO {
   public:
    PGO();
    ~PGO();

    /// 同步接收融合结果的回调类型。
    using GlobalOutputHandleFunction = std::function<void(const LocalizationResult& output_result)>;
    /// 启动输入前设置优化结果回调；优化帧时刻可能早于当前高频输出。
    /// 回调同步执行，异步使用参数须复制；不得重入获取同一数据锁的输入接口或 Reset()。
    void SetGlobalOutputHandleFunction(GlobalOutputHandleFunction handle);
    /// 启动输入前设置高频回调；输入处理可能持数据锁同步调用，上层仍会门控输出。
    /// 异步使用参数须复制；回调不得重入获取同一数据锁的输入接口或 Reset()。
    void SetHighFrequencyGlobalOutputHandleFunction(GlobalOutputHandleFunction handle);
    /// 控制输出位姿平滑。
    void SetDrSmoothingEnabled(bool enabled) { dr_smoothing_enabled_ = enabled; }
    /// 控制高频结果的 DR 外推。
    void SetDrExtrapolationEnabled(bool enabled) { dr_extrapolation_enabled_ = enabled; }
    /**
     * @brief 设置协方差传播噪声；仅赋值，不校验。
     * @param position 平移噪声（m/sqrt(s)）。
     * @param angle 旋转噪声（rad/sqrt(s)）。
     */
    void SetCovarianceProcessNoise(double position, double angle) {
        covariance_position_noise_ = position;
        covariance_angle_noise_ = angle;
    }

    /**
     * @brief 接收 DR 位姿并推进运动输出或静止保持。
     * @details 本入口缓存、外推或保持输出，不建图或求解；仅检查时间顺序，数值和位姿有效性由上游保证。
     * @param dr_result DR 状态，时间戳为传感器时间（s）。
     * @return 时间回退返回 false；接收成功不保证产生输出。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>接收有序 DR 并维护队列</b>：持锁拒绝时间回退并保存有界 DR 历史；相同时间戳仍可进入此函数。
     * 2. <b>更新静止保持切换</b>：仅有效融合结果可进入保持；退出保持时提升最后输出的时刻并建立旧地图结果拒绝水位。
     * 3. <b>按运动状态交付高频结果</b>：运动时调用 PubResult()；静止时保持位姿、传播有效协方差并按当前 DR 时刻交付。
     */
    bool ProcessDR(const NavState& dr_result);

    /**
     * @brief 接收 LIO 位姿，更新可靠性并按输入时刻尝试高频输出。
     * @details 本入口缓存并可能触发输出，不建图或求解；仅检查时间顺序，数值和位姿有效性由上游保证。
     * @param lio_result LO 状态，时间戳为传感器时间（s）。
     * @return 时间回退返回 false；接收成功不保证产生输出。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>接收有序 LO 并维护队列</b>：持锁拒绝时间回退并保存有界 LO 历史；相同时间戳不在此处拒绝。
     * 2. <b>更新 LO 与 DR 冲突状态</b>：不可靠 LO 立即置冲突；可靠输入通过计数逐步解除。
     * 3. <b>按时间关系尝试输出</b>：已有 DR 且 LO 时刻不早于最新 DR、车辆未静止时调用 PubResult()；不表示启用了 LO 位姿外推。
     */
    bool ProcessLidarOdom(const NavState& lio_result);

    /**
     * @brief 接收有效地图匹配帧，通过时间门控后触发 PGO 优化。
     * @param loc_result 地图定位结果，时间戳为传感器时间（s）。
     * @return 输入门控或优化失败返回 false；成功不保证产生高频输出。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>检查输入来源与时间水位</b>：持锁更新静止状态；要求相对位姿和有效地图匹配，拒绝保持释放前、重置前或时间回退的结果。
     * 2. <b>构造地图匹配帧</b>：复制地图位姿、时间和分数，并生成当前实现使用的六维匹配权重。
     * 3. <b>更新里程计可信度并优化</b>：更新平滑性和 LO 可信计数，再同步交给 ProcessPGOFrame()。
     */
    bool ProcessLidarLoc(const LocalizationResult& loc_result);

    /**
     * @brief 适用时比较静止保持位姿与新地图匹配，防止保持掩盖冲突。
     * @details 仅静止中、保持结果和地图匹配有效且输入不早于保持开始时检查。
     * @param loc_result 待检查的地图定位结果。
     * @return 不适用或平移不超过 0.5 m 且转角不超过 5 度时返回 true；
     * 冲突或数值非有限时将保持和融合结果置无效并返回 false。重置或重定位由上层处理。
     */
    bool CheckStationaryMapConsistency(const LocalizationResult& loc_result);

    /**
     * @brief 提交 PGO 帧，更新静止协方差并尝试输出。
     * @details 通常由已持有 data_mutex_ 的 ProcessLidarLoc() 调用；本函数不自行加锁。
     * @param frame 待优化帧，须为有效对象。
     * @return PGOImpl::AddPGOFrame() 失败返回 false，否则返回 true。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>执行优化并检查成功</b>：调用 PGOImpl::AddPGOFrame()；失败直接返回，成功路径中的全局回调此时已由实现层处理。
     * 2. <b>更新静止保持协方差</b>：保持输出位姿不变，将其与新融合均值的偏差、旋转轴变换及时间传播计入有效协方差。
     * 3. <b>维护匹配历史并尝试输出</b>：维护 50 帧历史；地图帧时刻不早于最新 DR 时调用 PubResult()，静止高频时钟仍由 DR 驱动。
     */
    bool ProcessPGOFrame(std::shared_ptr<PGOFrame> frame);

    /// 返回当前帧的共享可变对象，可能为空且不是快照；本函数无锁，调用方须协调访问。
    std::shared_ptr<PGOFrame> GetCurrentPGOFrame() const;

    /// 保留声明：当前仓库没有定义或调用，不能作为可用查询接口。
    bool GetCurrentLocalization() const;

    /**
     * @brief 重置窗口、平滑器及部分输入和保持状态；不清除 imu_interruption_tag_、回调或配置。
     * @return PGOImpl::Reset() 的返回值。
     */
    bool Reset();

    /**
     * @brief 在运动状态下构造高频结果并同步调用输出回调。
     * @details 由持有 data_mutex_ 的输入处理路径调用；本函数不自行加锁，静止输出由 ProcessDR() 推进。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>检查输出条件与匹配状态</b>：静止时直接返回；要求回调存在且融合结果有效，再更新匹配分数异常计数。
     * 2. <b>外推候选并记录对照</b>：复制融合结果、取得 DR 速度，按开关外推并记录对照；下方差异拒绝逻辑目前被注释。
     * 3. <b>平滑位姿并传播协方差</b>：保存外推状态与未平滑位姿；启用平滑时将修正量和旋转轴变换计入有效协方差。
     * 4. <b>调用输出回调并维护历史</b>：同步交付结果并维护有界输出历史；上层仍会执行时间门控。
     */
    void PubResult();

    /// 设置实现层调试开关；本函数不自行加锁。
    void SetDebug(bool debug = true);
    /// 持数据锁输出预留诊断表头；当前不打印实际帧数据。
    void LogWindowState();

   public:
    bool localization_unusual_tag_ = false;  ///< 输出路径中连续低分超过门限的提示，非最终定位失效判定。
    bool imu_interruption_tag_ = false;  ///< 外推时报告 DR 相对最新输入滞后，不能单独证明 IMU 健康。

   private:
    /// 检查 DR 和 LO 队列是否均为空；调用方须持数据锁。
    inline bool RelativePoseQueueEmpty() {
        return (impl_->dr_pose_queue_.empty() && impl_->lidar_odom_pose_queue_.empty());
    }

    /**
     * @brief 用 DR 增量外推融合结果，必要时由地图匹配历史补到最新输入时刻。
     * @details LO 当前只参与确定最新时刻；LO 位姿外推与 DR 增量拒绝条件均被注释。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>检查结果并确定最新时刻</b>：要求输入结果有效；由融合结果、最新 DR 与 LO 的时间戳确定目标时刻。
     * 2. <b>报告 DR 断流状态</b>：比较目标时刻与最新 DR 时刻并更新诊断；该标志本身不阻止后面的外推。
     * 3. <b>应用 DR 相对增量</b>：DR 比结果更新且旧结果时刻可插值时，右乘 DR 增量，并传播有效协方差。
     * 4. <b>使用地图历史补偿剩余延迟</b>：仍落后目标超过 0.05 s 时尝试地图位姿历史插值或外推；采用后将协方差标记无效。
     * @param[in,out] output_result 待外推的融合定位结果。
     * @return 输入无效返回 false；返回 true 不保证时间戳已推进。
     */
    bool ExtrapolateLocResult(LocalizationResult& output_result);

    std::unique_ptr<PGOImpl> impl_;  ///< 维护滑动窗口及优化约束的 PGO 实现。

    std::unique_ptr<PoseExtrapolator> pose_extrapolator_;  ///< 保留对象；当前仅构造，实际外推由 ExtrapolateLocResult() 完成。
    std::function<void(const LocalizationResult& output_result)> high_freq_output_func_;  ///< 高频融合结果的输出回调。
    LocalizationResult high_freq_result_;  ///< 保留状态；当前仅在 Reset() 清空，不保存实际高频输出。
    LocalizationResult parking_result_;  ///< 静止保持期间复用的定位结果。

    std::shared_ptr<PoseSmoother> smoother_;  ///< 用于抑制输出跳变的位姿平滑器。

    float imu_interruption_time_thd_ = 1.0;  ///< 最新 DR 相对最新输入的滞后阈值（s）。
    float lidar_loc_score_thd_ = 0.5;        ///< 激光定位分值阈值
    int localization_unusual_thd_ = 10;      ///< 连续低分计数阈值，超过时置 localization_unusual_tag_。
    int localization_unusual_count_ = 0;     ///< PubResult() 按融合结果时间推进统计的连续低分次数。
    double last_lidar_loc_time_ = 0.;        ///< 上次低分检查对应的融合结果时间戳（s）。
    double last_lidar_loc_input_time_ = -1.;  ///< 最近收到的地图匹配输入时间戳（s）。
    // 异步地图结果即使在退出静止后到达，释放时刻及以前的样本仍属于保持期间。
    double static_hold_release_watermark_ = -1.;  ///< 静止保持释放后拒绝旧结果的时间戳水位（s）。
    double static_hold_start_stamp_ = -1.;  ///< 本轮静止保持开始的传感器时间戳（s）。
    bool is_parking_ = false;  ///< PGO 综合判定的静止状态。
    bool dr_is_parking_ = false;  ///< DR 输入报告的静止状态。
    bool lidar_loc_is_parking_ = false;  ///< 地图匹配输入报告的静止状态。
    bool dr_smoothing_enabled_ = true;  ///< 是否使用 DR 约束位姿平滑。
    bool dr_extrapolation_enabled_ = true;  ///< 是否使用 DR 外推高频位姿。
    double covariance_position_noise_ = 0.20;  ///< 协方差传播平移噪声（m/sqrt(s)），需要台架标定。
    double covariance_angle_noise_ = 2.0 * M_PI / 180.0;  ///< 协方差传播旋转噪声（rad/sqrt(s)）。
};

}  // namespace lightning::loc
