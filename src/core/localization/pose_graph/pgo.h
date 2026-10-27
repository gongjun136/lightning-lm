#pragma once

#include "common/eigen_types.h"
#include "common/nav_state.h"
#include "core/localization/localization_result.h"

#include "pgo_impl.h"
#include "pose_extrapolator.h"
#include "smoother.h"

namespace lightning::loc {

/**
 * Pose Graph Optimization
 * 基于 Pose Graph 的多源信息融合。
 *
 * 输入：
 * -- DR，需要降低对此的依赖（轮速失效case较多），兼容没有DR
 * -- LiDAR Odom，频率可能不高，比如5Hz；lidarOdom应该提供自己的hessian矩阵的PCA分析，提供退化情况参考
 * -- LiDAR Loc，频率可能也不高，可动态配置
 *
 * 输出：
 * -- 融合定位结果，与IMU同频率发布
 * -- 判断各个信息源的outlier，并剔除之
 * -- 盲走一段时间误差控制在容许范围内
 *
 * 调用：
 * -- 外部可以调用所有接口塞数据和获取结果，是多线程安全的；
 * -- 但内部只有一个线程，不保证一定不阻塞，外部应当handle这种情况。
 */
class PGO {
   public:
    PGO();
    ~PGO();

    /// 向外输出全局定位结果
    using GlobalOutputHandleFunction = std::function<void(const LocalizationResult& output_result)>;
    void SetGlobalOutputHandleFunction(GlobalOutputHandleFunction handle);
    void SetHighFrequencyGlobalOutputHandleFunction(GlobalOutputHandleFunction handle);
    void SetDrSmoothingEnabled(bool enabled) { dr_smoothing_enabled_ = enabled; }
    void SetDrExtrapolationEnabled(bool enabled) { dr_extrapolation_enabled_ = enabled; }
    void SetCovarianceProcessNoise(double position, double angle) {
        covariance_position_noise_ = position;
        covariance_angle_noise_ = angle;
    }

    /**
     * @brief 接收 DR 位姿并推进运动输出或静止保持。
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
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>接收有序 LO 并维护队列</b>：持锁拒绝时间回退并保存有界 LO 历史；相同时间戳不在此处拒绝。
     * 2. <b>更新 LO 与 DR 冲突状态</b>：不可靠 LO 立即置冲突；可靠输入通过计数逐步解除。
     * 3. <b>按时间关系尝试输出</b>：已有 DR 且 LO 时刻不早于最新 DR、车辆未静止时调用 PubResult()；不表示启用了 LO 位姿外推。
     */
    bool ProcessLidarOdom(const NavState& lio_result);

    /**
     * @brief 接收有效地图匹配帧，通过时间门控后触发 PGO 优化。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>检查输入来源与时间水位</b>：持锁更新静止状态；要求相对位姿和有效地图匹配，拒绝保持释放前、重置前或时间回退的结果。
     * 2. <b>构造地图匹配帧</b>：复制地图位姿、时间和分数，并生成当前实现使用的六维匹配权重。
     * 3. <b>更新里程计可信度并优化</b>：更新平滑性和 LO 可信计数，再同步交给 ProcessPGOFrame()。
     */
    bool ProcessLidarLoc(const LocalizationResult& loc_result);

    /// A stationary hold may suppress jitter, but must not hide a map conflict.
    /// On failure invalidate held output until the caller resets/relocalizes.
    bool CheckStationaryMapConsistency(const LocalizationResult& loc_result);

    /**
     * @brief 提交 PGO 帧，更新静止协方差并尝试输出。
     * @details 通常由已持有 data_mutex_ 的 ProcessLidarLoc() 调用；外部直接调用仅用于单测。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>执行优化并检查成功</b>：调用 AddPGOFrame()；失败直接返回，成功路径中的全局回调此时已由实现层处理。
     * 2. <b>更新静止保持协方差</b>：保持输出位姿不变，将其与新融合均值的偏差、旋转轴变换及时间传播计入有效协方差。
     * 3. <b>维护匹配历史并尝试输出</b>：维护 50 帧历史；地图帧时刻不早于最新 DR 时调用 PubResult()，静止高频时钟仍由 DR 驱动。
     */
    bool ProcessPGOFrame(std::shared_ptr<PGOFrame> frame);

    std::shared_ptr<PGOFrame> GetCurrentPGOFrame() const;

    // 高频查询定位，使用DR做递推，多线程安全
    bool GetCurrentLocalization() const;

    /// 重置PGO
    bool Reset();

    /**
     * @brief 在运动状态下构造高频结果并同步调用输出回调。
     * @details 由持有 data_mutex_ 的输入处理路径调用；静止输出由 ProcessDR() 推进。
     * @par 主要流程
     * 以下阶段与函数体中的“阶段 N”使用相同名称，按实际执行顺序阅读；早退分支不会经过后续阶段。
     * 1. <b>检查输出条件与匹配状态</b>：静止时直接返回；要求回调存在且融合结果有效，再更新匹配分数异常计数。
     * 2. <b>外推候选并记录对照</b>：复制融合结果、取得 DR 速度，按开关外推并记录对照；下方差异拒绝逻辑目前被注释。
     * 3. <b>平滑位姿并传播协方差</b>：保存外推状态与未平滑位姿；启用平滑时将修正量和旋转轴变换计入有效协方差。
     * 4. <b>调用输出回调并维护历史</b>：同步交付结果并维护有界输出历史；上层仍会执行时间门控。
     */
    void PubResult();

    /// debug stuffs
    void SetDebug(bool debug = true);
    void LogWindowState();

   public:
    bool localization_unusual_tag_ = false;  ///< 定位失效标志位
    bool imu_interruption_tag_ = false;      ///< imu断流标志位

   private:
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

    // 在本层放一个extrapolator，被三个地方访问：DR，lidarOdom，LidarLoc，RTK
    std::unique_ptr<PoseExtrapolator> pose_extrapolator_;  ///< 根据近期 PGO 校正与里程计生成高频位姿。
    std::function<void(const LocalizationResult& output_result)> high_freq_output_func_;  ///< 高频融合结果的输出回调。
    LocalizationResult high_freq_result_;  ///< 最近生成的高频融合结果。
    LocalizationResult parking_result_;  ///< 静止保持期间复用的定位结果。

    std::shared_ptr<PoseSmoother> smoother_;  ///< 用于抑制输出跳变的位姿平滑器。

    float imu_interruption_time_thd_ = 1.0;  ///< imu数据断流时间阈值 s
    float lidar_loc_score_thd_ = 0.5;        ///< 激光定位分值阈值
    int localization_unusual_thd_ = 10;      ///< 激光定位失效次数阈值
    int localization_unusual_count_ = 0;     ///< 定位异常次数
    double last_lidar_loc_time_ = 0.;        ///< 上次激光定位时间戳
    double last_lidar_loc_input_time_ = -1.;  ///< 最近收到的地图匹配输入时间戳（s）。
    // LidarLoc is asynchronous. Results sampled no later than the most recent
    // DR hold release still belong to the parked epoch even if they arrive
    // after is_parking_ becomes false.
    double static_hold_release_watermark_ = -1.;  ///< 静止保持释放后拒绝旧结果的时间戳水位（s）。
    double static_hold_start_stamp_ = -1.;  ///< 本轮静止保持开始的传感器时间戳（s）。
    bool is_parking_ = false;  ///< PGO 综合判定的静止状态。
    bool dr_is_parking_ = false;  ///< DR 输入报告的静止状态。
    bool lidar_loc_is_parking_ = false;  ///< 地图匹配输入报告的静止状态。
    bool dr_smoothing_enabled_ = true;  ///< 是否使用 DR 约束位姿平滑。
    bool dr_extrapolation_enabled_ = true;  ///< 是否使用 DR 外推高频位姿。
    double covariance_position_noise_ = 0.20;  // m / sqrt(s), bench calibration required
    double covariance_angle_noise_ = 2.0 * M_PI / 180.0;  // rad / sqrt(s)
};

}  // namespace lightning::loc
