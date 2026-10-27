#ifndef LIGHTNING_BACKEND_BACKEND_PIPELINE_H
#define LIGHTNING_BACKEND_BACKEND_PIPELINE_H

#include <atomic>
#include <condition_variable>
#include <deque>
#include <functional>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "core/backend/btc_loop_detector.h"
#include "core/backend/hierarchical_bundle_adjustment.h"

namespace lightning::backend {

enum class BackendMode { kLegacy, kBaBtcHba, kDisabled };

BackendMode ReadBackendMode(const std::string& yaml_path);
const char* BackendModeName(BackendMode mode);

struct BackendPipelineOptions {
    bool online_mode = false;  ///< 是否通过后台队列异步处理关键帧。
    bool verbose = true;  ///< 是否输出详细后端日志。
    int pose_graph_iterations = 20;  ///< 位姿图优化的迭代次数。
    int pose_graph_outlier_iterations = 10;  ///< 剔除异常约束后再次优化的迭代次数。
    double motion_translation_noise = 0.10;  ///< 里程计相邻位姿约束的平移噪声（m）。
    double motion_rotation_noise_deg = 3.0;  ///< 里程计相邻位姿约束的旋转噪声（deg）。
    double loop_translation_noise = 0.20;  ///< 回环约束的平移噪声（m）。
    double loop_rotation_noise_deg = 3.0;  ///< 回环约束的旋转噪声（deg）。
    double robust_kernel_delta = 1.0;  ///< 位姿图鲁棒核的尺度参数。
    double loop_outlier_chi2 = 25.0;  ///< 回环约束的卡方异常判定阈值。
    BundleAdjustmentOptions local_ba;  ///< 局部体素 BA 配置。
    BtcLoopDetectorOptions btc;  ///< BTC 回环检测配置。
    HierarchicalBundleAdjustmentOptions hba;  ///< 分层体素 BA 配置。
};

struct BackendRuntimeSummary {
    std::size_t keyframes = 0;  ///< 已进入后端的关键帧数。
    std::size_t local_ba_attempts = 0;  ///< 局部 BA 尝试次数。
    std::size_t local_ba_accepted = 0;  ///< 局部 BA 接受结果的次数。
    std::size_t btc_descriptors = 0;  ///< 已生成的 BTC 描述子数量。
    std::size_t btc_candidates = 0;  ///< BTC 检索得到的回环候选数。
    std::size_t loops_accepted = 0;  ///< 通过几何与连续确认的回环数。
    std::size_t loops_applied = 0;  ///< 累计写入图约束容器的回环数，不代表经过 PGO 检验仍为内点。
    std::size_t loops_graph_inliers = 0;  ///< 位姿图保留的回环内点数。
    std::size_t hba_runs = 0;  ///< HBA 执行次数。
    std::size_t hba_accepted = 0;  ///< HBA 接受结果的次数。
    double local_ba_time_ms = 0.0;  ///< 局部 BA 耗时统计（ms）。
    double btc_time_ms = 0.0;  ///< BTC 检索与验证耗时统计（ms）。
    double pose_graph_time_ms = 0.0;  ///< 位姿图优化耗时统计（ms）。
    double hba_time_ms = 0.0;  ///< 分层优化耗时统计（ms）。
};

class BackendPipeline {
   public:
    BackendPipeline() = default;
    ~BackendPipeline();

    /**
     * @brief 初始化建图后端。
     *
     * @par 主要流程
     * 1. <b>读取配置与外参</b>：已初始化则直接返回，否则读取后端选项和雷达到 IMU 外参。
     * 2. <b>创建优化模块</b>：创建局部 BA、BTC、HBA 对象并清除停止标志。
     * 3. <b>按模式启动线程</b>：在线模式启动关键帧线程；启用 HBA 时启动其线程并标记初始化完成。
     * @note 离线模式也可能启动 HBA 线程；本函数未清空已有关键帧和运行统计，不能据此推断重复初始化会重置后端。
     */
    bool Init(const std::string& yaml_path, bool online_mode);
    /// @brief 在线入队，离线同步处理；未初始化或空关键帧直接返回。
    /// @see @ref backend_queue_ba_reading
    void AddKeyframe(const Keyframe::Ptr& keyframe);
    void RequestGlobalOptimization(const std::string& reason = "manual");
    /**
     * @brief 等待当前后端工作结束，可选触发全局优化。
     *
     * @par 主要流程
     * 1. <b>等待关键帧处理空闲</b>：未初始化则返回；在线模式等待队列为空且无正在处理的帧。
     * 2. <b>按需请求最终优化</b>：force_global_optimization 为 true 时请求一次全局优化。
     * 3. <b>等待 HBA 请求结束</b>：启用 HBA 时等待无待处理请求且未运行。
     * @note 不会禁止新的关键帧或优化请求进入；调用方应先停止输入，才能把返回视为数据结束后的排空。
     */
    void WaitUntilIdle(bool force_global_optimization);
    /**
     * @brief 停止并回收后端工作线程。
     *
     * @par 主要流程
     * 1. <b>等待现有任务结束</b>：未初始化则返回，否则等待空闲且不强制追加全局优化。
     * 2. <b>停止关键帧线程</b>：设置停止标志、唤醒等待并 join 在线工作线程。
     * 3. <b>停止 HBA 线程并清除状态</b>：唤醒并 join HBA 线程，最后清除初始化标记。
     * @note 调用方应先停止生产新任务；本函数不负责关闭上游传感器或阻止并发入队。
     */
    void Shutdown();

    using OptimizedCallback = std::function<void()>;
    /// @brief 设置关键帧处理或优化之后同步执行的通知函数。
    /// @note 调用时已释放 data_mutex_，但调用者可能仍持 optimization_mutex_；回调不应等待后端空闲或重入同步优化。
    /// @see @ref backend_map_callback_reading
    void SetOptimizedCallback(OptimizedCallback callback);

    SE3 GetMapToOdom() const;
    std::vector<BtcLoopResult> GetLoopResults() const;
    BackendRuntimeSummary GetRuntimeSummary() const;
    bool SaveDiagnostics(const std::string& directory,
                         const map_frame::Metadata* map_metadata = nullptr) const;
    bool SaveRelocalizationDatabase(
        const std::string& directory,
        const map_frame::Metadata* map_metadata = nullptr) const;
    const BackendPipelineOptions& GetOptions() const { return options_; }

   private:
    struct LoopConstraint {
        BtcLoopResult detection;  ///< 该图约束对应的回环检测结果。
        SE3 T_history_imu_current_imu;  ///< 将当前 IMU 坐标变换到历史 IMU 坐标的回环相对位姿。
        bool graph_inlier = true;  ///< 该回环约束是否被位姿图判为内点。
    };

    /**
     * @brief 按局部优化、回环检测和全局优化顺序处理关键帧。
     *
     * @par 主要流程
     * 1. <b>登记关键帧并选择窗口</b>：跳过与末帧相同的指针，按窗口大小和步长选取局部 BA 输入。
     * 2. <b>按条件执行局部 BA</b>：启用且窗口非空时串行优化，随后记录结果和耗时。
     * 3. <b>检测回环并登记约束</b>：BTC 给出检测结果后记录统计，仅 accepted 且 optimization_warranted 时写入图约束。
     * 4. <b>按入图结果触发全局优化</b>：实际写入约束后运行 PGO 并请求 HBA。
     * 5. <b>更新地图到里程计变换</b>：本帧处理结束时刷新变换并同步通知订阅回调。
     * @note loops_applied 表示约束曾入图，不等于 PGO 内点数；回调可能在持 optimization_mutex_ 时同步执行。
     */
    void HandleKeyframe(const Keyframe::Ptr& keyframe);
    /**
     * @brief 持续消费在线关键帧队列。
     *
     * @par 主要流程
     * 1. <b>等待并取出关键帧</b>：每轮持队列锁等待，取出一帧并标记处理中；停止且队列为空时退出。
     * 2. <b>释放队列锁后处理</b>：调用 HandleKeyframe 执行耗时工作。
     * 3. <b>结束本帧并通知空闲</b>：清除处理中标志，队列为空时唤醒等待者。
     * @note 停止请求不直接丢弃已入队关键帧；队列锁在算法处理前释放。
     */
    void WorkerLoop();
    /**
     * @brief 消费 HBA 请求并控制优化结果的提交。
     *
     * @par 主要流程
     * 1. <b>等待并消费优化请求</b>：每轮取出合并后的请求和前置 PGO 标志，设置 HBA 运行状态。
     * 2. <b>快照关键帧与入图历史</b>：读取当前关键帧和是否曾有回环约束入图。
     * 3. <b>串行优化并执行提交门控</b>：可选 PGO 后保存位姿快照，执行 HBA，按入图历史门控回退 HBA 修改并刷新变换。
     * 4. <b>记录结果并通知空闲</b>：更新统计和日志，清除运行标志；无后续请求时通知等待者。
     * @note 提交门控检查累计入图历史，不检查本次 PGO 的内点数；回退快照取自前置 PGO 之后，仅撤销本次 HBA 修改。
     */
    void HbaLoop();
    void RequestHba(bool rerun_pose_graph, const std::string& reason);
    /**
     * @brief 优化回环位姿图并写回关键帧。
     *
     * @par 主要流程
     * 1. <b>快照输入并检查规模</b>：复制关键帧和回环约束；少于两帧或无约束时返回 false。
     * 2. <b>建立顶点与里程计边</b>：用优化位姿初始化顶点、固定首帧，并以 LIO 相对运动建立相邻边。
     * 3. <b>加入回环边并首次优化</b>：跳过端点缺失的约束，给回环边加鲁棒核后运行优化。
     * 4. <b>筛选内点并写回位姿</b>：按卡方阈值筛选回环；存在内点时再次优化，随后写回顶点估计。
     * 5. <b>更新统计并返回内点结果</b>：回写约束的内点标记，累计耗时，按是否存在内点返回。
     * @note 零内点时也会写回首次优化的估计；返回 false 不等于关键帧位姿未改变。调用方负责串行化优化。
     */
    bool OptimizePoseGraph();
    void UpdateMapToOdomAndNotify();

    BackendPipelineOptions options_;  ///< 后端运行、图约束与各优化模块的配置。
    SE3 T_imu_lidar_;  ///< 将雷达坐标变换到 IMU 坐标的固定外参。
    std::unique_ptr<VoxelBundleAdjuster> local_ba_;  ///< 局部滑动窗口体素 BA 优化器。
    std::unique_ptr<BtcLoopDetector> btc_;  ///< 全局回环候选检索与几何验证器。
    std::unique_ptr<HierarchicalBundleAdjuster> hba_;  ///< 多层子地图的全局 BA 优化器。

    mutable std::mutex data_mutex_;  ///< 保护关键帧、回环约束、输出变换与统计。
    std::vector<Keyframe::Ptr> keyframes_;  ///< 后端持有的全部关键帧。
    std::vector<BtcLoopResult> loop_results_;  ///< 已记录的回环检测与验证结果。
    std::vector<LoopConstraint> loop_constraints_;  ///< 交给位姿图优化的回环约束。
    BackendRuntimeSummary runtime_;  ///< 各优化阶段的运行计数与耗时。
    SE3 T_map_odom_;  ///< 里程计坐标系到优化后地图坐标系的变换。
    OptimizedCallback optimized_callback_;  ///< 完成后端优化后的通知回调。

    std::mutex optimization_mutex_;  ///< 防止局部、图优化与分层优化并发修改位姿。

    std::mutex queue_mutex_;  ///< 保护关键帧工作队列与工作线程状态。
    std::condition_variable queue_cv_;  ///< 新关键帧或停止请求到达时唤醒工作线程。
    std::condition_variable queue_idle_cv_;  ///< 通知等待者关键帧处理队列已空闲。
    std::deque<Keyframe::Ptr> queue_;  ///< 在线模式等待后端处理的关键帧队列。
    std::thread worker_thread_;  ///< 执行关键帧处理的在线工作线程。
    bool worker_processing_ = false;  ///< 工作线程是否正在处理关键帧。
    bool worker_stop_ = false;  ///< 请求关键帧工作线程退出的标志。

    std::mutex hba_mutex_;  ///< 保护 HBA 请求、执行及退出状态。
    std::condition_variable hba_cv_;  ///< 唤醒 HBA 线程执行新的优化请求。
    std::condition_variable hba_idle_cv_;  ///< 通知等待者 HBA 已结束。
    std::thread hba_thread_;  ///< 执行全局分层优化的线程。
    bool hba_requested_ = false;  ///< 是否有尚未处理的 HBA 请求。
    bool hba_rerun_pose_graph_ = false;  ///< 本次 HBA 请求是否需要先重跑位姿图。
    bool hba_running_ = false;  ///< HBA 线程是否正在执行优化。
    bool hba_stop_ = false;  ///< 请求 HBA 线程退出的标志。
    std::string hba_reason_;  ///< 当前 HBA 请求的触发原因。

    bool initialized_ = false;  ///< 后端配置和模块是否初始化完成。
    std::size_t last_local_ba_keyframe_count_ = 0;  ///< 最近触发局部 BA 时的关键帧总数。
};

}  // namespace lightning::backend

#endif  // LIGHTNING_BACKEND_BACKEND_PIPELINE_H
