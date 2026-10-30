#pragma once

/**
 * @file
 * @brief 分阶段性能计时、窗口分布统计与日志格式化接口。
 * @see lightning::profiling::Stopwatch, lightning::profiling::MultiStageTimingWindow
 */

#include <chrono>
#include <cstddef>
#include <initializer_list>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

namespace lightning::profiling {

constexpr const char* kComputeProfilingEnv = "LIGHTNING_LM_COMPUTE_PROFILE";  ///< 性能采样开关的环境变量名。
/// 控制调用方裁剪非必要诊断开销的环境变量名。
constexpr const char* kReduceNonessentialOverheadEnv =
    "LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD";

/**
 * @brief 将环境变量文本解析为布尔开关，忽略字母大小写但不去除空白。
 * @param value 待解析的字符串，允许为 nullptr。
 * @return 仅 "1"、"true"、"yes"、"on" 返回 true；空指针和其他文本均返回 false。
 */
bool ParseBooleanFlag(const char* value);

/**
 * @brief 查询性能采样开关，首次调用后缓存结果。
 * @return LIGHTNING_LM_COMPUTE_PROFILE 经 ParseBooleanFlag() 解析后的值。
 * @note 后续修改环境变量不会刷新此结果；与因果事件追踪的开关分别控制。
 * @see Stopwatch, MultiStageTimingWindow
 */
bool ComputeProfilingEnabled();

/**
 * @brief 查询是否裁剪非必要诊断开销，首次调用后缓存结果。
 * @return LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD 经 ParseBooleanFlag() 解析后的值。
 * @details 调用方据此选择是否输出高频日志等；本函数只查询策略。
 */
bool ReduceNonessentialOverhead();

/**
 * @brief 单次处理阶段的耗时样本，所有时间字段均以毫秒为单位。
 * @note valid 为 false 时窗口会跳过该阶段；默认零值不代表实际测得耗时为零。
 */
struct TimingSample {
    bool valid = false;  ///< 是否将该阶段样本纳入窗口统计；禁用计时器时为 false。
    double wall_ms = 0.0;  ///< 单调时钟经过时间（ms），包含执行、阻塞和调度等待。
    double thread_cpu_ms = 0.0;  ///< 调用线程实际消耗的 CPU 时间（ms）。
    /**
     * @brief 阶段期间整个进程消耗的 CPU 时间（ms），可因并行执行而大于 wall_ms。
     * @details 包含并发工作线程及无关线程的开销，不能视为该阶段独占的 CPU 时间。
     */
    double process_cpu_ms = 0.0;
};

/**
 * @brief 构造时开始、显式 Stop() 时结束的一次性计时器。
 * @details 同时采集单调时钟、调用线程 CPU 和进程 CPU 时间；自身不输出日志。
 * @note 构造与 Stop() 应在同一线程执行，实例不提供并发同步。
 * CPU 时钟读取失败时以零替代该次读数，TimingSample::valid 不单独反映此类失败。
 */
class Stopwatch {
   public:
    /**
     * @brief 按开关决定是否立即读取起始时钟。
     * @param enabled 是否启用本次计时；默认读取 ComputeProfilingEnabled()。
     */
    explicit Stopwatch(bool enabled = ComputeProfilingEnabled());

    /**
     * @brief 结束计时并缓存样本，重复调用返回首次停止的结果。
     * @return 启用时返回 valid 为 true 的样本；禁用时返回默认的无效样本。
     */
    TimingSample Stop();

   private:
    bool enabled_ = false;  ///< 是否启用本次计时。
    bool stopped_ = false;  ///< 本次计时是否已结束，重复 Stop() 复用结果。
    std::chrono::steady_clock::time_point wall_start_;  ///< 本次计时开始的单调时钟时刻。
    double thread_cpu_start_ms_ = 0.0;  ///< 开始时调用线程的 CPU 累计耗时（ms）。
    double process_cpu_start_ms_ = 0.0;  ///< 开始时整个进程的 CPU 累计耗时（ms）。
    TimingSample result_;  ///< 停止计时后保存的耗时样本。
};

/**
 * @brief 一组有效样本的分布摘要，数值单位与输入样本一致。
 * @details 分位数按排序后的位置 p * (count - 1) 线性插值；无样本时各字段均为零，
 * 应通过 count 区分空分布与实际测得的零值。
 */
struct DistributionSummary {
    std::size_t count = 0;  ///< 参与该分布统计的有效样本数。
    double mean = 0.0;  ///< 样本的算术平均值。
    double p50 = 0.0;  ///< 样本的第 50 百分位数。
    double p95 = 0.0;  ///< 样本的第 95 百分位数。
    double p99 = 0.0;  ///< 样本的第 99 百分位数。
    double max = 0.0;  ///< 样本中的最大值。
};

/// @brief 单个命名阶段的经过时间、线程 CPU 和进程 CPU 分布，均以毫秒为单位。
struct StageTimingSummary {
    std::string name;  ///< 该组耗时统计对应的处理阶段名称。
    DistributionSummary wall_ms;  ///< 单调时钟经过时间的分布（ms），包含运行与等待。
    DistributionSummary thread_cpu_ms;  ///< 调用线程 CPU 时间的分布（ms）。
    DistributionSummary process_cpu_ms;  ///< 含并发线程的进程 CPU 时间分布（ms），不表示阶段独占开销。
};

/// @brief 一个实际统计窗口的汇总；各阶段的有效样本数可能少于样本组数。
struct TimingWindowSummary {
    double window_s = 0.0;  ///< 当前统计窗口实际覆盖的时间（s）。
    std::size_t samples = 0;  ///< 该窗口接收的样本组数，包含各阶段均无效的组。
    std::vector<StageTimingSummary> stages;  ///< 按固定阶段顺序排列的耗时统计结果。
};

/**
 * @brief 线程安全的固定阶段模式统计窗口，由 Add() 驱动收集和输出。
 * @details 首次 Add() 建立窗口起点；达到目标时长的那次调用先计入本组样本，
 * 再返回汇总、清空样本并以该次 now 作为下一窗口起点，边界样本不重复计入。
 * @note 没有后台定时器；输入停止时未到期的窗口不会自动输出。
 * @see TimingSample, TimingWindowSummary, FormatTimingWindow()
 */
class MultiStageTimingWindow {
   public:
    /**
     * @brief 固定阶段名称、排列顺序和目标统计时长。
     * @param stage_names 非空阶段列表；后续每组样本必须按此顺序提供。
     * @param window_seconds 目标窗口长度（s），必须大于零。
     * @throws std::invalid_argument 阶段列表为空或窗口长度不大于零。
     */
    explicit MultiStageTimingWindow(std::vector<std::string> stage_names,
                                    double window_seconds = 1.0);

    /**
     * @brief 按阶段顺序加入一组样本，到期时返回包含本组的窗口汇总。
     * @param samples 与阶段数相同的样本列表；valid 为 false 的阶段不加入分布。
     * @param now 本次统计的单调时钟时刻；显式传入时应保持时间不回退。
     * @return 未到期返回 std::nullopt；到期返回汇总并重置窗口。
     * @throws std::invalid_argument 样本数量与构造时的阶段数量不符。
     */
    std::optional<TimingWindowSummary> Add(
        std::initializer_list<TimingSample> samples,
        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now());
    /**
     * @brief 接收 vector 形式的一组样本，窗口与有效性规则与列表重载相同。
     * @param samples 按固定阶段顺序排列的样本，数量须与阶段数一致。
     * @param now 本次统计的单调时钟时刻；显式传入时应保持时间不回退。
     * @return 未到期返回 std::nullopt；到期返回包含本组的汇总并重置窗口。
     * @throws std::invalid_argument 样本数量与阶段数量不符。
     */
    std::optional<TimingWindowSummary> Add(
        const std::vector<TimingSample>& samples,
        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now());

   private:
    /// @brief 当前窗口中某一阶段的有效样本，三种计时指标按同一样本顺序保存。
    struct StageSamples {
        std::vector<double> wall_ms;  ///< 有效阶段的单调时钟经过时间样本（ms）。
        std::vector<double> thread_cpu_ms;  ///< 有效阶段的调用线程 CPU 时间样本（ms）。
        std::vector<double> process_cpu_ms;  ///< 有效阶段的进程 CPU 时间样本（ms），包含并发线程。
    };

    /**
     * @brief 两个 Add() 重载共用的窗口更新实现。
     * @pre 调用者已持有 mutex_。
     * @par 主要流程
     * 1. <b>检查样本模式</b>：拒绝阶段数不一致的输入。
     * 2. <b>收集本次样本</b>：必要时建立起点，按有效性收集并累计样本组数。
     * 3. <b>判断窗口期限</b>：尚未到期时返回空结果。
     * 4. <b>生成汇总并重置</b>：统计各阶段分布，清空样本并更新起点。
     * @param samples 本组连续存放的阶段样本。
     * @param sample_count 本组样本数量，必须等于阶段数。
     * @param now 本次统计的单调时钟时刻。
     * @return 未到期的空结果，或已结束窗口的完整汇总。
     * @throws std::invalid_argument 样本数量与阶段数量不符。
     */
    std::optional<TimingWindowSummary> AddLocked(
        const TimingSample* samples, std::size_t sample_count,
        std::chrono::steady_clock::time_point now);

    std::vector<std::string> stage_names_;  ///< 固定统计模式中各阶段的名称与顺序。
    std::vector<StageSamples> stage_samples_;  ///< 当前窗口按阶段积累的耗时样本。
    double window_seconds_ = 1.0;  ///< 目标统计窗口长度（s）。
    std::chrono::steady_clock::time_point window_start_;  ///< 当前窗口开始的单调时钟时刻。
    bool window_started_ = false;  ///< 当前统计窗口是否已接收到首个样本。
    std::size_t samples_ = 0;  ///< 当前窗口累积的样本组数。
    std::mutex mutex_;  ///< 保护窗口样本与起止状态。
};

/**
 * @brief 将一次计时格式化为带前缀的 key=value 字段，数值保留六位小数。
 * @param prefix 字段名前缀，用于标识处理阶段。
 * @param sample 待输出的耗时样本，时间单位为毫秒。
 * @return 包含经过时间、线程 CPU 和并发进程 CPU 时间的文本。
 * @note 不检查 sample.valid；调用者须区分无效样本与实际测得的零值。
 */
std::string FormatTimingSample(const std::string& prefix, const TimingSample& sample);

/**
 * @brief 将窗口长度、样本组数及各阶段分布格式化为 key=value 日志字段。
 * @param summary 已完成窗口的汇总，时间分别以秒和毫秒表示。
 * @return 含每种计时的 count、mean、p50、p95、p99、max 字段的文本。
 * @see MultiStageTimingWindow, TimingWindowSummary
 */
std::string FormatTimingWindow(const TimingWindowSummary& summary);

}  // namespace lightning::profiling
