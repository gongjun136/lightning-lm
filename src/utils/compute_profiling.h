#pragma once

#include <chrono>
#include <cstddef>
#include <initializer_list>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

namespace lightning::profiling {

constexpr const char* kComputeProfilingEnv = "LIGHTNING_LM_COMPUTE_PROFILE";
constexpr const char* kReduceNonessentialOverheadEnv =
    "LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD";

bool ParseBooleanFlag(const char* value);
bool ComputeProfilingEnabled();
bool ReduceNonessentialOverhead();

struct TimingSample {
    bool valid = false;  ///< 该耗时样本是否由实际计时得到。
    double wall_ms = 0.0;  ///< 经过的单调时钟耗时或其分布（ms），包含运行与等待。
    double thread_cpu_ms = 0.0;  ///< 调用线程实际消耗的 CPU 时间或其分布（ms）。
    // CPU consumed by the complete process while the stage is in flight. This
    // includes concurrent workers and unrelated process threads, so it is an
    // attribution upper bound rather than exclusive stage CPU.
    double process_cpu_ms = 0.0;  ///< 阶段期间整个进程消耗的 CPU 时间或其分布（ms），包含并发线程。
};

class Stopwatch {
   public:
    explicit Stopwatch(bool enabled = ComputeProfilingEnabled());
    TimingSample Stop();

   private:
    bool enabled_ = false;  ///< 是否启用本次计时。
    bool stopped_ = false;  ///< 本次计时是否已结束，重复 Stop() 复用结果。
    std::chrono::steady_clock::time_point wall_start_;  ///< 本次计时开始的单调时钟时刻。
    double thread_cpu_start_ms_ = 0.0;  ///< 开始时调用线程的 CPU 累计耗时（ms）。
    double process_cpu_start_ms_ = 0.0;  ///< 开始时整个进程的 CPU 累计耗时（ms）。
    TimingSample result_;  ///< 停止计时后保存的耗时样本。
};

struct DistributionSummary {
    std::size_t count = 0;  ///< 参与该分布统计的有效样本数。
    double mean = 0.0;  ///< 样本的算术平均值。
    double p50 = 0.0;  ///< 样本的第 50 百分位数。
    double p95 = 0.0;  ///< 样本的第 95 百分位数。
    double p99 = 0.0;  ///< 样本的第 99 百分位数。
    double max = 0.0;  ///< 样本中的最大值。
};

struct StageTimingSummary {
    std::string name;  ///< 该组耗时统计对应的处理阶段名称。
    DistributionSummary wall_ms;  ///< 经过的单调时钟耗时或其分布（ms），包含运行与等待。
    DistributionSummary thread_cpu_ms;  ///< 调用线程实际消耗的 CPU 时间或其分布（ms）。
    DistributionSummary process_cpu_ms;  ///< 阶段期间整个进程消耗的 CPU 时间或其分布（ms），包含并发线程。
};

struct TimingWindowSummary {
    double window_s = 0.0;  ///< 当前统计窗口实际覆盖的时间（s）。
    std::size_t samples = 0;  ///< 该窗口包含的样本组数。
    std::vector<StageTimingSummary> stages;  ///< 按固定阶段顺序排列的耗时统计结果。
};

// Thread-safe fixed-schema timing window. The first sample opens a window;
// the first sample at least window_seconds later closes and returns it.
class MultiStageTimingWindow {
   public:
    explicit MultiStageTimingWindow(std::vector<std::string> stage_names,
                                    double window_seconds = 1.0);

    std::optional<TimingWindowSummary> Add(
        std::initializer_list<TimingSample> samples,
        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now());
    std::optional<TimingWindowSummary> Add(
        const std::vector<TimingSample>& samples,
        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now());

   private:
    struct StageSamples {
        std::vector<double> wall_ms;  ///< 经过的单调时钟耗时或其分布（ms），包含运行与等待。
        std::vector<double> thread_cpu_ms;  ///< 调用线程实际消耗的 CPU 时间或其分布（ms）。
        std::vector<double> process_cpu_ms;  ///< 阶段期间整个进程消耗的 CPU 时间或其分布（ms），包含并发线程。
    };

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

std::string FormatTimingSample(const std::string& prefix, const TimingSample& sample);
std::string FormatTimingWindow(const TimingWindowSummary& summary);

}  // namespace lightning::profiling
