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
    bool valid = false;
    double wall_ms = 0.0;
    double thread_cpu_ms = 0.0;
    // CPU consumed by the complete process while the stage is in flight. This
    // includes concurrent workers and unrelated process threads, so it is an
    // attribution upper bound rather than exclusive stage CPU.
    double process_cpu_ms = 0.0;
};

class Stopwatch {
   public:
    explicit Stopwatch(bool enabled = ComputeProfilingEnabled());
    TimingSample Stop();

   private:
    bool enabled_ = false;
    bool stopped_ = false;
    std::chrono::steady_clock::time_point wall_start_;
    double thread_cpu_start_ms_ = 0.0;
    double process_cpu_start_ms_ = 0.0;
    TimingSample result_;
};

struct DistributionSummary {
    std::size_t count = 0;
    double mean = 0.0;
    double p50 = 0.0;
    double p95 = 0.0;
    double p99 = 0.0;
    double max = 0.0;
};

struct StageTimingSummary {
    std::string name;
    DistributionSummary wall_ms;
    DistributionSummary thread_cpu_ms;
    DistributionSummary process_cpu_ms;
};

struct TimingWindowSummary {
    double window_s = 0.0;
    std::size_t samples = 0;
    std::vector<StageTimingSummary> stages;
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
        std::vector<double> wall_ms;
        std::vector<double> thread_cpu_ms;
        std::vector<double> process_cpu_ms;
    };

    std::optional<TimingWindowSummary> AddLocked(
        const TimingSample* samples, std::size_t sample_count,
        std::chrono::steady_clock::time_point now);

    std::vector<std::string> stage_names_;
    std::vector<StageSamples> stage_samples_;
    double window_seconds_ = 1.0;
    std::chrono::steady_clock::time_point window_start_;
    bool window_started_ = false;
    std::size_t samples_ = 0;
    std::mutex mutex_;
};

std::string FormatTimingSample(const std::string& prefix, const TimingSample& sample);
std::string FormatTimingWindow(const TimingWindowSummary& summary);

}  // namespace lightning::profiling
