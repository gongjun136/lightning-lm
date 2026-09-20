#include "utils/compute_profiling.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <time.h>
#include <utility>

namespace lightning::profiling {
namespace {

double ClockMilliseconds(clockid_t clock_id) {
    timespec value{};
    if (clock_gettime(clock_id, &value) != 0) return 0.0;
    return static_cast<double>(value.tv_sec) * 1e3 +
           static_cast<double>(value.tv_nsec) * 1e-6;
}

double Quantile(const std::vector<double>& sorted, double probability) {
    if (sorted.empty()) return 0.0;
    const double position = probability * static_cast<double>(sorted.size() - 1);
    const auto lower = static_cast<std::size_t>(position);
    const auto upper = std::min(lower + 1, sorted.size() - 1);
    const double fraction = position - static_cast<double>(lower);
    return sorted[lower] * (1.0 - fraction) + sorted[upper] * fraction;
}

DistributionSummary Summarize(std::vector<double> values) {
    DistributionSummary summary;
    summary.count = values.size();
    if (values.empty()) return summary;
    summary.mean = std::accumulate(values.begin(), values.end(), 0.0) /
                   static_cast<double>(values.size());
    std::sort(values.begin(), values.end());
    summary.p50 = Quantile(values, 0.50);
    summary.p95 = Quantile(values, 0.95);
    summary.p99 = Quantile(values, 0.99);
    summary.max = values.back();
    return summary;
}

void AppendDistribution(std::ostringstream& stream, const std::string& prefix,
                        const DistributionSummary& summary) {
    stream << ' ' << prefix << "_count=" << summary.count
           << ' ' << prefix << "_mean_ms=" << summary.mean
           << ' ' << prefix << "_p50_ms=" << summary.p50
           << ' ' << prefix << "_p95_ms=" << summary.p95
           << ' ' << prefix << "_p99_ms=" << summary.p99
           << ' ' << prefix << "_max_ms=" << summary.max;
}

}  // namespace

bool ParseBooleanFlag(const char* value) {
    if (value == nullptr) return false;
    std::string normalized(value);
    std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                   [](unsigned char character) {
                       return static_cast<char>(std::tolower(character));
                   });
    return normalized == "1" || normalized == "true" || normalized == "yes" ||
           normalized == "on";
}

bool ComputeProfilingEnabled() {
    static const bool enabled = ParseBooleanFlag(std::getenv(kComputeProfilingEnv));
    return enabled;
}

bool ReduceNonessentialOverhead() {
    static const bool enabled =
        ParseBooleanFlag(std::getenv(kReduceNonessentialOverheadEnv));
    return enabled;
}

Stopwatch::Stopwatch(bool enabled) : enabled_(enabled) {
    if (!enabled_) return;
    wall_start_ = std::chrono::steady_clock::now();
    thread_cpu_start_ms_ = ClockMilliseconds(CLOCK_THREAD_CPUTIME_ID);
    process_cpu_start_ms_ = ClockMilliseconds(CLOCK_PROCESS_CPUTIME_ID);
}

TimingSample Stopwatch::Stop() {
    if (!enabled_) return result_;
    if (stopped_) return result_;
    const double thread_cpu_end_ms = ClockMilliseconds(CLOCK_THREAD_CPUTIME_ID);
    const double process_cpu_end_ms = ClockMilliseconds(CLOCK_PROCESS_CPUTIME_ID);
    const auto wall_end = std::chrono::steady_clock::now();
    result_.valid = true;
    result_.wall_ms =
        std::chrono::duration<double, std::milli>(wall_end - wall_start_).count();
    result_.thread_cpu_ms = std::max(0.0, thread_cpu_end_ms - thread_cpu_start_ms_);
    result_.process_cpu_ms = std::max(0.0, process_cpu_end_ms - process_cpu_start_ms_);
    stopped_ = true;
    return result_;
}

MultiStageTimingWindow::MultiStageTimingWindow(std::vector<std::string> stage_names,
                                               double window_seconds)
    : stage_names_(std::move(stage_names)),
      stage_samples_(stage_names_.size()),
      window_seconds_(window_seconds) {
    if (stage_names_.empty()) {
        throw std::invalid_argument("timing window requires at least one stage");
    }
    if (!(window_seconds_ > 0.0)) {
        throw std::invalid_argument("timing window duration must be positive");
    }
}

std::optional<TimingWindowSummary> MultiStageTimingWindow::Add(
    std::initializer_list<TimingSample> samples,
    std::chrono::steady_clock::time_point now) {
    std::lock_guard<std::mutex> lock(mutex_);
    return AddLocked(samples.begin(), samples.size(), now);
}

std::optional<TimingWindowSummary> MultiStageTimingWindow::Add(
    const std::vector<TimingSample>& samples, std::chrono::steady_clock::time_point now) {
    std::lock_guard<std::mutex> lock(mutex_);
    return AddLocked(samples.data(), samples.size(), now);
}

std::optional<TimingWindowSummary> MultiStageTimingWindow::AddLocked(
    const TimingSample* samples, std::size_t sample_count,
    std::chrono::steady_clock::time_point now) {
    if (sample_count != stage_names_.size()) {
        throw std::invalid_argument("timing sample count does not match stage schema");
    }
    if (!window_started_) {
        window_start_ = now;
        window_started_ = true;
    }
    for (std::size_t index = 0; index < sample_count; ++index) {
        if (!samples[index].valid) continue;
        stage_samples_[index].wall_ms.push_back(samples[index].wall_ms);
        stage_samples_[index].thread_cpu_ms.push_back(samples[index].thread_cpu_ms);
        stage_samples_[index].process_cpu_ms.push_back(samples[index].process_cpu_ms);
    }
    ++samples_;

    const double elapsed_seconds =
        std::chrono::duration<double>(now - window_start_).count();
    if (elapsed_seconds < window_seconds_) return std::nullopt;

    TimingWindowSummary summary;
    summary.window_s = elapsed_seconds;
    summary.samples = samples_;
    summary.stages.reserve(stage_names_.size());
    for (std::size_t index = 0; index < stage_names_.size(); ++index) {
        StageTimingSummary stage;
        stage.name = stage_names_[index];
        stage.wall_ms = Summarize(std::move(stage_samples_[index].wall_ms));
        stage.thread_cpu_ms = Summarize(std::move(stage_samples_[index].thread_cpu_ms));
        stage.process_cpu_ms = Summarize(std::move(stage_samples_[index].process_cpu_ms));
        summary.stages.push_back(std::move(stage));
        stage_samples_[index] = StageSamples{};
    }
    window_start_ = now;
    samples_ = 0;
    return summary;
}

std::string FormatTimingSample(const std::string& prefix, const TimingSample& sample) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(6)
           << prefix << "_wall_ms=" << sample.wall_ms
           << ' ' << prefix << "_thread_cpu_ms=" << sample.thread_cpu_ms
           << ' ' << prefix << "_process_cpu_ms_concurrent=" << sample.process_cpu_ms;
    return stream.str();
}

std::string FormatTimingWindow(const TimingWindowSummary& summary) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(6)
           << "window_s=" << summary.window_s << " samples=" << summary.samples;
    for (const auto& stage : summary.stages) {
        AppendDistribution(stream, stage.name + "_wall", stage.wall_ms);
        AppendDistribution(stream, stage.name + "_thread_cpu", stage.thread_cpu_ms);
        AppendDistribution(stream, stage.name + "_process_cpu_concurrent",
                           stage.process_cpu_ms);
    }
    return stream.str();
}

}  // namespace lightning::profiling
