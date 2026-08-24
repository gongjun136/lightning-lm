#include "utils/compute_profiling.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void Require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

bool Near(double lhs, double rhs, double tolerance = 1e-9) {
    return std::abs(lhs - rhs) <= tolerance;
}

}  // namespace

int main() {
    using lightning::profiling::MultiStageTimingWindow;
    using lightning::profiling::ParseBooleanFlag;
    using lightning::profiling::TimingSample;

    Require(!ParseBooleanFlag(nullptr), "null flag is false");
    Require(ParseBooleanFlag("1"), "numeric true flag");
    Require(ParseBooleanFlag("TRUE"), "case-insensitive true flag");
    Require(ParseBooleanFlag("yes"), "yes flag");
    Require(ParseBooleanFlag("On"), "on flag");
    Require(!ParseBooleanFlag("0"), "numeric false flag");
    Require(!ParseBooleanFlag("unexpected"), "unknown flag is false");

    MultiStageTimingWindow window({"first", "second"});
    const auto start = std::chrono::steady_clock::now();
    Require(!window.Add({TimingSample{true, 1.0, 0.5, 2.0},
                         TimingSample{true, 3.0, 1.5, 4.0}},
                        start),
            "first sample keeps the window open");
    const auto summary = window.Add(
        {TimingSample{true, 5.0, 2.5, 6.0}, TimingSample{true, 7.0, 3.5, 8.0}},
        start + std::chrono::milliseconds(1100));
    Require(summary.has_value(), "elapsed window emits a summary");
    Require(summary->samples == 2, "summary sample count");
    Require(summary->stages.size() == 2, "summary stage count");
    Require(Near(summary->stages[0].wall_ms.mean, 3.0), "wall mean");
    Require(Near(summary->stages[0].wall_ms.p50, 3.0), "wall median");
    Require(Near(summary->stages[1].thread_cpu_ms.max, 3.5), "thread CPU max");
    Require(Near(summary->stages[1].process_cpu_ms.p95, 7.8), "process CPU p95");
    Require(lightning::profiling::FormatTimingWindow(*summary).find(
                "first_wall_count=2") != std::string::npos,
            "formatted stage sample count");

    bool rejected = false;
    try {
        window.Add({TimingSample{true, 1.0, 1.0, 1.0}}, start);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    Require(rejected, "schema mismatch is rejected");

    std::cout << "compute profiling tests passed\n";
    return 0;
}
