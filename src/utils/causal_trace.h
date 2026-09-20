#pragma once

#include <chrono>
#include <cstdint>
#include <limits>

namespace lightning::profiling {

// Diagnostic-only bounded asynchronous trace. Literal kind/source/reason strings
// must outlive the process. Unmeasured numeric fields are NaN, never fake zeros.
struct CausalEvent {
    const char* kind = "unknown";
    const char* source = "unknown";
    const char* reason = "none";
    std::uint64_t frame_id = 0;
    double stamp = std::numeric_limits<double>::quiet_NaN();
    double input_stamp = std::numeric_limits<double>::quiet_NaN();
    double prior_stamp = std::numeric_limits<double>::quiet_NaN();
    double before_v = std::numeric_limits<double>::quiet_NaN();
    double after_v = std::numeric_limits<double>::quiet_NaN();
    double measured_v = std::numeric_limits<double>::quiet_NaN();
    double innovation = std::numeric_limits<double>::quiet_NaN();
    double nis = std::numeric_limits<double>::quiet_NaN();
    double stddev = std::numeric_limits<double>::quiet_NaN();
    double torque = std::numeric_limits<double>::quiet_NaN();
    double duration_ms = std::numeric_limits<double>::quiet_NaN();
    double queue_wait_ms = std::numeric_limits<double>::quiet_NaN();
    double thread_cpu_ms = std::numeric_limits<double>::quiet_NaN();
    double position_delta = std::numeric_limits<double>::quiet_NaN();
    double yaw_delta = std::numeric_limits<double>::quiet_NaN();
    int accepted = -1;
    int hold = -1;
};

bool CausalTraceEnabled();
void RecordCausalEvent(CausalEvent event);
void TraceLockWait(const char* source, double stamp,
                   std::chrono::steady_clock::time_point started);

}  // namespace lightning::profiling
