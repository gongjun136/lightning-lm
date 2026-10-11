#pragma once

#include <chrono>
#include <cstdint>
#include <limits>

namespace lightning::profiling {

// Diagnostic-only bounded asynchronous trace. Literal kind/source/reason strings
// must outlive the process. Unmeasured numeric fields are NaN, never fake zeros.
struct CausalEvent {
    const char* kind = "unknown";  ///< 事件类型标识，字符串须在进程生命周期内有效。
    const char* source = "unknown";  ///< 事件产生位置标识，字符串须在进程生命周期内有效。
    const char* reason = "none";  ///< 接受、拒绝或状态变化原因，字符串须在进程生命周期内有效。
    std::uint64_t frame_id = 0;  ///< 与事件关联的帧 ID。
    double stamp = std::numeric_limits<double>::quiet_NaN();  ///< 事件对应的传感器时间戳（s），未测量时为 NaN。
    double input_stamp = std::numeric_limits<double>::quiet_NaN();  ///< 本次输入的时间戳（s），未测量时为 NaN。
    double prior_stamp = std::numeric_limits<double>::quiet_NaN();  ///< 用于比较的此前状态时间戳（s）。
    double before_v = std::numeric_limits<double>::quiet_NaN();  ///< 状态更新前的速度值（m/s）。
    double after_v = std::numeric_limits<double>::quiet_NaN();  ///< 状态更新后的速度值（m/s）。
    double measured_v = std::numeric_limits<double>::quiet_NaN();  ///< 参与观测的量测速度（m/s）。
    double innovation = std::numeric_limits<double>::quiet_NaN();  ///< 速度量测与预测之差（m/s）。
    double nis = std::numeric_limits<double>::quiet_NaN();  ///< 归一化新息平方。
    double stddev = std::numeric_limits<double>::quiet_NaN();  ///< 该观测使用的速度标准差（m/s）。
    double torque = std::numeric_limits<double>::quiet_NaN();  ///< 与事件关联的电机扭矩值。
    double duration_ms = std::numeric_limits<double>::quiet_NaN();  ///< 事件对应处理阶段的总耗时（ms）。
    double queue_wait_ms = std::numeric_limits<double>::quiet_NaN();  ///< 入队到开始处理的等待耗时（ms），包含队列与唤醒等待。
    double thread_cpu_ms = std::numeric_limits<double>::quiet_NaN();  ///< 该阶段调用线程的 CPU 耗时（ms）。
    double position_delta = std::numeric_limits<double>::quiet_NaN();  ///< 事件前后的位置变化幅度（m）。
    double yaw_delta = std::numeric_limits<double>::quiet_NaN();  ///< 事件前后的航向变化量。
    int accepted = -1;  ///< 观测是否被接受，-1 表示该事件未记录此项。
    int hold = -1;  ///< 是否处于静止保持，-1 表示该事件未记录此项。
};

bool CausalTraceEnabled();
void RecordCausalEvent(CausalEvent event);
void TraceLockWait(const char* source, double stamp,
                   std::chrono::steady_clock::time_point started);

}  // namespace lightning::profiling
