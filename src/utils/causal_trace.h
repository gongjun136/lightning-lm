#pragma once

/**
 * @file
 * @brief 因果诊断事件定义、异步追踪入口和锁等待采样接口。
 * @see lightning::profiling::CausalEvent, @ref causal_diagnostics "因果链与时间关联"
 */

#include <chrono>
#include <cstdint>
#include <limits>

namespace lightning::profiling {

/**
 * @brief 供不同处理阶段复用的因果诊断事件，只解释该事件实际赋值的字段。
 * @details RecordCausalEvent() 按值复制事件后异步写入；kind、source、reason 仅复制指针，
 * 应使用进程生命周期内有效的非空字符串，内容不得包含 CSV 分隔逗号或换行。
 * @note 未测量的浮点字段默认是 NaN，accepted/hold 默认是 -1，均不能解释为零。
 * 个别调用点可用负的 queue_wait_ms 表示未取得入队时刻，见该字段说明。
 * @see RecordCausalEvent(), TraceLockWait()
 */
struct CausalEvent {
    const char* kind = "unknown";  ///< 事件类型标识，字符串须在进程生命周期内有效。
    const char* source = "unknown";  ///< 事件产生位置标识，字符串须在进程生命周期内有效。
    const char* reason = "none";  ///< 接受、拒绝或状态变化原因，字符串须在进程生命周期内有效。
    std::uint64_t frame_id = 0;  ///< 关联帧标识，编码由事件产生方指定；默认值 0 不表示 CSV 入队序号。
    double stamp = std::numeric_limits<double>::quiet_NaN();  ///< 事件关联的状态或传感器时间戳（s），未测量时为 NaN。
    double input_stamp = std::numeric_limits<double>::quiet_NaN();  ///< 本次输入的时间戳（s），未测量时为 NaN。
    double prior_stamp = std::numeric_limits<double>::quiet_NaN();  ///< 用于比较的此前状态时间戳（s）。
    double before_v = std::numeric_limits<double>::quiet_NaN();  ///< 状态更新前的速度值（m/s）。
    double after_v = std::numeric_limits<double>::quiet_NaN();  ///< 状态更新后的速度值（m/s）。
    double measured_v = std::numeric_limits<double>::quiet_NaN();  ///< 参与观测的量测速度（m/s）。
    double innovation = std::numeric_limits<double>::quiet_NaN();  ///< 速度量测与预测之差（m/s）。
    double nis = std::numeric_limits<double>::quiet_NaN();  ///< 归一化新息平方，无量纲。
    double stddev = std::numeric_limits<double>::quiet_NaN();  ///< 该观测使用的速度标准差（m/s）。
    double torque = std::numeric_limits<double>::quiet_NaN();  ///< 与事件关联的电机扭矩值。
    double duration_ms = std::numeric_limits<double>::quiet_NaN();  ///< 事件对应处理阶段的总耗时（ms）。
    double queue_wait_ms = std::numeric_limits<double>::quiet_NaN();  ///< 入队前采样到处理入口的耗时（ms），含入队锁、队列及唤醒等待；未取得时刻时调用方可填 -1。
    double thread_cpu_ms = std::numeric_limits<double>::quiet_NaN();  ///< 该阶段调用线程的 CPU 耗时（ms）。
    double position_delta = std::numeric_limits<double>::quiet_NaN();  ///< 事件前后的位置变化幅度（m）。
    double yaw_delta = std::numeric_limits<double>::quiet_NaN();  ///< 事件前后的航向变化量（rad）。
    int accepted = -1;  ///< 观测接受状态：1 接受、0 拒绝、-1 未记录或不适用。
    int hold = -1;  ///< 静止保持状态：1 开启、0 关闭、-1 未记录或不适用。
};

/**
 * @brief 获取共享追踪器并查询当前是否允许接收事件，首次调用会触发初始化。
 * @return 输出文件已成功打开且后台未报告 I/O 失败时为 true。
 * @note 首次初始化读取 LIGHTNING_LM_CAUSAL_TRACE_PATH，并以仅创建新文件的方式打开。
 * 路径未设置、为空或打开失败时关闭追踪，后续调用不重新读取路径或重试打开。
 * 在线入口应在启动传感器回调前调用，以提前完成文件创建和后台线程启动。
 * @see RecordCausalEvent()
 */
bool CausalTraceEnabled();

/**
 * @brief 尝试将事件交给共享追踪器，补充本机记录时刻和线程 ID 后异步写入 CSV。
 * @param event 待记录事件；字符串指针须满足 CausalEvent 的生命周期约束。
 * @note 未启用或发生 I/O 失败时直接返回；队列满或入队锁竞争时允许丢弃记录。
 * 此接口也可能触发首次初始化，调用方不能据其返回推断事件已经写入文件。
 * @see CausalTraceEnabled(), CausalEvent
 */
void RecordCausalEvent(CausalEvent event);

/**
 * @brief 记录从 started 到当前的经过时间达到 0.2 ms 的锁等待事件。
 * @param source 锁或调用位置标识，字符串须在进程生命周期内有效。
 * @param stamp 关联的状态或传感器时间戳（s），不参与等待时间计算。
 * @param started 等待开始的单调时钟时刻；默认构造的空时刻表示不采样。
 * @details 调用方在加锁前采样 started，获锁后调用本函数；空时刻或耗时不足阈值时返回。
 * 事件 kind 为 lock_wait，经过时间包含阻塞、调度等待及调用间开销。
 * @see RecordCausalEvent()
 */
void TraceLockWait(const char* source, double stamp,
                   std::chrono::steady_clock::time_point started);

}  // namespace lightning::profiling
