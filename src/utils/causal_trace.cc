/**
 * @file
 * @brief 共享因果追踪器及其有界队列、CSV 写入与锁等待采样实现。
 * @see causal_trace.h, lightning::BoundedAsyncQueue
 */

#include "utils/causal_trace.h"
#include "utils/bounded_async_queue.h"

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <sys/syscall.h>
#include <unistd.h>

namespace lightning::profiling {
namespace {
/// @brief 队列中的一条完整记录：业务事件加上 Record() 调用时的时间和线程信息。
struct Row {
    CausalEvent event;  ///< 生产者提供的业务事件，包含需异步读取的字符串指针。
    std::uint64_t sequence;  ///< 入队成功时由队列分配的序号，表示接收顺序。
    std::int64_t wall_ns;  ///< 系统时钟时间（ns），用于与外部日志关联。
    std::int64_t steady_ns;  ///< 单调时钟时间（ns），用于分析记录之间的时间间隔。
    long tid;  ///< 产生该记录的 Linux 线程 ID。
};

/**
 * @brief 接收因果诊断事件并异步写入 CSV 的共享输出器。
 * @details CausalTraceEnabled() 和 RecordCausalEvent() 均通过 TraceSink() 访问它。
 * 首次构造时读取 LIGHTNING_LM_CAUSAL_TRACE_PATH、创建文件并启动写线程；
 * Record() 在调用线程补充记录信息，
 * 后台线程从有界队列取出事件并调用 WriteBatch()，承担常规事件的文件 I/O。
 * @note 队列满或入队锁竞争时允许丢弃诊断事件；所有事件生产者结束后才可析构。
 * @see CausalEvent, lightning::BoundedAsyncQueue
 */
class Sink {
 public:
    /**
     * @brief 按环境变量初始化输出文件，并在文件可用时启动写线程。
     * @par 主要流程
     * 1. <b>检查输出路径</b>：仅在构造时读取路径，未配置或为空时返回。
     * 2. <b>创建追踪文件</b>：仅创建新文件，打开失败时返回，成功后写入 CSV 表头。
     * 3. <b>启动后台写入</b>：启动队列消费者，批次回调交给 WriteBatch()。
     */
    Sink() : queue_(kCapacity, [this](const auto& batch) { WriteBatch(batch); }) {
        // 阶段 1：检查输出路径。
        const char* path = std::getenv("LIGHTNING_LM_CAUSAL_TRACE_PATH");
        if (!path || !*path) return;
        // 阶段 2：创建追踪文件。
        // "wx" 只创建新文件，避免覆盖既有追踪；打开失败时保持追踪关闭。
        file_ = std::fopen(path, "wx");
        if (!file_) { std::perror("causal trace open failed"); return; }
        std::fprintf(file_, "sequence,wall_ns,steady_ns,tid,kind,source,reason,frame_id,stamp,input_stamp,prior_stamp,before_v,after_v,measured_v,innovation,nis,stddev,torque,duration_ms,position_delta,yaw_delta,accepted,hold,queue_wait_ms,thread_cpu_ms\n");
        // 阶段 3：启动后台写入。
        // 表头写完后再启动消费者；构造队列时注册的回调会调用 WriteBatch()。
        queue_.Start();
    }
    /**
     * @brief 等待后台线程处理剩余队列并退出，再输出计数、关闭文件。
     * @note 先停止所有生产者；join 完成后才读取后台线程维护的普通计数器。
     * I/O 已失败时剩余记录不再写入，末尾计数也不等同于数据持久化保证。
     */
    ~Sink() {
        if (!file_) return;
        queue_.Stop();
        std::fprintf(file_, "# received=%llu written=%llu dropped=%llu io_errors=%llu\n",
                     static_cast<unsigned long long>(queue_.Received()),
                     static_cast<unsigned long long>(written_),
                     static_cast<unsigned long long>(queue_.Dropped()),
                     static_cast<unsigned long long>(io_errors_));
        if (std::fclose(file_) != 0) std::perror("causal trace close failed");
    }
    /// @brief 查询文件已打开且后台未报告 I/O 失败的状态，不重新读取环境变量。
    bool Enabled() const { return file_ != nullptr && !failed_.load(); }

    /**
     * @brief 在事件生产线程补充本机记录时刻与线程 ID，再尝试入队。
     * @param event 业务事件，字符串生命周期由生产者保证。
     * @note TryPush() 不等待队列锁或磁盘，失败时由队列累计 dropped 计数。
     * 记录时刻与 event.stamp 的传感器/状态时间属于不同时间口径。
     * @see RecordCausalEvent(), WriteBatch()
     */
    void Record(CausalEvent event) {
        if (!Enabled()) return;
        const auto steady = std::chrono::steady_clock::now().time_since_epoch();
        const auto wall = std::chrono::system_clock::now().time_since_epoch();
        Row row{event, 0,
                std::chrono::duration_cast<std::chrono::nanoseconds>(wall).count(),
                std::chrono::duration_cast<std::chrono::nanoseconds>(steady).count(),
                syscall(SYS_gettid)};
        queue_.TryPush(row);
    }
 private:
    /**
     * @brief 在后台消费线程批量写入 CSV 并刷新 stdio 缓冲区。
     * @param batch 本批次记录，允许为空；由队列在未持有队列锁时交付。
     * @note 写入或刷新失败后设置 failed_，后续记录与未写批次会被跳过。
     * written_ 只累计 fprintf() 成功的记录，刷新失败时不能据此推断文件完整。
     */
    void WriteBatch(const std::vector<Row>& batch) {
        if (failed_) return;
        for (const auto& r : batch) {
            const auto& e = r.event;
            if (std::fprintf(file_, "%llu,%lld,%lld,%ld,%s,%s,%s,%llu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%.17g,%.17g\n",
                static_cast<unsigned long long>(r.sequence), static_cast<long long>(r.wall_ns),
                static_cast<long long>(r.steady_ns), r.tid, e.kind, e.source, e.reason,
                static_cast<unsigned long long>(e.frame_id), e.stamp, e.input_stamp, e.prior_stamp,
                e.before_v,e.after_v,e.measured_v,e.innovation,e.nis,e.stddev,e.torque,
                e.duration_ms,e.position_delta,e.yaw_delta,e.accepted,e.hold,
                e.queue_wait_ms,e.thread_cpu_ms) < 0) {
                ++io_errors_; failed_ = true; break;
            }
            ++written_;
        }
        if (std::fflush(file_) != 0) { ++io_errors_; failed_ = true; }
    }
    static constexpr std::size_t kCapacity = 8192;  ///< 待写缓冲区的事件数上限。
    FILE* file_ = nullptr;  ///< 构造时打开、析构时关闭，由后台线程写入事件。
    std::atomic<bool> failed_{false};  ///< 将后台 I/O 失败状态同步给事件生产线程。
    std::uint64_t written_ = 0 /**< 后台累计的 fprintf() 成功行数，join 后由析构读取。 */,
                  io_errors_ = 0 /**< 后台累计的事件写入及批次刷新错误数，join 后由析构读取。 */;
    BoundedAsyncQueue<Row> queue_;  ///< 连接事件生产线程与后台写线程的有界队列。
};
/**
 * @brief 取得本翻译单元的共享追踪器，函数内静态对象仅在首次调用时初始化。
 * @return 供所有追踪入口复用的 Sink 引用，存活至静态对象销毁阶段。
 * @note 在线入口应提前调用，以在传感器回调之前完成文件创建和写线程启动。
 * @see CausalTraceEnabled(), RecordCausalEvent()
 */
Sink& TraceSink() { static Sink sink; return sink; }
}  // namespace

bool CausalTraceEnabled() { return TraceSink().Enabled(); }
void RecordCausalEvent(CausalEvent event) { TraceSink().Record(event); }
void TraceLockWait(const char* source, double stamp,
                   std::chrono::steady_clock::time_point started) {
    // 空时刻用于调用方关闭采样；等待耗时始终由同一种单调时钟计算。
    if (started == std::chrono::steady_clock::time_point{}) return;
    const double ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    if (ms < 0.2) return;
    CausalEvent event;
    // 此处填写的是经过时间，不能将其当作线程实际消耗的 CPU 时间。
    event.kind = "lock_wait"; event.source = source;
    event.stamp = stamp; event.duration_ms = ms;
    RecordCausalEvent(event);
}
}  // namespace lightning::profiling
