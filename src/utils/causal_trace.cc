#include "utils/causal_trace.h"

#include <atomic>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <thread>
#include <vector>
#include <sys/syscall.h>
#include <unistd.h>

namespace lightning::profiling {
namespace {
struct Row {
    CausalEvent event;
    std::uint64_t sequence;
    std::int64_t wall_ns;
    std::int64_t steady_ns;
    long tid;
};

class Sink {
 public:
    Sink() {
        const char* path = std::getenv("LIGHTNING_LM_CAUSAL_TRACE_PATH");
        if (!path || !*path) return;
        // Never overwrite an earlier test. Failure disables tracing, not localization.
        file_ = std::fopen(path, "wx");
        if (!file_) { std::perror("causal trace open failed"); return; }
        pending_.reserve(kCapacity);
        std::fprintf(file_, "sequence,wall_ns,steady_ns,tid,kind,source,reason,frame_id,stamp,input_stamp,prior_stamp,before_v,after_v,measured_v,innovation,nis,stddev,torque,duration_ms,position_delta,yaw_delta,accepted,hold,queue_wait_ms,thread_cpu_ms\n");
        worker_ = std::thread([this] { WriteLoop(); });
    }
    ~Sink() {
        if (!file_) return;
        { std::lock_guard<std::mutex> lock(mutex_); stopping_ = true; }
        wake_.notify_one();
        worker_.join();
        std::fprintf(file_, "# received=%llu written=%llu dropped=%llu io_errors=%llu\n",
                     static_cast<unsigned long long>(received_.load()),
                     static_cast<unsigned long long>(written_),
                     static_cast<unsigned long long>(dropped_.load()),
                     static_cast<unsigned long long>(io_errors_));
        if (std::fclose(file_) != 0) std::perror("causal trace close failed");
    }
    bool Enabled() const { return file_ != nullptr && !failed_.load(); }
    void Record(CausalEvent event) {
        if (!Enabled()) return;
        const auto sequence = received_.fetch_add(1);
        const auto steady = std::chrono::steady_clock::now().time_since_epoch();
        const auto wall = std::chrono::system_clock::now().time_since_epoch();
        Row row{event, sequence,
                std::chrono::duration_cast<std::chrono::nanoseconds>(wall).count(),
                std::chrono::duration_cast<std::chrono::nanoseconds>(steady).count(),
                syscall(SYS_gettid)};
        std::unique_lock<std::mutex> lock(mutex_, std::try_to_lock);
        if (!lock.owns_lock() || pending_.size() >= kCapacity) { ++dropped_; return; }
        pending_.push_back(row);
    }
 private:
    void WriteLoop() {
        std::vector<Row> batch;
        batch.reserve(kCapacity);
        while (true) {
            bool done;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                wake_.wait_for(lock, std::chrono::milliseconds(100), [this] { return stopping_; });
                pending_.swap(batch);
                done = stopping_;
            }
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
            batch.clear();
            if (done || failed_) break;
        }
    }
    static constexpr std::size_t kCapacity = 8192;
    FILE* file_ = nullptr;
    std::mutex mutex_;
    std::condition_variable wake_;
    std::vector<Row> pending_;
    std::thread worker_;
    std::atomic<std::uint64_t> received_{0}, dropped_{0};
    std::atomic<bool> failed_{false};
    std::uint64_t written_ = 0, io_errors_ = 0;
    bool stopping_ = false;
};
Sink& TraceSink() { static Sink sink; return sink; }
}  // namespace

bool CausalTraceEnabled() { return TraceSink().Enabled(); }
void RecordCausalEvent(CausalEvent event) { TraceSink().Record(event); }
void TraceLockWait(const char* source, double stamp,
                   std::chrono::steady_clock::time_point started) {
    if (started == std::chrono::steady_clock::time_point{}) return;
    const double ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    if (ms < 0.2) return;
    CausalEvent event;
    event.kind = "lock_wait"; event.source = source;
    event.stamp = stamp; event.duration_ms = ms;
    RecordCausalEvent(event);
}
}  // namespace lightning::profiling
