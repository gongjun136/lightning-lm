#include "utils/online_logging.h"
#include "utils/bounded_async_queue.h"
#include "utils/causal_trace.h"

#include <glog/logging.h>
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <unistd.h>

namespace lightning {
namespace {
struct LogRow {
    std::uint64_t sequence = 0;
    std::size_t size = 0;
    bool truncated = false;
    std::array<char, 8192> text{};
};

}  // namespace

class AsyncOnlineLogger : public google::base::Logger {
 public:
    AsyncOnlineLogger() : queue_(2048, [this](const auto& batch) { WriteBatch(batch); }, 128) {
        // Own a separate stdio buffer; preserve the launcher's stderr path.
        const int fd = dup(STDERR_FILENO);
        if (fd >= 0) {
            file_ = fdopen(fd, "w");
            if (!file_) close(fd);
        }
        if (!file_) {
            // Initialization is before localization starts. Never silently fall
            // back to synchronous logging on a real-time thread.
            throw std::runtime_error("cannot initialize asynchronous stderr logger");
        }
        setvbuf(file_, nullptr, _IOFBF, 64 * 1024);
        queue_.Start();
    }
    ~AsyncOnlineLogger() override { Stop(); }
    void Stop() {
        if (!file_) return;
        queue_.Stop();
        std::fprintf(file_, "ASYNC_LOG_SUMMARY received=%llu written=%llu dropped=%llu critical_dropped=%llu io_errors=%llu truncated=%llu\n",
            static_cast<unsigned long long>(queue_.Received()),
            static_cast<unsigned long long>(written_),
            static_cast<unsigned long long>(queue_.Dropped()),
            static_cast<unsigned long long>(queue_.CriticalDropped()),
            static_cast<unsigned long long>(io_errors_),
            static_cast<unsigned long long>(truncated_));
        std::fclose(file_);
        file_ = nullptr;
    }
    void Write(bool, time_t, const char* message, int length) override {
        if (length <= 0) return;
        LogRow row;
        row.size = std::min(static_cast<std::size_t>(length), row.text.size());
        row.truncated = row.size < static_cast<std::size_t>(length);
        if (row.truncated) ++truncated_;
        std::memcpy(row.text.data(), message, row.size);
        const bool critical = message[0] == 'W' || message[0] == 'E' || message[0] == 'F';
        queue_.TryPush(std::move(row), critical);
    }
    // glog may Flush under its global lock, including from producer threads.
    // The writer flushes every batch; this callback must never perform I/O.
    void Flush() override {}
    google::uint32 LogSize() override { return 0; }

 private:
    void WriteBatch(const std::vector<LogRow>& batch) {
        for (const auto& row : batch) {
            auto size = row.size;
            if (size && row.text[size - 1] == '\n') --size;
            if (std::fwrite(row.text.data(), 1, size, file_) != size ||
                std::fprintf(file_, " async_seq=%llu%s\n",
                    static_cast<unsigned long long>(row.sequence),
                    row.truncated ? " async_truncated=1" : "") < 0) {
                ++io_errors_;
            } else {
                ++written_;
            }
        }
        if (!batch.empty() && std::fflush(file_) != 0) ++io_errors_;
    }
    FILE* file_ = nullptr;
    std::uint64_t written_ = 0, io_errors_ = 0;
    std::atomic<std::uint64_t> truncated_{0};
    BoundedAsyncQueue<LogRow> queue_;
};

namespace {
class DiscardDuplicateLogger : public google::base::Logger {
 public:
    void Write(bool, time_t, const char*, int) override {}
    void Flush() override {}
    google::uint32 LogSize() override { return 0; }
};
}  // namespace

OnlineLogging::OnlineLogging() {
    // glog sends every severity to the INFO destination as well. Only that
    // destination queues records, avoiding duplicate WARNING/ERROR messages.
    logger_ = new AsyncOnlineLogger;
    FLAGS_logtostderr = false;
    FLAGS_alsologtostderr = false;
    FLAGS_stderrthreshold = google::FATAL;
    FLAGS_colorlogtostderr = false;
    FLAGS_log_prefix = true;
    for (int severity = google::WARNING; severity <= google::FATAL; ++severity)
        google::base::SetLogger(severity, new DiscardDuplicateLogger);
    google::base::SetLogger(google::INFO, logger_);
    // Open/allocate the existing diagnostic sink before the first callback.
    profiling::CausalTraceEnabled();
}
OnlineLogging::~OnlineLogging() {
    logger_->Stop();
    google::ShutdownGoogleLogging();
}
}  // namespace lightning
