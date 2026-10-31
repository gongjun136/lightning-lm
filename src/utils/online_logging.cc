#include "utils/online_logging.h"
#include "utils/bounded_async_queue.h"
#include "utils/causal_trace.h"
#include "utils/compute_profiling.h"

#include <glog/logging.h>
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string_view>
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

// Each destination has its own queue and writer. A stalled profile file must
// not consume the operational log's capacity or reserved critical slots.
class AsyncLogWriter {
 public:
    explicit AsyncLogWriter(const std::string& path = {})
        : profile_(!path.empty()),
          queue_(2048, [this](const auto& batch) { WriteBatch(batch); }, profile_ ? 0 : 128) {
        if (profile_) {
            file_ = std::fopen(path.c_str(), "wx");
            if (!file_) throw std::runtime_error("cannot create compute profile: " + path);
            setvbuf(file_, nullptr, _IOFBF, 64 * 1024);
            if (std::fprintf(file_, "# COMPUTE_PROFILE schema=1 sampling=all sequence=profile_seq\n") < 0 ||
                std::fflush(file_) != 0) {
                std::fclose(file_);
                file_ = nullptr;
                throw std::runtime_error("cannot write compute profile header: " + path);
            }
        } else {
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
        }
        queue_.Start();
    }
    ~AsyncLogWriter() { Stop(); }
    std::string Summary() const {
        std::ostringstream stream;
        stream << (profile_ ? "COMPUTE_PROFILE_SUMMARY" : "ASYNC_LOG_SUMMARY")
               << " received=" << queue_.Received() << " written=" << written_
               << " dropped=" << queue_.Dropped() << " critical_dropped=" << queue_.CriticalDropped()
               << " io_errors=" << io_errors_ << " truncated=" << truncated_.load();
        return stream.str();
    }
    void Stop() {
        if (!file_) return;
        queue_.Stop();
        if (std::fprintf(file_, "%s\n", Summary().c_str()) < 0 || std::fflush(file_) != 0) ++io_errors_;
        if (std::fclose(file_) != 0) ++io_errors_;
        file_ = nullptr;
    }
    void Write(const char* message, int length) {
        if (length <= 0) return;
        LogRow row;
        row.size = std::min(static_cast<std::size_t>(length), row.text.size());
        row.truncated = row.size < static_cast<std::size_t>(length);
        if (row.truncated) ++truncated_;
        std::memcpy(row.text.data(), message, row.size);
        const bool critical = message[0] == 'W' || message[0] == 'E' || message[0] == 'F';
        queue_.TryPush(std::move(row), critical);
    }
 private:
    void WriteBatch(const std::vector<LogRow>& batch) {
        for (const auto& row : batch) {
            auto size = row.size;
            if (size && row.text[size - 1] == '\n') --size;
            if (std::fwrite(row.text.data(), 1, size, file_) != size ||
                (profile_ && std::fprintf(file_, " profile_seq=%llu",
                    static_cast<unsigned long long>(row.sequence)) < 0) ||
                std::fprintf(file_, "%s\n",
                    row.truncated ? " async_truncated=1" : "") < 0) {
                ++io_errors_;
            } else {
                ++written_;
            }
        }
        if (!batch.empty() && std::fflush(file_) != 0) ++io_errors_;
        // This runs on the profile writer, outside glog's producer lock. Never
        // log recursively from Write(), which glog calls under that lock.
        if (profile_ && !reported_incomplete_ && (io_errors_ || queue_.Dropped() || truncated_.load())) {
            reported_incomplete_ = true;
            LOG(ERROR) << "COMPUTE_PROFILE_INCOMPLETE " << Summary();
        }
    }
    const bool profile_;
    bool reported_incomplete_ = false;
    FILE* file_ = nullptr;
    std::uint64_t written_ = 0, io_errors_ = 0;
    std::atomic<std::uint64_t> truncated_{0};
    BoundedAsyncQueue<LogRow> queue_;
};

class AsyncOnlineLogger : public google::base::Logger {
 public:
    AsyncOnlineLogger() {
        if (!profiling::ComputeProfilingEnabled()) return;
        const char* configured = std::getenv("LIGHTNING_LM_COMPUTE_PROFILE_PATH");
        profile_path_ = std::filesystem::absolute(
            configured && *configured ? configured : "compute_profile.log").string();
        profile_ = std::make_unique<AsyncLogWriter>(profile_path_);
    }
    ~AsyncOnlineLogger() override { Stop(); }
    void Stop() {
        if (stopped_) return;
        stopped_ = true;
        if (profile_) {
            profile_->Stop();
            LOG(INFO) << "COMPUTE_PROFILE_STATUS path=" << profile_path_ << ' ' << profile_->Summary();
        }
        operational_.Stop();
    }
    void AnnounceProfile() const {
        LOG(INFO) << "COMPUTE_PROFILE_STATUS enabled=" << (profile_ ? 1 : 0)
                  << " storage=dedicated-v1 path=" << (profile_ ? profile_path_ : "<disabled>");
    }
    void Write(bool, time_t, const char* message, int length) override {
        if (length <= 0) return;
        const std::string_view line(message, length);
        const auto prefix = line.find("] ");
        const auto body = prefix == std::string_view::npos ? line : line.substr(prefix + 2);
        const bool timing = body.rfind("COMPUTE_BENCH_", 0) == 0 || body.rfind("LIO_BENCH_FRAME ", 0) == 0;
        const bool operational_summary =
            body.rfind("COMPUTE_BENCH_SUMMARY module=lidar_end_to_end ", 0) == 0;
        if (profile_ && timing && !operational_summary) {
            profile_->Write(message, length);
        } else {
            operational_.Write(message, length);
        }
    }
    // glog may Flush under its global lock, including from producer threads.
    // Writers flush every batch; this callback must never perform I/O.
    void Flush() override {}
    google::uint32 LogSize() override { return 0; }

 private:
    AsyncLogWriter operational_;
    std::unique_ptr<AsyncLogWriter> profile_;
    std::string profile_path_;
    bool stopped_ = false;
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
    logger_->AnnounceProfile();
    // Open/allocate the existing diagnostic sink before the first callback.
    profiling::CausalTraceEnabled();
}
OnlineLogging::~OnlineLogging() {
    logger_->Stop();
    google::ShutdownGoogleLogging();
}
}  // namespace lightning
