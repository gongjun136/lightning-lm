#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <mutex>
#include <stdexcept>

namespace lightning::functional_safety {
// No ROS subscriptions or worker thread: the existing Safety timer evaluates
// this short, independent lock even when the sensor/assembler thread is idle.
class LidarFrameWatchdog {
 public:
    using Clock = std::chrono::steady_clock;
    struct Config {
        double warning_sec = 0.2, timeout_sec = 0.3, startup_grace_sec = 8.0;
        unsigned recovery_frames = 3;
    };
    struct Snapshot {
        bool started = false, seen = false, warning = false, fault = false;
        double age_sec = std::numeric_limits<double>::infinity();
        unsigned recovery_count = 0;
        std::uint64_t accepted = 0, rejected = 0;
    };
    explicit LidarFrameWatchdog(Config config) : config_(config) {
        if (!std::isfinite(config.warning_sec) || !std::isfinite(config.timeout_sec) ||
            !std::isfinite(config.startup_grace_sec) || config.warning_sec <= 0 ||
            config.timeout_sec <= config.warning_sec || config.startup_grace_sec < config.timeout_sec ||
            config.recovery_frames == 0)
            throw std::invalid_argument("invalid lidar frame watchdog thresholds");
    }
    // Idempotent: a pose/relocalization request cannot restart the grace period.
    void Start(Clock::time_point now = Clock::now()) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (!started_) { started_ = true; start_ = now; }
    }
    bool Observe(double stamp, Clock::time_point received_at, bool usable,
                 Clock::time_point now = Clock::now()) {
        std::lock_guard<std::mutex> lock(mutex_);
        Expire(now);
        const double age = std::chrono::duration<double>(now - received_at).count();
        if (!started_ || !usable || !std::isfinite(stamp) || stamp <= 0 ||
            (seen_ && (stamp <= last_stamp_ || received_at <= received_)) ||
            received_at < start_ || received_at > now || age >= config_.timeout_sec) {
            ++rejected_;
            recovery_count_ = 0;
            return false;
        }
        if (fault_) {
            const bool continuous = recovery_count_ != 0 &&
                received_at >= observed_at_ &&
                std::chrono::duration<double>(received_at - received_).count() <= config_.warning_sec;
            recovery_count_ = age <= config_.warning_sec ? (continuous ? recovery_count_ + 1 : 1) : 0;
            if (recovery_count_ >= config_.recovery_frames) { fault_ = false; recovery_count_ = 0; }
        }
        seen_ = true;
        last_stamp_ = stamp;
        received_ = received_at;  // Raw-input grouping progress; never an algorithm completion time.
        observed_at_ = now;  // A cached batch cannot count as consecutive live recovery arrivals.
        ++accepted_;
        return true;
    }
    Snapshot Evaluate(Clock::time_point now = Clock::now()) {
        std::lock_guard<std::mutex> lock(mutex_);
        Expire(now);
        Snapshot out;
        out.started = started_; out.seen = seen_; out.fault = fault_;
        if (seen_) out.age_sec = std::max(0.0, std::chrono::duration<double>(now - received_).count());
        out.warning = started_ && (fault_ || (seen_ && out.age_sec >= config_.warning_sec));
        out.recovery_count = recovery_count_; out.accepted = accepted_; out.rejected = rejected_;
        return out;
    }
 private:
    void Expire(Clock::time_point now) {
        if (started_ && std::chrono::duration<double>(now - (seen_ ? received_ : start_)).count() >=
                            (seen_ ? config_.timeout_sec : config_.startup_grace_sec)) {
            fault_ = true;
            recovery_count_ = 0;
        }
    }
    const Config config_;
    std::mutex mutex_;
    bool started_ = false, seen_ = false, fault_ = false;
    Clock::time_point start_{}, received_{}, observed_at_{};
    double last_stamp_ = 0;
    unsigned recovery_count_ = 0;
    std::uint64_t accepted_ = 0, rejected_ = 0;
};
}
