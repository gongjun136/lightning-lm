#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <mutex>

namespace lightning {

struct TimestampGateDecision {
    bool accepted = false;  ///< 本次输入是否通过门控。
    bool valid = false;  ///< 本次时间戳是否有限且为正值。
    double timestamp = 0.0;  ///< 本次待检查的时间戳（s）。
    double reference_timestamp = 0.0;  ///< 门控检查使用的此前时间戳水位（s）。
    double lag_sec = 0.0;  ///< 本次输入相对水位的滞后或回退量（s）。
};

/// Thread-safe strict monotonic timestamp gate for externally visible streams.
class MonotonicTimestampGate {
   public:
    TimestampGateDecision Observe(double timestamp) {
        std::lock_guard<std::mutex> lock(mutex_);
        TimestampGateDecision decision;
        decision.timestamp = timestamp;
        decision.reference_timestamp = last_accepted_timestamp_;
        decision.valid = std::isfinite(timestamp) && timestamp > 0.0;
        if (!decision.valid ||
            (last_accepted_timestamp_ > 0.0 && timestamp <= last_accepted_timestamp_)) {
            ++rejected_count_;
            if (decision.valid && last_accepted_timestamp_ > timestamp) {
                decision.lag_sec = last_accepted_timestamp_ - timestamp;
                worst_rollback_sec_ = std::max(worst_rollback_sec_, decision.lag_sec);
            }
            return decision;
        }
        last_accepted_timestamp_ = timestamp;
        decision.accepted = true;
        return decision;
    }

    void Reset() {
        std::lock_guard<std::mutex> lock(mutex_);
        last_accepted_timestamp_ = 0.0;
        rejected_count_ = 0;
        worst_rollback_sec_ = 0.0;
    }

    double LastAcceptedTimestamp() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return last_accepted_timestamp_;
    }

    std::uint64_t RejectedCount() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return rejected_count_;
    }

    double WorstRollbackSec() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return worst_rollback_sec_;
    }

   private:
    mutable std::mutex mutex_;  ///< 保护时间戳水位与拒绝统计。
    double last_accepted_timestamp_ = 0.0;  ///< 严格递增门控最近接受的时间戳（s）。
    std::uint64_t rejected_count_ = 0;  ///< 因无效时间戳或不满足门控条件而拒绝的次数。
    double worst_rollback_sec_ = 0.0;  ///< 被拒输入相对上次接受时间戳的最大回退（s）。
};

/// Thread-safe gate for dropping sensor packets that arrive too far behind the
/// newest timestamp already observed across an input group.
class MaximumLagTimestampGate {
   public:
    void SetMaximumLag(double maximum_lag_sec) {
        std::lock_guard<std::mutex> lock(mutex_);
        maximum_lag_sec_ = maximum_lag_sec;
    }

    TimestampGateDecision Observe(double timestamp) {
        std::lock_guard<std::mutex> lock(mutex_);
        TimestampGateDecision decision;
        decision.timestamp = timestamp;
        decision.reference_timestamp = latest_timestamp_;
        decision.valid = std::isfinite(timestamp) && timestamp > 0.0;
        if (!decision.valid) {
            ++rejected_count_;
            return decision;
        }
        if (maximum_lag_sec_ > 0.0 && latest_timestamp_ > timestamp) {
            decision.lag_sec = latest_timestamp_ - timestamp;
            if (decision.lag_sec > maximum_lag_sec_) {
                ++rejected_count_;
                worst_lag_sec_ = std::max(worst_lag_sec_, decision.lag_sec);
                return decision;
            }
        }
        latest_timestamp_ = std::max(latest_timestamp_, timestamp);
        decision.accepted = true;
        return decision;
    }

    void Reset() {
        std::lock_guard<std::mutex> lock(mutex_);
        latest_timestamp_ = 0.0;
        rejected_count_ = 0;
        worst_lag_sec_ = 0.0;
    }

    std::uint64_t RejectedCount() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return rejected_count_;
    }

    double WorstLagSec() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return worst_lag_sec_;
    }

   private:
    mutable std::mutex mutex_;  ///< 保护时间戳水位与拒绝统计。
    double maximum_lag_sec_ = 0.0;  ///< 允许落后最新输入的时间（s），0 关闭滞后限制。
    double latest_timestamp_ = 0.0;  ///< 输入组中已观察到的最大有效时间戳（s）。
    std::uint64_t rejected_count_ = 0;  ///< 因无效时间戳或不满足门控条件而拒绝的次数。
    double worst_lag_sec_ = 0.0;  ///< 被拒输入相对最新时间戳的最大滞后（s）。
};

}  // namespace lightning
