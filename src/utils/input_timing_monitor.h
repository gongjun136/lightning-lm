#pragma once
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>

namespace lightning::diagnostics {
// 调用方复用输入统计锁；每条消息只更新固定内存，不格式化日志、不启动线程。
class InputTimingMonitor {
 public:
    using Clock = std::chrono::steady_clock;
    enum class Kind { kStart, kInterval, kSilence, kSummary };
    struct Record {
        Kind kind = Kind::kStart;
        std::uint64_t sample = 0, event = 0, lost = 0;
        double stamp = 0, previous_stamp = 0, source_gap = -1, arrival_gap = -1;
        double arrival_steady = 0, threshold = 0;
        double max_source_gap = 0, max_arrival_gap = 0;
        std::uint64_t source_gaps = 0, arrival_gaps = 0, duplicates = 0, rollbacks = 0, invalid = 0;
    };
    explicit InputTimingMonitor(double threshold = 0.1) : threshold_(threshold) {}
    void Observe(double stamp, Clock::time_point now) {
        Record r;
        r.kind = seen_ ? Kind::kInterval : Kind::kStart;
        r.sample = ++samples_; r.stamp = stamp; r.previous_stamp = last_stamp_;
        r.threshold = threshold_; r.arrival_steady = Seconds(now.time_since_epoch());
        const bool valid = std::isfinite(stamp) && stamp > 0;
        if (!valid) ++invalid_;
        if (seen_) {
            r.arrival_gap = Seconds(now - last_arrival_);
            if (valid && std::isfinite(last_stamp_) && last_stamp_ > 0) {
                r.source_gap = stamp - last_stamp_;
                if (r.source_gap == 0) ++duplicates_;
                if (r.source_gap < 0) ++rollbacks_;
                if (r.source_gap > threshold_) ++source_gaps_;
                max_source_gap_ = std::max(max_source_gap_, r.source_gap);
            }
            if (r.arrival_gap > threshold_) ++arrival_gaps_;
            max_arrival_gap_ = std::max(max_arrival_gap_, r.arrival_gap);
        }
        // 只缓存异常相邻样本；源时间差和到达时间差始终属于同一对消息。
        if (!seen_ || !valid || r.source_gap <= 0 ||
            r.source_gap > threshold_ || r.arrival_gap > threshold_) Push(r);
        seen_ = true; silence_reported_ = false;
        previous_stamp_ = last_stamp_; last_source_gap_ = r.source_gap; last_arrival_gap_ = r.arrival_gap;
        last_stamp_ = stamp; last_arrival_ = now;
    }
    void Poll(Clock::time_point now) {
        if (!seen_ || silence_reported_ || Seconds(now - last_arrival_) <= threshold_) return;
        Record r;
        r.kind = Kind::kSilence; r.sample = samples_; r.stamp = last_stamp_;
        r.arrival_gap = Seconds(now - last_arrival_);
        r.arrival_steady = Seconds(now.time_since_epoch()); r.threshold = threshold_;
        Push(r); silence_reported_ = true;
    }
    bool Pop(Record& r) {
        if (!size_) return false;
        r = pending_[head_]; r.lost = lost_;
        head_ = (head_ + 1) % pending_.size(); --size_; return true;
    }
    Record Summary() const {
        Record r;
        r.kind = Kind::kSummary; r.sample = samples_; r.stamp = last_stamp_;
        r.event = event_sequence_;
        r.previous_stamp = previous_stamp_; r.source_gap = last_source_gap_; r.arrival_gap = last_arrival_gap_;
        if (seen_) r.arrival_steady = Seconds(last_arrival_.time_since_epoch());
        r.threshold = threshold_; r.lost = lost_;
        r.max_source_gap = max_source_gap_; r.max_arrival_gap = max_arrival_gap_;
        r.source_gaps = source_gaps_; r.arrival_gaps = arrival_gaps_;
        r.duplicates = duplicates_; r.rollbacks = rollbacks_; r.invalid = invalid_;
        return r;
    }
    static const char* Name(Kind kind) {
        switch (kind) {
            case Kind::kStart: return "start";
            case Kind::kInterval: return "interval";
            case Kind::kSilence: return "silence";
            default: return "summary";
        }
    }
 private:
    template<class Duration> static double Seconds(Duration d) { return std::chrono::duration<double>(d).count(); }
    void Push(const Record& r) {
        ++event_sequence_;
        if (size_ == pending_.size()) { ++lost_; return; }
        auto& target = pending_[(head_ + size_) % pending_.size()];
        target = r; target.event = event_sequence_; ++size_;
    }
    double threshold_, last_stamp_ = 0, max_source_gap_ = 0, max_arrival_gap_ = 0;
    double previous_stamp_ = 0, last_source_gap_ = -1, last_arrival_gap_ = -1;
    Clock::time_point last_arrival_;
    bool seen_ = false, silence_reported_ = false;
    std::uint64_t samples_ = 0, event_sequence_ = 0, source_gaps_ = 0, arrival_gaps_ = 0;
    std::uint64_t duplicates_ = 0, rollbacks_ = 0, invalid_ = 0, lost_ = 0;
    std::array<Record, 16> pending_{};
    std::size_t head_ = 0, size_ = 0;
};
}  // namespace lightning::diagnostics
