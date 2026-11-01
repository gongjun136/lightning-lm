#pragma once
#include <iterator>
#include <utility>
#include "core/lio/multi_lidar_fusion.h"
#include "utils/lidar_frame_watchdog.h"

namespace lightning::loc {
// Observe raw subscription input before preprocessing, throttling or algorithm
// queues. No algorithm progress, latency, quality or device DTC enters this check.
class RadarInputFrameMonitor {
 public:
    using Watchdog = functional_safety::LidarFrameWatchdog;
    using Clock = Watchdog::Clock;
    RadarInputFrameMonitor(MultiLidarConfig config, Watchdog::Config timing)
        : config_(std::move(config)), max_input_age_sec_(timing.timeout_sec), watchdog_(timing) {
        assembler_.Reset(config_, MultiLidarFrameAssembler::Mode::kInputStamps);
    }
    void Start(Clock::time_point now = Clock::now()) { watchdog_.Start(now); }
    Watchdog::Snapshot Evaluate(Clock::time_point now = Clock::now()) { return watchdog_.Evaluate(now); }
    void ObserveInput(int id, double stamp, bool nonempty, Clock::time_point now = Clock::now()) {
        if (!nonempty || !std::isfinite(stamp) || stamp <= 0) return;
        std::lock_guard<std::mutex> lock(input_mutex_);
        if (!config_.enabled) {
            if (id == 0) watchdog_.Observe(stamp, now, true, now);
            return;
        }
        if (id == config_.primary_lidar_id) {
            primary_arrivals_.emplace(stamp, now);  // Keep first arrival for duplicates.
            while (primary_arrivals_.size() > 256) primary_arrivals_.erase(primary_arrivals_.begin());
        }
        if (!assembler_.AddStamp(id, stamp)) return;
        FusedLidarFrame frame;
        while (assembler_.PopReady(frame)) {
            const auto& stats = frame.stats;
            const bool has_primary = std::find(stats.present_lidar_ids.begin(), stats.present_lidar_ids.end(),
                config_.primary_lidar_id) != stats.present_lidar_ids.end();
            const auto arrival = primary_arrivals_.lower_bound(stats.begin_time - 5e-7);
            const bool fresh_primary = arrival != primary_arrivals_.end() &&
                arrival->first <= stats.begin_time + config_.match_tolerance + 5e-7 &&
                arrival->second <= now &&
                std::chrono::duration<double>(now - arrival->second).count() < max_input_age_sec_;
            // Completion is at raw-input grouping, never LIO processing. An old
            // pending primary cannot be revived by a later auxiliary timestamp.
            watchdog_.Observe(stats.begin_time, now, has_primary && fresh_primary, now);
            if (has_primary && arrival != primary_arrivals_.end() &&
                arrival->first <= stats.begin_time + config_.match_tolerance + 5e-7)
                primary_arrivals_.erase(primary_arrivals_.begin(), std::next(arrival));
        }
    }
 private:
    const MultiLidarConfig config_;
    const double max_input_age_sec_;
    std::mutex input_mutex_;  // Never shared with an algorithm/processing lock.
    MultiLidarFrameAssembler assembler_;
    std::map<double, Clock::time_point> primary_arrivals_;
    Watchdog watchdog_;
};
}
