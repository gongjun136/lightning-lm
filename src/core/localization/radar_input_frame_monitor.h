#pragma once
#include <iterator>
#include <iomanip>
#include <ostream>
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
    struct InputEvidence {
        int id = 0;
        std::string topic;
        std::uint64_t samples = 0, window_samples = 0;
        double stamp = 0, callback_ros_sec = 0;
        double lag_sec = NAN, window_max_lag_sec = NAN;
        double source_gap_sec = NAN, arrival_gap_sec = NAN, window_max_arrival_gap_sec = NAN;
        double callback_age_sec = NAN;
        bool nonempty = false;
        const char* stamp_state = "unseen";
        const char* assembly = "unseen";
        Clock::time_point arrival;
    };
    struct Evidence {
        bool available = false, consistent = false, pending = false;
        double steady_sec = 0, after_safety_sec = 0, sensor_progress_sec = 0;
        std::uint64_t accepted = 0, safety_accepted = 0;
        MultiLidarFrameStats oldest;
        std::vector<InputEvidence> inputs;
        const char* group_result = "unseen";
        double group_age_sec = NAN;
    };
    RadarInputFrameMonitor(MultiLidarConfig config, Watchdog::Config timing)
        : config_(std::move(config)), max_input_age_sec_(timing.timeout_sec), watchdog_(timing) {
        assembler_.Reset(config_, MultiLidarFrameAssembler::Mode::kInputStamps);
        if (config_.enabled) {
            for (const auto& sensor : config_.lidars) {
                inputs_[sensor.id].id = sensor.id;
                inputs_[sensor.id].topic = sensor.lidar_topic;
            }
        } else {
            inputs_[0].id = 0;
        }
    }
    // 单雷达入口也保留实际订阅话题；只在创建订阅前调用。
    void SetInputTopic(int id, const std::string& topic) {
        std::lock_guard<std::mutex> lock(input_mutex_);
        if (inputs_.count(id)) inputs_.at(id).topic = topic;
    }
    void Start(Clock::time_point now = Clock::now()) { watchdog_.Start(now); }
    Watchdog::Snapshot Evaluate(Clock::time_point now = Clock::now()) { return watchdog_.Evaluate(now); }
    void ObserveInput(int id, double stamp, bool nonempty, Clock::time_point now = Clock::now(),
                      double callback_ros_sec = NAN) {
        std::lock_guard<std::mutex> lock(input_mutex_);
        auto input = inputs_.find(id);
        if (input == inputs_.end()) return;
        auto& r = input->second;
        r.source_gap_sec = r.samples && std::isfinite(stamp) && stamp > 0 && r.stamp > 0 ? stamp - r.stamp : NAN;
        r.arrival_gap_sec = r.samples ? Seconds(now - r.arrival) : NAN;
        r.stamp_state = !std::isfinite(stamp) || stamp <= 0 ? "invalid" :
            r.source_gap_sec == 0 ? "duplicate" : r.source_gap_sec < 0 ? "rollback" : "advancing";
        r.stamp = stamp; r.callback_ros_sec = callback_ros_sec; r.arrival = now; r.nonempty = nonempty;
        r.lag_sec = std::isfinite(stamp) && stamp > 0 && std::isfinite(callback_ros_sec) && callback_ros_sec > 0 ?
            callback_ros_sec - stamp : NAN;
        UpdateMax(r.window_max_lag_sec, r.lag_sec);
        UpdateMax(r.window_max_arrival_gap_sec, r.arrival_gap_sec);
        ++r.samples; ++r.window_samples;
        if (!nonempty) { r.assembly = "empty"; return; }
        if (!std::isfinite(stamp) || stamp <= 0) { r.assembly = "invalid_stamp"; return; }
        if (!config_.enabled) {
            const bool accepted = watchdog_.Observe(stamp, now, true, now);
            if (accepted) ++accepted_groups_;
            r.assembly = group_result_ = accepted ? "accepted" : "watchdog_rejected";
            group_time_ = now; group_seen_ = true;
            return;
        }
        if (id == config_.primary_lidar_id) {
            primary_arrivals_.emplace(stamp, now);  // Keep first arrival for duplicates.
            while (primary_arrivals_.size() > 256) primary_arrivals_.erase(primary_arrivals_.begin());
        }
        const auto late = assembler_.LateDropCount(), duplicate = assembler_.DuplicateDropCount();
        const auto tolerance = assembler_.ToleranceDropCount();
        if (!assembler_.AddStamp(id, stamp)) {
            r.assembly = assembler_.LateDropCount() != late ? "late_bucket" :
                assembler_.DuplicateDropCount() != duplicate ? "duplicate_bucket" :
                assembler_.ToleranceDropCount() != tolerance ? "tolerance" : "invalid";
            return;
        }
        r.assembly = "accepted";  // 仅表示进入组帧器，不声称已形成合格组帧。
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
            const bool accepted = watchdog_.Observe(stats.begin_time, now, has_primary && fresh_primary, now);
            if (accepted) ++accepted_groups_;
            group_result_ = accepted ? "accepted" : !has_primary ? "missing_primary" :
                !fresh_primary ? "stale_primary" : "watchdog_rejected";
            group_time_ = now; group_seen_ = true;
            if (has_primary && arrival != primary_arrivals_.end() &&
                arrival->first <= stats.begin_time + config_.match_tolerance + 5e-7)
                primary_arrivals_.erase(primary_arrivals_.begin(), std::next(arrival));
        }
    }
    Evidence SnapshotEvidence(std::uint64_t safety_accepted, Clock::time_point evaluated,
                              Clock::time_point now = Clock::now()) {
        Evidence e;
        e.safety_accepted = safety_accepted;
        e.steady_sec = Seconds(now.time_since_epoch()); e.after_safety_sec = Seconds(now - evaluated);
        // 诊断不能等待输入组帧锁；忙时明确缺证，不拖住 Safety 发布。
        std::unique_lock<std::mutex> lock(input_mutex_, std::try_to_lock);
        if (!lock.owns_lock()) return e;
        e.available = true;
        e.accepted = accepted_groups_;
        e.consistent = safety_accepted == e.accepted;
        e.pending = assembler_.DescribeOldestPending(e.oldest, e.sensor_progress_sec);
        e.group_result = group_result_;
        if (group_seen_) e.group_age_sec = std::max(0.0, Seconds(now - group_time_));
        for (auto& [id, r] : inputs_) {
            auto copy = r;
            if (r.samples) copy.callback_age_sec = std::max(0.0, Seconds(now - r.arrival));
            e.inputs.push_back(std::move(copy));
            r.window_samples = 0; r.window_max_lag_sec = NAN; r.window_max_arrival_gap_sec = NAN;
        }
        return e;
    }
    void WriteEvidence(std::ostream& out, const Evidence& e) const {
        const auto precision = out.precision();
        out << std::setprecision(17) << " radar_rx_schema=1 radar_rx_available=" << e.available
            << " radar_rx_consistent=" << e.consistent << " radar_rx_steady_sec=" << e.steady_sec
            << " radar_rx_after_safety_sec=" << e.after_safety_sec
            << " radar_rx_safety_accepted=" << e.safety_accepted << " radar_rx_snapshot_accepted=" << e.accepted;
        if (!e.available) { out.precision(precision); return; }
        out << " radar_rx_clock=ros_header_unverified radar_rx_primary=" << (config_.enabled ? config_.primary_lidar_id : 0)
            << " radar_rx_timeout_sec=" << max_input_age_sec_
            << " radar_rx_min_lidars=" << (config_.enabled ? (config_.adaptive_load.enabled ?
                config_.adaptive_load.tracking_min_lidars : config_.min_lidars) : 1)
            << " radar_rx_tolerance_sec=" << config_.match_tolerance
            << " radar_rx_reorder_sec=" << config_.reorder_window
            << " radar_rx_pending=" << e.pending << " radar_rx_pending_progress_sec=" << e.sensor_progress_sec
            << " radar_rx_pending_stamp_sec=" << e.oldest.begin_time
            << " radar_rx_group_result=" << e.group_result << " radar_rx_group_age_sec=" << e.group_age_sec;
        auto ids = [&](const auto& list) {
            if (list.empty()) out << "none";
            for (std::size_t i = 0; i < list.size(); ++i) out << (i ? "|" : "") << list[i];
        };
        out << " radar_rx_pending_present="; ids(e.oldest.present_lidar_ids);
        out << " radar_rx_pending_missing="; ids(e.oldest.missing_lidar_ids);
        std::vector<int> configured;
        for (const auto& r : e.inputs) configured.push_back(r.id);
        out << " radar_rx_ids="; ids(configured);
        for (const auto& r : e.inputs) {
            const std::string p = " radar_rx_" + std::to_string(r.id) + "_";
            out << p << "topic=" << (r.topic.empty() ? "unknown" : r.topic)
                << p << "samples=" << r.samples << p << "window_samples=" << r.window_samples
                << p << "stamp_sec=" << r.stamp << p << "callback_ros_sec=" << r.callback_ros_sec
                << p << "lag_sec=" << r.lag_sec << p << "max_lag_sec=" << r.window_max_lag_sec
                << p << "age_sec=" << r.callback_age_sec << p << "source_gap_sec=" << r.source_gap_sec
                << p << "arrival_gap_sec=" << r.arrival_gap_sec
                << p << "max_arrival_gap_sec=" << r.window_max_arrival_gap_sec
                << p << "nonempty=" << r.nonempty << p << "stamp_state=" << r.stamp_state
                << p << "assembly=" << r.assembly;
        }
        out.precision(precision);
    }
 private:
    template <class Duration> static double Seconds(Duration d) { return std::chrono::duration<double>(d).count(); }
    static void UpdateMax(double& value, double sample) {
        if (std::isfinite(sample) && (!std::isfinite(value) || sample > value)) value = sample;
    }
    const MultiLidarConfig config_;
    const double max_input_age_sec_;
    std::mutex input_mutex_;  // Never shared with an algorithm/processing lock.
    MultiLidarFrameAssembler assembler_;
    std::map<double, Clock::time_point> primary_arrivals_;
    Watchdog watchdog_;
    std::map<int, InputEvidence> inputs_;
    const char* group_result_ = "unseen";
    Clock::time_point group_time_;
    bool group_seen_ = false;
    std::uint64_t accepted_groups_ = 0;
};
}
