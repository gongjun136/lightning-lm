#pragma once
#include <algorithm>
#include <array>
#include <chrono>
#include <mutex>
#include <set>
#include <stdexcept>
#include <livox_ros_driver2/msg/lidar_diagnostic_state.hpp>

namespace lightning::functional_safety {
/// Safety consumes the configured devices directly, preserving all four role IDs.
class RadarDtcState {
 public:
    using Clock = std::chrono::steady_clock;
    struct Status {
        bool all_failed = false, fault = false, all_pass = false;
        bool fault_without_time_sync = false, all_pass_without_time_sync = false;
        uint8_t time_sync_fault_mask = 0;  // front, rear, left, right bits
        bool bounded_fault = false, bounded_all_pass = false;
        uint8_t bounded_sync_warning_mask = 0;
    };
    RadarDtcState(double timeout, double hms_timeout,
                  std::array<uint16_t, 4> ids = {{7, 12, 8, 9}}, double startup_grace = 8.0,
                  Clock::time_point start = Clock::now(),
                  std::array<bool, 4> required = {{true, true, true, true}})
        : timeout_(timeout), hms_timeout_(hms_timeout), ids_(ids),
          startup_grace_(startup_grace), start_(start), required_(required) {
        if (std::none_of(required_.begin(), required_.end(), [](bool enabled) { return enabled; }))
            throw std::invalid_argument("at least one radar must be required");
    }

    bool Observe(const livox_ros_driver2::msg::LidarDiagnosticState& msg,
                 Clock::time_point now = Clock::now()) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (msg.role >= 4 || !required_[msg.role] || msg.comm_header.source_id != ids_[msg.role] ||
            msg.boot_id == 0 || msg.config_version != "livox-dtc-v3-bounded-sync" ||
            msg.usability > 2 || msg.observations.size() != offsets_.size()) return false;
        for (size_t i = 0; i < offsets_.size(); ++i)
            if (msg.observations[i].dtc != bases_[msg.role] + offsets_[i] ||
                msg.observations[i].physical_condition > 2) return false;
        auto& role = sources_[msg.role];
        if (role.retired.count(msg.boot_id) || (role.seen && role.state.boot_id == msg.boot_id &&
            static_cast<int32_t>(msg.comm_header.seq - role.state.comm_header.seq) <= 0)) return false;
        if (role.seen && role.state.boot_id != msg.boot_id) role.retired.insert(role.state.boot_id);
        role.seen = true;
        role.state = msg;
        role.received = now;
        return true;
    }

    Status Evaluate(Clock::time_point now = Clock::now()) const {
        std::lock_guard<std::mutex> lock(mutex_);
        Status result;
        result.all_pass = result.all_pass_without_time_sync = result.bounded_all_pass = true;
        size_t unusable = 0, usable = 0, required_count = 0;
        for (size_t i = 0; i < sources_.size(); ++i) {
            if (!required_[i]) continue;
            ++required_count;
            const auto& source = sources_[i];
            const double age = std::chrono::duration<double>(now - source.received).count();
            const bool fresh = source.seen && age <= timeout_;
            const uint8_t usability = fresh ? source.state.usability :
                (source.seen || std::chrono::duration<double>(now - start_).count() >= startup_grace_ ? 2 : 0);
            unusable += usability == 2;
            usable += usability == 1;
            result.all_pass &= fresh && usability == 1;
            result.all_pass_without_time_sync &= fresh && usability == 1;
            result.bounded_all_pass &= fresh && usability == 1;
            result.fault |= usability == 2;
            result.fault_without_time_sync |= usability == 2;
            result.bounded_fault |= usability == 2;
            if (!fresh) continue;
            for (size_t category = 0; category < offsets_.size(); ++category) {
                const auto& obs = source.state.observations[category];
                const bool evidence_fresh = obs.evaluation_age_us <= hms_timeout_ * 1e6 &&
                    obs.evaluation_age_us * 1e-6 + age <= hms_timeout_;
                const bool history = category >= offsets_.size() - 2;
                const bool sync_fail = category == 1 && evidence_fresh && obs.physical_condition == 2;
                const auto bit = uint8_t(1u << i);
                if (sync_fail) result.time_sync_fault_mask |= bit;
                const auto& state = source.state;
                const auto proof_age = std::max(state.time_sync_evidence_age_us, obs.evaluation_age_us);
                const bool advisory = sync_fail && usability == 1 && state.bounded_time_sync_warning &&
                    state.time_sync_offset_ns >= -1000000 && state.time_sync_offset_ns <= 1000000 &&
                    proof_age <= 2500000 && proof_age * 1e-6 + age <= 2.5 &&
                    proof_age * 1e-6 + age <= hms_timeout_;
                if (advisory) result.bounded_sync_warning_mask |= bit;
                const bool pass = obs.physical_condition == 1 || (history && obs.physical_condition == 2);
                result.fault |= evidence_fresh && obs.physical_condition == 2 && !history;
                result.all_pass &= evidence_fresh && pass;
                result.fault_without_time_sync |= evidence_fresh && obs.physical_condition == 2 && !history && category != 1;
                result.all_pass_without_time_sync &= evidence_fresh && (pass || sync_fail);
                result.bounded_fault |= evidence_fresh && obs.physical_condition == 2 && !history && !advisory;
                result.bounded_all_pass &= evidence_fresh && (pass || advisory);
            }
        }
        // A restart/UNKNOWN cannot recover a previously established all-failed
        // state. At least one fresh usable device proves that it is no longer all.
        if (unusable == required_count) all_failed_ = true;
        if (usable != 0) all_failed_ = false;
        result.all_failed = all_failed_;
        return result;
    }

 private:
    inline static constexpr std::array<uint32_t, 4> bases_{{0x110476C4u, 0x11037660u, 0x11017598u, 0x110275FCu}};
    inline static constexpr std::array<uint32_t, 11> offsets_{{0, 1, 3, 4, 7, 8, 9, 10, 11, 12, 13}};
    struct Source {
        livox_ros_driver2::msg::LidarDiagnosticState state;
        Clock::time_point received;
        std::set<uint64_t> retired;
        bool seen = false;
    };
    double timeout_, hms_timeout_;  ///< Receiver silence and source evidence limits, in seconds.
    std::array<uint16_t, 4> ids_;  ///< Registered front/rear/left/right source IDs.
    double startup_grace_;  ///< Missing devices become unavailable after this many seconds.
    Clock::time_point start_;  ///< Steady-clock start of the receiver's grace period.
    const std::array<bool, 4> required_;  ///< Immutable front/rear/left/right deployment topology.
    mutable std::mutex mutex_;  ///< Protects source snapshots and the severe-state latch.
    mutable bool all_failed_ = false;  ///< Held across UNKNOWN until a fresh usable device arrives.
    std::array<Source, 4> sources_;  ///< Latest validated per-device evidence and retired boots.
};
}
