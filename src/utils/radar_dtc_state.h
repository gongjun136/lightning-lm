#pragma once
#include <algorithm>
#include <chrono>
#include <mutex>
#include <set>
#include <array>
#include <livox_ros_driver2/msg/lidar_diagnostic_state.hpp>
#include <diagnostic_monitor_interfaces/msg/dtc_observation_array.hpp>

namespace lightning::functional_safety {
class RadarDtcState {
 public:
    using Clock = std::chrono::steady_clock;
    struct Status {
        bool all_failed = false, fault = false, all_pass = false;
        // Counterfactual evidence only: never use this to authorize motion until
        // a bounded timing-error policy has been defined and validated.
        bool fault_without_time_sync = false, all_pass_without_time_sync = false;
        uint8_t time_sync_fault_mask = 0;  // front, rear, left, right bits
        bool bounded_fault = false, bounded_all_pass = false;
        uint8_t bounded_sync_warning_mask = 0;
    };
    RadarDtcState(uint16_t source, double timeout, double hms_timeout, bool check_availability = false,
                  std::array<uint16_t, 4> ids = {{7, 12, 8, 9}})
        : source_(source), timeout_(timeout), hms_timeout_(hms_timeout),
          check_availability_(check_availability), ids_(ids) {}
    bool ObserveAvailability(const livox_ros_driver2::msg::LidarDiagnosticState& msg,
                             Clock::time_point now = Clock::now()) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (msg.role >= 4 || msg.comm_header.source_id != ids_[msg.role] ||
            msg.boot_id == 0 || msg.config_version.empty() || msg.usability > 2) return false;
        auto& role = availability_[msg.role];
        if (role.retired.count(msg.boot_id) || (role.seen && role.boot == msg.boot_id &&
            static_cast<int32_t>(msg.comm_header.seq - role.seq) <= 0)) return false;
        if (role.seen && role.boot != msg.boot_id) role.retired.insert(role.boot);
        role.seen = true; role.boot = msg.boot_id; role.seq = msg.comm_header.seq;
        role.usability = msg.usability; role.received = now;
        role.bounded_warning = msg.config_version == "livox-dtc-v3-bounded-sync" &&
            msg.bounded_time_sync_warning && msg.time_sync_offset_ns >= -1000000 &&
            msg.time_sync_offset_ns <= 1000000 && msg.time_sync_evidence_age_us <= 2500000 &&
            msg.observations.size() == 11 &&
            TimeSyncRoleMask(msg.observations[1].dtc) == (1u << msg.role) &&
            msg.observations[1].physical_condition == 2;
        role.proof_age_us = msg.time_sync_evidence_age_us;
        if (msg.observations.size() == 11)
            role.proof_age_us = std::max(role.proof_age_us, msg.observations[1].evaluation_age_us);
        return true;
    }
    bool Observe(const diagnostic_monitor_interfaces::msg::DtcObservationArray& msg,
                 Clock::time_point now = Clock::now()) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (msg.comm_header.source_id != source_ || msg.boot_id == 0 || msg.interface_version != 1 ||
            msg.config_version.empty() || msg.observations.size() != 46) return false;
        std::set<uint32_t> codes;
        for (const auto& obs : msg.observations)
            if (!Known(obs.dtc) || obs.physical_condition > 2 || !codes.insert(obs.dtc).second) return false;
        if (seen_ && msg.boot_id == state_.boot_id &&
            static_cast<int32_t>(msg.comm_header.seq - state_.comm_header.seq) <= 0) return false;
        if (retired_.count(msg.boot_id)) return false;
        if (seen_ && msg.boot_id != state_.boot_id) retired_.insert(state_.boot_id);
        for (const auto& obs : msg.observations) {
            if (obs.dtc == 0x10F00002 && obs.evaluation_age_us <= timeout_ * 1e6) {
                if (obs.physical_condition == 2) all_failed_ = true;
                if (obs.physical_condition == 1) all_failed_ = false;
            }
        }
        state_ = msg; received_ = now; seen_ = true;
        return true;
    }
    Status Evaluate(Clock::time_point now = Clock::now()) const {
        std::lock_guard<std::mutex> lock(mutex_);
        Status result;
        result.all_failed = all_failed_;
        const double age = std::chrono::duration<double>(now - received_).count();
        result.all_pass = seen_ && age <= timeout_;
        result.all_pass_without_time_sync = result.all_pass;
        result.bounded_all_pass = result.all_pass;
        if (result.all_pass) for (const auto& obs : state_.observations) {
            const double ttl = obs.dtc == 0x10F00001 || obs.dtc == 0x10F00002 ? timeout_ : hms_timeout_;
            const bool fresh = obs.evaluation_age_us <= ttl * 1e6 && obs.evaluation_age_us * 1e-6 + age <= ttl;
            // Recovered-link/synchronization-history DTCs remain visible to
            // diagnostics, but are not evidence of a current safety fault.
            // UNKNOWN/stale history still cannot stand in for valid evidence.
            const bool history = Historical(obs.dtc);
            const uint8_t sync_role = TimeSyncRoleMask(obs.dtc);
            const bool sync_fail = fresh && obs.physical_condition == 2 && sync_role != 0;
            if (sync_fail) result.time_sync_fault_mask |= sync_role;
            bool advisory = false;
            if (sync_fail && check_availability_) {
                for (size_t i = 0; i < availability_.size(); ++i) if (sync_role == (1u << i)) {
                    const auto& device = availability_[i];
                    const double arrival_age = std::chrono::duration<double>(now - device.received).count();
                    advisory = device.seen && arrival_age <= timeout_ && device.usability == 1 &&
                        device.bounded_warning && device.proof_age_us * 1e-6 + arrival_age <= 2.5 &&
                        device.proof_age_us * 1e-6 + arrival_age <= hms_timeout_;
                }
            }
            if (advisory) result.bounded_sync_warning_mask |= sync_role;
            result.fault |= fresh && obs.physical_condition == 2 && !history;
            result.all_pass &= fresh && (obs.physical_condition == 1 || (history && obs.physical_condition == 2));
            result.fault_without_time_sync |= fresh && obs.physical_condition == 2 && !history && !sync_role;
            // UNKNOWN, stale evidence and non-sync faults stay blocking even in shadow.
            result.all_pass_without_time_sync &= fresh &&
                (obs.physical_condition == 1 || (history && obs.physical_condition == 2) || sync_fail);
            result.bounded_fault |= fresh && obs.physical_condition == 2 && !history && !advisory;
            result.bounded_all_pass &= fresh &&
                (obs.physical_condition == 1 || (history && obs.physical_condition == 2) || advisory);
        }
        if (check_availability_) {
            bool all_unusable = true;
            for (const auto& role : availability_) {
                const bool fresh = role.seen && std::chrono::duration<double>(now - role.received).count() <= timeout_;
                result.all_pass &= fresh && role.usability == 1;
                result.fault |= fresh && role.usability == 2;
                result.all_pass_without_time_sync &= fresh && role.usability == 1;
                result.fault_without_time_sync |= fresh && role.usability == 2;
                result.bounded_all_pass &= fresh && role.usability == 1;
                result.bounded_fault |= fresh && role.usability == 2;
                all_unusable &= fresh && role.usability == 2;
            }
            result.all_failed |= all_unusable;
        }
        return result;
    }
    static bool Known(uint32_t code) {
        if (code == 0x10F00001 || code == 0x10F00002) return true;
        for (uint32_t base : {0x110476C4u, 0x11037660u, 0x11017598u, 0x110275FCu})
            for (uint32_t offset : {0u, 1u, 3u, 4u, 7u, 8u, 9u, 10u, 11u, 12u, 13u})
                if (code == base + offset) return true;
        return false;
    }
    static bool Historical(uint32_t code) {
        for (uint32_t base : {0x110476C4u, 0x11037660u, 0x11017598u, 0x110275FCu})
            if (code == base + 12 || code == base + 13) return true;
        return false;
    }
    static uint8_t TimeSyncRoleMask(uint32_t code) {
        const std::array<uint32_t, 4> codes{{0x110476C5u, 0x11037661u, 0x11017599u, 0x110275FDu}};
        for (size_t i = 0; i < codes.size(); ++i) if (code == codes[i]) return uint8_t(1u << i);
        return 0;
    }
 private:
    uint16_t source_;
    double timeout_, hms_timeout_;
    mutable std::mutex mutex_;
    diagnostic_monitor_interfaces::msg::DtcObservationArray state_;
    Clock::time_point received_;
    std::set<uint64_t> retired_;
    bool seen_ = false, all_failed_ = false;
    struct Availability {
        bool seen = false;
        uint64_t boot = 0;
        uint32_t seq = 0;
        uint8_t usability = 0;
        bool bounded_warning = false;
        uint64_t proof_age_us = 0;
        Clock::time_point received;
        std::set<uint64_t> retired;
    };
    bool check_availability_;
    std::array<uint16_t, 4> ids_;
    std::array<Availability, 4> availability_;
};
}
