#pragma once
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
    struct Status { bool all_failed = false, fault = false, all_pass = false; };
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
        if (result.all_pass) for (const auto& obs : state_.observations) {
            const double ttl = obs.dtc == 0x10F00001 || obs.dtc == 0x10F00002 ? timeout_ : hms_timeout_;
            const bool fresh = obs.evaluation_age_us <= ttl * 1e6 && obs.evaluation_age_us * 1e-6 + age <= ttl;
            result.fault |= fresh && obs.physical_condition == 2;
            result.all_pass &= fresh && obs.physical_condition == 1;
        }
        if (check_availability_) {
            bool all_unusable = true;
            for (const auto& role : availability_) {
                const bool fresh = role.seen && std::chrono::duration<double>(now - role.received).count() <= timeout_;
                result.all_pass &= fresh && role.usability == 1;
                result.fault |= fresh && role.usability == 2;
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
        Clock::time_point received;
        std::set<uint64_t> retired;
    };
    bool check_availability_;
    std::array<uint16_t, 4> ids_;
    std::array<Availability, 4> availability_;
};
}
