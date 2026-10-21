#include "utils/radar_dtc_state.h"
#include "utils/sany_radar_topology.h"
#include <iostream>
#include <stdexcept>

using lightning::functional_safety::RadarDtcState;
using Message = livox_ros_driver2::msg::LidarDiagnosticState;
namespace {
void Require(bool value, const char* reason) { if (!value) throw std::runtime_error(reason); }
auto At(int ms) { return RadarDtcState::Clock::time_point{} + std::chrono::milliseconds(ms); }
constexpr std::array<uint16_t, 4> ids{{7, 12, 8, 9}};
Message Healthy(uint8_t role) {
    constexpr std::array<uint32_t, 4> bases{{0x110476C4, 0x11037660, 0x11017598, 0x110275FC}};
    Message msg;
    msg.role = role; msg.comm_header.source_id = ids[role]; msg.comm_header.seq = 1;
    msg.boot_id = 1; msg.config_version = "livox-dtc-v3-bounded-sync"; msg.usability = 1;
    for (uint32_t offset : {0u, 1u, 3u, 4u, 7u, 8u, 9u, 10u, 11u, 12u, 13u}) {
        diagnostic_monitor_interfaces::msg::DtcObservation obs;
        obs.dtc = bases[role] + offset; obs.physical_condition = 1;
        msg.observations.push_back(obs);
    }
    return msg;
}
}
int main() {
    using lightning::functional_safety::RequiredSanyRadarRoles;
    const std::array<bool, 4> three_roles{{true, false, true, true}};
    auto config = YAML::Load(R"(
vehicle_configuration: {LIDAR: MASTER AND SLAVE, lidar_count: 3}
multi_lidar:
  enabled: true
  primary_lidar_id: 0
  topics: {lidar_0: /front, imu_0: /imu0, lidar_1: /left, imu_1: /imu1, lidar_2: /right, imu_2: /imu2}
)");
    Require(RequiredSanyRadarRoles(config) == three_roles, "three-lidar snapshot preserves physical DTC roles");
    auto four = YAML::Clone(config);
    four["vehicle_configuration"]["LIDAR"] = "MASTER";
    four["vehicle_configuration"]["lidar_count"] = 4;
    four["multi_lidar"]["topics"]["lidar_3"] = "/rear";
    four["multi_lidar"]["topics"]["imu_3"] = "/imu3";
    Require(RequiredSanyRadarRoles(four) == std::array<bool, 4>{{true, true, true, true}},
            "four-lidar snapshot requires rear");
    Require(RequiredSanyRadarRoles(YAML::Load("{}")) == RequiredSanyRadarRoles(four),
            "legacy non-SANY configuration retains four roles");
    for (int failure = 0; failure < 7; ++failure) {
        auto broken = YAML::Clone(config);
        switch (failure) {
            case 0: broken["vehicle_configuration"]["lidar_count"] = 4; break;
            case 1: broken["vehicle_configuration"]["LIDAR"] = "SLAVE"; break;
            case 2: broken["multi_lidar"]["enabled"] = false; break;
            case 3: broken["multi_lidar"]["topics"].remove("lidar_2"); break;
            case 4: broken["multi_lidar"]["topics"]["lidar_3"] = "/extra"; break;
            case 5: broken["sany_vehicle_config"]["calibrations"] = "unresolved"; break;
            case 6: broken["multi_lidar"]["topics"]["lidar_2"] = "/left"; break;
        }
        bool rejected = false;
        try { RequiredSanyRadarRoles(broken); } catch (const std::exception&) { rejected = true; }
        Require(rejected, "inconsistent topology must reject startup rather than omit a required radar");
    }
    bool empty_rejected = false;
    try { RadarDtcState empty(0.5, 2.5, ids, 8.0, At(0), {{false, false, false, false}}); }
    catch (const std::invalid_argument&) { empty_rejected = true; }
    Require(empty_rejected, "empty required set cannot fabricate healthy Safety");

    RadarDtcState three(0.5, 2.5, ids, 8.0, At(0), RequiredSanyRadarRoles(config));
    Require(!three.Evaluate(At(7999)).all_failed && three.Evaluate(At(8000)).all_failed,
            "three missing required devices retain startup grace and severe timeout");
    for (uint8_t role : {0, 2, 3}) Require(three.Observe(Healthy(role), At(8010)), "three active roles accepted");
    Require(three.Evaluate(At(8010)).all_pass && !three.Evaluate(At(8010)).fault,
            "healthy front left right pass without rear");
    auto rear = Healthy(1); rear.usability = 2;
    Require(!three.Observe(rear, At(8010)) && three.Evaluate(At(8010)).all_pass,
            "inactive rear cannot introduce fault");
    auto left = Healthy(2); left.comm_header.seq = 2; left.usability = 2;
    Require(three.Observe(left, At(8020)) && three.Evaluate(At(8020)).fault &&
            !three.Evaluate(At(8020)).all_failed, "one active radar fault remains blocking");
    for (uint8_t role : {0, 2, 3}) {
        auto failed = Healthy(role); failed.comm_header.seq = 3; failed.usability = 2;
        Require(three.Observe(failed, At(8030)), "three failed sources accepted");
    }
    Require(three.Evaluate(At(8030)).all_failed, "all configured radars failed remains severe");
    Require(!three.Observe(Healthy(1), At(8040)) && three.Evaluate(At(8040)).all_failed,
            "inactive rear cannot clear severe state");
    auto restarted = Healthy(0); restarted.boot_id = 2; restarted.usability = 0;
    Require(three.Observe(restarted, At(8050)) && three.Evaluate(At(8050)).all_failed,
            "active restart UNKNOWN cannot clear severe state");
    restarted.usability = 1; ++restarted.comm_header.seq;
    Require(three.Observe(restarted, At(8060)) && !three.Evaluate(At(8060)).all_failed &&
            !three.Evaluate(At(8060)).all_pass, "active usable recovery disproves all failed only");
    for (uint8_t role : {0, 2, 3}) {
        auto healthy = Healthy(role); healthy.boot_id = 2; healthy.comm_header.seq = 4;
        Require(three.Observe(healthy, At(8070)), "three recovered sources accepted");
    }
    Require(three.Evaluate(At(8070)).all_pass && three.Evaluate(At(8571)).all_failed,
            "three-source outage still becomes severe");

    RadarDtcState startup(0.5, 2.5, ids, 8.0, At(0));
    Require(!startup.Evaluate(At(0)).all_pass && !startup.Evaluate(At(7999)).all_failed,
            "startup missing evidence is UNKNOWN during grace");
    Require(startup.Evaluate(At(8000)).all_failed, "missing four devices after grace remains severe");
    for (uint8_t role = 0; role < 4; ++role) Require(startup.Observe(Healthy(role), At(8010)), "fresh role after startup");
    Require(startup.Evaluate(At(8010)).all_pass && !startup.Evaluate(At(8010)).all_failed, "fresh recovery after startup timeout");

    RadarDtcState radar(0.5, 2.5, ids, 8.0, At(0));
    std::array<Message, 4> states;
    for (uint8_t role = 0; role < 4; ++role) {
        states[role] = Healthy(role);
        Require(radar.Observe(states[role], At(0)), "four complete source snapshots accepted");
    }
    Require(radar.Evaluate(At(0)).all_pass && !radar.Evaluate(At(0)).fault, "all healthy without aggregator");
    Require(!radar.Observe(states[0], At(300)), "replay must not refresh source age");
    Require(radar.Evaluate(At(501)).all_failed, "four timed-out sources count as unavailable");
    for (auto& state : states) {
        state.boot_id = 2; state.comm_header.seq = 0; state.usability = 0;
        for (auto& obs : state.observations) obs.physical_condition = 0;
        Require(radar.Observe(state, At(510)), "source restart accepted");
    }
    Require(radar.Evaluate(At(510)).all_failed, "restart UNKNOWN cannot clear previous severe state");
    auto retired = Healthy(0); retired.comm_header.seq = 2;
    Require(!radar.Observe(retired, At(520)), "retired boot rejected");
    states[0] = Healthy(0); states[0].boot_id = 2;
    Require(radar.Observe(states[0], At(530)) && !radar.Evaluate(At(530)).all_failed &&
            !radar.Evaluate(At(530)).all_pass, "one fresh usable device disproves all-failed but not all-healthy");
    for (uint8_t role = 0; role < 4; ++role) {
        states[role] = Healthy(role); states[role].boot_id = 2; states[role].comm_header.seq = 2;
        Require(radar.Observe(states[role], At(540)), "all sources recover");
    }
    Require(radar.Evaluate(At(540)).all_pass, "fresh complete recovery");
    auto invalid = states[0]; ++invalid.comm_header.seq;
    invalid.observations[0].dtc = invalid.observations[1].dtc;
    Require(!radar.Observe(invalid, At(550)), "duplicate/foreign code rejected");
    invalid = states[0]; ++invalid.comm_header.seq; invalid.comm_header.source_id = 8;
    Require(!radar.Observe(invalid, At(550)), "wrong source ID rejected");
    invalid = states[0]; ++invalid.comm_header.seq; invalid.observations.pop_back();
    Require(!radar.Observe(invalid, At(550)), "incomplete source snapshot rejected");
    invalid = states[0]; ++invalid.comm_header.seq; invalid.observations[0].physical_condition = 3;
    Require(!radar.Observe(invalid, At(550)), "invalid tri-state rejected");
    invalid = states[0]; ++invalid.comm_header.seq; invalid.config_version = "livox-dtc-v2-sync-source";
    Require(!radar.Observe(invalid, At(550)), "mixed driver semantics rejected");
    for (auto& state : states) {
        ++state.comm_header.seq;
        state.observations[9].physical_condition = state.observations[10].physical_condition = 2;
        Require(radar.Observe(state, At(560)), "history retained");
    }
    Require(radar.Evaluate(At(560)).all_pass && !radar.Evaluate(At(560)).fault, "history-only FAIL does not stop Safety");
    auto& front = states[0]; ++front.comm_header.seq; front.observations[9].physical_condition = 0;
    Require(radar.Observe(front, At(570)) && !radar.Evaluate(At(570)).all_pass, "unknown history cannot fabricate healthy evidence");
    ++front.comm_header.seq; front.observations[9].physical_condition = 1; front.usability = 2;
    Require(radar.Observe(front, At(580)) && radar.Evaluate(At(580)).fault &&
            !radar.Evaluate(At(580)).all_failed, "one unmapped native fatal still stops Safety");
    ++states[1].comm_header.seq; states[1].usability = 2;
    Require(radar.Observe(states[1], At(585)) && radar.Evaluate(At(585)).fault &&
            !radar.Evaluate(At(585)).all_failed, "two failures remain a radar fault without claiming all four failed");
    for (auto& state : states) {
        ++state.comm_header.seq; state.usability = 2;
        Require(radar.Observe(state, At(590)), "unusable devices");
    }
    Require(radar.Evaluate(At(590)).all_failed, "four physical failures remain severe");

    for (uint8_t role = 0; role < 4; ++role) {
        RadarDtcState audit(0.5, 2.5, ids, 8.0, At(0));
        for (uint8_t other = 0; other < 4; ++other) {
            auto state = Healthy(other);
            if (role == other) state.observations[1].physical_condition = 2;
            Require(audit.Observe(state, At(2000)), "sync fault fixture");
        }
        auto r = audit.Evaluate(At(2000));
        Require(r.time_sync_fault_mask == (1u << role) && r.fault && !r.all_pass &&
                r.all_pass_without_time_sync && !r.bounded_all_pass, "strict and shadow isolate only current sync fault");
        auto proof = Healthy(role); proof.comm_header.seq = 2;
        proof.observations[1].physical_condition = 2;
        proof.bounded_time_sync_warning = true; proof.time_sync_offset_ns = 1000000;
        Require(audit.Observe(proof, At(2000)), "independent bounded proof");
        r = audit.Evaluate(At(2000));
        Require(r.fault && !r.all_pass && !r.bounded_fault && r.bounded_all_pass &&
                r.bounded_sync_warning_mask == (1u << role), "bounded policy does not clear original DTC");
        ++proof.comm_header.seq; proof.time_sync_offset_ns = -1000001;
        Require(audit.Observe(proof, At(2000)) && !audit.Evaluate(At(2000)).bounded_all_pass, "negative bound enforced");
        ++proof.comm_header.seq; proof.time_sync_offset_ns = 1000001;
        Require(audit.Observe(proof, At(2000)) && !audit.Evaluate(At(2000)).bounded_all_pass, "positive bound enforced");
        ++proof.comm_header.seq; proof.time_sync_offset_ns = -1000000;
        Require(audit.Observe(proof, At(2000)) && audit.Evaluate(At(2000)).bounded_all_pass, "negative boundary inclusive");
        ++proof.comm_header.seq; proof.time_sync_offset_ns = 0; proof.time_sync_evidence_age_us = 2490000;
        Require(audit.Observe(proof, At(2000)) && !audit.Evaluate(At(2011)).bounded_all_pass, "transport residence included in proof age");
        ++proof.comm_header.seq; proof.time_sync_evidence_age_us = 0; proof.bounded_time_sync_warning = false;
        Require(audit.Observe(proof, At(2000)) && !audit.Evaluate(At(2000)).bounded_all_pass, "revoked proof blocks immediately");
        ++proof.comm_header.seq; proof.bounded_time_sync_warning = true; proof.observations[0].physical_condition = 2;
        Require(audit.Observe(proof, At(2010)) && audit.Evaluate(At(2010)).bounded_fault &&
                !audit.Evaluate(At(2010)).bounded_all_pass, "additional non-sync failure still blocks");
        ++proof.comm_header.seq; proof.observations[0].physical_condition = 0;
        Require(audit.Observe(proof, At(2020)) && !audit.Evaluate(At(2020)).all_pass_without_time_sync, "UNKNOWN non-sync stays blocking");
        ++proof.comm_header.seq; proof.observations[0].physical_condition = 1; proof.observations[1].physical_condition = 0;
        Require(audit.Observe(proof, At(2030)) && !audit.Evaluate(At(2030)).bounded_all_pass, "UNKNOWN sync cannot be waived");
        ++proof.comm_header.seq; proof.observations[1].physical_condition = 2; proof.observations[1].evaluation_age_us = 2500001;
        Require(audit.Observe(proof, At(2040)) && !audit.Evaluate(At(2040)).all_pass_without_time_sync, "stale HMS blocks even with fresh transport");
        Require(!audit.Evaluate(At(2541)).bounded_all_pass, "expired transport blocks bounded policy");
    }
    std::cout << "radar_dtc_state_test passed: three/four-source Safety, topology, startup/outage/restart, identity, history and bounded PTP\n";
}
