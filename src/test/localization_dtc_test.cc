#include "core/localization/localization_dtc_monitor.h"
#include "core/localization/pose_covariance.h"
#include "utils/radar_dtc_state.h"
#include <cstdlib>
#include <iostream>
#include <cmath>

using namespace lightning;
using namespace lightning::loc::dtc;
namespace {
void Require(bool value, const char* message) {
    if (!value) { std::cerr << message << std::endl; std::exit(1); }
}
auto At(int ms) { return Clock::time_point{} + std::chrono::milliseconds(ms); }
Mat6d SmallCov() { return 0.0001 * Mat6d::Identity(); }
SE3 Pose(double x, double y = 0, double yaw = 0) { return SE3(SO3::exp(Vec3d(0, 0, yaw)), Vec3d(x, y, 0)); }
void Feed(Monitor& monitor, int ms, const SE3& pose, double speed = 1, double gyro = 0,
          uint64_t epoch = 1, const SE3* raw = nullptr, const Mat6d* covariance = nullptr) {
    const double stamp = 100 + ms * 1e-3;
    monitor.ObserveSpeed(stamp, speed, At(ms));
    monitor.ObserveGyro(stamp, gyro, At(ms));
    monitor.ObservePose(stamp, pose, raw ? *raw : pose, covariance ? *covariance : SmallCov(), true, false, epoch, At(ms));
}
}
int main() {
    Config cfg;
    Monitor normal(cfg, At(0));
    Require(SafetyValue(normal.Evaluate(false, At(0)), false, false, true) == 4, "startup must not publish normal");
    for (int ms = 0; ms <= 1000; ms += 20) Feed(normal, ms, Pose(ms * 0.001));
    auto snapshot = normal.Evaluate(false, At(1000));
    Require(snapshot.observations[0].condition == kPass && snapshot.observations[1].condition == kPass, "straight motion and valid covariance pass");
    Require(SafetyValue(snapshot, false, false, true) == 0, "all confirmed healthy gives zero");
    Require(SafetyValue(snapshot, false, true, true) == 2 && SafetyValue(snapshot, false, false, false) == 2, "fault or stale radar gives two");
    snapshot.observations[0].condition = kFail;
    Require(SafetyValue(snapshot, false, true, false) == 4, "localization wins over radar");
    Require(SafetyValue(snapshot, true, true, false) == -9, "all-radar failure is severe");
    snapshot.severe = true;
    Require(SafetyValue(snapshot, false, false, true) == -9, "localization timeout is severe");
    Require(normal.Evaluate(false, At(1301)).severe, "real output timeout stops normal");
    Require(normal.Evaluate(true, At(1000)).severe, "output gate loss is severe immediately");

    Monitor reverse(cfg, At(0)), curve(cfg, At(0));
    for (int ms = 0; ms <= 1000; ms += 20) {
        const double t = ms * 0.001;
        Feed(reverse, ms, Pose(-t), -1);
        Feed(curve, ms, Pose(std::sin(0.2 * t) / 0.2, (1 - std::cos(0.2 * t)) / 0.2, 0.2 * t), 1, 0.2);
    }
    Require(reverse.Evaluate(false, At(1000)).observations[0].condition == kPass, "signed reverse wheel speed");
    Require(curve.Evaluate(false, At(1000)).observations[0].condition == kPass, "independent curved integration");
    Feed(normal, 1020, Pose(3.02));
    Require(normal.Evaluate(false, At(1020)).observations[0].condition == kFail, "position jump fails without debounce");
    const SE3 hidden = Pose(4.04);
    Feed(curve, 1020, Pose(std::sin(0.204) / 0.2, (1 - std::cos(0.204)) / 0.2, 0.204), 1, 0.2, 1, &hidden);
    Require(curve.Evaluate(false, At(1020)).observations[0].condition == kFail, "smoothing cannot hide a raw jump");
    Feed(normal, 1040, Pose(20), 1, 0, 2);
    Require(normal.Evaluate(false, At(1040)).observations[0].condition == kUnknown, "relocalization starts a new warmup epoch");
    for (int ms = 1060; ms <= 1600; ms += 20) Feed(normal, ms, Pose(20 + (ms - 1040) * 0.001), 1, 0, 2);
    Require(normal.Evaluate(false, At(1600)).observations[0].condition == kPass, "fresh epoch recovers with complete evidence");
    for (int ms = 1620; ms <= 1740; ms += 20) {
        normal.ObserveSpeed(101.6, 1, At(ms));  // frozen sensor clock
        normal.ObserveGyro(100 + ms * 0.001, 0, At(ms));
        normal.ObservePose(100 + ms * 0.001, Pose(20 + (ms - 1040) * 0.001), Pose(20 + (ms - 1040) * 0.001), SmallCov(), true, false, 2, At(ms));
    }
    Require(normal.Evaluate(false, At(1740)).observations[0].condition == kUnknown, "frozen CAN must not refresh evidence");
    Require(!normal.Evaluate(false, At(1740)).severe, "CAN unknown is separate from pose loss");
    Require(normal.Evaluate(false, At(1740)).speed_stamp_lag_sec > 0.13,
            "diagnostics expose the frozen CAN sensor-time lag");
    Mat6d large = SmallCov(); large(0, 0) = 1;
    Feed(reverse, 1020, Pose(-1.02), -1, 0, 1, nullptr, &large);
    Require(reverse.Evaluate(false, At(1020)).observations[1].condition == kFail, "major covariance ellipse exceeds limit");
    large(0, 0) = -1;
    Feed(reverse, 1040, Pose(-1.04), -1, 0, 1, nullptr, &large);
    Require(reverse.Evaluate(false, At(1040)).observations[1].condition == kUnknown, "negative covariance cannot pass");
    reverse.ObservePose(101.06, Pose(-1.06), Pose(-1.06), SmallCov(), false, true, 1, At(1060));
    Require(reverse.Evaluate(false, At(1060)).observations[1].condition == kFail, "unobservable graph is an uncertainty fault");
    reverse.ObservePose(101.08, Pose(-1.08), Pose(-1.08), SmallCov(), false, false, 1, At(1080));
    Require(reverse.Evaluate(false, At(1080)).observations[1].condition == kUnknown, "missing covariance cannot pass");
    reverse.ObservePose(101.10, Pose(-1.10), Pose(-1.10), Mat6d::Zero(), true, false, 1, At(1100));
    Require(reverse.Evaluate(false, At(1100)).observations[1].condition == kUnknown, "placeholder zero covariance cannot pass");

    // Physical 1-sigma limits remain active after removing the uncalibrated
    // default inflation. These cases separately exercise lateral and yaw limits.
    Mat6d lateral_fault = SmallCov(); lateral_fault(1, 1) = 0.31 * 0.31;
    Feed(reverse, 1120, Pose(-1.12), -1, 0, 1, nullptr, &lateral_fault);
    Require(reverse.Evaluate(false, At(1120)).observations[1].condition == kFail,
            "31 cm lateral uncertainty still fails the 30 cm limit");
    Mat6d yaw_fault = SmallCov(); yaw_fault(5, 5) = std::pow(2.1 * M_PI / 180.0, 2);
    Feed(reverse, 1140, Pose(-1.14), -1, 0, 1, nullptr, &yaw_fault);
    Require(reverse.Evaluate(false, At(1140)).observations[1].condition == kFail,
            "2.1 degree yaw uncertainty still fails the 2 degree limit");

    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(12, 12);
    h.topLeftCorner<6, 6>() = 15 * Mat6d::Identity();
    h.bottomRightCorner<6, 6>() = 5 * Mat6d::Identity();
    h.topRightCorner<6, 6>() = h.bottomLeftCorner<6, 6>() = -5 * Mat6d::Identity();
    Mat6d marginal;
    Require(loc::MarginalPoseCovariance(h, marginal) && (marginal - 0.3 * Mat6d::Identity()).norm() < 1e-9,
            "Schur marginal includes eliminated-state uncertainty");
    h.topLeftCorner<6, 6>() = 5 * Mat6d::Identity();
    Require(!loc::MarginalPoseCovariance(h, marginal), "unanchored graph cannot gain covariance from damping");
    Mat6d p = SmallCov(); p(5, 5) = 0.04;
    const auto rear = loc::TransformPoseCovariance(p, SE3(), Pose(-2), SO3());
    Require(rear(1, 1) > 0.15 && std::abs(rear(1, 5)) > 0.07, "rear axle lever arm and cross covariance retained");
    const auto propagated = loc::PropagatePoseCovariance(p, SE3(), Pose(1), 1, 0.2, 0.02);
    Require(propagated(1, 1) > p(1, 1) && propagated(5, 5) > p(5, 5), "dead reckoning grows uncertainty");

    diagnostic_monitor_interfaces::msg::DtcObservationArray msg;
    msg.comm_header.source_id = 13; msg.comm_header.seq = 1; msg.boot_id = 1;
    msg.interface_version = 1; msg.config_version = "test";
    for (uint32_t base : {0x110476C4u, 0x11037660u, 0x11017598u, 0x110275FCu})
        for (uint32_t offset : {0u, 1u, 3u, 4u, 7u, 8u, 9u, 10u, 11u, 12u, 13u}) {
            diagnostic_monitor_interfaces::msg::DtcObservation obs; obs.dtc = base + offset; obs.physical_condition = 1;
            msg.observations.push_back(obs);
        }
    for (uint32_t code : {0x10F00001u, 0x10F00002u}) {
        diagnostic_monitor_interfaces::msg::DtcObservation obs; obs.dtc = code; obs.physical_condition = 1;
        msg.observations.push_back(obs);
    }
    functional_safety::RadarDtcState radar(13, 0.5, 2.5);
    Require(radar.Observe(msg, At(0)) && radar.Evaluate(At(0)).all_pass, "complete radar snapshot accepted");
    Require(!radar.Observe(msg, At(300)) && !radar.Evaluate(At(501)).all_pass, "cached replay does not refresh heartbeat");
    msg.comm_header.seq = 2; msg.observations.back().physical_condition = 2;
    Require(radar.Observe(msg, At(510)) && radar.Evaluate(At(510)).all_failed, "all-failed flag accepted");
    Require(radar.Evaluate(At(1100)).all_failed, "timeout cannot silently clear previous severe state");
    msg.boot_id = 2; msg.comm_header.seq = 0; msg.observations.back().physical_condition = 0;
    Require(radar.Observe(msg, At(1110)) && radar.Evaluate(At(1110)).all_failed, "restart UNKNOWN cannot clear severe");
    msg.boot_id = 1; msg.comm_header.seq = 3;
    Require(!radar.Observe(msg, At(1120)), "retired boot cannot return");
    msg.boot_id = 2; msg.comm_header.seq = 1; msg.observations.back().physical_condition = 1;
    Require(radar.Observe(msg, At(1130)) && !radar.Evaluate(At(1130)).all_failed, "fresh PASS recovers severe state");
    msg.comm_header.seq = 2; msg.observations[0].dtc = msg.observations[1].dtc;
    Require(!radar.Observe(msg, At(1140)), "duplicate DTC invalidates snapshot");
    msg.observations[0].dtc = 0x110476C4;
    auto history_msg = msg;
    for (auto& obs : history_msg.observations) {
        obs.physical_condition = 1;
        for (uint32_t base : {0x110476C4u, 0x11037660u, 0x11017598u, 0x110275FCu})
            if (obs.dtc == base + 12 || obs.dtc == base + 13) obs.physical_condition = 2;
    }
    functional_safety::RadarDtcState history(13, 0.5, 2.5);
    Require(history.Observe(history_msg, At(1140)) && history.Evaluate(At(1140)).all_pass &&
            !history.Evaluate(At(1140)).fault, "fresh recovered/history warnings alone are not a current safety fault");
    Require(!history.Evaluate(At(1641)).all_pass, "historical warning exception does not waive heartbeat expiry");
    ++history_msg.comm_header.seq;
    history_msg.observations[1].physical_condition = 2;  // front current time-sync DTC
    Require(history.Observe(history_msg, At(1150)) && history.Evaluate(At(1150)).fault &&
            !history.Evaluate(At(1150)).all_pass, "current synchronization fault still requests safety action");
    ++history_msg.comm_header.seq;
    history_msg.observations[1].physical_condition = 1;
    history_msg.observations[9].physical_condition = 0;  // front recovered-link history
    Require(history.Observe(history_msg, At(1160)) && !history.Evaluate(At(1160)).all_pass,
            "unknown history is not fabricated as valid evidence");
    functional_safety::RadarDtcState availability(13, 0.5, 2.5, true);
    Require(availability.Observe(msg, At(1140)) && !availability.Evaluate(At(1140)).all_pass,
            "DTC PASS alone does not prove four devices usable");
    const std::array<uint16_t, 4> ids{{7, 12, 8, 9}};
    for (uint8_t role = 0; role < 4; ++role) {
        livox_ros_driver2::msg::LidarDiagnosticState state;
        state.role = role; state.comm_header.source_id = ids[role];
        state.comm_header.seq = 1; state.boot_id = 1; state.config_version = "test"; state.usability = 1;
        Require(availability.ObserveAvailability(state, At(1140)), "role availability accepted");
    }
    Require(availability.Evaluate(At(1140)).all_pass, "normal requires both valid DTCs and device usability");
    history_msg.observations[9].physical_condition = 2;
    ++history_msg.comm_header.seq;
    Require(availability.Observe(history_msg, At(1140)) && availability.Evaluate(At(1140)).all_pass,
            "healthy physical devices plus fresh history-only diagnostics allow normal");
    livox_ros_driver2::msg::LidarDiagnosticState fatal;
    fatal.role = 0; fatal.comm_header.source_id = 7; fatal.comm_header.seq = 2;
    fatal.boot_id = 1; fatal.config_version = "test"; fatal.usability = 2;
    Require(availability.ObserveAvailability(fatal, At(1150)) && availability.Evaluate(At(1150)).fault &&
            !availability.Evaluate(At(1150)).all_pass, "single unmapped native fatal gives radar fault");
    std::cout << "localization_dtc_test passed" << std::endl;
}
