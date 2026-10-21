#include "core/localization/localization_dtc_monitor.h"
#include "core/localization/pose_covariance.h"
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
    Require(SafetyValue(snapshot, false, true, false, false) == 2, "radar mask defaults off");
    Require(SafetyValue(snapshot, false, true, false, true) == 0, "explicit mask suppresses radar-only fault");
    Require(SafetyValue(snapshot, false, false, false, true) == 0, "explicit mask also suppresses unknown radar branch");
    Require(SafetyValue(snapshot, true, true, false, true) == -9, "mask preserves all-lidar-failed");
    auto masked_severe = snapshot;
    masked_severe.severe = true;
    Require(SafetyValue(masked_severe, false, true, false, true) == -9, "mask preserves localization severe");
    auto masked_unknown = snapshot;
    masked_unknown.observations[0].condition = kUnknown;
    Require(SafetyValue(masked_unknown, false, true, false, true) == 4, "mask preserves localization unknown");
    masked_unknown.observations[0].condition = kFail;
    Require(SafetyValue(masked_unknown, false, true, false, true) == 4, "mask preserves localization failure");
    const auto healthy_loc = snapshot;
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
    Require(normal.Evaluate(false, At(1140)).observations[0].condition == kFail,
            "input expiry must not overwrite a definite motion failure");
    Require(SafetyValue(normal.Evaluate(false, At(1140)), false, false, true) == 4,
            "definite failure remains a stop while output is fresh");
    const SE3 hidden = Pose(4.04);
    Feed(curve, 1020, Pose(std::sin(0.204) / 0.2, (1 - std::cos(0.204)) / 0.2, 0.204), 1, 0.2, 1, &hidden);
    Require(curve.Evaluate(false, At(1020)).observations[0].condition == kFail, "smoothing cannot hide a raw jump");
    Feed(normal, 1040, Pose(20), 1, 0, 2);
    Require(normal.Evaluate(false, At(1040)).observations[0].condition == kUnknown, "relocalization starts a new warmup epoch");
    Require(normal.Evaluate(false, At(1040)).confirmed_failure,
            "new epoch warmup cannot clear a confirmed fault");
    for (int ms = 1060; ms <= 1600; ms += 20) Feed(normal, ms, Pose(20 + (ms - 1040) * 0.001), 1, 0, 2);
    Require(normal.Evaluate(false, At(1600)).observations[0].condition == kPass, "fresh epoch recovers with complete evidence");
    Require(!normal.Evaluate(false, At(1600)).confirmed_failure,
            "fresh passing evidence clears a confirmed fault");
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
    Require(SafetyValue(reverse.Evaluate(false, At(1080)), false, false, true) == 4,
            "missing covariance cannot clear the prior unobservable failure");

    Config grace_cfg = cfg;
    grace_cfg.consistency_unknown_grace_sec = 0.1;
    Monitor grace(grace_cfg, At(0));
    Require(!grace.Evaluate(false, At(0)).consistency_unknown_tolerated,
            "startup without evidence never receives grace");
    for (int ms = 0; ms <= 1000; ms += 20) Feed(grace, ms, Pose(ms * 0.001));
    // A small hole in independent wheel-speed samples makes integration
    // unknown, while newest input, output, and prior PASS remain fresh.
    Monitor hole(grace_cfg, At(0));
    for (int ms = 0; ms <= 1000; ms += 20) Feed(hole, ms, Pose(ms * 0.001));
    hole.ObserveSpeed(101.11, 1, At(1020));
    hole.ObserveGyro(101.02, 0, At(1020));
    hole.ObservePose(101.02, Pose(1.02), Pose(1.02), SmallCov(), true, false, 1, At(1020));
    const auto short_hole = hole.Evaluate(false, At(1020));
    Require(short_hole.observations[0].condition == kUnknown &&
            short_hole.consistency_unknown_tolerated &&
            SafetyValue(short_hole, false, false, true) == 0,
            "bounded missing integration evidence remains UNKNOWN without immediate stop");
    Require(!hole.Evaluate(false, At(1101)).consistency_unknown_tolerated,
            "UNKNOWN cannot extend its own evidence deadline");
    hole.ObserveSpeed(104.6, 1, At(1040));
    hole.ObserveGyro(104.6, 0, At(1040));
    hole.ObservePose(104.6, Pose(1.04), Pose(1.04), SmallCov(), true, false, 1, At(1040));
    Require(!hole.Evaluate(false, At(1040)).consistency_unknown_tolerated,
            "simultaneous sensor clock jump cannot reuse pre-jump passing evidence");
    auto tolerated = grace.Evaluate(false, At(1000));
    tolerated.observations[0].condition = kUnknown;
    tolerated.consistency_unknown_tolerated = true;
    Require(SafetyValue(tolerated, false, false, true) == 0,
            "eligible UNKNOWN keeps its raw state while allowing bounded operation");
    tolerated.confirmed_failure = true;
    Require(SafetyValue(tolerated, false, false, true) == 4,
            "grace never overrides a confirmed fault");
    tolerated.severe = true;
    Require(SafetyValue(tolerated, false, false, true) == -9,
            "grace never overrides severe conditions");
    for (int ms = 1020; ms <= 1160; ms += 20) {
        grace.ObserveGyro(100 + ms * 0.001, 0, At(ms));
        grace.ObservePose(100 + ms * 0.001, Pose(ms * 0.001), Pose(ms * 0.001),
                          SmallCov(), true, false, 1, At(ms));
    }
    Require(!grace.Evaluate(false, At(1160)).consistency_unknown_tolerated &&
            SafetyValue(grace.Evaluate(false, At(1160)), false, false, true) == 4,
            "frozen CAN cannot be hidden by continuous fresh poses");
    Require(grace.Evaluate(false, At(1461)).severe,
            "grace does not extend the 300 ms output timeout");
    Feed(grace, 1180, Pose(1.18), 1, 0, 2);
    Require(!grace.Evaluate(false, At(1180)).consistency_unknown_tolerated,
            "epoch change invalidates prior passing evidence");
    bool invalid_grace_rejected = false;
    try {
        auto invalid = cfg;
        invalid.consistency_unknown_grace_sec = invalid.output_timeout_sec;
        Monitor rejected(invalid, At(0));
    } catch (const std::invalid_argument&) { invalid_grace_rejected = true; }
    Require(invalid_grace_rejected, "grace must remain below the unchanged severe output timeout");

    // Reproduce the stationary 2026-10-21 11:28 event: both independent
    // arrivals are 107 ms old, a pose just arrived, covariance still passes.
    Config bounded_cfg = cfg;
    bounded_cfg.consistency_unknown_grace_sec = 0.15;
    Monitor strict_gap(cfg, At(0)), bounded_gap(bounded_cfg, At(0));
    for (int ms = 0; ms <= 1000; ms += 20) {
        Feed(strict_gap, ms, Pose(0), 0);
        Feed(bounded_gap, ms, Pose(0), 0);
    }
    for (auto* monitor : {&strict_gap, &bounded_gap})
        monitor->ObservePose(101.107, Pose(0), Pose(0), SmallCov(), true, false, 1, At(1107));
    const auto delayed = bounded_gap.Evaluate(false, At(1107));
    Require(SafetyValue(strict_gap.Evaluate(false, At(1107)), false, false, true) == 4,
            "zero grace reproduces the original stop at 107 ms");
    Require(delayed.observations[0].condition == kUnknown &&
            delayed.observations[1].condition == kPass && delayed.consistency_unknown_tolerated &&
            SafetyValue(delayed, false, false, true) == 0,
            "107 ms arrival gap may use bounded prior evidence without pretending consistency passed");
    Require(SafetyValue(delayed, false, true, false) == 2 &&
            SafetyValue(delayed, true, true, false) == -9,
            "motion grace never bypasses radar faults or all-radar failure");
    Monitor delayed_delivery(bounded_cfg, At(0));
    for (int ms = 0; ms <= 1000; ms += 20) Feed(delayed_delivery, ms, Pose(0), 0);
    delayed_delivery.ObserveSpeed(101.108674, 0, At(1035));
    delayed_delivery.ObserveGyro(101.105777, 0, At(1035));
    delayed_delivery.ObservePose(101.105777, Pose(0), Pose(0), SmallCov(), true, false, 1, At(1142));
    Require(delayed_delivery.Evaluate(false, At(1142)).consistency_unknown_tolerated &&
            !delayed_delivery.Evaluate(false, At(1151)).consistency_unknown_tolerated,
            "a delayed pose may share the last IMU stamp, but the deadline remains anchored to the older PASS");
    Require(bounded_gap.Evaluate(false, At(1150)).consistency_unknown_tolerated,
            "150 ms total evidence-age boundary is inclusive");
    bounded_gap.ObservePose(101.151, Pose(0), Pose(0), SmallCov(), true, false, 1, At(1151));
    Require(!bounded_gap.Evaluate(false, At(1151)).consistency_unknown_tolerated &&
            SafetyValue(bounded_gap.Evaluate(false, At(1151)), false, false, true) == 4,
            "fresh poses and repeated UNKNOWN cannot extend the 150 ms total deadline");
    Require(bounded_gap.Evaluate(false, At(1452)).severe && delayed.confirmed_failure == false,
            "output timeout remains 300 ms and missing evidence is not a confirmed failure");
    Feed(bounded_gap, 1160, Pose(0), 0);
    Require(bounded_gap.Evaluate(false, At(1160)).observations[0].condition == kUnknown,
            "receiving one new input does not prove a missing integration window is complete");
    for (int ms = 1180; ms <= 1800; ms += 20) Feed(bounded_gap, ms, Pose(0), 0);
    Require(SafetyValue(bounded_gap.Evaluate(false, At(1800)), false, false, true) == 0,
            "complete fresh evidence recovers after the missing interval leaves the window");

    Monitor clock_jump(bounded_cfg, At(0));
    for (int ms = 0; ms <= 1000; ms += 20) Feed(clock_jump, ms, Pose(0), 0);
    clock_jump.ObserveSpeed(101.14, 0, At(1020));
    clock_jump.ObserveGyro(101.14, 0, At(1020));
    clock_jump.ObservePose(101.14, Pose(0), Pose(0), SmallCov(), true, false, 1, At(1020));
    Require(!clock_jump.Evaluate(false, At(1020)).consistency_unknown_tolerated,
            "a source-clock jump cannot borrow pre-jump PASS even inside the grace interval");

    Monitor moving_gap(bounded_cfg, At(0));
    for (int ms = 0; ms <= 1000; ms += 20) Feed(moving_gap, ms, Pose(ms * 0.005), 5);
    moving_gap.ObservePose(101.107, Pose(5.535), Pose(5.535), SmallCov(), true, false, 1, At(1107));
    Require(moving_gap.Evaluate(false, At(1107)).consistency_unknown_tolerated,
            "the same bounded timing rule applies at the stated 5 m/s operating limit");
    // A new definite mismatch must still stop immediately, even with grace enabled.
    for (int ms = 1120; ms <= 1800; ms += 20) Feed(moving_gap, ms, Pose(ms * 0.005), 5);
    Feed(moving_gap, 1820, Pose(12.1), 5);
    Require(moving_gap.Evaluate(false, At(1820)).confirmed_failure &&
            SafetyValue(moving_gap.Evaluate(false, At(1820)), false, false, true) == 4,
            "a definite moving-vehicle jump remains an immediate stop");
    moving_gap.ObservePose(101.927, Pose(12.635), Pose(12.635), SmallCov(), true, false, 1, At(1927));
    Require(moving_gap.Evaluate(false, At(1927)).confirmed_failure &&
            !moving_gap.Evaluate(false, At(1927)).consistency_unknown_tolerated &&
            SafetyValue(moving_gap.Evaluate(false, At(1927)), false, false, true) == 4,
            "a later input gap cannot hide the confirmed moving-vehicle jump");
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

    std::cout << "localization_dtc_test passed" << std::endl;
}
