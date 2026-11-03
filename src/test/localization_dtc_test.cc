#include "core/localization/localization_dtc_monitor.h"
#include "core/localization/pose_covariance.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <yaml-cpp/yaml.h>

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
void FeedWithoutCan(Monitor& monitor, int ms, const SE3& pose, double gyro = 0,
                    const SE3* raw = nullptr) {
    const double stamp = 100 + ms * 1e-3;
    monitor.ObserveGyro(stamp, gyro, At(ms));
    monitor.ObservePose(stamp, pose, raw ? *raw : pose, SmallCov(), true, false, 1, At(ms));
}
void CheckNativeYawIntegration() {
    Config cfg;  // Deliberately keep the original 2-degree base threshold.
    Require(std::abs(cfg.yaw_rad * 180 / M_PI - 2) < 1e-12, "generic yaw threshold stays two degrees");
    Monitor vibration(cfg, At(0));
    double legacy_angle = 0, legacy_turn = 0;
    for (int ms = 0; ms <= 500; ms += 5) {
        const double omega = .7 * std::sin(2 * M_PI * 47 * ms * .001);
        vibration.ObserveGyro(100 + ms * .001, omega, At(ms));
        if (ms % 20 == 0) vibration.ObservePose(100 + ms * .001, Pose(0), Pose(0),
            SmallCov(), true, false, 1, At(ms));
        if (ms % 20 == 10) { legacy_angle += omega * .02; legacy_turn += std::abs(omega) * .02; }
    }
    Require(std::abs(legacy_angle) > cfg.yaw_rad + cfg.angle_ratio * legacy_turn,
            "47 Hz / 200 Hz fixture reproduces the old 20 ms midpoint false failure");
    const auto state = vibration.Evaluate(false, At(500));
    const auto& e = state.yaw_evidence;
    Require(state.yaw_condition == kPass && SafetyValue(state, false, false, true) == 0 &&
            e.available && e.window && e.segments == 100 && e.held_sec == 0 &&
            std::abs(e.gyro_delta_rad) < .005 && e.published_error_rad < .005,
            "all 100 native segments reject the alias without relaxing the original threshold");

    // Native source jitter and a pose cadence asynchronous to the IMU must not
    // reintroduce pose-rate quadrature or skip the partial boundary segments.
    for (int cadence : {17, 23}) {
        Monitor asynchronous(cfg, At(0));
        for (int ms = 0, sample = 0; ms <= 900; ++ms) {
            if (ms == sample * 5 + sample % 2) {
                asynchronous.ObserveGyro(100 + ms * .001,
                    .7 * std::sin(2 * M_PI * 47 * ms * .001), At(ms));
                ++sample;
            }
            if (ms % cadence == 0) {
                asynchronous.ObservePose(100 + ms * .001, Pose(0), Pose(0), SmallCov(), true, false, 1, At(ms));
                if (ms > 550) Require(asynchronous.Evaluate(false, At(ms)).yaw_condition == kPass,
                    "jittered native IMU and asynchronous pose boundaries still pass");
            }
        }
    }
    Config short_cfg = cfg; short_cfg.window_sec = .02;
    for (double x : {0.0, 1.0}) {
        Monitor partial(short_cfg, At(0));
        partial.ObserveSpeed(100, 0, At(0));
        partial.ObservePose(100, Pose(0), Pose(0), SmallCov(), true, false, 1, At(0));
        partial.ObserveGyro(100.001, 0, At(1));
        partial.ObserveGyro(100.02, 0, At(20));
        partial.ObserveSpeed(100.02, 0, At(20));
        partial.ObservePose(100.02, Pose(x), Pose(x), SmallCov(), true, false, 1, At(20));
        const auto unavailable = partial.Evaluate(false, At(20));
        Require(unavailable.yaw_condition == kUnknown && !unavailable.yaw_evidence.available &&
                unavailable.translation_condition == (x == 0 ? kPass : kFail) &&
                SafetyValue(unavailable, false, false, true) == (x == 0 ? 0 : 4),
                "missing native yaw endpoint must preserve the original independent translation decision");
    }
    Monitor crossing(short_cfg, At(0));
    crossing.ObserveGyro(100, 1, At(0));
    crossing.ObservePose(100, Pose(0), Pose(0), SmallCov(), true, false, 1, At(0));
    crossing.ObserveGyro(100.02, -1, At(20));
    crossing.ObservePose(100.02, Pose(0), Pose(0), SmallCov(), true, false, 1, At(20));
    const auto triangle = crossing.Evaluate(false, At(20)).yaw_evidence;
    Require(triangle.available && std::abs(triangle.gyro_delta_rad) < 1e-12 &&
            std::abs(triangle.gyro_turn_rad - .01) < 1e-12 &&
            std::abs(triangle.limit_rad - cfg.yaw_rad - .001) < 1e-12,
            "zero crossing uses the exact absolute linear area for the dynamic threshold");

    Config gap_cfg = cfg; gap_cfg.window_sec = .2;
    Monitor gap(gap_cfg, At(0));
    gap.ObserveGyro(100, 0, At(0));
    gap.ObservePose(100, Pose(0), Pose(0), SmallCov(), true, false, 1, At(0));
    gap.ObserveGyro(100.2, 0, At(200));
    gap.ObservePose(100.2, Pose(0), Pose(0), SmallCov(), true, false, 1, At(200));
    Require(gap.Evaluate(false, At(200)).yaw_condition == kUnknown &&
            !gap.Evaluate(false, At(200)).yaw_evidence.available,
            "exact endpoint samples cannot conceal a source gap longer than 100 ms");

    Config sany = cfg; sany.yaw_rad = 3 * M_PI / 180;
    for (bool raw_only : {false, true}) {
        Monitor real_jump(sany, At(0));
        for (int ms = 0; ms <= 500; ms += 20) FeedWithoutCan(real_jump, ms, Pose(0));
        const SE3 small = Pose(0, 0, 2.5 * M_PI / 180);
        FeedWithoutCan(real_jump, 520, small);
        Require(SafetyValue(real_jump.Evaluate(false, At(520)), false, false, true) == 0,
                "SANY three-degree base tolerates a 2.5-degree zero-turn residual");
        const SE3 jump = Pose(0, 0, 9 * M_PI / 180);
        FeedWithoutCan(real_jump, 540, raw_only ? small : jump, 0, raw_only ? &jump : &small);
        const auto failure = real_jump.Evaluate(false, At(540));
        Require(failure.yaw_condition == kFail && failure.confirmed_failure &&
                SafetyValue(failure, false, false, true) == 4 && !failure.yaw_evidence.window,
                "published-only and raw-only jumps still stop immediately under the relaxed threshold");
        Require((raw_only ? failure.yaw_evidence.raw_error_rad : failure.yaw_evidence.published_error_rad) >
                    failure.yaw_evidence.limit_rad &&
                (raw_only ? failure.yaw_evidence.published_error_rad : failure.yaw_evidence.raw_error_rad) < 1e-10,
                "split evidence identifies the actual failed branch and comparison");
    }
    Monitor fallback(cfg, At(0));
    fallback.ObserveGyro(100, 0, At(0));
    fallback.ObservePose(100, Pose(0), Pose(0), SmallCov(), true, false, 1, At(0), false);
    fallback.ObserveGyro(100.02, 0, At(20));
    fallback.ObservePose(100.02, Pose(0), Pose(0), SmallCov(), true, false, 1, At(20));
    Require(!fallback.Evaluate(false, At(20)).yaw_evidence.raw_independent,
            "a fallback at either endpoint is not independent raw-branch evidence");
}
void CheckConfirmedLossPolicy() {
    Config cfg;
    Require(cfg.output_timeout_sec == .5 && cfg.evaluation_timeout_sec == .5,
            "output and output-quality freshness use the same 500 ms boundary");
    Monitor startup(cfg, At(0));
    TrackingEvidence tracking;
    for (int ms : {0, 7999, 8000}) {
        const auto state = startup.Evaluate(tracking, At(ms));
        Require(state.startup_pending && !state.startup_expired && !state.localization_lost &&
                SafetyValue(state, false, false, true) == 4,
                "initial relocalization stays 4 through the inclusive eight-second boundary");
    }
    const auto expired = startup.Evaluate(tracking, At(8001));
    Require(expired.startup_expired && !expired.localization_lost &&
            SafetyValue(expired, false, false, true) == -9,
            "initial timeout is severe with a distinct startup cause");
    tracking.initialized = tracking.output_gate_open = tracking.latest_match_valid = true;
    tracking.valid_match_age_sec = 0;
    Require(SafetyValue(startup.Evaluate(tracking, At(8002)), false, false, true) == 4 &&
            !startup.Evaluate(tracking, At(8002)).startup_expired,
            "an accepted first map match ends startup grace even before business output exists");

    Monitor running(cfg, At(0));
    for (int ms = 0; ms <= 1000; ms += 20) Feed(running, ms, Pose(ms * .001));
    const auto healthy = running.Evaluate(tracking, At(1000));
    Require(healthy.output_continuous && HealthyLocalizationForDegradedRadar(healthy),
            "sustained fresh output and accepted map/quality evidence authorize degraded operation");
    functional_safety::LidarFrameWatchdog::Snapshot frames;
    frames.started = frames.seen = true;
    Require(EvaluateFrameSafetyDecision(healthy, frames, true).published == 0,
            "v6 healthy frame supply has no device-DTC stop input");
    frames.fault = true;
    Require(EvaluateFrameSafetyDecision(healthy, frames, true).published == 2,
            "v6 frame timeout stops despite healthy static/extrapolated pose");
    auto consequence = healthy;
    consequence.output_stale = consequence.severe = consequence.startup_expired = true;
    Require(EvaluateFrameSafetyDecision(consequence, frames, true).published == 2,
            "frame supply root cause precedes consequent output/startup timeout");
    consequence.localization_lost = true;
    Require(EvaluateFrameSafetyDecision(consequence, frames, true).published == -9,
            "explicit confirmed tracking loss keeps severe priority");
    const auto disabled_frame_fault = EvaluateFrameSafetyDecision(healthy, frames, false);
    Require(disabled_frame_fault.raw == 2 && disabled_frame_fault.published == 0,
            "reporting override preserves frame fault evidence");
    frames.fault = false;
    consequence = healthy; consequence.confirmed_failure = true;
    Require(EvaluateFrameSafetyDecision(consequence, frames, true, {true, true}).published == 4,
            "frame health does not hide localization quality failure");
    consequence = healthy;
    consequence.startup_pending = consequence.startup_expired = consequence.severe = true;
    frames.seen = false;
    Require(EvaluateFrameSafetyDecision(consequence, frames, true).published == 4,
            "an earlier monitor startup cannot bypass the frame watchdog's no-input grace");
    frames.fault = true;
    Require(EvaluateFrameSafetyDecision(consequence, frames, true).published == 2,
            "no first frame at the watchdog startup deadline is a supply fault");
    frames.fault = false; frames.seen = true;
    Require(EvaluateFrameSafetyDecision(consequence, frames, true).published == -9,
            "startup expiry still applies when valid fused input cannot establish localization");
    const auto degraded = EvaluateSafetyDecision(healthy, false, true, false, true, false);
    Require(degraded.raw == 2 && degraded.selected == 0 && degraded.published == 0 && degraded.partial_radar_suppressed,
            "partial radar failure preserves raw 2 and suppresses only its stop decision with healthy localization");
    const auto reporting_off = EvaluateSafetyDecision(healthy, false, true, false, false, false);
    Require(reporting_off.raw == 2 && reporting_off.partial_radar_suppressed && reporting_off.published == 0,
            "disabling Safety reporting preserves the partial radar cause and suppression evidence");
    Require(EvaluateSafetyDecision(healthy, true, true, false, true, false).selected == -9,
            "healthy localization never overrides all required radars unavailable");
    for (int missing = 0; missing < 5; ++missing) {
        auto incomplete = healthy;
        if (missing == 0) incomplete.output_continuous = false;
        if (missing == 1) incomplete.match_evidence_available = false;
        if (missing == 2) incomplete.latest_match_valid = false;
        if (missing == 3) incomplete.covariance_valid = false;
        if (missing == 4) incomplete.match_stale = true;
        const auto decision = EvaluateSafetyDecision(incomplete, false, true, false, true, false);
        Require(!decision.partial_radar_suppressed && decision.selected != 0,
                "missing positive evidence never enables partial-radar suppression");
    }
    auto unknown = healthy;
    unknown.observations[0].condition = kUnknown;
    unknown.yaw_condition = kUnknown;
    unknown.consistency_unknown_tolerated = true;
    Require(EvaluateSafetyDecision(unknown, false, true, false, true, false).selected == 2,
            "bounded UNKNOWN grace is insufficient evidence for radar-degraded operation");
    unknown.consistency_unknown_tolerated = false;
    Require(EvaluateSafetyDecision(unknown, false, true, false, true, false).selected == 2,
            "motion UNKNOWN is advisory but does not authorize suppressing a radar fault");
    auto failed = healthy;
    failed.confirmed_failure = true;
    Require(EvaluateSafetyDecision(failed, false, true, false, true, false).selected == 4,
            "positive sensor evidence cannot erase a confirmed localization-quality failure");
    Monitor optional(cfg, At(0));
    for (int ms = 0; ms <= 1000; ms += 20) FeedWithoutCan(optional, ms, Pose(ms * .001));
    Require(EvaluateSafetyDecision(optional.Evaluate(tracking, At(1000)), false, true, false, true, false).partial_radar_suppressed,
            "legal optional CAN absence still requires explicit passing yaw, map and uncertainty evidence");

    tracking.latest_match_valid = tracking.output_gate_open = false;
    const auto rejected = running.Evaluate(tracking, At(1000));
    Require(!rejected.localization_lost && !rejected.output_stale &&
            SafetyValue(rejected, false, false, true) == 4,
            "one rejected match closes outputs and reports 4 without inventing confirmed loss");
    tracking.latest_match_valid = true;
    tracking.match_stale = true;
    const auto stale_match = running.Evaluate(tracking, At(1000));
    Require(!stale_match.localization_lost && SafetyValue(stale_match, false, false, true) == 4,
            "map-match freshness loss alone remains 4 while a recent output exists");
    tracking.confirmed_lost = true;
    Require(SafetyValue(running.Evaluate(tracking, At(1000)), false, false, true) == -9,
            "confirmed runtime loss has no new startup grace");
    tracking.confirmed_lost = tracking.match_stale = false;
    tracking.output_gate_open = true;
    // Only outputs pause. Independent input callbacks continue, avoiding a separate input-gap fault.
    for (int ms = 1020; ms <= 1500; ms += 20) {
        running.ObserveSpeed(100 + ms * .001, 1, At(ms));
        running.ObserveGyro(100 + ms * .001, 0, At(ms));
    }
    for (int ms : {1490, 1500}) {
        const auto state = running.Evaluate(tracking, At(ms));
        Require(!state.output_stale && !state.severe && SafetyValue(state, false, false, true) == 0,
                "fresh inputs do not turn an output-only pause into an early 300 ms quality timeout");
    }
    const auto timed_out = running.Evaluate(tracking, At(1501));
    Require(timed_out.output_stale && !timed_out.severe && !timed_out.localization_lost &&
            SafetyValue(timed_out, false, false, true) == 4,
            "strictly over 500 ms output silence is 4, independently of confirmed loss");
    Feed(running, 1520, Pose(1.52));
    Require(!running.Evaluate(tracking, At(1520)).output_continuous,
            "one recovered output cannot immediately certify continuous degraded operation");
    for (int ms = 1540; ms <= 2020; ms += 20) Feed(running, ms, Pose(ms * .001));
    Require(running.Evaluate(tracking, At(2020)).output_continuous,
            "one complete motion window of recovered output requalifies continuity");
    Feed(running, 2040, Pose(2.04), 1, 0, 2);
    Require(!running.Evaluate(tracking, At(2040)).output_continuous,
            "a new pose epoch cannot borrow continuous output evidence from before relocalization");
}
void CheckQualityFaultSwitches() {
    const QualityFaultReporting defaults;
    Require(!defaults.consistency && !defaults.uncertainty, "quality stop switches default off");
    Monitor monitor(Config{}, At(0));
    for (int ms = 0; ms <= 1000; ms += 20) Feed(monitor, ms, Pose(ms * .001));
    TrackingEvidence tracking;
    tracking.initialized = tracking.output_gate_open = tracking.latest_match_valid = true;
    const auto healthy = monitor.Evaluate(tracking, At(1000));
    functional_safety::LidarFrameWatchdog::Snapshot frames;
    frames.started = frames.seen = true;
    for (bool consistency : {false, true}) for (bool uncertainty : {false, true}) {
        const QualityFaultReporting quality{consistency, uncertainty};
        for (int kind = 0; kind < 7; ++kind) {
            auto state = healthy;
            const bool consistency_failure = kind == 0 || kind == 1 || kind == 4 || kind == 6;
            const bool uncertainty_failure = kind == 2 || kind == 5 || kind == 6;
            state.confirmed_consistency_failure = consistency_failure;
            state.confirmed_uncertainty_failure = uncertainty_failure;
            state.confirmed_failure = consistency_failure || uncertainty_failure;
            if (kind == 0 || kind == 6) state.translation_condition = state.observations[0].condition = kFail;
            if (kind == 1) state.yaw_condition = state.observations[0].condition = kFail;
            if (kind == 2 || kind == 6) state.observations[1].condition = kFail;
            if (kind == 3 || kind == 5) state.observations[1].condition = kUnknown;
            if (kind == 4) state.translation_condition = state.observations[0].condition = kUnknown;
            const int expected = ((consistency && consistency_failure) ||
                (uncertainty && (uncertainty_failure || kind == 3))) ? 4 : 0;
            const auto decision = EvaluateFrameSafetyDecision(state, frames, true, quality);
            Require(decision.raw == 4 && decision.selected == expected && decision.published == expected,
                "independent quality switches select the correct current and latched category");
            Require(state.confirmed_failure == (consistency_failure || uncertainty_failure),
                "quality switches preserve original failure evidence");
            const auto disabled = EvaluateFrameSafetyDecision(state, frames, false, quality);
            Require(disabled.raw == 4 && disabled.selected == expected && disabled.published == 0,
                "global reporting override does not erase either decision");
        }
        auto state = healthy;
        state.output_stale = true;
        Require(EvaluateFrameSafetyDecision(state, frames, true, quality).published == 4,
            "quality switches never mask output timeout");
        state = healthy; state.output_gate_open = false;
        Require(EvaluateFrameSafetyDecision(state, frames, true, quality).published == 4,
            "quality switches never mask closed match output gate");
        state = healthy; state.startup_pending = true;
        Require(EvaluateFrameSafetyDecision(state, frames, true, quality).published == 4,
            "quality switches never mask initial localization wait");
        state = healthy; state.severe = state.localization_lost = true;
        frames.fault = true;
        Require(EvaluateFrameSafetyDecision(state, frames, true, quality).published == -9,
            "confirmed loss keeps priority over frame fault and quality switches");
        Require(EvaluateFrameSafetyDecision(healthy, frames, true, quality).published == 2,
            "quality switches never mask raw frame supply fault");
        frames.fault = false;
    }
    Feed(monitor, 1020, Pose(3.02));
    const auto failed = monitor.Evaluate(tracking, At(1020));
    Require(failed.confirmed_consistency_failure && !failed.confirmed_uncertainty_failure,
        "monitor exposes the actual consistency failure category");
    for (int ms = 1040; ms <= 1400; ms += 20) FeedWithoutCan(monitor, ms, Pose(2 + ms * .001));
    const auto latched = monitor.Evaluate(tracking, At(1400));
    Require(latched.confirmed_consistency_failure && !latched.confirmed_uncertainty_failure &&
        latched.observations[0].condition == kUnknown &&
        EvaluateFrameSafetyDecision(latched, frames, true, {false, true}).published == 0 &&
        EvaluateFrameSafetyDecision(latched, frames, true, {true, false}).published == 4,
        "an actual consistency latch cannot leak into the uncertainty-only stop branch");
    Monitor covariance(Config{}, At(0));
    for (int ms = 0; ms <= 1000; ms += 20) Feed(covariance, ms, Pose(ms * .001));
    const Mat6d large = Mat6d::Identity(), invalid = Mat6d::Zero();
    Feed(covariance, 1020, Pose(1.02), 1, 0, 1, nullptr, &large);
    Feed(covariance, 1040, Pose(1.04), 1, 0, 1, nullptr, &invalid);
    const auto uncertain = covariance.Evaluate(tracking, At(1040));
    Require(uncertain.confirmed_uncertainty_failure && !uncertain.confirmed_consistency_failure &&
        EvaluateFrameSafetyDecision(uncertain, frames, true, {true, false}).published == 0 &&
        EvaluateFrameSafetyDecision(uncertain, frames, true, {false, true}).published == 4,
        "an actual covariance failure remains independently selectable after UNKNOWN");
}
void CheckMotionUnknownAdvisory() {
    for (double grace : {0.0, 0.15}) {
        for (double speed : {0.0, 5.0}) {
            for (int gap : {150, 151, 208, 400}) {
                Config cfg;
                cfg.consistency_unknown_grace_sec = grace;
                Monitor monitor(cfg, At(0));
                for (int ms = 0; ms <= 1000; ms += 20)
                    Feed(monitor, ms, Pose(speed * ms * .001), speed);
                const int now = 1000 + gap;
                const double stamp = 100 + now * .001;
                const SE3 pose = Pose(speed * now * .001);
                monitor.ObserveSpeed(stamp, speed, At(now));
                monitor.ObservePose(stamp, pose, pose, SmallCov(), true, false, 1, At(now));
                TrackingEvidence tracking;
                tracking.initialized = tracking.output_gate_open = tracking.latest_match_valid = true;
                tracking.valid_match_age_sec = 0;
                const auto state = monitor.Evaluate(tracking, At(now));
                Require(state.observations[0].condition == kUnknown && state.yaw_condition == kUnknown &&
                        state.observations[1].condition == kPass && !state.confirmed_failure,
                        "an IMU evidence gap stays UNKNOWN rather than PASS or confirmed failure");
                Require(std::abs(state.gyro_arrival_age_sec - gap * .001) < 1e-9,
                        "the real input evidence age is retained for diagnosis");
                functional_safety::LidarFrameWatchdog::Snapshot frames;
                frames.started = frames.seen = true;
                const auto decision = EvaluateFrameSafetyDecision(state, frames, true);
                Require(decision.raw == 0 && decision.selected == 0 && decision.published == 0,
                        "150-208 ms and longer motion UNKNOWN does not independently publish 4 with fresh accepted output");
                frames.fault = true;
                Require(EvaluateFrameSafetyDecision(state, frames, true).published == 2,
                        "motion UNKNOWN advisory cannot bypass a raw radar frame-supply fault");
                Require(SafetyValue(monitor.Evaluate(tracking, At(now + 501)), false, false, true) == 4,
                        "output timeout remains independent of advisory motion UNKNOWN");
            }
        }
    }
}
void CheckSanyUncertaintyConfig(const char* path) {
    const auto diagnostics = YAML::LoadFile(path)["diagnostics"];
    Config cfg;
    cfg.major_sigma_m = diagnostics["major_sigma_limit_m"].as<double>();
    cfg.lateral_sigma_m = diagnostics["lateral_sigma_limit_m"].as<double>();
    cfg.yaw_sigma_rad = diagnostics["yaw_sigma_limit_deg"].as<double>() * M_PI / 180.0;
    cfg.covariance_scale = diagnostics["covariance_variance_scale"].as<double>();
    cfg.position_sigma_floor_m = diagnostics["position_sigma_floor_m"].as<double>();
    cfg.angle_sigma_floor_rad = diagnostics["angle_sigma_floor_deg"].as<double>() * M_PI / 180.0;
    cfg.output_timeout_sec = diagnostics["output_timeout_sec"].as<double>();
    cfg.evaluation_timeout_sec = diagnostics["evaluation_timeout_sec"].as<double>();
    cfg.longitudinal_m = diagnostics["longitudinal_tolerance_m"].as<double>();
    cfg.lateral_m = diagnostics["lateral_tolerance_m"].as<double>();
    cfg.yaw_rad = diagnostics["yaw_tolerance_deg"].as<double>() * M_PI / 180.0;
    Require(!diagnostics["consistency_fault_reporting_enabled"].as<bool>() &&
            !diagnostics["uncertainty_fault_reporting_enabled"].as<bool>(),
            "SANY quality stop switches are explicitly disabled by default");
    Require(std::abs(cfg.yaw_rad * 180 / M_PI - 3) < 1e-12 &&
            diagnostics["angle_tolerance_ratio"].as<double>() == .1,
            "SANY yaw base is three degrees with the existing ten-percent turn term");
    // Logged sigma already includes the configured additive variance floors.
    // Construct input covariance for the desired final sigma at identity pose.
    const auto covariance_for = [&](double x_m, double y_m, double yaw_deg) {
        Mat6d covariance = SmallCov();
        const double floor_position = std::pow(cfg.position_sigma_floor_m, 2);
        covariance(0, 0) = (x_m * x_m - floor_position) / cfg.covariance_scale;
        covariance(1, 1) = (y_m * y_m - floor_position) / cfg.covariance_scale;
        covariance(5, 5) = (std::pow(yaw_deg * M_PI / 180.0, 2) -
                            std::pow(cfg.angle_sigma_floor_rad, 2)) / cfg.covariance_scale;
        return covariance;
    };
    // Conservative envelope of the logged Oct 26-27 finite sigma peaks.
    // This rechecks the configured monitor; it is not a raw sensor replay.
    const Mat6d recorded = covariance_for(0.369499, 0.343525, 4.324468);
    Monitor monitor(cfg, At(0));
    for (int ms = 0; ms <= 1000; ms += 20)
        Feed(monitor, ms, Pose(0), 0, 0, 1, nullptr, &recorded);
    Require(SafetyValue(monitor.Evaluate(false, At(1000)), false, false, true) == 0,
            "SANY relaxed uncertainty limits accept the logged peak envelope");
    const auto below = covariance_for(0.369499, 0.49, 4.99);
    Feed(monitor, 1020, Pose(0), 0, 0, 1, nullptr, &below);
    Require(SafetyValue(monitor.Evaluate(false, At(1020)), false, false, true) == 0,
            "SANY 49 cm lateral and 4.99 degree yaw uncertainty pass");
    auto above = covariance_for(0.369499, 0.343525, 5.01);
    Feed(monitor, 1040, Pose(0), 0, 0, 1, nullptr, &above);
    Require(SafetyValue(monitor.Evaluate(false, At(1040)), false, false, true) == 4,
            "SANY 5.01 degree yaw uncertainty still stops");
    above = covariance_for(0.369499, 0.51, 4.324468);
    Feed(monitor, 1060, Pose(0), 0, 0, 1, nullptr, &above);
    Require(SafetyValue(monitor.Evaluate(false, At(1060)), false, false, true) == 4,
            "SANY 51 cm lateral uncertainty still stops");
    above = covariance_for(0.51, 0.343525, 4.324468);
    Feed(monitor, 1080, Pose(0), 0, 0, 1, nullptr, &above);
    Require(SafetyValue(monitor.Evaluate(false, At(1080)), false, false, true) == 4,
            "SANY major-axis uncertainty limit remains 50 cm");
    Feed(monitor, 1100, Pose(0), 0, 0, 1, nullptr, &recorded);
    Require(SafetyValue(monitor.Evaluate(false, At(1100)), false, false, true) == 0,
            "passing uncertainty clears its own confirmed failure");
    Require(SafetyValue(monitor.Evaluate(false, At(1601)), false, false, true) == 4,
            "SANY running output timeout remains a stop at 500 ms without claiming tracking loss");
    Require(SafetyValue(monitor.Evaluate(true, At(1100)), false, false, true) == -9,
            "SANY relaxed uncertainty does not mask confirmed tracking loss");
    Feed(monitor, 1120, Pose(2), 0, 0, 1, nullptr, &recorded);
    Require(SafetyValue(monitor.Evaluate(false, At(1120)), false, false, true) == 4,
            "SANY relaxed uncertainty does not mask a confirmed position jump");
}
}
int main(int argc, char** argv) {
    CheckNativeYawIntegration();
    if (argc == 2) {
        CheckSanyUncertaintyConfig(argv[1]);
        return 0;
    }
    Config cfg;
    CheckConfirmedLossPolicy();
    CheckQualityFaultSwitches();
    CheckMotionUnknownAdvisory();
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
    Require(SafetyValue(masked_unknown, false, true, false, true) == 0, "motion UNKNOWN alone does not stop");
    masked_unknown.observations[0].condition = kFail;
    Require(SafetyValue(masked_unknown, false, true, false, true) == 4, "mask preserves localization failure");
    const auto healthy_loc = snapshot;
    Require(snapshot.observations[0].condition == kPass && snapshot.observations[1].condition == kPass, "straight motion and valid covariance pass");
    Require(SafetyValue(snapshot, false, false, true) == 0, "all confirmed healthy gives zero");
    for (int expected : {0, 2, 4, -9}) {
        auto input = healthy_loc;
        const bool radar_fault = expected == 2;
        if (expected == 4) input.observations[1].condition = kFail;
        if (expected == -9) input.severe = true;
        const auto disabled = EvaluateSafetyDecision(input, false, radar_fault, !radar_fault, false, false);
        Require(disabled.raw == expected && disabled.selected == expected && disabled.published == 0,
                "disabled reporting preserves each original 2/4/-9 diagnosis for the log");
        const auto enabled = EvaluateSafetyDecision(input, false, radar_fault, !radar_fault, true, false);
        Require(enabled.raw == expected && enabled.published == expected,
                "reenabled reporting publishes the unchanged diagnosis");
    }
    const auto masked_radar = EvaluateSafetyDecision(healthy_loc, false, true, false, true, true);
    Require(masked_radar.raw == 2 && masked_radar.selected == 0 && masked_radar.published == 0,
            "radar-only masking preserves the original radar fault for the log");
    const auto disabled_severe = EvaluateSafetyDecision(healthy_loc, true, true, false, false, true);
    Require(disabled_severe.raw == -9 && disabled_severe.selected == -9 && disabled_severe.published == 0,
            "all-radar failure survives radar masking and remains auditable with global reporting off");
    Require(SafetyValue(snapshot, false, true, true) == 2 && SafetyValue(snapshot, false, false, false) == 2, "fault or stale radar gives two");
    snapshot.observations[0].condition = kFail;
    Require(SafetyValue(snapshot, false, true, false) == 4, "localization wins over radar");
    Require(SafetyValue(snapshot, true, true, false) == -9, "all-radar failure is severe");
    snapshot.severe = true;
    Require(SafetyValue(snapshot, false, false, true) == -9, "confirmed localization loss is severe");
    Require(normal.Evaluate(false, At(1501)).output_stale && !normal.Evaluate(false, At(1501)).severe,
            "running output timeout stops without inventing confirmed loss");
    Require(normal.Evaluate(true, At(1000)).severe, "confirmed tracking loss is severe immediately");

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
    Require(!reverse.Evaluate(false, At(1040)).covariance_valid &&
            std::isnan(reverse.Evaluate(false, At(1040)).major_sigma_m),
            "invalid covariance must not log sigma from the previous failing pose");
    reverse.ObservePose(101.06, Pose(-1.06), Pose(-1.06), SmallCov(), false, true, 1, At(1060));
    Require(reverse.Evaluate(false, At(1060)).observations[1].condition == kFail, "unobservable graph is an uncertainty fault");
    Require(reverse.Evaluate(false, At(1060)).pose_unobservable &&
            std::isnan(reverse.Evaluate(false, At(1060)).yaw_sigma_rad),
            "unobservable graph logs its explicit cause rather than stale sigma");
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
            short_hole.optional_can_unavailable &&
            SafetyValue(short_hole, false, false, true) == 0,
            "CAN-only integration hole remains UNKNOWN without stopping required healthy checks");
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
            grace.Evaluate(false, At(1160)).optional_can_unavailable &&
            SafetyValue(grace.Evaluate(false, At(1160)), false, false, true) == 0,
            "frozen optional CAN does not stop while independent yaw and uncertainty pass");
    Require(grace.Evaluate(false, At(1661)).output_stale && !grace.Evaluate(false, At(1661)).severe,
            "grace does not extend the 500 ms output timeout or turn it into confirmed loss");
    Feed(grace, 1180, Pose(1.18), 1, 0, 2);
    Require(!grace.Evaluate(false, At(1180)).consistency_unknown_tolerated,
            "epoch change invalidates prior passing evidence");
    bool invalid_grace_rejected = false;
    try {
        auto invalid = cfg;
        invalid.consistency_unknown_grace_sec = invalid.output_timeout_sec;
        Monitor rejected(invalid, At(0));
    } catch (const std::invalid_argument&) { invalid_grace_rejected = true; }
    Require(invalid_grace_rejected, "grace must remain below the output timeout");

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
    Require(SafetyValue(strict_gap.Evaluate(false, At(1107)), false, false, true) == 0,
            "motion UNKNOWN does not require a configured grace to avoid an independent stop");
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
            SafetyValue(bounded_gap.Evaluate(false, At(1151)), false, false, true) == 0,
            "the historical evidence deadline remains visible without independently stopping fresh outputs");
    Require(bounded_gap.Evaluate(false, At(1652)).output_stale &&
            !bounded_gap.Evaluate(false, At(1652)).severe && delayed.confirmed_failure == false,
            "output timeout is 500 ms and missing evidence is not confirmed loss");
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

    // CAN is optional even from startup or after a prolonged loss; the missing
    // independent translation check must stay UNKNOWN rather than become PASS.
    Monitor no_can(bounded_cfg, At(0));
    FeedWithoutCan(no_can, 0, Pose(0));
    Require(SafetyValue(no_can.Evaluate(false, At(0)), false, false, true) == 0,
            "motion-window warmup alone is advisory once valid localization output exists");
    for (int ms = 20; ms <= 60000; ms += 20) FeedWithoutCan(no_can, ms, Pose(ms * 0.005));
    const auto optional = no_can.Evaluate(false, At(60000));
    Require(optional.translation_condition == kUnknown && optional.yaw_condition == kPass &&
            optional.observations[0].condition == kUnknown && optional.optional_can_unavailable &&
            SafetyValue(optional, false, false, true) == 0,
            "long-term absent CAN allows 5 m/s motion with explicit unavailable translation evidence");
    Require(SafetyValue(optional, false, true, false) == 2 &&
            SafetyValue(optional, true, true, false) == -9,
            "optional CAN cannot mask radar faults");
    no_can.ObservePose(160.107, Pose(300.535), Pose(300.535), SmallCov(), true, false, 1, At(60107));
    Require(no_can.Evaluate(false, At(60107)).consistency_unknown_tolerated,
            "required IMU still uses bounded prior yaw evidence without CAN");
    no_can.ObservePose(160.151, Pose(300.755), Pose(300.755), SmallCov(), true, false, 1, At(60151));
    Require(SafetyValue(no_can.Evaluate(false, At(60151)), false, false, true) == 0,
            "expired IMU evidence remains UNKNOWN without an independent stop while output stays valid");
    Require(SafetyValue(no_can.Evaluate(false, At(60652)), false, false, true) == 4 &&
            SafetyValue(no_can.Evaluate(true, At(60151)), false, false, true) == -9,
            "optional CAN never masks output timeout or confirmed tracking loss");

    Monitor no_can_yaw(cfg, At(0));
    for (int ms = 0; ms <= 1000; ms += 20) FeedWithoutCan(no_can_yaw, ms, Pose(ms * .001));
    const SE3 raw_yaw_jump = Pose(1.02, 0, 0.2);
    FeedWithoutCan(no_can_yaw, 1020, Pose(1.02), 0, &raw_yaw_jump);
    Require(no_can_yaw.Evaluate(false, At(1020)).yaw_condition == kFail &&
            SafetyValue(no_can_yaw.Evaluate(false, At(1020)), false, false, true) == 4,
            "raw yaw mismatch still stops with CAN absent");
    for (int ms = 1040; ms <= 1800; ms += 20) FeedWithoutCan(no_can_yaw, ms, Pose(ms * .001));
    Require(SafetyValue(no_can_yaw.Evaluate(false, At(1800)), false, false, true) == 0,
            "independently verified yaw recovery can clear its own failure without CAN");
    Mat6d bad_optional_cov = SmallCov(); bad_optional_cov(0, 0) = 1;
    no_can_yaw.ObservePose(101.82, Pose(1.82), Pose(1.82), bad_optional_cov, true, false, 1, At(1820));
    Require(SafetyValue(no_can_yaw.Evaluate(false, At(1820)), false, false, true) == 4,
            "uncertainty failure still stops without CAN");

    Monitor dropout(cfg, At(0));
    for (int ms = 0; ms <= 1000; ms += 20) Feed(dropout, ms, Pose(ms * .001));
    Feed(dropout, 1020, Pose(3.02));
    for (int ms = 1040; ms <= 2200; ms += 20) FeedWithoutCan(dropout, ms, Pose(2 + ms * .001));
    Require(dropout.Evaluate(false, At(2200)).confirmed_failure &&
            SafetyValue(dropout.Evaluate(false, At(2200)), false, false, true) == 4,
            "CAN loss and fresh yaw cannot clear a confirmed translation failure");
    Feed(dropout, 2220, Pose(4.22));
    Require(dropout.Evaluate(false, At(2220)).confirmed_failure,
            "one recovered CAN sample cannot certify a missing translation window");
    for (int ms = 2240; ms <= 2900; ms += 20) Feed(dropout, ms, Pose(2 + ms * .001));
    Require(SafetyValue(dropout.Evaluate(false, At(2900)), false, false, true) == 0 &&
            dropout.Evaluate(false, At(2900)).translation_condition == kPass,
            "complete recovered CAN evidence clears only the verified translation failure");
    Feed(dropout, 2920, Pose(8.92));
    Require(SafetyValue(dropout.Evaluate(false, At(2920)), false, false, true) == 4,
            "translation checking resumes after CAN recovers");

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
