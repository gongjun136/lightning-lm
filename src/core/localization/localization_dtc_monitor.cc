#include "core/localization/localization_dtc_monitor.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lightning::loc::dtc {
namespace {
double Yaw(const SE3& pose) { return std::atan2(pose.rotationMatrix()(1, 0), pose.rotationMatrix()(0, 0)); }
double Wrap(double a) { return std::atan2(std::sin(a), std::cos(a)); }
double Age(Clock::time_point now, Clock::time_point time) { return std::chrono::duration<double>(now - time).count(); }
Condition Combined(Condition translation, Condition yaw) {
    if (translation == kFail || yaw == kFail) return kFail;
    return translation == kPass && yaw == kPass ? kPass : kUnknown;
}
}
Monitor::Monitor(Config config, Clock::time_point now) : config_(config), start_(now) {
    for (double v : {config.window_sec, config.max_input_gap_sec, config.longitudinal_m,
                    config.lateral_m, config.yaw_rad, config.major_sigma_m, config.lateral_sigma_m,
                    config.yaw_sigma_rad, config.output_timeout_sec, config.evaluation_timeout_sec,
                    config.startup_grace_sec, config.max_speed_mps, config.max_yaw_rate_radps})
        if (!std::isfinite(v) || v <= 0) throw std::invalid_argument("DTC thresholds must be finite and positive");
    if (!std::isfinite(config.covariance_scale) || config.covariance_scale < 1.0 ||
        !std::isfinite(config.consistency_unknown_grace_sec) ||
        config.consistency_unknown_grace_sec < 0 ||
        config.consistency_unknown_grace_sec >= config.output_timeout_sec ||
        !std::isfinite(config.distance_ratio) || config.distance_ratio < 0 ||
        !std::isfinite(config.angle_ratio) || config.angle_ratio < 0 ||
        !std::isfinite(config.position_sigma_floor_m) || config.position_sigma_floor_m < 0 ||
        !std::isfinite(config.angle_sigma_floor_rad) || config.angle_sigma_floor_rad < 0)
        throw std::invalid_argument("invalid DTC uncertainty scale or tolerance");
    snapshot_.observations[0].code = kConsistency;
    snapshot_.observations[1].code = kUncertainty;
}
void Monitor::Observe(std::deque<Sample>& q, double stamp, double value, double limit, Clock::time_point now) {
    if (!std::isfinite(stamp) || stamp <= 0 || !std::isfinite(value) || std::abs(value) > limit) {
        if (&q == &gyros_) yaw_pass_seen_ = false;
        q.clear(); consistency_pass_seen_ = false; return;
    }
    if (!q.empty() && stamp <= q.back().stamp) return;  // frozen/replayed input cannot refresh validity
    q.push_back({stamp, value, now});
    while (!q.empty() && (stamp - q.front().stamp > 2 * config_.window_sec + 1.0 || q.size() > 10000)) q.pop_front();
}
void Monitor::ObserveSpeed(double stamp, double value, Clock::time_point now) {
    std::lock_guard<std::mutex> lock(mutex_); Observe(speeds_, stamp, value, config_.max_speed_mps, now);
}
void Monitor::ObserveGyro(double stamp, double value, Clock::time_point now) {
    std::lock_guard<std::mutex> lock(mutex_); Observe(gyros_, stamp, value, config_.max_yaw_rate_radps, now);
}
bool Monitor::Interpolate(const std::deque<Sample>& q, double stamp, double& value, Clock::time_point now) const {
    if (q.empty() || Age(now, q.back().arrival) > config_.max_input_gap_sec || stamp < q.front().stamp) return false;
    auto upper = std::lower_bound(q.begin(), q.end(), stamp, [](const Sample& s, double t) { return s.stamp < t; });
    if (upper == q.end()) {
        if (stamp - q.back().stamp > config_.max_input_gap_sec) return false;
        value = q.back().value; return true;  // bounded hold at the live input boundary
    }
    if (upper->stamp == stamp) { value = upper->value; return true; }
    if (upper == q.begin()) return false;
    auto lower = std::prev(upper);
    if (upper->stamp - lower->stamp > config_.max_input_gap_sec) return false;
    const double ratio = (stamp - lower->stamp) / (upper->stamp - lower->stamp);
    value = lower->value + ratio * (upper->value - lower->value);
    return std::isfinite(value);
}
Monitor::MotionResult Monitor::Motion(const Pose& from, const Pose& to, Clock::time_point now) {
    const double dt = to.stamp - from.stamp;
    if (dt <= 0 || dt > 2 * config_.window_sec + config_.max_input_gap_sec) return {};
    double dx = 0, dy = 0, theta = 0, distance = 0, turn = 0;
    bool translation_available = true;
    const int steps = std::max(1, static_cast<int>(std::ceil(dt / std::min(0.02, config_.max_input_gap_sec / 4))));
    const double step = dt / steps;
    for (int i = 0; i < steps; ++i) {
        double speed = 0, omega;
        const double stamp = from.stamp + (i + 0.5) * step;
        if (!Interpolate(gyros_, stamp, omega, now)) return {};
        translation_available &= Interpolate(speeds_, stamp, speed, now);
        dx += speed * step * std::cos(theta + omega * step / 2);
        dy += speed * step * std::sin(theta + omega * step / 2);
        theta += omega * step;
        distance += std::abs(speed) * step;
        turn += std::abs(omega) * step;
    }
    bool translation_fail = false, yaw_fail = false;
    for (bool raw : {false, true}) {
        const SE3& base = raw ? from.raw : from.final;
        const SE3& current = raw ? to.raw : to.final;
        const double yaw = Yaw(base);
        const Vec3d delta = current.translation() - base.translation();
        const double x = std::cos(yaw) * delta.x() + std::sin(yaw) * delta.y();
        const double y = -std::sin(yaw) * delta.x() + std::cos(yaw) * delta.y();
        const double ex = std::abs(x - dx), ey = std::abs(y - dy);
        const double ea = std::abs(Wrap(Yaw(current) - yaw - theta));
        if (translation_available) {
            snapshot_.error_longitudinal_m = std::max(snapshot_.error_longitudinal_m, ex);
            snapshot_.error_lateral_m = std::max(snapshot_.error_lateral_m, ey);
            translation_fail |= ex > config_.longitudinal_m + config_.distance_ratio * distance ||
                                ey > config_.lateral_m + config_.distance_ratio * distance;
        }
        snapshot_.error_yaw_rad = std::max(snapshot_.error_yaw_rad, ea);
        yaw_fail |= ea > config_.yaw_rad + config_.angle_ratio * turn;
    }
    return {translation_available ? (translation_fail ? kFail : kPass) : kUnknown,
            yaw_fail ? kFail : kPass};
}
void Monitor::Uncertainty(const SE3& pose, const Mat6d& p, bool valid, bool unobservable) {
    auto& condition = snapshot_.observations[1].condition;
    condition = kUnknown;
    snapshot_.covariance_valid = false;
    snapshot_.pose_unobservable = unobservable;
    // An unavailable covariance must not display sigma values from an older pose.
    snapshot_.major_sigma_m = snapshot_.lateral_sigma_m = snapshot_.yaw_sigma_rad =
        std::numeric_limits<double>::quiet_NaN();
    if (unobservable) { condition = kFail; return; }
    if (!valid || !p.allFinite() || (p - p.transpose()).norm() > 1e-6 * std::max(1.0, p.norm())) return;
    Eigen::SelfAdjointEigenSolver<Mat6d> check(p);
    if (check.info() != Eigen::Success || check.eigenvalues().minCoeff() <= 1e-12) return;
    Mat6d covariance = config_.covariance_scale * p;
    covariance.diagonal().head<3>().array() += config_.position_sigma_floor_m * config_.position_sigma_floor_m;
    covariance.diagonal().tail<3>().array() += config_.angle_sigma_floor_rad * config_.angle_sigma_floor_rad;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> ellipse(covariance.topLeftCorner<2, 2>());
    if (ellipse.info() != Eigen::Success) return;
    const double yaw = Yaw(pose);
    Eigen::Vector2d lateral(-std::sin(yaw), std::cos(yaw));
    // Numerical yaw Jacobian in the RIGHT rotation tangent; works with roll/pitch.
    Vec3d yaw_j;
    if (std::hypot(pose.rotationMatrix()(0, 0), pose.rotationMatrix()(1, 0)) < 1e-6) return;
    for (int i = 0; i < 3; ++i) {
        Vec3d d = Vec3d::Zero(); d[i] = 1e-6;
        const SE3 plus(pose.so3() * SO3::exp(d), pose.translation());
        const SE3 minus(pose.so3() * SO3::exp(-d), pose.translation());
        yaw_j[i] = Wrap(Yaw(plus) - Yaw(minus)) / 2e-6;
    }
    snapshot_.major_sigma_m = std::sqrt(std::max(0.0, ellipse.eigenvalues().maxCoeff()));
    snapshot_.lateral_sigma_m = std::sqrt(std::max(0.0, lateral.dot(covariance.topLeftCorner<2, 2>() * lateral)));
    snapshot_.yaw_sigma_rad = std::sqrt(std::max(0.0, yaw_j.dot(covariance.bottomRightCorner<3, 3>() * yaw_j)));
    snapshot_.covariance_valid = true;
    condition = snapshot_.major_sigma_m > config_.major_sigma_m || snapshot_.lateral_sigma_m > config_.lateral_sigma_m ||
                snapshot_.yaw_sigma_rad > config_.yaw_sigma_rad ? kFail : kPass;
}
void Monitor::ObservePose(double stamp, const SE3& pose, const SE3& raw, const Mat6d& covariance,
                          bool valid, bool unobservable, uint64_t epoch, Clock::time_point now) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!std::isfinite(stamp) || stamp <= 0 || !pose.matrix().allFinite() || !raw.matrix().allFinite()) return;
    if (epoch != epoch_) {
        poses_.clear(); epoch_ = epoch;
        consistency_pass_seen_ = false;
        yaw_pass_seen_ = false;
    }
    if (!poses_.empty() && stamp <= poses_.back().stamp) return;
    if (!poses_.empty()) {
        // A delivery gap alone is not a clock jump. Compare source progress
        // with steady elapsed time before reusing any earlier passing evidence.
        const double stamp_step = stamp - poses_.back().stamp;
        const double arrival_step = Age(now, output_time_);
        if (std::abs(stamp_step - arrival_step) > config_.max_input_gap_sec) {
            consistency_pass_seen_ = false;
            yaw_pass_seen_ = false;
        }
    }
    output_time_ = evaluation_time_ = now;
    output_seen_ = true;
    for (auto& obs : snapshot_.observations) ++obs.sequence;
    snapshot_.error_longitudinal_m = snapshot_.error_lateral_m = snapshot_.error_yaw_rad = 0;
    const Pose current{stamp, pose, raw};
    MotionResult motion;
    // Check adjacent increments as well, so a jump is detected during window warmup.
    if (!poses_.empty()) motion = Motion(poses_.back(), current, now);
    auto anchor = poses_.end();
    for (auto it = poses_.begin(); it != poses_.end() && it->stamp <= stamp - config_.window_sec; ++it) anchor = it;
    if (anchor != poses_.end()) {
        const auto window = Motion(*anchor, current, now);
        if (motion.translation != kFail) motion.translation = window.translation;
        if (motion.yaw != kFail) motion.yaw = window.yaw;
    } else {
        if (motion.translation != kFail) motion.translation = kUnknown;
        if (motion.yaw != kFail) motion.yaw = kUnknown;
    }
    snapshot_.translation_condition = motion.translation;
    snapshot_.yaw_condition = motion.yaw;
    const Condition condition = Combined(motion.translation, motion.yaw);
    snapshot_.observations[0].condition = condition;
    Uncertainty(pose, covariance, valid, unobservable);
    // [dtc-confirmed-failure-latch]
    const std::array<Condition, 3> checks{{motion.translation, motion.yaw, snapshot_.observations[1].condition}};
    for (size_t i = 0; i < failed_checks_.size(); ++i) {
        if (checks[i] == kFail) failed_checks_[i] = true;
        if (checks[i] == kPass) failed_checks_[i] = false;
    }
    // [dtc-confirmed-failure-latch]
    if (condition == kPass) {
        consistency_pass_seen_ = true;
        consistency_pass_time_ = now;
    }
    if (motion.yaw == kPass) {
        yaw_pass_seen_ = true;
        yaw_pass_time_ = now;
    }
    poses_.push_back(current);
    while (poses_.size() > 1 && poses_[1].stamp < stamp - 2 * config_.window_sec) poses_.pop_front();
}
Snapshot Monitor::Evaluate(bool lost, Clock::time_point now) const {
    std::lock_guard<std::mutex> lock(mutex_);
    Snapshot out = snapshot_;
    if (!poses_.empty()) {
        if (!speeds_.empty()) out.speed_stamp_lag_sec = poses_.back().stamp - speeds_.back().stamp;
        if (!gyros_.empty()) out.gyro_stamp_lag_sec = poses_.back().stamp - gyros_.back().stamp;
    }
    // [dtc-output-severity]
    const bool stale = output_seen_ ? Age(now, output_time_) > config_.output_timeout_sec : Age(now, start_) > config_.startup_grace_sec;
    out.severe = lost || stale;
    out.output_stale = stale;
    out.localization_lost = lost;
    // [dtc-output-severity]
    const double missing = std::numeric_limits<double>::infinity();
    out.output_age_sec = output_seen_ ? Age(now, output_time_) : missing;
    out.speed_arrival_age_sec = speeds_.empty() ? missing : Age(now, speeds_.back().arrival);
    out.gyro_arrival_age_sec = gyros_.empty() ? missing : Age(now, gyros_.back().arrival);
    out.consistency_pass_age_sec = consistency_pass_seen_ ? Age(now, consistency_pass_time_) : missing;
    out.yaw_pass_age_sec = yaw_pass_seen_ ? Age(now, yaw_pass_time_) : missing;
    out.confirmed_failure = failed_checks_[0] || failed_checks_[1] || failed_checks_[2];
    for (auto& obs : out.observations) {
        obs.age_us = output_seen_ ? static_cast<uint64_t>(std::max(0.0, Age(now, evaluation_time_)) * 1e6) : std::numeric_limits<uint64_t>::max();
        if (obs.condition != kFail &&
            (!output_seen_ || Age(now, evaluation_time_) > config_.evaluation_timeout_sec || lost))
            obs.condition = kUnknown;
    }
    const auto expire = [](Condition& condition) { if (condition != kFail) condition = kUnknown; };
    if (!output_seen_ || Age(now, evaluation_time_) > config_.evaluation_timeout_sec || lost) {
        expire(out.translation_condition);
        expire(out.yaw_condition);
    }
    // Missing CAN disables translation evidence only. IMU remains required.
    if (out.speed_arrival_age_sec > config_.max_input_gap_sec) expire(out.translation_condition);
    if (out.gyro_arrival_age_sec > config_.max_input_gap_sec) {
        expire(out.translation_condition);
        expire(out.yaw_condition);
    }
    out.observations[0].condition = Combined(out.translation_condition, out.yaw_condition);
    // [dtc-optional-can-evidence]
    out.optional_can_unavailable = !out.severe && !out.confirmed_failure &&
        out.observations[1].condition == kPass && out.translation_condition == kUnknown &&
        out.yaw_condition == kPass;
    // [dtc-optional-can-evidence]
    // Do not report UNKNOWN as PASS. Grace affects only the stop decision and
    // cannot refresh itself, survive an epoch change, hide a failure, stale
    // output, missing uncertainty, or a detected clock discontinuity. The
    // integration gap limit decides whether evidence is complete; the separate
    // grace is a TOTAL evidence-age limit, not extra time after that gap limit.
    // [dtc-unknown-evidence-age]
    const double grace = config_.consistency_unknown_grace_sec;
    out.consistency_unknown_tolerated =
        grace > 0 && !out.severe && !out.confirmed_failure &&
        out.observations[0].condition == kUnknown &&
        out.observations[1].condition == kPass &&
        out.yaw_pass_age_sec <= grace &&
        out.output_age_sec <= grace &&
        out.gyro_arrival_age_sec <= grace &&
        out.yaw_condition == kUnknown &&
        std::abs(out.gyro_stamp_lag_sec) <= grace;
    // [dtc-unknown-evidence-age]
    return out;
}
// [dtc-safety-priority]
int8_t SafetyValue(const Snapshot& loc, bool all_failed, bool radar_fault, bool radar_all_pass, bool mask_radar_fault) {
    if (loc.severe || all_failed) return -9;
    if (loc.confirmed_failure || loc.observations[0].condition == kFail ||
        loc.observations[1].condition != kPass) return 4;
    if (loc.observations[0].condition != kPass && !loc.consistency_unknown_tolerated &&
        !loc.optional_can_unavailable) return 4;
    if (!mask_radar_fault && (radar_fault || !radar_all_pass)) return 2;
    return 0;
}
// [dtc-safety-priority]
SafetyDecision EvaluateSafetyDecision(const Snapshot& localization, bool all_failed,
                                     bool radar_fault, bool radar_all_pass,
                                     bool reporting_enabled, bool mask_radar_fault) {
    SafetyDecision decision;
    decision.raw = SafetyValue(localization, all_failed, radar_fault, radar_all_pass);
    decision.selected = SafetyValue(localization, all_failed, radar_fault, radar_all_pass, mask_radar_fault);
    decision.published = reporting_enabled ? decision.selected : 0;
    return decision;
}
}
