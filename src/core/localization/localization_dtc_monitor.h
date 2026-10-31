#pragma once
#include <array>
#include <chrono>
#include <cstdint>
#include <deque>
#include <limits>
#include <mutex>
#include "common/eigen_types.h"

namespace lightning::loc::dtc {
using Clock = std::chrono::steady_clock;
enum Condition : uint8_t { kUnknown = 0, kPass = 1, kFail = 2 };
constexpr uint32_t kConsistency = 0x30811813;
constexpr uint32_t kUncertainty = 0x40450198;
constexpr const char* kSafetyPolicyVersion = "localization-dtc-v5-confirmed-loss";
struct Config {
    double window_sec = 0.5, max_input_gap_sec = 0.1;
    double longitudinal_m = 0.5, lateral_m = 0.3, yaw_rad = 2.0 * M_PI / 180.0;
    double distance_ratio = 0.1, angle_ratio = 0.1;
    double major_sigma_m = 0.5, lateral_sigma_m = 0.3, yaw_sigma_rad = 2.0 * M_PI / 180.0;
    double covariance_scale = 1.0;  ///< Additional variance inflation; calibrate against independent error evidence.
    double position_sigma_floor_m = 0.10, angle_sigma_floor_rad = 0.5 * M_PI / 180.0;
    double output_timeout_sec = 0.5, evaluation_timeout_sec = 0.5, startup_grace_sec = 8.0;
    double consistency_unknown_grace_sec = 0.0;  ///< Disabled by default; maximum age of required passing yaw evidence, steady seconds.
    double max_speed_mps = 30.0, max_yaw_rate_radps = 3.0;
};
struct Observation {
    uint32_t code = 0, sequence = 0;
    Condition condition = kUnknown;
    uint64_t age_us = std::numeric_limits<uint64_t>::max();
};
struct TrackingEvidence {
    bool available = true;
    bool initialized = false;  ///< At least one accepted map match, independent of pose publication.
    bool confirmed_lost = false;  ///< Explicit tracking failure/relocalization, never freshness alone.
    bool output_gate_open = false, match_stale = false, latest_match_valid = false;
    double valid_match_age_sec = std::numeric_limits<double>::infinity();
};
struct Snapshot {
    std::array<Observation, 2> observations;
    bool severe = false;
    bool confirmed_failure = false;  ///< A definite failure remains active until the corresponding check passes again.
    bool consistency_unknown_tolerated = false;  ///< Only motion UNKNOWN may use bounded prior passing evidence.
    bool optional_can_unavailable = false;  ///< Translation UNKNOWN alone does not stop when required yaw evidence passes.
    Condition translation_condition = kUnknown, yaw_condition = kUnknown;
    double output_age_sec = 0.0;  ///< Steady time since the last accepted output.
    double speed_arrival_age_sec = 0.0, gyro_arrival_age_sec = 0.0;
    double consistency_pass_age_sec = 0.0;  ///< Steady age of last PASS; infinity before first PASS in this epoch.
    double yaw_pass_age_sec = 0.0;  ///< Age of required IMU/pose evidence, independent of optional CAN.
    bool output_stale = false, localization_lost = false;  ///< Only confirmed loss is severe during operation.
    bool startup_pending = true, startup_expired = false;
    double startup_age_sec = 0.0;
    bool match_evidence_available = false, output_gate_open = false, match_stale = false, latest_match_valid = false;
    double valid_match_age_sec = std::numeric_limits<double>::infinity();
    bool output_continuous = false;  ///< Advancing output over the configured motion window without a timeout/epoch jump.
    double error_longitudinal_m = 0, error_lateral_m = 0, error_yaw_rad = 0;
    bool covariance_valid = false, pose_unobservable = false;
    double major_sigma_m = std::numeric_limits<double>::quiet_NaN();
    double lateral_sigma_m = std::numeric_limits<double>::quiet_NaN();
    double yaw_sigma_rad = std::numeric_limits<double>::quiet_NaN();
    double speed_stamp_lag_sec = std::numeric_limits<double>::quiet_NaN();  ///< Latest pose minus latest CAN stamp; NaN if absent.
    double gyro_stamp_lag_sec = std::numeric_limits<double>::quiet_NaN();  ///< Latest pose minus latest IMU stamp; NaN if absent.
};
// Raw CAN speed and body-frame IMU angular rate are independent of fused pose velocity.
class Monitor {
 public:
    explicit Monitor(Config config, Clock::time_point now = Clock::now());
    void ObserveSpeed(double stamp, double signed_speed, Clock::time_point now = Clock::now());
    void ObserveGyro(double stamp, double body_yaw_rate, Clock::time_point now = Clock::now());
    void ObservePose(double stamp, const SE3& rear_pose, const SE3& raw_rear_pose,
                     const Mat6d& covariance, bool covariance_valid, bool unobservable,
                     uint64_t epoch, Clock::time_point now = Clock::now());
    Snapshot Evaluate(bool lost, Clock::time_point now = Clock::now()) const;
    Snapshot Evaluate(const TrackingEvidence& tracking, Clock::time_point now = Clock::now()) const;
 private:
    struct Sample { double stamp, value; Clock::time_point arrival; };
    struct Pose { double stamp; SE3 final, raw; };
    struct MotionResult { Condition translation = kUnknown, yaw = kUnknown; };
    void Observe(std::deque<Sample>& queue, double stamp, double value, double limit, Clock::time_point now);
    bool Interpolate(const std::deque<Sample>& queue, double stamp, double& value, Clock::time_point now) const;
    MotionResult Motion(const Pose& from, const Pose& to, Clock::time_point now);
    void Uncertainty(const SE3& pose, const Mat6d& covariance, bool valid, bool unobservable);
    Snapshot EvaluateLocked(const TrackingEvidence& tracking, Clock::time_point now) const;
    Config config_;
    mutable std::mutex mutex_;
    std::deque<Sample> speeds_, gyros_;
    std::deque<Pose> poses_;
    Clock::time_point start_, output_time_, evaluation_time_, continuous_output_start_;
    uint64_t epoch_ = 0;
    bool output_seen_ = false;
    std::array<bool, 3> failed_checks_{{false, false, false}};  ///< Translation, yaw, uncertainty; UNKNOWN never clears failure.
    bool consistency_pass_seen_ = false;
    Clock::time_point consistency_pass_time_{};
    bool yaw_pass_seen_ = false;
    Clock::time_point yaw_pass_time_{};
    Snapshot snapshot_;
};
// Severe > localization > lidar > normal. An explicit radar-only override
// suppresses the lidar branch (including UNKNOWN), never the first two branches.
int8_t SafetyValue(const Snapshot& localization, bool all_lidar_failed,
                   bool radar_fault, bool radar_all_pass, bool mask_radar_fault = false);
struct SafetyDecision {
    int8_t raw = 0, selected = 0, published = 0;
    bool partial_radar_suppressed = false;
};
bool HealthyLocalizationForDegradedRadar(const Snapshot& localization);
// Preserve the diagnosis before either commissioning override is applied.
SafetyDecision EvaluateSafetyDecision(const Snapshot& localization, bool all_lidar_failed,
                                     bool radar_fault, bool radar_all_pass,
                                     bool reporting_enabled, bool mask_radar_fault);
}
