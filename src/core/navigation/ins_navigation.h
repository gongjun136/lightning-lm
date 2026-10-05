#pragma once

#include <array>
#include <deque>
#include <string>
#include <GeographicLib/LocalCartesian.hpp>
#include "common/eigen_types.h"

namespace lightning::ins {

// All device outputs refer to the configured rear axle. Never apply the
// GNSS/IMU lever arm again in this layer.
enum Field : std::size_t { Latitude, Longitude, Altitude, Attitude, Velocity,
                          PositionSigma, AttitudeSigma, Status, VelocitySigma, FieldCount };
struct Sample {
    std::array<double, FieldCount> time{};       // corrected ROS-domain receive stamps
    std::array<double, FieldCount> arrival{};    // monotonic receive times
    std::array<bool, FieldCount> valid{};
    double latitude = 0, longitude = 0, altitude = 0;
    double heading = 0, pitch = 0, roll = 0;
    Vec3d velocity = Vec3d::Zero();              // local ENU, m/s
    Vec3d position_sigma = Vec3d::Zero();        // E/N/U, m
    Vec3d attitude_sigma = Vec3d::Zero();        // heading/pitch/roll, deg
    Vec3d velocity_sigma = Vec3d::Zero();
    int system_state = 0, satellite_status = 0;
    double differential_age = 0;
    double Stamp() const { return time[Latitude]; }
};

struct QualityPolicy {
    double max_age = 0.15;
    double max_skew = 0.02;
    double max_future = 0.02;
    double max_differential_age = 2.0;
    double max_horizontal_sigma = 0.10;
    double max_vertical_sigma = 0.20;
    double max_attitude_sigma = 0.5;
    double max_velocity_sigma = 0.20;
    int recovery_samples = 10;
};

std::string CheckSample(const Sample& sample, const QualityPolicy& policy,
                        double ros_now, double steady_now);

class RecoveryGate {
 public:
    explicit RecoveryGate(int required) : required_(required) {}
    void Reject() { count_ = 0; good_ = false; }
    bool Observe(double stamp, bool valid);
    bool Good() const { return good_; }
 private:
    int required_, count_ = 0;
    double last_stamp_ = 0;
    bool good_ = false;
};

class GeoReference {
 public:
    GeoReference(double latitude, double longitude, double ellipsoid_height);
    Vec3d Forward(double latitude, double longitude, double ellipsoid_height,
                  Mat3d* fixed_from_local_enu = nullptr) const;
    Vec3d Reverse(const Vec3d& enu) const;
 private:
    GeographicLib::LocalCartesian local_;
};

// CHC: north-clockwise heading, nose-up pitch, right-down roll.
// Result maps rear-body FLU axes to the local ENU axes at the measurement.
SO3 VehicleAttitude(double heading_deg, double pitch_deg, double roll_deg);

struct Pose {
    double stamp = 0;
    SE3 map_from_rear;
    Vec3d map_velocity = Vec3d::Zero();
    bool global_valid = false;
    double SignedSpeed() const { return (map_from_rear.so3().inverse() * map_velocity).x(); }
};
Pose Convert(const Sample& sample, const GeoReference& reference);

// No extrapolation; an interval must be covered by finite, monotonic samples.
// A generation break prevents interpolation over a device state/quality reset.
class PoseBuffer {
 public:
    PoseBuffer(double horizon, double max_gap) : horizon_(horizon), max_gap_(max_gap) {}
    bool Add(const Pose& pose);
    void Clear() { poses_.clear(); }
    bool At(double stamp, Pose& result, bool require_global = true) const;
    bool Covers(double begin, double end, bool require_global = true) const;
    double Latest() const { return poses_.empty() ? 0 : poses_.back().stamp; }
 private:
    std::deque<Pose> poses_;
    double horizon_, max_gap_;
};

}  // namespace lightning::ins
