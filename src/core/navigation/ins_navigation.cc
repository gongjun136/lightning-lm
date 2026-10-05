#include "core/navigation/ins_navigation.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace lightning::ins {
namespace {
bool PositiveFinite(double v) { return std::isfinite(v) && v > 0; }
bool Sigma(const Vec3d& v) { return v.allFinite() && (v.array() >= 0).all(); }
}

std::string CheckSample(const Sample& s, const QualityPolicy& p, double now, double steady) {
    if (!std::isfinite(now) || !std::isfinite(steady)) return "invalid_clock";
    double earliest = s.time[0], latest = s.time[0];
    for (std::size_t i = 0; i < FieldCount; ++i) {
        if (!s.valid[i] || !PositiveFinite(s.time[i]) || !PositiveFinite(s.arrival[i]))
            return "missing_or_invalid_field_" + std::to_string(i);
        if (now - s.time[i] > p.max_age || steady - s.arrival[i] > p.max_age)
            return "stale_field_" + std::to_string(i);
        if (s.time[i] - now > p.max_future || steady < s.arrival[i])
            return "clock_mismatch_" + std::to_string(i);
        earliest = std::min(earliest, s.time[i]); latest = std::max(latest, s.time[i]);
    }
    if (latest - earliest > p.max_skew) return "field_time_skew";
    if (s.system_state != 2) return "not_integrated_navigation";
    if (s.satellite_status != 4) return "not_rtk_fixed_and_heading_valid";
    if (!std::isfinite(s.latitude) || std::abs(s.latitude) > 90 ||
        !std::isfinite(s.longitude) || std::abs(s.longitude) > 180 ||
        !std::isfinite(s.altitude) || !std::isfinite(s.heading) || s.heading < 0 || s.heading >= 360 ||
        !std::isfinite(s.pitch) || std::abs(s.pitch) > 90 ||
        !std::isfinite(s.roll) || std::abs(s.roll) > 180 || !s.velocity.allFinite())
        return "invalid_navigation_value";
    if (!std::isfinite(s.differential_age) || s.differential_age < 0 ||
        s.differential_age > p.max_differential_age) return "differential_age";
    if (!Sigma(s.position_sigma) || !Sigma(s.attitude_sigma) || !Sigma(s.velocity_sigma))
        return "invalid_sigma";
    if (s.position_sigma.head<2>().maxCoeff() > p.max_horizontal_sigma ||
        s.position_sigma.z() > p.max_vertical_sigma) return "position_sigma";
    if (s.attitude_sigma.maxCoeff() > p.max_attitude_sigma) return "attitude_sigma";
    if (s.velocity_sigma.maxCoeff() > p.max_velocity_sigma) return "velocity_sigma";
    return {};
}

bool RecoveryGate::Observe(double stamp, bool valid) {
    if (!valid || !PositiveFinite(stamp) || stamp <= last_stamp_) { Reject(); return false; }
    last_stamp_ = stamp;
    count_ = std::min(required_, count_ + 1);
    return good_ = count_ >= required_;
}

GeoReference::GeoReference(double lat, double lon, double height) : local_(lat, lon, height) {
    if (!std::isfinite(lat) || std::abs(lat) >= 90 || !std::isfinite(lon) ||
        std::abs(lon) > 180 || !std::isfinite(height))
        throw std::invalid_argument("invalid fixed geodetic origin (poles unsupported)");
}

Vec3d GeoReference::Forward(double lat, double lon, double height, Mat3d* rotation) const {
    if (!std::isfinite(lat) || std::abs(lat) >= 90 || !std::isfinite(lon) ||
        std::abs(lon) > 180 || !std::isfinite(height)) throw std::invalid_argument("invalid LLH");
    Vec3d p;
    std::vector<double> m(9);
    local_.Forward(lat, lon, height, p.x(), p.y(), p.z(), m);
    if (rotation) for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) (*rotation)(r,c) = m[3*r+c];
    return p;
}
Vec3d GeoReference::Reverse(const Vec3d& p) const {
    if (!p.allFinite()) throw std::invalid_argument("invalid ENU");
    Vec3d llh;
    local_.Reverse(p.x(), p.y(), p.z(), llh.x(), llh.y(), llh.z());
    return llh;
}

SO3 VehicleAttitude(double heading, double pitch, double roll) {
    const double rad = M_PI / 180;
    // FLU positive pitch is nose-down, positive roll is left-up/right-down.
    return SO3(Quatd(Eigen::AngleAxisd((90 - heading)*rad, Vec3d::UnitZ()) *
                     Eigen::AngleAxisd(-pitch*rad, Vec3d::UnitY()) *
                     Eigen::AngleAxisd(roll*rad, Vec3d::UnitX())));
}
Pose Convert(const Sample& s, const GeoReference& ref) {
    Mat3d rotation;
    const Vec3d p = ref.Forward(s.latitude, s.longitude, s.altitude, &rotation);
    Pose result;
    result.stamp = s.Stamp();
    result.map_from_rear = SE3(SO3(rotation) * VehicleAttitude(s.heading, s.pitch, s.roll), p);
    result.map_velocity = rotation * s.velocity;
    return result;
}

bool PoseBuffer::Add(const Pose& p) {
    if (!PositiveFinite(p.stamp) || !p.map_from_rear.matrix().allFinite() ||
        (!poses_.empty() && p.stamp <= poses_.back().stamp)) return false;
    poses_.push_back(p);
    while (poses_.size() > 2 && poses_[1].stamp < p.stamp - horizon_) poses_.pop_front();
    // A hard bound also protects against malformed stamps with tiny increments.
    while (poses_.size() > 10000) poses_.pop_front();
    return true;
}
bool PoseBuffer::At(double t, Pose& result, bool global) const {
    if (!std::isfinite(t) || poses_.empty() || t < poses_.front().stamp || t > poses_.back().stamp)
        return false;
    auto hi = std::lower_bound(poses_.begin(), poses_.end(), t,
                               [](const Pose& p, double v) { return p.stamp < v; });
    if (hi == poses_.end()) return false;
    if (hi->stamp == t) { result = *hi; return !global || result.global_valid; }
    if (hi == poses_.begin()) return false;
    const auto& lo = *std::prev(hi);
    if (hi->stamp - lo.stamp > max_gap_ || (global && (!lo.global_valid || !hi->global_valid))) return false;
    const double alpha = (t - lo.stamp) / (hi->stamp - lo.stamp);
    result.stamp = t;
    result.map_from_rear = SE3(lo.map_from_rear.unit_quaternion().slerp(alpha, hi->map_from_rear.unit_quaternion()),
                               (1-alpha)*lo.map_from_rear.translation() + alpha*hi->map_from_rear.translation());
    result.map_velocity = (1-alpha)*lo.map_velocity + alpha*hi->map_velocity;
    result.global_valid = lo.global_valid && hi->global_valid;
    return true;
}
bool PoseBuffer::Covers(double begin, double end, bool global) const {
    Pose p;
    if (begin > end || !At(begin,p,global) || !At(end,p,global)) return false;
    for (std::size_t i=1; i<poses_.size(); ++i) {
        if (poses_[i].stamp <= begin || poses_[i-1].stamp >= end) continue;
        if (poses_[i].stamp - poses_[i-1].stamp > max_gap_ ||
            (global && (!poses_[i].global_valid || !poses_[i-1].global_valid))) return false;
    }
    return true;
}
}  // namespace lightning::ins
