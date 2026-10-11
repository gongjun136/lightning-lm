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
    std::array<double, FieldCount> time{};       ///< corrected ROS-domain receive stamps
    std::array<double, FieldCount> arrival{};    ///< monotonic receive times
    std::array<bool, FieldCount> valid{};  ///< 各导航字段的设备有效性与时间戳检查结果。
    double latitude = 0 /**< 设备报告的纬度（deg）。 */, longitude = 0 /**< 设备报告的经度（deg）。 */, altitude = 0 /**< 设备报告的高度（m），高程基准由配置确定。 */;
    double heading = 0 /**< 设备报告的北起顺时针航向角（deg）。 */, pitch = 0 /**< 设备报告的抬头为正俯仰角（deg）。 */, roll = 0 /**< 设备报告的右侧下倾为正滚转角（deg）。 */;
    Vec3d velocity = Vec3d::Zero();              ///< local ENU, m/s
    Vec3d position_sigma = Vec3d::Zero();        ///< E/N/U, m
    Vec3d attitude_sigma = Vec3d::Zero();        ///< heading/pitch/roll, deg
    Vec3d velocity_sigma = Vec3d::Zero();  ///< 设备报告的 ENU 速度标准差（m/s）。
    int system_state = 0 /**< 设备导航解算状态码。 */, satellite_status = 0 /**< 设备卫星定位状态码。 */;
    double differential_age = 0;  ///< 设备报告的差分龄期（s）。
    double Stamp() const { return time[Latitude]; }
};

struct QualityPolicy {
    double max_age = 0.15;  ///< 各导航字段允许的最大年龄（s）。
    double max_skew = 0.02;  ///< 同组导航字段允许的最大时间偏差（s）。
    double max_future = 0.02;  ///< 允许导航字段超前当前时间的容差（s）。
    double max_differential_age = 2.0;  ///< 允许的最大差分龄期（s）。
    double max_horizontal_sigma = 0.10;  ///< 允许的最大水平位置标准差（m）。
    double max_vertical_sigma = 0.20;  ///< 允许的最大垂直位置标准差（m）。
    double max_attitude_sigma = 0.5;  ///< 允许的最大姿态标准差（deg）。
    double max_velocity_sigma = 0.20;  ///< 允许的最大速度标准差（m/s）。
    int recovery_samples = 10;  ///< 质量恢复后重新发布所需的连续有效样本数。
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
    int required_ /**< 恢复有效输出所需的连续有效样本数。 */, count_ = 0 /**< 当前已连续满足恢复条件的样本数。 */;
    double last_stamp_ = 0;  ///< 恢复门控最近处理的传感器时间戳（s）。
    bool good_ = false;  ///< 恢复门控是否允许输出有效导航结果。
};

class GeoReference {
 public:
    GeoReference(double latitude, double longitude, double ellipsoid_height);
    Vec3d Forward(double latitude, double longitude, double ellipsoid_height,
                  Mat3d* fixed_from_local_enu = nullptr) const;
    Vec3d Reverse(const Vec3d& enu) const;
 private:
    GeographicLib::LocalCartesian local_;  ///< 以配置参考点为原点的大地坐标到局部 ENU 转换器。
};

// CHC: north-clockwise heading, nose-up pitch, right-down roll.
// Result maps rear-body FLU axes to the local ENU axes at the measurement.
SO3 VehicleAttitude(double heading_deg, double pitch_deg, double roll_deg);

struct Pose {
    double stamp = 0;  ///< 后轴位姿对应的传感器时间戳（s）。
    SE3 map_from_rear;  ///< 将后轴车体坐标变换到局部地图坐标的位姿。
    Vec3d map_velocity = Vec3d::Zero();  ///< 局部地图坐标系中的速度（m/s）。
    bool global_valid = false;  ///< 该位姿是否通过设备质量与恢复门控。
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
    std::deque<Pose> poses_;  ///< 按时间保存的有效 INS 后轴位姿历史。
    double horizon_ /**< 位姿历史保留时间窗口（s）。 */, max_gap_ /**< 允许进行位姿插值的最大样本时间间隔（s）。 */;
};

}  // namespace lightning::ins
