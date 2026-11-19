// 跟踪阶段独立重力一致性门限策略。
//
// gravity_cos 定义为 (candidate_pose * odometry_pose.inverse()).rotationMatrix()(2,2)，
// 表示候选位姿与里程计位姿之间相对旋转的 z 轴对齐程度，并不是车辆本身的坡度。
//
// 本文件只负责决定“本次判定使用哪个余弦阈值”：
//   - 跟踪阶段：使用 lidar_loc.tracking_map_max_gravity_angle_deg 在配置加载时
//     一次性转换得到的余弦值（默认 20 度，精确为 cos(20*pi/180)，不是 0.94）；
//   - 初始化 / 重定位阶段：原样使用调用者传入的共享可配置余弦值（通常 0.95），
//     不校验、不改写，保持原有允许范围和语义。
// 比较语义不变：通过条件仍是 gravity_cos >= min_cos，相等通过，NaN 自然失败。
// 与重合率开关完全独立：本文件中没有任何函数接收重合率开关参数。
// 不改算法、评分、报警时序或其它门限。
#pragma once

#include <cmath>
#include <limits>
#include <string>

#include "utils/tracking_overlap_policy.h"  // MatchPhase

namespace lightning::diagnostics {

// 圆周率常量，与常见 M_PI 的 double 值一致；不依赖非标准宏。
inline constexpr double kGravityPolicyPi = 3.14159265358979323846;

// 跟踪阶段默认最大重力夹角（度）。
inline constexpr double kDefaultTrackingMaxGravityAngleDeg = 20.0;

// 跟踪阶段夹角合法范围为 [0, 90)，此为开区间上界（度）。
inline constexpr double kTrackingMaxGravityAngleExclusiveUpperDeg = 90.0;

// 配置键名（机器字段）。
inline constexpr const char* kTrackingMaxGravityAngleKey =
    "lidar_loc.tracking_map_max_gravity_angle_deg";

// 跟踪夹角配置是否合法：有限值且位于 [0, 90)。不做任何 clamp。
inline bool IsValidTrackingMaxGravityAngleDeg(double angle_deg) {
  return std::isfinite(angle_deg) && angle_deg >= 0.0 &&
         angle_deg < kTrackingMaxGravityAngleExclusiveUpperDeg;
}

// 角度（度）转余弦。调用前应先通过 IsValidTrackingMaxGravityAngleDeg 校验。
inline double GravityAngleDegToMinCos(double angle_deg) {
  return std::cos(angle_deg * kGravityPolicyPi / 180.0);
}

// 默认跟踪余弦阈值，精确等于 cos(20*pi/180)。
inline double DefaultTrackingMinGravityCos() {
  return GravityAngleDegToMinCos(kDefaultTrackingMaxGravityAngleDeg);
}

// 未知 / 未执行判定时的阈值占位（统计、快照中使用）。
inline double UnknownMinGravityCos() {
  return std::numeric_limits<double>::quiet_NaN();
}

// 配置加载时调用一次：校验角度并转换为余弦。
// 非法时返回 false，*min_cos 保持不变，*error（若非空）写入原因；调用方应拒绝启动。
inline bool ConvertTrackingMaxGravityAngleDeg(double angle_deg, double* min_cos,
                                              std::string* error) {
  if (min_cos == nullptr) {
    if (error != nullptr) {
      *error = std::string(kTrackingMaxGravityAngleKey) + "：输出指针为空";
    }
    return false;
  }
  if (!IsValidTrackingMaxGravityAngleDeg(angle_deg)) {
    if (error != nullptr) {
      *error = std::string(kTrackingMaxGravityAngleKey) +
               " 非法：必须是有限值且位于 [0, 90) 度，实际为 " +
               std::to_string(angle_deg);
    }
    return false;
  }
  *min_cos = GravityAngleDegToMinCos(angle_deg);
  return true;
}

// 按调用者明确给出的阶段选择本次实际使用的余弦阈值。
// shared_min_cos：初始化/重定位共享的原有可配置值，原样透传。
// tracking_min_cos：配置加载时已转换好的跟踪阈值。
inline double SelectMinGravityCos(MatchPhase phase, double shared_min_cos,
                                  double tracking_min_cos) {
  switch (phase) {
    case MatchPhase::kTracking:
      return tracking_min_cos;
    case MatchPhase::kInitialization:
    case MatchPhase::kRelocalization:
      return shared_min_cos;
  }
  // 非法枚举值：回到共享阈值（与引入本策略前的行为一致）。
  return shared_min_cos;
}

// 原有比较语义：相等通过；任一操作数为 NaN 时 >= 自然为 false。
// 不引入任何允许 NaN 的分支。
inline bool PassesMinGravityCos(double gravity_cos, double min_cos) {
  return gravity_cos >= min_cos;
}

}  // namespace lightning::diagnostics
