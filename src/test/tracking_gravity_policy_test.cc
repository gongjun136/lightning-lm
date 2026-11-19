// tracking_gravity_policy.h 的独立 assert 测试，不依赖测试框架。
// 依赖工程中已有的 utils/tracking_overlap_policy.h（提供 MatchPhase）。
// 构建示例（在工程根目录）：
//   g++ -std=c++17 -Wall -Wextra -I. utils/tracking_gravity_policy_test.cc -o tracking_gravity_policy_test
#undef NDEBUG
#include <cassert>
#include <cmath>
#include <cstdio>
#include <initializer_list>
#include <limits>
#include <string>

#include "utils/tracking_gravity_policy.h"

namespace {

namespace dg = lightning::diagnostics;
using dg::MatchPhase;

constexpr double kSharedMinCos = 0.95;  // 初始化/重定位常见配置值
const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kInf = std::numeric_limits<double>::infinity();

const MatchPhase kAllPhases[] = {MatchPhase::kTracking, MatchPhase::kInitialization,
                                 MatchPhase::kRelocalization};

double CosDeg(double deg) { return std::cos(deg * dg::kGravityPolicyPi / 180.0); }

// 模拟调用者的重力门：先按阶段选阈值，再做原有 >= 比较。
bool GravityGate(MatchPhase phase, double gravity_cos, double shared, double tracking) {
  return dg::PassesMinGravityCos(gravity_cos,
                                 dg::SelectMinGravityCos(phase, shared, tracking));
}

// 模拟原有地图判定中的门限组合：重合率开关只影响重合率门。
struct SimResult {
  bool overlap_ok;
  bool gravity_ok;
  double min_cos;
};

SimResult SimulateCheck(MatchPhase phase, bool overlap_check_enabled, double overlap,
                        double min_overlap, double gravity_cos, double shared,
                        double tracking) {
  SimResult r{};
  r.min_cos = dg::SelectMinGravityCos(phase, shared, tracking);
  r.overlap_ok = !overlap_check_enabled || overlap >= min_overlap;
  r.gravity_ok = dg::PassesMinGravityCos(gravity_cos, r.min_cos);
  return r;
}

void TestDefaultAngle() {
  assert(dg::kDefaultTrackingMaxGravityAngleDeg == 20.0);
  const double expected = std::cos(20.0 * dg::kGravityPolicyPi / 180.0);
  const double t = dg::DefaultTrackingMinGravityCos();
  assert(t == expected);
  assert(t != 0.94);
  assert(std::fabs(t - 0.9396926207859084) < 1e-15);
  double converted = 0.0;
  std::string err;
  assert(dg::ConvertTrackingMaxGravityAngleDeg(20.0, &converted, &err));
  assert(converted == expected);
  assert(err.empty());
  assert(std::isnan(dg::UnknownMinGravityCos()));
}

void TestEqualityAndBothSides() {
  const double t = dg::DefaultTrackingMinGravityCos();
  const double below = std::nextafter(t, -kInf);
  const double above = std::nextafter(t, kInf);
  assert(dg::PassesMinGravityCos(t, t));
  assert(!dg::PassesMinGravityCos(below, t));
  assert(dg::PassesMinGravityCos(above, t));
  assert(GravityGate(MatchPhase::kTracking, t, kSharedMinCos, t));
  assert(!GravityGate(MatchPhase::kTracking, below, kSharedMinCos, t));
  // 共享阈值同样相等通过、略低失败。
  assert(GravityGate(MatchPhase::kRelocalization, kSharedMinCos, kSharedMinCos, t));
  assert(GravityGate(MatchPhase::kInitialization, kSharedMinCos, kSharedMinCos, t));
  assert(!GravityGate(MatchPhase::kRelocalization, std::nextafter(kSharedMinCos, 0.0),
                      kSharedMinCos, t));
}

void Test19DegTrackingOnly() {
  const double t = dg::DefaultTrackingMinGravityCos();
  const double g = CosDeg(19.0);
  assert(GravityGate(MatchPhase::kTracking, g, kSharedMinCos, t));
  assert(!GravityGate(MatchPhase::kInitialization, g, kSharedMinCos, t));
  assert(!GravityGate(MatchPhase::kRelocalization, g, kSharedMinCos, t));
}

void Test21DegFlipAndNaNFail() {
  const double t = dg::DefaultTrackingMinGravityCos();
  for (MatchPhase phase : kAllPhases) {
    assert(!GravityGate(phase, CosDeg(21.0), kSharedMinCos, t));
    assert(!GravityGate(phase, CosDeg(90.0), kSharedMinCos, t));
    assert(!GravityGate(phase, CosDeg(180.0), kSharedMinCos, t));  // 翻转
    assert(!GravityGate(phase, -1.0, kSharedMinCos, t));
    assert(!GravityGate(phase, kNaN, kSharedMinCos, t));  // NaN 自然失败
  }
}

void TestOtherPhasesUseCallerThreshold() {
  const double t = dg::DefaultTrackingMinGravityCos();
  for (double shared : {0.95, 0.97, 0.5, 0.99}) {
    assert(dg::SelectMinGravityCos(MatchPhase::kInitialization, shared, t) == shared);
    assert(dg::SelectMinGravityCos(MatchPhase::kRelocalization, shared, t) == shared);
    assert(dg::SelectMinGravityCos(MatchPhase::kTracking, shared, t) == t);
  }
  // 跟踪角度配置变化不影响其它阶段。
  const double t5 = CosDeg(5.0);
  assert(dg::SelectMinGravityCos(MatchPhase::kInitialization, kSharedMinCos, t5) ==
         kSharedMinCos);
  assert(dg::SelectMinGravityCos(MatchPhase::kRelocalization, kSharedMinCos, t5) ==
         kSharedMinCos);
  assert(dg::SelectMinGravityCos(MatchPhase::kTracking, kSharedMinCos, t5) == t5);
}

void TestOverlapSwitchIndependent() {
  const double t = dg::DefaultTrackingMinGravityCos();
  for (MatchPhase phase : kAllPhases) {
    for (double deg : {0.0, 19.0, 20.0, 21.0, 180.0}) {
      const double g = CosDeg(deg);
      const SimResult on = SimulateCheck(phase, true, 0.1, 0.5, g, kSharedMinCos, t);
      const SimResult off = SimulateCheck(phase, false, 0.1, 0.5, g, kSharedMinCos, t);
      assert(on.min_cos == off.min_cos);
      assert(on.gravity_ok == off.gravity_ok);
      assert(!on.overlap_ok);
      assert(off.overlap_ok);
    }
  }
}

void TestInvalidConfig() {
  const double bad[] = {kNaN, kInf, -kInf, -1e-9, -1.0, 90.0,
                        std::nextafter(90.0, kInf), 180.0, 1e300};
  for (double b : bad) {
    assert(!dg::IsValidTrackingMaxGravityAngleDeg(b));
    double out = 123.0;
    std::string err;
    assert(!dg::ConvertTrackingMaxGravityAngleDeg(b, &out, &err));
    assert(out == 123.0);  // 不 clamp、不写出
    assert(!err.empty());
  }
  const double good[] = {0.0, 1e-9, 19.0, 20.0, 45.0, std::nextafter(90.0, 0.0)};
  for (double g : good) {
    assert(dg::IsValidTrackingMaxGravityAngleDeg(g));
    double out = 123.0;
    std::string err;
    assert(dg::ConvertTrackingMaxGravityAngleDeg(g, &out, &err));
    assert(out == CosDeg(g));
    assert(err.empty());
  }
  double zero_out = 0.0;
  assert(dg::ConvertTrackingMaxGravityAngleDeg(0.0, &zero_out, nullptr));
  assert(zero_out == 1.0);
  std::string err;
  assert(!dg::ConvertTrackingMaxGravityAngleDeg(20.0, nullptr, &err));
  assert(!err.empty());
}

}  // namespace

int main() {
  TestDefaultAngle();
  TestEqualityAndBothSides();
  Test19DegTrackingOnly();
  Test21DegFlipAndNaNFail();
  TestOtherPhasesUseCallerThreshold();
  TestOverlapSwitchIndependent();
  TestInvalidConfig();
  std::printf("tracking_gravity_policy_test: all passed\n");
  return 0;
}
