// tracking_overlap_policy.h 的独立确定性测试。
// 使用自定义检查函数而非 assert，Release（NDEBUG）下同样生效；不休眠，仅使用标准库与合成数值。
#include "utils/tracking_overlap_policy.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <initializer_list>
#include <iostream>
#include <limits>

namespace {

using lightning::diagnostics::MatchPhase;
using lightning::diagnostics::OverlapGate;
using lightning::diagnostics::OverlapGateEnabled;
using lightning::diagnostics::OverlapObservation;
using lightning::diagnostics::OverlapResultName;
using lightning::diagnostics::PassOverlapGate;
using lightning::diagnostics::SelectOverlapGate;
using lightning::diagnostics::TrackingOverlapObserver;
using lightning::diagnostics::ValidOverlap;

int g_failures = 0;
int g_checks = 0;

void Check(bool condition, const char* expression, int line) {
  ++g_checks;
  if (!condition) {
    ++g_failures;
    std::cerr << "检查失败（第 " << line << " 行）：" << expression << std::endl;
  }
}

#define CHECK(condition) Check(static_cast<bool>(condition), #condition, __LINE__)

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kThreshold = 0.60;
constexpr double kLow = 0.59;
constexpr std::int64_t kSecondNs = 1000000000;
constexpr std::int64_t kStep = kSecondNs / 10;
constexpr std::int64_t kMaxStamp = std::numeric_limits<std::int64_t>::max();
constexpr std::int64_t kMinStamp = std::numeric_limits<std::int64_t>::min();

bool Near(double a, double b) { return std::fabs(a - b) <= 1e-12; }

bool Named(const char* actual, const char* expected) {
  return actual != nullptr && std::strcmp(actual, expected) == 0;
}

// 合成的“已接受、已评估、仅观测、低于阈值”帧。
OverlapObservation Low(TrackingOverlapObserver& observer, std::int64_t stamp_ns, double steady_seconds) {
  return observer.Observe(true, true, OverlapGate::kObserveOnly, kLow, kThreshold, stamp_ns, steady_seconds);
}

static_assert(SelectOverlapGate(MatchPhase::kTracking, true) == OverlapGate::kEnforce, "跟踪开启应强制执行");
static_assert(SelectOverlapGate(MatchPhase::kTracking, false) == OverlapGate::kObserveOnly, "跟踪关闭应仅观测");
static_assert(SelectOverlapGate(MatchPhase::kInitialization, false) == OverlapGate::kEnforce, "初始化始终强制执行");
static_assert(SelectOverlapGate(MatchPhase::kRelocalization, false) == OverlapGate::kEnforce, "重定位始终强制执行");
static_assert(OverlapGateEnabled(OverlapGate::kEnforce), "强制门控应为启用");
static_assert(!OverlapGateEnabled(OverlapGate::kObserveOnly), "仅观测门控应为关闭");

void TestGateSelection() {
  for (bool enabled : {true, false}) {
    CHECK(SelectOverlapGate(MatchPhase::kInitialization, enabled) == OverlapGate::kEnforce);
    CHECK(SelectOverlapGate(MatchPhase::kRelocalization, enabled) == OverlapGate::kEnforce);
  }
  CHECK(SelectOverlapGate(MatchPhase::kTracking, true) == OverlapGate::kEnforce);
  CHECK(SelectOverlapGate(MatchPhase::kTracking, false) == OverlapGate::kObserveOnly);
  CHECK(OverlapGateEnabled(SelectOverlapGate(MatchPhase::kTracking, true)));
  CHECK(!OverlapGateEnabled(SelectOverlapGate(MatchPhase::kTracking, false)));
}

void TestPhaseTransitions() {
  // 初始化 -> 跟踪 -> 丢失（全局重定位预检）-> 重定位最终验证 -> 恢复跟踪
  const MatchPhase sequence[] = {MatchPhase::kInitialization, MatchPhase::kTracking, MatchPhase::kRelocalization,
                                 MatchPhase::kRelocalization, MatchPhase::kTracking};
  const OverlapGate expected_when_off[] = {OverlapGate::kEnforce, OverlapGate::kObserveOnly, OverlapGate::kEnforce,
                                           OverlapGate::kEnforce, OverlapGate::kObserveOnly};
  for (int i = 0; i < 5; ++i) {
    CHECK(SelectOverlapGate(sequence[i], false) == expected_when_off[i]);
    CHECK(SelectOverlapGate(sequence[i], true) == OverlapGate::kEnforce);
  }

  TrackingOverlapObserver observer;
  const std::int64_t t0 = 1 * kSecondNs;
  CHECK(Low(observer, t0, 10.0).streak_frames == 1);
  CHECK(Low(observer, t0 + kStep, 10.1).streak_frames == 2);
  // 丢失：阶段切换显式重置连续性
  observer.ResetContinuity();
  const OverlapObservation reloc =
      observer.Observe(false, true, SelectOverlapGate(MatchPhase::kRelocalization, false), kLow, kThreshold,
                       t0 + 2 * kStep, 10.2);
  CHECK(!reloc.low);
  CHECK(reloc.streak_frames == 0);
  CHECK(reloc.total_low_frames == 2);
  // 即使在强制门控下帧被接受，也不是观测低重叠事件
  const OverlapObservation reloc_accepted =
      observer.Observe(true, true, SelectOverlapGate(MatchPhase::kRelocalization, false), kLow, kThreshold,
                       t0 + 3 * kStep, 10.3);
  CHECK(!reloc_accepted.low);
  CHECK(reloc_accepted.total_low_frames == 2);
  // 恢复跟踪：重新使用仅观测门控，片段从 1 开始
  observer.ResetContinuity();
  const OverlapObservation back =
      observer.Observe(true, true, SelectOverlapGate(MatchPhase::kTracking, false), kLow, kThreshold,
                       t0 + 4 * kStep, 10.4);
  CHECK(back.low);
  CHECK(back.streak_frames == 1);
  CHECK(Near(back.streak_seconds, 0.0));
  CHECK(back.total_low_frames == 3);
  // 跟踪开启时的强制门控同样不产生观测事件
  const OverlapObservation enforced =
      observer.Observe(true, true, SelectOverlapGate(MatchPhase::kTracking, true), kLow, kThreshold,
                       t0 + 5 * kStep, 10.5);
  CHECK(!enforced.low);
  CHECK(enforced.total_low_frames == 3);
}

void TestThresholds() {
  for (OverlapGate gate : {OverlapGate::kEnforce, OverlapGate::kObserveOnly}) {
    CHECK(PassOverlapGate(0.60, kThreshold, gate));  // 与阈值相等即通过
    CHECK(PassOverlapGate(1.0, kThreshold, gate));
    CHECK(Named(OverlapResultName(true, 0.60, kThreshold, gate), "pass"));
  }
  CHECK(!PassOverlapGate(kLow, kThreshold, OverlapGate::kEnforce));
  CHECK(PassOverlapGate(kLow, kThreshold, OverlapGate::kObserveOnly));
  CHECK(Named(OverlapResultName(true, kLow, kThreshold, OverlapGate::kEnforce), "fail"));
  CHECK(Named(OverlapResultName(true, kLow, kThreshold, OverlapGate::kObserveOnly), "bypassed_low"));
  CHECK(ValidOverlap(0.0, 0.0));
  CHECK(ValidOverlap(1.0, 1.0));
  CHECK(ValidOverlap(kLow, kThreshold));
}

void TestInvalidEvidence() {
  const double bad_values[] = {kNaN, kInf, -kInf, -0.01, 1.01};
  for (OverlapGate gate : {OverlapGate::kEnforce, OverlapGate::kObserveOnly}) {
    for (double bad : bad_values) {
      CHECK(!ValidOverlap(bad, kThreshold));
      CHECK(!ValidOverlap(0.7, bad));
      CHECK(!PassOverlapGate(bad, kThreshold, gate));  // 仅观测模式也不放行无效证据
      CHECK(!PassOverlapGate(0.7, bad, gate));
      CHECK(Named(OverlapResultName(true, bad, kThreshold, gate), "invalid"));
      CHECK(Named(OverlapResultName(true, 0.7, bad, gate), "invalid"));
      // 未执行检查：无论残留数值如何都是 unavailable
      CHECK(Named(OverlapResultName(false, bad, kThreshold, gate), "unavailable"));
    }
    CHECK(Named(OverlapResultName(false, 0.9, kThreshold, gate), "unavailable"));
    CHECK(Named(OverlapResultName(false, kLow, kThreshold, gate), "unavailable"));
  }
}

void TestNonLowFrames() {
  TrackingOverlapObserver observer;
  const OverlapObservation cases[] = {
      observer.Observe(false, true, OverlapGate::kObserveOnly, kLow, kThreshold, kSecondNs, 1.0),
      observer.Observe(true, false, OverlapGate::kObserveOnly, kLow, kThreshold, kSecondNs, 1.0),
      observer.Observe(true, true, OverlapGate::kEnforce, kLow, kThreshold, kSecondNs, 1.0),
      observer.Observe(true, true, OverlapGate::kObserveOnly, kNaN, kThreshold, kSecondNs, 1.0),
      observer.Observe(true, true, OverlapGate::kObserveOnly, kLow, kInf, kSecondNs, 1.0),
      observer.Observe(true, true, OverlapGate::kObserveOnly, -0.1, kThreshold, kSecondNs, 1.0),
      observer.Observe(true, true, OverlapGate::kObserveOnly, 0.60, kThreshold, kSecondNs, 1.0),
  };
  for (const OverlapObservation& obs : cases) {
    CHECK(!obs.low);
    CHECK(!obs.log_due);
    CHECK(!obs.capture_due);
    CHECK(obs.streak_frames == 0);
    CHECK(obs.total_low_frames == 0);
    CHECK(obs.capture_suppressed == 0);
    CHECK(Near(obs.streak_seconds, 0.0));
  }
  // 非低重叠帧（例如未接受）打断片段，但不清零累计值
  const OverlapObservation first = Low(observer, 2 * kSecondNs, 2.0);
  CHECK(first.streak_frames == 1);
  CHECK(first.log_due);
  CHECK(first.capture_due);
  const OverlapObservation second = Low(observer, 2 * kSecondNs + kStep, 2.1);
  CHECK(second.streak_frames == 2);
  CHECK(Near(second.streak_seconds, 0.1));
  CHECK(second.capture_suppressed == 1);
  const OverlapObservation rejected =
      observer.Observe(false, true, OverlapGate::kObserveOnly, kLow, kThreshold, 2 * kSecondNs + 2 * kStep, 2.2);
  CHECK(!rejected.low);
  CHECK(rejected.streak_frames == 0);
  CHECK(rejected.total_low_frames == 2);
  CHECK(rejected.capture_suppressed == 1);
  const OverlapObservation third = Low(observer, 2 * kSecondNs + 3 * kStep, 2.3);
  CHECK(third.streak_frames == 1);
  CHECK(Near(third.streak_seconds, 0.0));
  CHECK(third.total_low_frames == 3);
  CHECK(third.capture_suppressed == 2);
}

void TestStreakContinuity() {
  TrackingOverlapObserver observer;
  const std::int64_t t0 = 1 * kSecondNs;
  const std::int64_t tc = t0 + kStep + kSecondNs / 2;  // 与上一帧间隔恰为 0.5 s
  const OverlapObservation results[] = {
      Low(observer, t0, 10.0),                       // 0：首帧，时长 0
      Low(observer, t0 + kStep, 10.1),               // 1：延续
      Low(observer, tc, 10.2),                       // 2：间隔 0.5 s，延续
      Low(observer, tc + kSecondNs / 2 + 1, 10.3),   // 3：间隔超过 0.5 s，新片段
      Low(observer, tc + kSecondNs / 2 + 1, 10.4),   // 4：重复时间戳，新片段
      Low(observer, tc, 10.5),                       // 5：时间戳回退，新片段
      Low(observer, 0, 10.6),                        // 6：零时间戳，新片段
      Low(observer, -kSecondNs, 10.7),               // 7：负时间戳，新片段
      Low(observer, 5 * kSecondNs, 10.8),            // 8：上一帧无锚点，新片段
      Low(observer, 5 * kSecondNs + 1, 10.9),        // 9：延续 1 ns
  };
  CHECK(results[0].streak_frames == 1);
  CHECK(Near(results[0].streak_seconds, 0.0));
  CHECK(results[1].streak_frames == 2);
  CHECK(Near(results[1].streak_seconds, 0.1));
  CHECK(results[2].streak_frames == 3);
  CHECK(Near(results[2].streak_seconds, 0.6));
  for (int i = 3; i <= 8; ++i) {
    CHECK(results[i].streak_frames == 1);
    CHECK(Near(results[i].streak_seconds, 0.0));
  }
  CHECK(results[9].streak_frames == 2);
  CHECK(Near(results[9].streak_seconds, 1e-9));
  CHECK(results[9].total_low_frames == 10);
  for (const OverlapObservation& obs : results) {
    CHECK(obs.low);
    CHECK(std::isfinite(obs.streak_seconds));
    CHECK(obs.streak_seconds >= 0.0);
  }
}

void TestTimestampOverflow() {
  TrackingOverlapObserver observer;
  CHECK(Low(observer, 1, 0.0).streak_frames == 1);
  const OverlapObservation jump = Low(observer, kMaxStamp, 0.0);  // 巨大间隔：新片段，不溢出
  CHECK(jump.streak_frames == 1);
  CHECK(Near(jump.streak_seconds, 0.0));
  const OverlapObservation back = Low(observer, kMaxStamp - 1, 0.0);  // 回退：新片段
  CHECK(back.streak_frames == 1);
  const OverlapObservation edge = Low(observer, kMaxStamp, 0.0);  // 最大值处延续
  CHECK(edge.streak_frames == 2);
  CHECK(Near(edge.streak_seconds, 1e-9));
  const OverlapObservation minimum = Low(observer, kMinStamp, 0.0);
  CHECK(minimum.streak_frames == 1);
  CHECK(Near(minimum.streak_seconds, 0.0));
  CHECK(std::isfinite(minimum.streak_seconds));
}

void TestRateLimitsAcrossEpisodes() {
  TrackingOverlapObserver observer;
  std::int64_t stamp = 100 * kSecondNs;
  const OverlapObservation r1 = Low(observer, stamp, 100.0);  // 首个合格帧立即到期
  CHECK(r1.log_due);
  CHECK(r1.capture_due);
  CHECK(r1.capture_suppressed == 0);
  stamp += kStep;
  const OverlapObservation r2 = Low(observer, stamp, 101.0);
  CHECK(!r2.log_due);
  CHECK(!r2.capture_due);
  CHECK(r2.capture_suppressed == 1);
  stamp += kStep;
  const OverlapObservation r3 = Low(observer, stamp, 105.0);  // 日志间隔相等即到期
  CHECK(r3.log_due);
  CHECK(!r3.capture_due);
  CHECK(r3.capture_suppressed == 2);
  // 片段被未接受帧打断并显式重置（模拟丢失 -> 重定位 -> 恢复跟踪），全局限频不得被重置
  stamp += kStep;
  CHECK(!observer.Observe(false, true, OverlapGate::kObserveOnly, kLow, kThreshold, stamp, 106.0).low);
  observer.ResetContinuity();
  stamp += 60 * kSecondNs;
  const OverlapObservation r4 = Low(observer, stamp, 107.0);
  CHECK(r4.streak_frames == 1);
  CHECK(!r4.log_due);
  CHECK(!r4.capture_due);
  CHECK(r4.capture_suppressed == 3);
  stamp += kStep;
  const OverlapObservation r5 = Low(observer, stamp, 110.5);
  CHECK(r5.log_due);
  CHECK(!r5.capture_due);
  CHECK(r5.capture_suppressed == 4);
  stamp += kStep;
  const OverlapObservation r6 = Low(observer, stamp, 130.0);  // 采集间隔相等即到期
  CHECK(r6.capture_due);
  CHECK(r6.log_due);
  CHECK(r6.capture_suppressed == 4);
  stamp += kStep;
  const OverlapObservation r7 = Low(observer, stamp, 130.5);
  CHECK(!r7.log_due);
  CHECK(!r7.capture_due);
  CHECK(r7.capture_suppressed == 5);
  // 长时间间隔后只到期一次，不补发
  stamp += kStep;
  const OverlapObservation r8 = Low(observer, stamp, 10000.0);
  CHECK(r8.log_due);
  CHECK(r8.capture_due);
  CHECK(r8.capture_suppressed == 5);
  stamp += kStep;
  const OverlapObservation r9 = Low(observer, stamp, 10000.1);
  CHECK(!r9.log_due);
  CHECK(!r9.capture_due);
  stamp += kStep;
  const OverlapObservation r10 = Low(observer, stamp, 10001.0);
  CHECK(!r10.log_due);
  CHECK(!r10.capture_due);
  CHECK(r10.capture_suppressed == 7);
  CHECK(r10.total_low_frames == 10);
}

void TestBackwardSteadyClock() {
  TrackingOverlapObserver observer;
  std::int64_t stamp = kSecondNs;
  const OverlapObservation b1 = Low(observer, stamp, 100.0);
  CHECK(b1.log_due);
  CHECK(b1.capture_due);
  stamp += kStep;
  const OverlapObservation b2 = Low(observer, stamp, 50.0);  // 回退：重新锚定，不突发
  CHECK(!b2.log_due);
  CHECK(!b2.capture_due);
  stamp += kStep;
  const OverlapObservation b3 = Low(observer, stamp, 54.9);
  CHECK(!b3.log_due);
  CHECK(!b3.capture_due);
  stamp += kStep;
  const OverlapObservation b4 = Low(observer, stamp, 55.0);
  CHECK(b4.log_due);
  CHECK(!b4.capture_due);
  stamp += kStep;
  const OverlapObservation b5 = Low(observer, stamp, 79.9);
  CHECK(b5.log_due);
  CHECK(!b5.capture_due);
  stamp += kStep;
  const OverlapObservation b6 = Low(observer, stamp, 80.0);
  CHECK(!b6.log_due);
  CHECK(b6.capture_due);
  stamp += kStep;
  const OverlapObservation b7 = Low(observer, stamp, 0.0);  // 再次回退
  CHECK(!b7.log_due);
  CHECK(!b7.capture_due);
  CHECK(b7.capture_suppressed == 5);
  CHECK(b7.total_low_frames == 7);
}

void TestInvalidSteadyTime() {
  TrackingOverlapObserver observer;
  std::int64_t stamp = kSecondNs;
  const double invalid_times[] = {kNaN, kInf, -kInf, -1.0};
  std::uint64_t expected = 0;
  for (double t : invalid_times) {
    const OverlapObservation obs = Low(observer, stamp, t);
    ++expected;
    CHECK(obs.low);
    CHECK(!obs.log_due);
    CHECK(!obs.capture_due);
    CHECK(obs.total_low_frames == expected);
    CHECK(obs.capture_suppressed == expected);
    stamp += kStep;
  }
  const OverlapObservation usable = Low(observer, stamp, 200.0);  // 首个可用稳态时间立即到期
  CHECK(usable.log_due);
  CHECK(usable.capture_due);
  CHECK(usable.capture_suppressed == 4);
  CHECK(usable.total_low_frames == 5);
  stamp += kStep;
  const OverlapObservation again_invalid = Low(observer, stamp, kNaN);
  CHECK(!again_invalid.log_due);
  CHECK(!again_invalid.capture_due);
  CHECK(again_invalid.capture_suppressed == 5);
  stamp += kStep;
  const OverlapObservation limited = Low(observer, stamp, 201.0);  // 限频状态保留
  CHECK(!limited.log_due);
  CHECK(!limited.capture_due);
  CHECK(limited.capture_suppressed == 6);
  CHECK(limited.streak_frames == 7);  // 稳态时间不影响传感器时间戳连续性
  CHECK(Near(limited.streak_seconds, 0.6));
}

}  // namespace

int main() {
  TestGateSelection();
  TestPhaseTransitions();
  TestThresholds();
  TestInvalidEvidence();
  TestNonLowFrames();
  TestStreakContinuity();
  TestTimestampOverflow();
  TestRateLimitsAcrossEpisodes();
  TestBackwardSteadyClock();
  TestInvalidSteadyTime();
  if (g_failures != 0) {
    std::cerr << "共 " << g_failures << " 项检查失败（总计 " << g_checks << " 项）" << std::endl;
    return 1;
  }
  std::cout << "全部 " << g_checks << " 项检查通过" << std::endl;
  return 0;
}
