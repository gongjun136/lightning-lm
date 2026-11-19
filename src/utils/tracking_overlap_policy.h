// 跟踪阶段地图重叠门控策略与低重叠观测器（C++17，仅头文件，仅依赖标准库）。
// 本文件只提供判定与观测辅助：不决定帧是否被接受，也不修改调用方的任何计数器。
#pragma once

#include <cmath>
#include <cstdint>
#include <limits>

namespace lightning::diagnostics {

// 重叠门控：kEnforce 表示低于阈值即拒绝；kObserveOnly 表示仅观测、不拒绝。
enum class OverlapGate { kEnforce, kObserveOnly };

// 匹配阶段：只有正常跟踪受配置开关影响。
enum class MatchPhase { kTracking, kInitialization, kRelocalization };

// 正常跟踪按 map_overlap_check_enabled 选择门控；初始化与重定位始终强制执行。
constexpr OverlapGate SelectOverlapGate(MatchPhase phase, bool tracking_enabled) noexcept {
  return (phase == MatchPhase::kTracking && !tracking_enabled) ? OverlapGate::kObserveOnly
                                                                : OverlapGate::kEnforce;
}

// 门控是否会拒绝低重叠帧。
constexpr bool OverlapGateEnabled(OverlapGate gate) noexcept { return gate == OverlapGate::kEnforce; }

// 重叠证据有效：重叠率与阈值均为有限值且位于 [0,1]。
inline bool ValidOverlap(double value, double threshold) noexcept {
  return std::isfinite(value) && std::isfinite(threshold) && value >= 0.0 && value <= 1.0 &&
         threshold >= 0.0 && threshold <= 1.0;
}

// 重叠谓词：无效证据在任何门控下都不通过（仅观测模式也不会隐藏无效证据）；与阈值相等视为通过。
inline bool PassOverlapGate(double value, double threshold, OverlapGate gate) noexcept {
  if (!ValidOverlap(value, threshold)) {
    return false;
  }
  if (gate == OverlapGate::kObserveOnly) {
    return true;
  }
  return value >= threshold;  // 任何非 kObserveOnly 的取值均按强制执行处理
}

// 结果名称：未检查 -> unavailable；无效 -> invalid；>= 阈值 -> pass；否则 fail 或 bypassed_low。
inline const char* OverlapResultName(bool evaluated, double value, double threshold, OverlapGate gate) noexcept {
  if (!evaluated) {
    return "unavailable";
  }
  if (!ValidOverlap(value, threshold)) {
    return "invalid";
  }
  if (value >= threshold) {
    return "pass";
  }
  return gate == OverlapGate::kObserveOnly ? "bypassed_low" : "fail";
}

// 单帧观测结果（无堆分配）。
struct OverlapObservation {
  bool low = false;          // 本帧是否为“已接受、已评估、仅观测门控下低于阈值”的帧
  bool log_due = false;      // 是否应输出人类可读日志（全局 5 s 间隔）
  bool capture_due = false;  // 是否应请求点云采集（全局 30 s 间隔）
  std::uint64_t streak_frames = 0;
  std::uint64_t total_low_frames = 0;
  std::uint64_t capture_suppressed = 0;
  double streak_seconds = 0.0;
};

// 低重叠观测器：只观测，不含隐藏阈值，不升级为失败。
// 每个尝试的跟踪帧在最终接受判定之后调用一次 Observe。
class TrackingOverlapObserver {
 public:
  static constexpr std::int64_t kMaxStampGapNs = 500000000;  // 连续性允许的最大时间戳间隔 0.5 s（含相等）
  static constexpr double kLogIntervalSeconds = 5.0;
  static constexpr double kCaptureIntervalSeconds = 30.0;

  OverlapObservation Observe(bool accepted, bool evaluated, OverlapGate gate, double ratio, double threshold,
                             std::int64_t stamp_ns, double steady_seconds) noexcept {
    OverlapObservation out;
    const bool low = accepted && evaluated && gate == OverlapGate::kObserveOnly &&
                     ValidOverlap(ratio, threshold) && ratio < threshold;
    if (!low) {
      // 非低重叠帧只重置连续片段，保留累计计数与全局限频状态。
      ResetContinuity();
      Fill(out);
      return out;
    }
    out.low = true;
    total_low_frames_ = SaturatingIncrement(total_low_frames_);
    ExtendStreak(stamp_ns);
    // 稳态时间无效时安全地不发出日志/采集，但仍计数。
    if (std::isfinite(steady_seconds) && steady_seconds >= 0.0) {
      out.log_due = Due(log_clock_, steady_seconds, kLogIntervalSeconds);
      out.capture_due = Due(capture_clock_, steady_seconds, kCaptureIntervalSeconds);
    }
    if (!out.capture_due) {
      capture_suppressed_ = SaturatingIncrement(capture_suppressed_);
    }
    Fill(out);
    return out;
  }

  // 仅重置连续片段与时间戳历史（例如阶段切换）；累计计数与全局限频保持不变。
  void ResetContinuity() noexcept {
    streak_frames_ = 0;
    has_anchor_ = false;
    first_stamp_ns_ = 0;
    last_stamp_ns_ = 0;
  }

 private:
  struct RateClock {
    bool armed = false;
    double last_seconds = 0.0;
  };

  static std::uint64_t SaturatingIncrement(std::uint64_t value) noexcept {
    return value == std::numeric_limits<std::uint64_t>::max() ? value : value + 1;
  }

  // 首个可用时刻立即到期；之后按间隔到期（相等即到期），长时间间隔后只到期一次，不补发。
  // 稳态时钟回退时重新锚定到当前时刻且不到期，避免突发。
  static bool Due(RateClock& clock, double now, double interval) noexcept {
    if (!clock.armed) {
      clock.armed = true;
      clock.last_seconds = now;
      return true;
    }
    if (now < clock.last_seconds) {
      clock.last_seconds = now;
      return false;
    }
    if (now - clock.last_seconds >= interval) {
      clock.last_seconds = now;
      return true;
    }
    return false;
  }

  void ExtendStreak(std::int64_t stamp_ns) noexcept {
    const bool stamp_ok = stamp_ns > 0;
    // 此分支中 stamp_ns > last_stamp_ns_ > 0，差值不会溢出。
    if (stamp_ok && has_anchor_ && streak_frames_ > 0 && stamp_ns > last_stamp_ns_ &&
        stamp_ns - last_stamp_ns_ <= kMaxStampGapNs) {
      last_stamp_ns_ = stamp_ns;
      streak_frames_ = SaturatingIncrement(streak_frames_);
      return;
    }
    // 重复、回退、间隔过大或无效时间戳：开始新的片段，时长为 0。
    streak_frames_ = 1;
    has_anchor_ = stamp_ok;
    first_stamp_ns_ = stamp_ok ? stamp_ns : 0;
    last_stamp_ns_ = first_stamp_ns_;
  }

  void Fill(OverlapObservation& out) const noexcept {
    out.streak_frames = streak_frames_;
    out.total_low_frames = total_low_frames_;
    out.capture_suppressed = capture_suppressed_;
    // last >= first > 0，差值非负且不溢出，因此时长不会为负。
    out.streak_seconds = (has_anchor_ && streak_frames_ > 0)
                             ? static_cast<double>(last_stamp_ns_ - first_stamp_ns_) * 1e-9
                             : 0.0;
  }

  std::uint64_t streak_frames_ = 0;
  std::uint64_t total_low_frames_ = 0;
  std::uint64_t capture_suppressed_ = 0;
  bool has_anchor_ = false;
  std::int64_t first_stamp_ns_ = 0;
  std::int64_t last_stamp_ns_ = 0;
  RateClock log_clock_;
  RateClock capture_clock_;
};

}  // namespace lightning::diagnostics
