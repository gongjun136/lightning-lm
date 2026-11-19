// translation_evidence.h
// 平移一致性检查的诊断证据(仅诊断用途;C++17,仅依赖标准库)。
//
// 本头文件只记录调用方已完成的计算与判定:
//   * 不计算阈值,不重新比较误差与阈值,不选择安全判定;
//   * 不改变积分、调度、锁存、恢复或任何防抖行为;
//   * RecordTranslationComparison 不分配内存、不做 I/O、不抛异常。
// 格式化函数只应在现有有界安全日志调度或后台事件写入器中、在监视器
// 互斥锁之外调用,不得在每次比较时调用。
#ifndef MOTION_DIAGNOSTICS_TRANSLATION_EVIDENCE_H_
#define MOTION_DIAGNOSTICS_TRANSLATION_EVIDENCE_H_

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <ios>
#include <limits>
#include <locale>
#include <ostream>

namespace motion_diagnostics {

inline constexpr int kTranslationEvidenceVersion = 1;

// selected_window:-1 无判定区间,0 相邻区间,1 窗口(较长)区间。
inline constexpr int kSelectedWindowNone = -1;
inline constexpr int kSelectedWindowAdjacent = 0;
inline constexpr int kSelectedWindowLong = 1;

inline constexpr std::size_t kTranslationSlotCount = 4;

// 缺失数值证据的默认值:NaN,而不是可能被误读的 0。
inline constexpr double kTranslationMissing = std::numeric_limits<double>::quiet_NaN();

// 一次已执行比较的原样记录;不包含安全判定的选择。
struct TranslationComparison {
  bool available = false;
  bool window = false;  // false:相邻区间;true:窗口区间
  bool raw = false;     // false:发布(最终平滑)分支;true:未平滑分支
  // 仅对未平滑分支有意义:两个端点都有真实的未平滑输出。
  // 这不表示统计独立;为 false 表示至少一个端点回退为发布位姿。
  bool raw_distinct = false;
  double from_stamp = kTranslationMissing;     // 区间起点源时间戳 [s]
  double to_stamp = kTranslationMissing;       // 区间终点源时间戳 [s]
  double distance_m = kTranslationMissing;     // 本区间行驶长度 D(非净位移)
  double predicted_x_m = kTranslationMissing;  // 轮速/陀螺积分 dx(起点机体平面)
  double predicted_y_m = kTranslationMissing;  // 轮速/陀螺积分 dy
  double measured_x_m = kTranslationMissing;   // 位姿位移按该分支起点航向反旋转后的 x
  double measured_y_m = kTranslationMissing;   // 同上 y
  double residual_x_m = kTranslationMissing;   // x - dx,带符号,仅诊断
  double residual_y_m = kTranslationMissing;   // y - dy,带符号,仅诊断
  double error_x_m = kTranslationMissing;      // 调用方计算的 ex(绝对值)
  double error_y_m = kTranslationMissing;      // 调用方计算的 ey(绝对值)
  double limit_x_m = kTranslationMissing;      // 调用方实际使用的 X 阈值
  double limit_y_m = kTranslationMissing;      // 调用方实际使用的 Y 阈值
  bool fail_x = false;  // 调用方原有比较结果,此处不重新计算
  bool fail_y = false;
};

// 槽位固定:0 相邻/发布,1 相邻/未平滑,2 窗口/发布,3 窗口/未平滑。
// 保存全部可用比较,即使相邻失败已决定判定而窗口结果不同。
struct TranslationEvidence {
  std::array<TranslationComparison, kTranslationSlotCount> comparisons{};
  std::uint64_t sequence = 0;
  double evaluation_stamp = kTranslationMissing;
  int selected_window = kSelectedWindowNone;
  bool latched = false;  // true:展示的是保留的最近一次失败证据,而非新比较
};

inline constexpr std::size_t TranslationSlotIndex(bool window, bool raw) noexcept {
  return (window ? std::size_t{2} : std::size_t{0}) + (raw ? std::size_t{1} : std::size_t{0});
}

// 开始一次新的位姿评估:清空证据、序号加一、记录评估时刻。
inline void BeginTranslationEvaluation(TranslationEvidence& evidence,
                                       double evaluation_stamp) noexcept {
  const std::uint64_t next_sequence = evidence.sequence + 1u;
  evidence = TranslationEvidence{};
  evidence.sequence = next_sequence;
  evidence.evaluation_stamp = evaluation_stamp;
}

// 原样记录调用方已完成的一次比较。不重新计算或修改调用方判定;
// 带符号残差 x-dx、y-dy 仅用于诊断。无动态分配、无 I/O、不抛异常。
inline void RecordTranslationComparison(TranslationEvidence& evidence, bool window, bool raw,
                                        bool raw_distinct, double from, double to,
                                        double distance, double dx, double dy, double x,
                                        double y, double ex, double ey, double limit_x,
                                        double limit_y, bool fail_x, bool fail_y) noexcept {
  TranslationComparison& c = evidence.comparisons[TranslationSlotIndex(window, raw)];
  c.available = true;
  c.window = window;
  c.raw = raw;
  c.raw_distinct = raw_distinct;
  c.from_stamp = from;
  c.to_stamp = to;
  c.distance_m = distance;
  c.predicted_x_m = dx;
  c.predicted_y_m = dy;
  c.measured_x_m = x;
  c.measured_y_m = y;
  c.residual_x_m = x - dx;
  c.residual_y_m = y - dy;
  c.error_x_m = ex;
  c.error_y_m = ey;
  c.limit_x_m = limit_x;
  c.limit_y_m = limit_y;
  c.fail_x = fail_x;
  c.fail_y = fail_y;
}

// 所选区间至少有一个可用比较时为 true(与其他区间是否可用无关)。
inline bool TranslationEvidenceAvailable(const TranslationEvidence& evidence) noexcept {
  if (evidence.selected_window != kSelectedWindowAdjacent &&
      evidence.selected_window != kSelectedWindowLong) {
    return false;
  }
  const bool window = evidence.selected_window == kSelectedWindowLong;
  return evidence.comparisons[TranslationSlotIndex(window, false)].available ||
         evidence.comparisons[TranslationSlotIndex(window, true)].available;
}

inline const char* TranslationSelectedWindowName(int selected_window) noexcept {
  switch (selected_window) {
    case kSelectedWindowAdjacent:
      return "adjacent";
    case kSelectedWindowLong:
      return "window";
    default:
      return "none";
  }
}

inline const char* TranslationEvidenceOriginName(const TranslationEvidence& evidence) noexcept {
  return evidence.latched ? "latched_failure" : "current";
}

// 选择发布快照中的平移证据(仅诊断,不改变锁存状态):
// 平移锁存有效且当前判定不是 FAIL 时,返回保留的最近失败证据并标记 latched=true;
// 其余情况返回当前证据并标记 latched=false。保留证据维持原序号与评估时刻。
inline TranslationEvidence SelectTranslationSnapshot(const TranslationEvidence& current,
                                                     const TranslationEvidence& last_failure,
                                                     bool translation_latch_active,
                                                     bool current_is_fail) noexcept {
  const bool use_retained =
      translation_latch_active && !current_is_fail && last_failure.sequence != 0u;
  TranslationEvidence out = use_retained ? last_failure : current;
  out.latched = use_retained;
  return out;
}

namespace detail {

inline constexpr const char* const kTranslationSlotPrefixes[kTranslationSlotCount] = {
    "translation_adjacent_final_", "translation_adjacent_raw_", "translation_window_final_",
    "translation_window_raw_"};

enum class TranslationOutputMode { kKeyValue, kCsvHeader, kCsvValues };

// 在作用域内使用固定格式(十进制、17 位有效数字、经典区域设置),
// 退出时恢复调用方的标志、精度、宽度与区域设置。
class TranslationStreamGuard {
 public:
  explicit TranslationStreamGuard(std::ostream& os)
      : os_(os),
        flags_(os.flags()),
        precision_(os.precision()),
        width_(os.width()),
        locale_(os.getloc()) {
    os_.flags(std::ios_base::dec);
    os_.precision(std::numeric_limits<double>::max_digits10);
    os_.width(0);
    os_.imbue(std::locale::classic());
  }
  ~TranslationStreamGuard() {
    os_.imbue(locale_);
    os_.width(width_);
    os_.precision(precision_);
    os_.flags(flags_);
  }
  TranslationStreamGuard(const TranslationStreamGuard&) = delete;
  TranslationStreamGuard& operator=(const TranslationStreamGuard&) = delete;

 private:
  std::ostream& os_;
  std::ios_base::fmtflags flags_;
  std::streamsize precision_;
  std::streamsize width_;
  std::locale locale_;
};

class TranslationFieldEmitter {
 public:
  TranslationFieldEmitter(std::ostream& os, TranslationOutputMode mode) : os_(os), mode_(mode) {}

  void Text(const char* prefix, const char* name, const char* value) {
    if (Begin(prefix, name)) os_ << value;
  }
  void Flag(const char* prefix, const char* name, bool value) {
    if (Begin(prefix, name)) os_ << (value ? '1' : '0');
  }
  void Unsigned(const char* prefix, const char* name, std::uint64_t value) {
    if (Begin(prefix, name)) os_ << value;
  }
  void Integer(const char* prefix, const char* name, int value) {
    if (Begin(prefix, name)) os_ << value;
  }
  // 非有限值写为固定文本 nan / inf / -inf,不依赖平台格式。
  void Number(const char* prefix, const char* name, double value) {
    if (!Begin(prefix, name)) return;
    if (std::isnan(value)) {
      os_ << "nan";
    } else if (std::isinf(value)) {
      os_ << (value > 0.0 ? "inf" : "-inf");
    } else {
      os_ << value;
    }
  }

 private:
  // 写出字段前缀;返回 true 表示接着需要写值。
  bool Begin(const char* prefix, const char* name) {
    switch (mode_) {
      case TranslationOutputMode::kKeyValue:
        os_ << ' ' << prefix << name << '=';
        return true;
      case TranslationOutputMode::kCsvHeader:
        os_ << ',' << prefix << name;
        return false;
      case TranslationOutputMode::kCsvValues:
        os_ << ',';
        return true;
    }
    return false;
  }

  std::ostream& os_;
  TranslationOutputMode mode_;
};

// 三种输出共用同一字段序列,保证表头、数值与 key=value 一一对应。
inline void EmitTranslationEvidence(std::ostream& os, const TranslationEvidence& evidence,
                                    TranslationOutputMode mode) {
  TranslationStreamGuard guard(os);
  TranslationFieldEmitter out(os, mode);
  out.Integer("", "translation_evidence_version", kTranslationEvidenceVersion);
  out.Unsigned("", "translation_evidence_sequence", evidence.sequence);
  out.Number("", "translation_evaluation_stamp", evidence.evaluation_stamp);
  out.Text("", "translation_selected_window",
           TranslationSelectedWindowName(evidence.selected_window));
  out.Text("", "translation_evidence_origin", TranslationEvidenceOriginName(evidence));
  out.Flag("", "translation_evidence_available", TranslationEvidenceAvailable(evidence));
  for (std::size_t i = 0; i < kTranslationSlotCount; ++i) {
    const char* p = kTranslationSlotPrefixes[i];
    const TranslationComparison& c = evidence.comparisons[i];
    out.Flag(p, "available", c.available);
    out.Flag(p, "raw_distinct", c.raw_distinct);
    out.Number(p, "from_stamp", c.from_stamp);
    out.Number(p, "to_stamp", c.to_stamp);
    out.Number(p, "dt_sec", c.available ? c.to_stamp - c.from_stamp : kTranslationMissing);
    out.Number(p, "distance_m", c.distance_m);
    out.Number(p, "predicted_x_m", c.predicted_x_m);
    out.Number(p, "predicted_y_m", c.predicted_y_m);
    out.Number(p, "measured_x_m", c.measured_x_m);
    out.Number(p, "measured_y_m", c.measured_y_m);
    out.Number(p, "residual_x_m", c.residual_x_m);
    out.Number(p, "residual_y_m", c.residual_y_m);
    out.Number(p, "error_x_m", c.error_x_m);
    out.Number(p, "error_y_m", c.error_y_m);
    out.Number(p, "limit_x_m", c.limit_x_m);
    out.Number(p, "limit_y_m", c.limit_y_m);
    out.Flag(p, "fail_x", c.fail_x);
    out.Flag(p, "fail_y", c.fail_y);
  }
}

}  // namespace detail

// 追加以空格开头的 key=value 字段(不换行)。
inline void WriteTranslationFields(std::ostream& os, const TranslationEvidence& evidence) {
  detail::EmitTranslationEvidence(os, evidence, detail::TranslationOutputMode::kKeyValue);
}

// 追加以逗号开头的 CSV 列名(不换行),与 WriteTranslationCsvValues 一一对应。
inline void WriteTranslationCsvHeader(std::ostream& os) {
  const TranslationEvidence empty{};
  detail::EmitTranslationEvidence(os, empty, detail::TranslationOutputMode::kCsvHeader);
}

// 追加以逗号开头的 CSV 数值(不换行)。
inline void WriteTranslationCsvValues(std::ostream& os, const TranslationEvidence& evidence) {
  detail::EmitTranslationEvidence(os, evidence, detail::TranslationOutputMode::kCsvValues);
}

}  // namespace motion_diagnostics

#endif  // MOTION_DIAGNOSTICS_TRANSLATION_EVIDENCE_H_
