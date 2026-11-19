// translation_evidence_test.cc
// 平移诊断证据的确定性单元测试:独立可执行程序,不依赖测试框架。
// 只验证证据记录与序列化契约,不能替代整车回放验证。
#include "utils/translation_evidence.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <sstream>
#include <string>
#include <vector>

namespace md = motion_diagnostics;

namespace {

int g_failures = 0;

void Check(bool ok, const char* what, int line) {
  if (!ok) {
    ++g_failures;
    std::cerr << "检查失败(第 " << line << " 行):" << what << std::endl;
  }
}

#define CHECK(cond) Check((cond), #cond, __LINE__)

constexpr double kT0 = 1767225600.1234567;  // 合成的类纪元时间戳

// 以逗号作小数点、点作千位分隔符的区域设置,用于验证输出不受调用方区域设置影响。
struct CommaDecimal : std::numpunct<char> {
  char do_decimal_point() const override { return ','; }
  char do_thousands_sep() const override { return '.'; }
  std::string do_grouping() const override { return std::string(1, static_cast<char>(3)); }
};

std::string KeyValues(const md::TranslationEvidence& e) {
  std::ostringstream os;
  md::WriteTranslationFields(os, e);
  return os.str();
}

std::string CsvHeader() {
  std::ostringstream os;
  md::WriteTranslationCsvHeader(os);
  return os.str();
}

std::string CsvValues(const md::TranslationEvidence& e) {
  std::ostringstream os;
  md::WriteTranslationCsvValues(os, e);
  return os.str();
}

bool HasNewline(const std::string& text) { return text.find('\n') != std::string::npos; }

// 输出以分隔符开头,因此丢弃首个空片段。
std::vector<std::string> SplitPrefixed(const std::string& text, char sep) {
  std::vector<std::string> parts;
  if (text.empty() || text[0] != sep) return parts;
  std::size_t start = 1;
  while (true) {
    const std::size_t next = text.find(sep, start);
    parts.push_back(
        text.substr(start, next == std::string::npos ? std::string::npos : next - start));
    if (next == std::string::npos) break;
    start = next + 1;
  }
  return parts;
}

std::string ValueOf(const std::string& kv, const std::string& key) {
  const std::string needle = " " + key + "=";
  const std::size_t pos = kv.find(needle);
  if (pos == std::string::npos) return "<missing>";
  const std::size_t start = pos + needle.size();
  const std::size_t end = kv.find(' ', start);
  return kv.substr(start, end == std::string::npos ? std::string::npos : end - start);
}

double ParseDouble(const std::string& text) { return std::strtod(text.c_str(), nullptr); }

// 模拟调用方:误差为绝对残差,判定为严格大于阈值(相等通过)。仅测试用。
void RecordLikeCaller(md::TranslationEvidence& e, bool window, bool raw, bool raw_distinct,
                      double from, double to, double distance, double dx, double dy, double x,
                      double y, double limit_x, double limit_y) {
  const double ex = std::fabs(x - dx);
  const double ey = std::fabs(y - dy);
  md::RecordTranslationComparison(e, window, raw, raw_distinct, from, to, distance, dx, dy, x, y,
                                  ex, ey, limit_x, limit_y, ex > limit_x, ey > limit_y);
}

md::TranslationEvidence FullEvidence() {
  md::TranslationEvidence e;
  md::BeginTranslationEvaluation(e, kT0 + 0.5);
  RecordLikeCaller(e, false, false, false, kT0 + 0.4, kT0 + 0.5, 0.1, 0.1, 0.0, 0.75, 0.0, 0.51,
                   0.31);
  RecordLikeCaller(e, false, true, true, kT0 + 0.4, kT0 + 0.5, 0.1, 0.1, 0.0, 0.12, 0.0, 0.51,
                   0.31);
  RecordLikeCaller(e, true, false, false, kT0, kT0 + 0.5, 0.5, 0.5, 0.0, 0.52, 0.01, 0.55, 0.35);
  RecordLikeCaller(e, true, true, false, kT0, kT0 + 0.5, 0.5, 0.5, 0.0, 0.53, 0.01, 0.55, 0.35);
  e.selected_window = md::kSelectedWindowAdjacent;
  return e;
}

void TestDefaultsAreMissingNotZero() {
  const md::TranslationEvidence e{};
  CHECK(e.sequence == 0u);
  CHECK(std::isnan(e.evaluation_stamp));
  CHECK(e.selected_window == md::kSelectedWindowNone);
  CHECK(!e.latched);
  for (const md::TranslationComparison& c : e.comparisons) {
    CHECK(!c.available);
    CHECK(std::isnan(c.from_stamp));
    CHECK(std::isnan(c.to_stamp));
    CHECK(std::isnan(c.distance_m));
    CHECK(std::isnan(c.error_x_m));
    CHECK(std::isnan(c.limit_x_m));
    CHECK(std::isnan(c.limit_y_m));
    CHECK(!c.fail_x && !c.fail_y);
  }
  const std::string kv = KeyValues(e);
  CHECK(ValueOf(kv, "translation_evidence_version") == "1");
  CHECK(ValueOf(kv, "translation_evidence_available") == "0");
  CHECK(ValueOf(kv, "translation_selected_window") == "none");
  CHECK(ValueOf(kv, "translation_evidence_origin") == "current");
  CHECK(ValueOf(kv, "translation_adjacent_final_limit_x_m") == "nan");
  CHECK(ValueOf(kv, "translation_window_raw_dt_sec") == "nan");
  CHECK(kv.find("limit_x_m=0") == std::string::npos);
  CHECK(kv.find("limit_y_m=0") == std::string::npos);
}

void TestRecordIsNoexcept() {
  md::TranslationEvidence e;
  static_assert(noexcept(md::RecordTranslationComparison(e, false, false, false, 0.0, 0.0, 0.0,
                                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                         0.0, false, false)),
                "RecordTranslationComparison 必须为 noexcept");
  static_assert(noexcept(md::BeginTranslationEvaluation(e, 0.0)),
                "BeginTranslationEvaluation 必须为 noexcept");
  CHECK(!e.comparisons[0].available);
}

void TestRecordStoresCallerValues() {
  md::TranslationEvidence e;
  md::BeginTranslationEvaluation(e, kT0 + 0.5);
  md::RecordTranslationComparison(e, true, true, true, kT0, kT0 + 0.5, 0.8, 0.7, 0.05, 1.4, -0.3,
                                  0.7, 0.35, 0.58, 0.38, true, false);
  const md::TranslationComparison& c = e.comparisons[3];
  CHECK(c.available && c.window && c.raw && c.raw_distinct);
  CHECK(c.from_stamp == kT0);
  CHECK(c.to_stamp == kT0 + 0.5);
  CHECK(c.distance_m == 0.8);
  CHECK(c.predicted_x_m == 0.7);
  CHECK(c.predicted_y_m == 0.05);
  CHECK(c.measured_x_m == 1.4);
  CHECK(c.measured_y_m == -0.3);
  CHECK(c.residual_x_m == 1.4 - 0.7);
  CHECK(c.residual_y_m == -0.3 - 0.05);
  CHECK(c.error_x_m == 0.7);
  CHECK(c.error_y_m == 0.35);
  CHECK(c.limit_x_m == 0.58);
  CHECK(c.limit_y_m == 0.38);
  CHECK(c.fail_x);
  CHECK(!c.fail_y);
  for (std::size_t i = 0; i < 3; ++i) CHECK(!e.comparisons[i].available);
}

// 辅助函数只存储调用方判定;刻意构造的矛盾值也原样保存,由离线解析报告矛盾。
void TestCallerVerdictIsNotRecomputed() {
  md::TranslationEvidence e;
  md::BeginTranslationEvaluation(e, kT0);
  md::RecordTranslationComparison(e, false, false, false, kT0, kT0 + 0.1, 0.1, 0.0, 0.0, 0.6, 0.1,
                                  0.6, 0.1, 0.5, 0.3, false, true);
  CHECK(!e.comparisons[0].fail_x);
  CHECK(e.comparisons[0].fail_y);
  const std::string kv = KeyValues(e);
  CHECK(ValueOf(kv, "translation_adjacent_final_fail_x") == "0");
  CHECK(ValueOf(kv, "translation_adjacent_final_fail_y") == "1");
}

void TestSlotMapping() {
  CHECK(md::TranslationSlotIndex(false, false) == 0u);
  CHECK(md::TranslationSlotIndex(false, true) == 1u);
  CHECK(md::TranslationSlotIndex(true, false) == 2u);
  CHECK(md::TranslationSlotIndex(true, true) == 3u);
  md::TranslationEvidence e;
  md::BeginTranslationEvaluation(e, kT0);
  RecordLikeCaller(e, false, true, false, kT0, kT0 + 0.1, 0.1, 0.1, 0.0, 0.1, 0.0, 0.51, 0.31);
  CHECK(!e.comparisons[0].available);
  CHECK(e.comparisons[1].available && e.comparisons[1].raw && !e.comparisons[1].window);
  CHECK(!e.comparisons[1].raw_distinct);  // 回退:至少一个端点使用发布位姿
  const std::string kv = KeyValues(e);
  CHECK(ValueOf(kv, "translation_adjacent_raw_available") == "1");
  CHECK(ValueOf(kv, "translation_adjacent_raw_raw_distinct") == "0");
  CHECK(ValueOf(kv, "translation_adjacent_final_available") == "0");
}

void TestEqualityPassesAsRecorded() {
  md::TranslationEvidence e;
  md::BeginTranslationEvaluation(e, kT0);
  RecordLikeCaller(e, false, false, false, kT0, kT0 + 0.1, 0.0, 0.0, 0.0, 0.5, 0.3, 0.5, 0.3);
  CHECK(!e.comparisons[0].fail_x);
  CHECK(!e.comparisons[0].fail_y);
  const std::string kv = KeyValues(e);
  CHECK(ValueOf(kv, "translation_adjacent_final_error_x_m") ==
        ValueOf(kv, "translation_adjacent_final_limit_x_m"));
  CHECK(ValueOf(kv, "translation_adjacent_final_fail_x") == "0");
  CHECK(ValueOf(kv, "translation_adjacent_final_fail_y") == "0");
}

void TestReverseMotionSigned() {
  md::TranslationEvidence e;
  md::BeginTranslationEvaluation(e, kT0);
  RecordLikeCaller(e, false, false, false, kT0, kT0 + 0.1, 0.4, -0.4, 0.0, -0.45, 0.0, 0.54, 0.34);
  const md::TranslationComparison& c = e.comparisons[0];
  CHECK(c.predicted_x_m < 0.0);
  CHECK(c.residual_x_m < 0.0);
  CHECK(c.error_x_m > 0.0);
  const std::string kv = KeyValues(e);
  CHECK(ValueOf(kv, "translation_adjacent_final_predicted_x_m").rfind("-", 0) == 0);
  CHECK(ValueOf(kv, "translation_adjacent_final_residual_x_m").rfind("-", 0) == 0);
}

void TestBeginEvaluation() {
  md::TranslationEvidence e;
  md::BeginTranslationEvaluation(e, 10.0);
  CHECK(e.sequence == 1u);
  RecordLikeCaller(e, true, false, false, 9.5, 10.0, 0.5, 0.5, 0.0, 0.5, 0.0, 0.55, 0.35);
  e.selected_window = md::kSelectedWindowLong;
  md::BeginTranslationEvaluation(e, 10.1);
  CHECK(e.sequence == 2u);
  CHECK(e.evaluation_stamp == 10.1);
  CHECK(e.selected_window == md::kSelectedWindowNone);
  CHECK(!e.latched);
  for (const md::TranslationComparison& c : e.comparisons) CHECK(!c.available);
}

void TestAvailabilityFollowsSelectedInterval() {
  md::TranslationEvidence e;
  md::BeginTranslationEvaluation(e, kT0);
  RecordLikeCaller(e, true, false, false, kT0 - 0.5, kT0, 0.5, 0.5, 0.0, 0.5, 0.0, 0.55, 0.35);
  e.selected_window = md::kSelectedWindowAdjacent;  // 所选区间没有比较
  CHECK(!md::TranslationEvidenceAvailable(e));
  CHECK(ValueOf(KeyValues(e), "translation_evidence_available") == "0");
  e.selected_window = md::kSelectedWindowLong;
  CHECK(md::TranslationEvidenceAvailable(e));
  CHECK(ValueOf(KeyValues(e), "translation_selected_window") == "window");
  CHECK(ValueOf(KeyValues(e), "translation_evidence_available") == "1");
  e.selected_window = md::kSelectedWindowNone;
  CHECK(!md::TranslationEvidenceAvailable(e));
  e.selected_window = 5;  // 非法值按 none 输出
  CHECK(ValueOf(KeyValues(e), "translation_selected_window") == "none");
  CHECK(ValueOf(KeyValues(e), "translation_evidence_available") == "0");
}

void TestSnapshotRetainsLatchedFailure() {
  md::TranslationEvidence current;
  for (int i = 0; i < 5; ++i) md::BeginTranslationEvaluation(current, kT0 + 0.1 * i);
  RecordLikeCaller(current, false, false, false, kT0 + 0.3, kT0 + 0.4, 0.1, 0.1, 0.0, 0.75, 0.0,
                   0.51, 0.31);
  current.selected_window = md::kSelectedWindowAdjacent;
  const md::TranslationEvidence last_failure = current;  // 序号 5
  md::BeginTranslationEvaluation(current, kT0 + 0.6);    // 序号 6:无比较(UNKNOWN)

  const md::TranslationEvidence latched =
      md::SelectTranslationSnapshot(current, last_failure, true, false);
  CHECK(latched.latched);
  CHECK(latched.sequence == 5u);
  CHECK(latched.evaluation_stamp == last_failure.evaluation_stamp);
  const std::string kv = KeyValues(latched);
  CHECK(ValueOf(kv, "translation_evidence_origin") == "latched_failure");
  CHECK(ValueOf(kv, "translation_evidence_sequence") == "5");
  CHECK(ValueOf(kv, "translation_adjacent_final_fail_x") == "1");

  const md::TranslationEvidence cleared =
      md::SelectTranslationSnapshot(current, last_failure, false, false);
  CHECK(!cleared.latched);
  CHECK(cleared.sequence == 6u);
  CHECK(ValueOf(KeyValues(cleared), "translation_evidence_origin") == "current");

  const md::TranslationEvidence failing_now =
      md::SelectTranslationSnapshot(current, last_failure, true, true);
  CHECK(!failing_now.latched);
  CHECK(failing_now.sequence == 6u);
}

void TestCsvParityAndDeterminism() {
  const md::TranslationEvidence e = FullEvidence();
  const std::string header_text = CsvHeader();
  const std::string values_text = CsvValues(e);
  const std::string kv_text = KeyValues(e);
  CHECK(!HasNewline(header_text));
  CHECK(!HasNewline(values_text));
  CHECK(!HasNewline(kv_text));
  CHECK(header_text == CsvHeader());
  CHECK(values_text == CsvValues(e));
  const std::vector<std::string> header = SplitPrefixed(header_text, ',');
  const std::vector<std::string> values = SplitPrefixed(values_text, ',');
  const std::vector<std::string> pairs = SplitPrefixed(kv_text, ' ');
  CHECK(header.size() == 6u + 4u * 18u);
  CHECK(values.size() == header.size());
  CHECK(pairs.size() == header.size());
  if (values.size() == header.size() && pairs.size() == header.size()) {
    for (std::size_t i = 0; i < header.size(); ++i) {
      const std::size_t eq = pairs[i].find('=');
      CHECK(eq != std::string::npos);
      if (eq == std::string::npos) continue;
      CHECK(pairs[i].substr(0, eq) == header[i]);
      CHECK(pairs[i].substr(eq + 1) == values[i]);
    }
  }
  CHECK(!header.empty() && header.front() == "translation_evidence_version");
  CHECK(!header.empty() && header.back() == "translation_window_raw_fail_y");
}

void TestStreamStateRestored() {
  std::ostringstream os;
  const std::locale comma(std::locale::classic(), new CommaDecimal);
  os.imbue(comma);
  os << std::fixed << std::showpos << std::uppercase << std::setprecision(3);
  os.width(9);
  const std::ios_base::fmtflags flags = os.flags();
  const std::streamsize precision = os.precision();
  const std::streamsize width = os.width();
  const md::TranslationEvidence e = FullEvidence();

  md::WriteTranslationFields(os, e);
  CHECK(os.flags() == flags);
  CHECK(os.precision() == precision);
  CHECK(os.width() == width);
  CHECK(os.getloc() == comma);
  const std::string kv = os.str();
  CHECK(ValueOf(kv, "translation_evidence_sequence") == "1");  // 无 '+',无千位分隔
  CHECK(ParseDouble(ValueOf(kv, "translation_evaluation_stamp")) == kT0 + 0.5);
  CHECK(ValueOf(kv, "translation_adjacent_final_limit_x_m").find(',') == std::string::npos);

  md::WriteTranslationCsvHeader(os);
  md::WriteTranslationCsvValues(os, e);
  CHECK(os.flags() == flags);
  CHECK(os.precision() == precision);
  CHECK(os.width() == width);
  CHECK(os.getloc() == comma);
}

void TestTimestampRoundTrip() {
  md::TranslationEvidence e;
  const double stamp = std::nextafter(kT0 + 0.5, std::numeric_limits<double>::infinity());
  md::BeginTranslationEvaluation(e, stamp);
  const double from = kT0;
  const double to = std::nextafter(kT0 + 0.1, std::numeric_limits<double>::infinity());
  RecordLikeCaller(e, false, false, false, from, to, 0.1, 0.1, 0.0, 0.1, 0.0, 0.51, 0.31);
  const std::string kv = KeyValues(e);
  CHECK(ParseDouble(ValueOf(kv, "translation_adjacent_final_from_stamp")) == from);
  CHECK(ParseDouble(ValueOf(kv, "translation_adjacent_final_to_stamp")) == to);
  CHECK(ParseDouble(ValueOf(kv, "translation_evaluation_stamp")) == stamp);
  CHECK(ParseDouble(ValueOf(kv, "translation_adjacent_final_dt_sec")) == to - from);
}

void TestNonfiniteInputs() {
  md::TranslationEvidence e;
  md::BeginTranslationEvaluation(e, kT0);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  md::RecordTranslationComparison(e, false, false, false, kT0, kT0 + 0.1, 0.1, 0.1, 0.0, 0.1, 0.0,
                                  nan, -inf, inf, 0.31, false, false);
  const std::string kv = KeyValues(e);
  CHECK(ValueOf(kv, "translation_adjacent_final_error_x_m") == "nan");
  CHECK(ValueOf(kv, "translation_adjacent_final_error_y_m") == "-inf");
  CHECK(ValueOf(kv, "translation_adjacent_final_limit_x_m") == "inf");
}

}  // namespace

int main() {
  TestDefaultsAreMissingNotZero();
  TestRecordIsNoexcept();
  TestRecordStoresCallerValues();
  TestCallerVerdictIsNotRecomputed();
  TestSlotMapping();
  TestEqualityPassesAsRecorded();
  TestReverseMotionSigned();
  TestBeginEvaluation();
  TestAvailabilityFollowsSelectedInterval();
  TestSnapshotRetainsLatchedFailure();
  TestCsvParityAndDeterminism();
  TestStreamStateRestored();
  TestTimestampRoundTrip();
  TestNonfiniteInputs();
  if (g_failures == 0) {
    std::cout << "translation_evidence_test:全部检查通过" << std::endl;
    return EXIT_SUCCESS;
  }
  std::cout << "translation_evidence_test:失败 " << g_failures << " 项" << std::endl;
  return EXIT_FAILURE;
}
