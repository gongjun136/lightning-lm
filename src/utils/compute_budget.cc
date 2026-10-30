/**
 * @file
 * @brief 统一计算预算的配置转换、环境覆盖与提交前校验实现。
 * @see compute_budget.h
 */

#include "utils/compute_budget.h"

#include "utils/thread_scheduling.h"

#include <charconv>
#include <cstdlib>
#include <string_view>

#include <yaml-cpp/yaml.h>

namespace lightning::compute {
namespace {

constexpr int kMinThreads = 1;  ///< 工作线程数量的最小允许值。
constexpr int kMaxThreads = 128;  ///< 工作线程数量的最大允许值。
constexpr int kMinNice = 0;  ///< SOLiD nice 配置的最小允许值。
constexpr int kMaxNice = 19;  ///< SOLiD nice 配置的最大允许值。

/**
 * @brief 读取可选的 YAML 整数字段；缺少节点或键时保留传入值。
 * @param node compute_budget 配置节点，存在时应为可索引的映射。
 * @param key 字段名。
 * @param value 转换成功后接收整数值。
 * @param error 可选的转换错误说明输出，包含 compute_budget 字段路径。
 * @return 缺少字段或转换成功返回 true，捕获到 YAML 转换异常返回 false。
 */
bool ReadYamlValue(const YAML::Node& node, const char* key, int& value,
                   std::string* error) {
    if (!node || !node[key]) return true;
    try {
        value = node[key].as<int>();
        return true;
    } catch (const YAML::Exception& exception) {
        if (error) {
            *error = std::string("compute_budget.") + key + ": " + exception.what();
        }
        return false;
    }
}

/**
 * @brief 读取可选的整数环境覆盖，要求完整字符串符合 from_chars() 的整数格式。
 * @param name 环境变量名。
 * @param value 解析成功后接收覆盖值；未设置、空字符串或解析失败时保留原值。
 * @param error 可选的解析错误说明输出。
 * @return 未设置、空字符串或解析成功返回 true；格式错误或整数溢出返回 false。
 */
bool ReadEnvironmentValue(const char* name, int& value, std::string* error) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') return true;
    const std::string_view text(raw);
    int parsed = 0;
    const auto result = std::from_chars(text.data(), text.data() + text.size(), parsed);
    // 不接受前后空白或尾随字符；范围约束由最终预算校验统一处理。
    if (result.ec != std::errc() || result.ptr != text.data() + text.size()) {
        if (error) *error = std::string(name) + " must be an integer";
        return false;
    }
    value = parsed;
    return true;
}

/**
 * @brief 读取可选的 YAML 字符串字段，缺少字段时保留原值。
 * @param node compute_budget 配置节点，存在时应为可索引的映射。
 * @param key 字段名。
 * @param value 转换成功后接收字符串。
 * @param error 可选的 YAML 转换错误说明输出。
 * @return 缺少字段或转换成功返回 true，捕获到转换异常返回 false。
 */
bool ReadYamlString(const YAML::Node& node, const char* key, std::string& value,
                    std::string* error) {
    if (!node || !node[key]) return true;
    try {
        value = node[key].as<std::string>();
        return true;
    } catch (const YAML::Exception& exception) {
        if (error) {
            *error = std::string("compute_budget.") + key + ": " + exception.what();
        }
        return false;
    }
}

/**
 * @brief 读取字符串环境覆盖；显式空字符串也会覆盖原值，用于清除 CPU 列表。
 * @param name 环境变量名。
 * @param value 变量存在时被替换，未设置时保留原值。
 */
void ReadEnvironmentString(const char* name, std::string& value) {
    const char* raw = std::getenv(name);
    if (raw != nullptr) value = raw;
}

/**
 * @brief 检查线程数是否在 [1, 128] 内，失败时填写带配置名的可选错误说明。
 * @param name 配置字段名。
 * @param value 待检查的线程数。
 * @param error 可选的错误说明输出。
 * @return 线程数在允许范围内返回 true。
 */
bool ValidateValue(const char* name, int value, std::string* error) {
    if (value >= kMinThreads && value <= kMaxThreads) return true;
    if (error) {
        *error = std::string(name) + " must be in [" + std::to_string(kMinThreads) +
                 ", " + std::to_string(kMaxThreads) + "]";
    }
    return false;
}

/**
 * @brief 检查 SOLiD 工作线程的 nice 配置是否在 [0, 19] 内。
 * @param value 待检查的 nice 下限配置。
 * @param error 可选的错误说明输出。
 * @return 配置在允许范围内返回 true。
 */
bool ValidateNice(int value, std::string* error) {
    if (value >= kMinNice && value <= kMaxNice) return true;
    if (error) {
        *error = "compute_budget.solid_worker_nice must be in [" +
                 std::to_string(kMinNice) + ", " + std::to_string(kMaxNice) + "]";
    }
    return false;
}

/**
 * @brief 复用 CPU 列表解析器检查语法与编号范围，空列表表示不额外限制。
 * @param value 逗号分隔的 CPU 编号或闭区间列表。
 * @param error 可选的错误说明输出，补充 compute_budget 配置路径。
 * @return 列表可解析返回 true；此处不实际应用亲和性。
 * @see lightning::threading::ParseCpuList()
 */
bool ValidateCpuAffinity(const std::string& value, std::string* error) {
    std::vector<int> cpus;
    std::string parse_error;
    if (threading::ParseCpuList(value, cpus, &parse_error)) return true;
    if (error) {
        *error = "compute_budget.solid_cpu_affinity: " + parse_error;
    }
    return false;
}

}  // namespace

bool LoadComputeBudget(const YAML::Node& root, ComputeBudget& budget,
                       std::string* error) {
    // 阶段 1：加载候选配置。
    // 所有读取均落在副本上，任一转换失败或异常都不会留下部分更新的 budget。
    ComputeBudget candidate = budget;
    const YAML::Node config = root["compute_budget"];
    if (!ReadYamlValue(config, "lio_threads", candidate.lio_threads, error) ||
        !ReadYamlValue(config, "ndt_threads", candidate.ndt_threads, error) ||
        !ReadYamlValue(config, "ndt_max_points", candidate.ndt_max_points, error) ||
        !ReadYamlValue(config, "solid_icp_workers", candidate.solid_icp_workers, error) ||
        !ReadYamlValue(config, "solid_worker_nice", candidate.solid_worker_nice,
                       error) ||
        !ReadYamlString(config, "solid_cpu_affinity", candidate.solid_cpu_affinity,
                        error) ||
        !ReadEnvironmentValue(kLioThreadsEnv, candidate.lio_threads, error) ||
        !ReadEnvironmentValue(kNdtThreadsEnv, candidate.ndt_threads, error) ||
        !ReadEnvironmentValue(kNdtMaxPointsEnv, candidate.ndt_max_points, error) ||
        !ReadEnvironmentValue(kSolidIcpWorkersEnv, candidate.solid_icp_workers, error) ||
        !ReadEnvironmentValue(kSolidWorkerNiceEnv, candidate.solid_worker_nice,
                              error)) {
        return false;
    }
    // 阶段 2：校验最终预算。
    // 字符串环境变量允许显式空值覆盖；整数环境变量的空值则表示不覆盖。
    ReadEnvironmentString(kSolidCpuAffinityEnv, candidate.solid_cpu_affinity);
    if (candidate.ndt_max_points != 0 &&
        (candidate.ndt_max_points < 100 || candidate.ndt_max_points > 1000000)) {
        if (error) *error = "compute_budget.ndt_max_points must be 0 or in [100, 1000000]";
        return false;
    }
    if (!ValidateValue("compute_budget.lio_threads", candidate.lio_threads, error) ||
        !ValidateValue("compute_budget.ndt_threads", candidate.ndt_threads, error) ||
        !ValidateValue("compute_budget.solid_icp_workers",
                       candidate.solid_icp_workers, error) ||
        !ValidateNice(candidate.solid_worker_nice, error) ||
        !ValidateCpuAffinity(candidate.solid_cpu_affinity, error)) {
        return false;
    }
    // 阶段 3：提交预算。
    budget = candidate;
    return true;
}

}  // namespace lightning::compute
