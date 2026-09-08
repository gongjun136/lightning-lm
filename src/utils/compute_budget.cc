#include "utils/compute_budget.h"

#include "utils/thread_scheduling.h"

#include <charconv>
#include <cstdlib>
#include <string_view>

#include <yaml-cpp/yaml.h>

namespace lightning::compute {
namespace {

constexpr int kMinThreads = 1;
constexpr int kMaxThreads = 128;
constexpr int kMinNice = 0;
constexpr int kMaxNice = 19;

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

bool ReadEnvironmentValue(const char* name, int& value, std::string* error) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') return true;
    const std::string_view text(raw);
    int parsed = 0;
    const auto result = std::from_chars(text.data(), text.data() + text.size(), parsed);
    if (result.ec != std::errc() || result.ptr != text.data() + text.size()) {
        if (error) *error = std::string(name) + " must be an integer";
        return false;
    }
    value = parsed;
    return true;
}

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

void ReadEnvironmentString(const char* name, std::string& value) {
    const char* raw = std::getenv(name);
    if (raw != nullptr) value = raw;
}

bool ValidateValue(const char* name, int value, std::string* error) {
    if (value >= kMinThreads && value <= kMaxThreads) return true;
    if (error) {
        *error = std::string(name) + " must be in [" + std::to_string(kMinThreads) +
                 ", " + std::to_string(kMaxThreads) + "]";
    }
    return false;
}

bool ValidateNice(int value, std::string* error) {
    if (value >= kMinNice && value <= kMaxNice) return true;
    if (error) {
        *error = "compute_budget.solid_worker_nice must be in [" +
                 std::to_string(kMinNice) + ", " + std::to_string(kMaxNice) + "]";
    }
    return false;
}

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
    budget = candidate;
    return true;
}

}  // namespace lightning::compute
