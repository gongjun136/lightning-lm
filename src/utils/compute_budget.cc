#include "utils/compute_budget.h"

#include <charconv>
#include <cstdlib>
#include <string_view>

#include <yaml-cpp/yaml.h>

namespace lightning::compute {
namespace {

constexpr int kMinThreads = 1;
constexpr int kMaxThreads = 128;

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

bool ValidateValue(const char* name, int value, std::string* error) {
    if (value >= kMinThreads && value <= kMaxThreads) return true;
    if (error) {
        *error = std::string(name) + " must be in [" + std::to_string(kMinThreads) +
                 ", " + std::to_string(kMaxThreads) + "]";
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
        !ReadYamlValue(config, "solid_icp_workers", candidate.solid_icp_workers, error) ||
        !ReadEnvironmentValue(kLioThreadsEnv, candidate.lio_threads, error) ||
        !ReadEnvironmentValue(kNdtThreadsEnv, candidate.ndt_threads, error) ||
        !ReadEnvironmentValue(kSolidIcpWorkersEnv, candidate.solid_icp_workers, error) ||
        !ValidateValue("compute_budget.lio_threads", candidate.lio_threads, error) ||
        !ValidateValue("compute_budget.ndt_threads", candidate.ndt_threads, error) ||
        !ValidateValue("compute_budget.solid_icp_workers",
                       candidate.solid_icp_workers, error)) {
        return false;
    }
    budget = candidate;
    return true;
}

}  // namespace lightning::compute
