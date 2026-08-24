#pragma once

#include <string>

namespace YAML {
class Node;
}

namespace lightning::compute {

constexpr const char* kLioThreadsEnv = "LIGHTNING_LM_LIO_THREADS";
constexpr const char* kNdtThreadsEnv = "LIGHTNING_LM_NDT_THREADS";
constexpr const char* kSolidIcpWorkersEnv = "LIGHTNING_LM_SOLID_ICP_WORKERS";

// A single configuration object owns every high-cost online localization
// worker pool. Callers may seed legacy-compatible defaults before loading it.
struct ComputeBudget {
    int lio_threads = 1;
    int ndt_threads = 4;
    int solid_icp_workers = 8;
};

// Read optional compute_budget.* values and then apply environment overrides.
// Every configured value must be in [1, 128]. On failure, budget is unchanged.
bool LoadComputeBudget(const YAML::Node& root, ComputeBudget& budget,
                       std::string* error = nullptr);

}  // namespace lightning::compute
