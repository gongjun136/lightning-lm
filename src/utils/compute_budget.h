#pragma once

#include <string>

namespace YAML {
class Node;
}

namespace lightning::compute {

constexpr const char* kLioThreadsEnv = "LIGHTNING_LM_LIO_THREADS";
constexpr const char* kNdtThreadsEnv = "LIGHTNING_LM_NDT_THREADS";
constexpr const char* kNdtMaxPointsEnv = "LIGHTNING_LM_NDT_MAX_POINTS";
constexpr const char* kSolidIcpWorkersEnv = "LIGHTNING_LM_SOLID_ICP_WORKERS";
constexpr const char* kSolidWorkerNiceEnv = "LIGHTNING_LM_SOLID_WORKER_NICE";
constexpr const char* kSolidCpuAffinityEnv = "LIGHTNING_LM_SOLID_CPU_AFFINITY";

// A single configuration object owns every high-cost online localization
// worker pool. Callers may seed legacy-compatible defaults before loading it.
struct ComputeBudget {
    int lio_threads = 1;
    int ndt_threads = 4;
    int solid_icp_workers = 8;
    int solid_worker_nice = 0;
    std::string solid_cpu_affinity;
    int ndt_max_points = 0;  // 0: uncapped; tracking only, never global initialization
};

// Read optional compute_budget.* values and then apply environment overrides.
// Thread counts must be in [1, 128], SOLiD nice in [0, 19], and affinity must
// use Linux taskset CPU-list syntax. On failure, budget is unchanged.
bool LoadComputeBudget(const YAML::Node& root, ComputeBudget& budget,
                       std::string* error = nullptr);

}  // namespace lightning::compute
