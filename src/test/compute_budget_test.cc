#include "utils/compute_budget.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include <yaml-cpp/yaml.h>

namespace {

void Require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

class ScopedEnvironment {
   public:
    ScopedEnvironment(const char* name, const char* value) : name_(name) {
        const char* previous = std::getenv(name);
        if (previous) {
            had_previous_ = true;
            previous_ = previous;
        }
        if (value) {
            setenv(name, value, 1);
        } else {
            unsetenv(name);
        }
    }

    ~ScopedEnvironment() {
        if (had_previous_) {
            setenv(name_.c_str(), previous_.c_str(), 1);
        } else {
            unsetenv(name_.c_str());
        }
    }

   private:
    std::string name_;
    std::string previous_;
    bool had_previous_ = false;
};

}  // namespace

int main() {
    using lightning::compute::ComputeBudget;
    using lightning::compute::LoadComputeBudget;

    ScopedEnvironment clear_lio(lightning::compute::kLioThreadsEnv, nullptr);
    ScopedEnvironment clear_ndt(lightning::compute::kNdtThreadsEnv, nullptr);
    ScopedEnvironment clear_solid(lightning::compute::kSolidIcpWorkersEnv, nullptr);
    ScopedEnvironment clear_solid_nice(lightning::compute::kSolidWorkerNiceEnv,
                                       nullptr);
    ScopedEnvironment clear_solid_affinity(
        lightning::compute::kSolidCpuAffinityEnv, nullptr);

    ComputeBudget budget{12, 4, 8, 0, {}};
    std::string error;
    const YAML::Node root = YAML::Load(R"(
compute_budget:
  lio_threads: 6
  ndt_threads: 3
  solid_icp_workers: 2
  solid_worker_nice: 5
  solid_cpu_affinity: 1-3,5
)");
    Require(LoadComputeBudget(root, budget, &error), "valid YAML budget loads");
    Require(budget.lio_threads == 6 && budget.ndt_threads == 3 &&
                budget.solid_icp_workers == 2 && budget.solid_worker_nice == 5 &&
                budget.solid_cpu_affinity == "1-3,5",
            "all YAML worker pools share one budget");

    {
        ScopedEnvironment override_ndt(lightning::compute::kNdtThreadsEnv, "5");
        ComputeBudget overridden{12, 4, 8, 0, {}};
        Require(LoadComputeBudget(root, overridden, &error), "environment override loads");
        Require(overridden.lio_threads == 6 && overridden.ndt_threads == 5 &&
                    overridden.solid_icp_workers == 2 &&
                    overridden.solid_worker_nice == 5,
                "environment overrides only its worker pool");
    }

    {
        ScopedEnvironment override_nice(lightning::compute::kSolidWorkerNiceEnv,
                                        "7");
        ScopedEnvironment override_affinity(
            lightning::compute::kSolidCpuAffinityEnv, "2-4,6");
        ComputeBudget overridden{12, 4, 8, 0, {}};
        Require(LoadComputeBudget(root, overridden, &error),
                "SOLiD scheduling overrides load");
        Require(overridden.solid_worker_nice == 7 &&
                    overridden.solid_cpu_affinity == "2-4,6",
                "SOLiD scheduling environment overrides apply");
    }

    {
        ScopedEnvironment invalid_lio(lightning::compute::kLioThreadsEnv, "0");
        ComputeBudget unchanged{7, 4, 3, 0, {}};
        Require(!LoadComputeBudget(root, unchanged, &error), "zero threads rejected");
        Require(unchanged.lio_threads == 7 && unchanged.ndt_threads == 4 &&
                    unchanged.solid_icp_workers == 3,
                "failed load is transactional");
    }

    {
        ScopedEnvironment invalid_nice(lightning::compute::kSolidWorkerNiceEnv,
                                       "20");
        ComputeBudget unchanged{7, 4, 3, 0, {}};
        Require(!LoadComputeBudget(root, unchanged, &error),
                "out-of-range SOLiD nice is rejected");
        Require(unchanged.solid_worker_nice == 0,
                "invalid nice leaves budget unchanged");
    }

    {
        ScopedEnvironment invalid_affinity(
            lightning::compute::kSolidCpuAffinityEnv, "7-3");
        ComputeBudget unchanged{7, 4, 3, 0, {}};
        Require(!LoadComputeBudget(root, unchanged, &error),
                "invalid SOLiD affinity is rejected");
        Require(unchanged.solid_cpu_affinity.empty(),
                "invalid affinity leaves budget unchanged");
    }

    std::cout << "compute budget tests passed\n";
    return 0;
}
