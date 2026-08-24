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

    ComputeBudget budget{12, 4, 8};
    std::string error;
    const YAML::Node root = YAML::Load(R"(
compute_budget:
  lio_threads: 6
  ndt_threads: 3
  solid_icp_workers: 2
)");
    Require(LoadComputeBudget(root, budget, &error), "valid YAML budget loads");
    Require(budget.lio_threads == 6 && budget.ndt_threads == 3 &&
                budget.solid_icp_workers == 2,
            "all YAML worker pools share one budget");

    {
        ScopedEnvironment override_ndt(lightning::compute::kNdtThreadsEnv, "5");
        ComputeBudget overridden{12, 4, 8};
        Require(LoadComputeBudget(root, overridden, &error), "environment override loads");
        Require(overridden.lio_threads == 6 && overridden.ndt_threads == 5 &&
                    overridden.solid_icp_workers == 2,
                "environment overrides only its worker pool");
    }

    {
        ScopedEnvironment invalid_lio(lightning::compute::kLioThreadsEnv, "0");
        ComputeBudget unchanged{7, 4, 3};
        Require(!LoadComputeBudget(root, unchanged, &error), "zero threads rejected");
        Require(unchanged.lio_threads == 7 && unchanged.ndt_threads == 4 &&
                    unchanged.solid_icp_workers == 3,
                "failed load is transactional");
    }

    std::cout << "compute budget tests passed\n";
    return 0;
}
