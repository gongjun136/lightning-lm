//
// Created by xiang on 24-4-26.
//
#include "opti_algo_gauss_newton.h"
#include "core/graph/optimizer.h"
#include "core/solver/solver.h"

#include <glog/logging.h>

#include <utility>

namespace lightning::miao {

OptimizationAlgorithmGaussNewton::OptimizationAlgorithmGaussNewton(std::shared_ptr<Solver> solver)
    : OptimizationAlgorithm(std::move(solver)) {}

OptimizationAlgorithm::SolverResult OptimizationAlgorithmGaussNewton::Solve(int iteration) {
    bool ok = true;

    // here so that correct component for max-mixtures can be computed before the
    // build structure
    optimizer_->ComputeActiveErrors();

    if (iteration == 0) {  // built up the CCS structure, here due to easy time measure
        ok = solver_->BuildStructure();
        if (!ok) {
            LOG(WARNING) << "Failure while building CCS structure";
            return OptimizationAlgorithm::SolverResult::Fail;
        }
    }

    if (!solver_->BuildSystem()) {
        return SolverResult::Fail;
    }
    ok = solver_->Solve();
    if (!ok) {
        return SolverResult::Fail;
    }
    optimizer_->Update(solver_->GetX());

    return SolverResult::OK;
}

}  // namespace lightning::miao
