//
// Created by xiang on 24-5-13.
//

#ifndef MIAO_OPTI_ALGO_LM_H
#define MIAO_OPTI_ALGO_LM_H

#include "optimization_algorithm.h"

namespace lightning::miao {

/**
 * \brief Implementation of the Levenberg Algorithm
 */
class OptimizationAlgorithmLevenberg : public OptimizationAlgorithm {
   public:
    /**
     * construct the Levenberg Algorithm, which will use the given Solver for
     * solving the linearized system.
     */
    explicit OptimizationAlgorithmLevenberg(std::shared_ptr<Solver> solver);

    virtual ~OptimizationAlgorithmLevenberg() = default;

    virtual SolverResult Solve(int iteration);

   protected:
    /**
     * helper for Levenberg, this function computes the initial damping factor, if
     * the user did not specify an own value, see setUserLambdaInit()
     */
    double ComputeLambdaInit() const;

    double ComputeScale() const;

    // Levenberg parameters
    int max_trails_after_failure_ = 10;  ///< 一次迭代失败后的试探次数上限。
    double user_lambda_limit_ = 0;  ///< 用户设置的阻尼参数限制。
    double current_lambda_ = -1;  ///< 当前 LM 阻尼系数。
    double tau_ = 1e-5;  ///< 自动选择初始阻尼时使用的比例系数。
    double good_step_lower_ = 1.0 / 3.0;  ///< lower bound for lambda decrease if a good  LM step
    double good_step_upper_ = 2.0 / 3.0;  ///< upper bound for lambda decrease if a good  LM step
    double ni_ = 2.0;  ///< 失败试探后增大阻尼的倍率状态。
    int lm_iter_ = 0;  ///< the number of levenberg iterations performed to accept the last step
};

}  // namespace lightning::miao

#endif  // MIAO_OPTI_ALGO_LM_H
