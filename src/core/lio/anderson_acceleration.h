#ifndef ANDERSONACCELERATION_H_
#define ANDERSONACCELERATION_H_

#include <omp.h>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <vector>

#include "common/eigen_types.h"

namespace lightning {
/**
 * AA 加速器
 * @tparam S Scalar type
 * @tparam D dimension
 * @tparam m 允许前多少次迭代, 最多10次，不允许动态设置，只能编译期设置
 */
template <typename S, int D, int m>
class AndersonAcceleration {
   public:
    using Scalar = S;
    using Vec = Eigen::Matrix<S, D, 1>;    // 长度为 D 的列向量
    using MatDM = Eigen::Matrix<S, D, m>;  // D 行、m 列的矩阵，用来存历史量
    using MatDD = Eigen::Matrix<S, D, D>;  // D×D 矩阵（法方程矩阵等）

    /**
     * 为 g 计算 Anderson 加速之后的结果
     * @param g 当前迭代给出的“原始更新量”/新解 G_k
     * @return  加速之后的更新量（新的 u_{k+1}）
     */
    Vec compute(const Vec& g) {
        assert(iter_ >= 0);
        Vec G = g;                    // 为了语义清晰，记作 G_k
        current_F_ = G - current_u_;  // 当前残差 F_k = G_k - u_k

        if (iter_ == 0) {
            // 第一次迭代：只初始化历史差分，不做线性组合
            prev_dF_.col(0) = -current_F_;  // dF_0 = F_0 - F_-1，这里等价于 -F_0（假定 F_-1 = 0）
            prev_dG_.col(0) = -G;           // dG_0 = G_0 - G_-1，这里等价于 -G_0（假定 G_-1 = 0）
            current_u_ = G;                 // 直接采用 G 作为新的 u_1
        } else {
            // 非第一次迭代，更新最新一列的 dF、dG
            prev_dF_.col(col_idx_) += current_F_;  // dF_j = F_k - F_{k-1}
            prev_dG_.col(col_idx_) += G;           // dG_j = G_k - G_{k-1}

            // 对该列 dF 做归一化，提升数值稳定性
            Scalar eps = 1e-14;
            Scalar scale = std::max(eps, prev_dF_.col(col_idx_).norm());
            dF_scale_(col_idx_) = scale;      // 记录缩放因子
            prev_dF_.col(col_idx_) /= scale;  // 归一化后的 dF_j

            // 实际使用的历史长度 m_k = min(设定的 m, 已迭代次数)
            int m_k = std::min(m, iter_);

            if (m_k == 1) {
                // 只用一列历史时，可以写成简单的一维情形
                theta_(0) = 0;
                Scalar dF_sqrnorm = prev_dF_.col(col_idx_).squaredNorm();
                M_(0, 0) = dF_sqrnorm;  // M = dFᵀ dF（标量）
                Scalar dF_norm = std::sqrt(dF_sqrnorm);

                if (dF_norm > eps) {
                    // theta = (dF / ||dF||)ᵀ (F / ||dF||)
                    theta_(0) = (prev_dF_.col(col_idx_) / dF_norm).dot(current_F_ / dF_norm);
                }
            } else {
                // m_k > 1：更新法方程矩阵 M = dFᵀ dF 的新行/列

                // 计算新列 dF 与已有 m_k 列 dF 的内积：dF_newᵀ * dF_i
                VecXd new_inner_prod = (prev_dF_.col(col_idx_).transpose() * prev_dF_.block(0, 0, D, m_k)).transpose();

                // 写入 M 的第 col_idx_ 行和第 col_idx_ 列（对称）
                M_.block(col_idx_, 0, 1, m_k) = new_inner_prod.transpose();
                M_.block(0, col_idx_, m_k, 1) = new_inner_prod;

                // 解法方程 M theta = dFᵀ F，得到 theta
                cod_.compute(M_.block(0, 0, m_k, m_k));
                theta_.head(m_k) = cod_.solve(prev_dF_.block(0, 0, D, m_k).transpose() * current_F_);
            }

            // 用 rescaled theta 计算新的 u：
            // u_{k+1} = G_k - sum_j (theta_j / scale_j) * dG_j
            current_u_ =
                G - prev_dG_.block(0, 0, D, m_k) * ((theta_.head(m_k).array() / dF_scale_.head(m_k).array()).matrix());

            // 循环选择下一次要写入的历史列索引
            col_idx_ = (col_idx_ + 1) % D;
            // 为下一轮迭代准备初始 dF、dG（等同于 -F_k, -G_k，后面再加上 F_{k+1}, G_{k+1}）
            prev_dF_.col(col_idx_) = -current_F_;
            prev_dG_.col(col_idx_) = -G;
        }

        iter_++;            // 迭代计数 +1
        return current_u_;  // 返回加速后的更新量
    }

    /**
     * 重置 AA 加速器，仅重置迭代状态，不清空历史矩阵内容
     * @param u 重置后的当前解 u
     */
    void reset(const Vec& u) {
        iter_ = 0;       // 从第 0 次迭代重新开始
        col_idx_ = 0;    // 历史列索引归零
        current_u_ = u;  // 当前解赋值
    }

    /**
     * 初始化 Anderson Acceleration（完全清零）
     * @param u0: 初始变量值 u_0
     */
    void init(const Vec& u0) {
        // 当前解与残差清零
        current_u_.setZero();
        current_F_.setZero();
        // 历史差分清零
        prev_dG_.setZero();
        prev_dF_.setZero();

        // 法方程矩阵和相关向量清零
        M_.setZero();
        theta_.setZero();
        dF_scale_.setZero();

        // 设置初始解
        current_u_ = u0;

        // 重置迭代计数与列索引
        iter_ = 0;
        col_idx_ = 0;
    }

   private:
    Vec current_u_ = Vec::Zero();  // 当前解 u_k（或当前加速后的更新）
    Vec current_F_ = Vec::Zero();  // 当前残差 F_k = G_k - u_k

    MatDM prev_dG_ = MatDM::Zero();  // 历史 dG_j 列向量：G_j - G_{j-1}
    MatDM prev_dF_ = MatDM::Zero();  // 历史 dF_j 列向量：F_j - F_{j-1}

    MatDD M_ = MatDD::Zero();     // 法方程矩阵 M = dFᵀ dF（只用前 m_k×m_k 块）
    Vec theta_ = Vec::Zero();     // 通过法方程求得的系数 theta（前 m_k 元素有效）
    Vec dF_scale_ = Vec::Zero();  // 每一列 dF 的缩放因子，用于反缩放 theta

    // 使用动态大小 MatXd 来做分解，因为 m_k 在迭代中会变化
    Eigen::CompleteOrthogonalDecomposition<MatXd> cod_;

    int iter_ = 0;      // 自初始化以来的迭代次数
    int col_idx_ = -1;  // 下一次写入历史矩阵的列索引（循环使用）
};

}  // namespace lightning

#endif /* ANDERSONACCELERATION_H_ */
