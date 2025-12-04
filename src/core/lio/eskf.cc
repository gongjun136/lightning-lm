//
// Created by xiang on 2022/2/15.
//

#include "core/lio/eskf.hpp"

namespace lightning {

void ESKF::Predict(const double& dt, const ESKF::ProcessNoiseType& Q, const Vec3d& gyro, const Vec3d& acce) {
    Eigen::Matrix<double, 24, 1> f_ = x_.get_f(gyro, acce);  // 状态雅可比: f(x,u) - 刚体运动方程
    Eigen::Matrix<double, 24, 23> f_x_ = x_.df_dx(acce);     // 部分误差状态雅可比: ∂f/∂x - 用于线性化

    Eigen::Matrix<double, 24, 12> f_w_ = x_.df_dw();          // 部分噪声雅可比: ∂f/∂w - 噪声传播矩阵
    Eigen::Matrix<double, 23, process_noise_dim_> f_w_final;  // 23*12，噪声雅可比，需要*dt的项

    NavState x_before = x_;  // 保存前一时刻状态（用于S2流形计算）
    x_.oplus(f_, dt);        // 名义状态积分: x_{k+1} = x_k ⊕ f(x_k,u_k)*dt

    F_x1_ = CovType::Identity();  // 存储部分误差状态雅可比，不需要*dt部分

    // 构建完整的状态雅可比和噪声雅可比
    CovType f_x_final;  // 23x23，误差状态雅可比，需要*dt的项
    for (auto st : x_.vect_states_) {
        int idx = st.idx_;
        int dim = st.dim_;
        int dof = st.dof_;

        for (int i = 0; i < 23; i++) {
            for (int j = 0; j < dof; j++) {
                f_x_final(idx + j, i) = f_x_(dim + j, i);  // 填充向量状态的雅可比
            }
        }

        for (int i = 0; i < process_noise_dim_; i++) {
            for (int j = 0; j < dof; j++) {
                f_w_final(idx + j, i) = f_w_(dim + j, i);  // 填充向量状态的噪声传播
            }
        }
    }
    // 1.处理SO3旋转状态
    Mat3d res_temp_SO3;
    Vec3d seg_SO3;
    for (auto st : x_.SO3_states_) {
        int idx = st.idx_;
        int dim = st.dim_;
        for (int i = 0; i < 3; i++) {
            seg_SO3(i) = -1 * f_(dim + i) * dt;  // 角轴向量 = -ω*dt
        }

        F_x1_.block<3, 3>(idx, idx) = math::exp(seg_SO3, 0.5).matrix();  // SO3状态转移矩阵

        res_temp_SO3 = math::A_matrix(seg_SO3);  // -v代入= 李代数右雅可比
        for (int i = 0; i < state_dim_; i++) {
            // [gj-2025-11-26] 为啥没有添加负号？？？
            f_x_final.template block<3, 1>(idx, i) = res_temp_SO3 * (f_x_.block<3, 1>(dim, i));
        }

        for (int i = 0; i < process_noise_dim_; i++) {
            // [gj-2025-11-26] 为啥没有添加负号？？？
            f_w_final.template block<3, 1>(idx, i) = res_temp_SO3 * (f_w_.block<3, 1>(dim, i));
        }
    }
    // 2.处理S2球面状态（重力向量）      这是什么玩意？？
    Eigen::Matrix<double, 2, 3> res_temp_S2;
    Vec3d seg_S2;
    for (auto st : x_.S2_states_) {
        int idx = st.idx_;
        int dim = st.dim_;
        for (int i = 0; i < 3; i++) {
            seg_S2(i) = f_(dim + i) * dt;  // S2状态的李代数增量
        }

        SO3 res = math::exp(seg_S2, 0.5f);  // 指数映射到SO3

        Vec2d vec = Vec2d::Zero();
        Eigen::Matrix<double, 2, 3> Nx = x_.grav_.S2_Nx_yy();  // S2流形投影矩阵
        Eigen::Matrix<double, 3, 2> Mx = x_before.grav_.S2_Mx(vec);

        F_x1_.block<2, 2>(idx, idx) = Nx * res.matrix() * Mx;  // S2状态转移矩阵

        Eigen::Matrix<double, 3, 3> x_before_hat = x_before.grav_.S2_hat();  // 反对称矩阵
        res_temp_S2 = -Nx * res.matrix() * x_before_hat * math::A_matrix(seg_S2).transpose();

        for (int i = 0; i < state_dim_; i++) {
            f_x_final.block<2, 1>(idx, i) = res_temp_S2 * (f_x_.block<3, 1>(dim, i));
        }
        for (int i = 0; i < process_noise_dim_; i++) {
            f_w_final.block<2, 1>(idx, i) = res_temp_S2 * (f_w_.block<3, 1>(dim, i));
        }
    }

    F_x1_ += f_x_final * dt;  // 完整状态转移矩阵: F = I + J*dt
    P_ = (F_x1_)*P_ * (F_x1_).transpose() + (dt * f_w_final) * Q * (dt * f_w_final).transpose();  // 卡尔曼协方差传播
}

/**
 * 原版的迭代过程中，收敛次数大于1才会结果，所以需要两次收敛。
 * 在未收敛时，实际上不会计算最近邻，也就回避了一次ObsModel的计算
 * 如果这边对每次迭代都计算最近邻的话，时间明显会变长一些，并不是非常合理。。
 *
 * @param obs 观测类型，包括LIDAR、WHEEL_SPEED、GPS等，决定使用哪种观测模型
 * @param R   观测噪声方差，控制观测对状态估计的权重（越小权重越高，典型值1e-3）
 */
void ESKF::Update(ESKF::ObsType obs, const double& R) {
    custom_obs_model_.valid_ = true;
    custom_obs_model_.converge_ = true;

    CovType P_propagated = P_;  // 保存预测阶段的协方差矩阵

    Eigen::Matrix<double, 23, 1> K_r;   // 卡尔曼增益与残差的乘积项 K*r
    Eigen::Matrix<double, 23, 23> K_H;  // 卡尔曼增益与雅可比的乘积项 K*H

    StateVecType dx_current = StateVecType::Zero();  // 本轮迭代的dx

    NavState start_x = x_;  // 迭代的起点
    NavState last_x = x_;   // last_x 不是上一次迭代的状态，而是AA算法优化后的状态

    int converged_times = 0;    // 连续收敛次数统计
    double last_lidar_res = 0;  // 上一轮迭代的激光雷达残差

    double init_res = 0.0;           // 初始残差（用于相对收敛判断）
    static double iterated_num = 0;  // 总迭代次数统计
    static double update_num = 0;    // 总更新次数统计
    update_num += 1;
    for (int i = -1; i < maximum_iter_; i++) {
        custom_obs_model_.valid_ = true;

        /// 计算observation function，主要是residual_, h_x_, s_
        /// x_ 在每次迭代中都是更新的，线性化点也会更新
        if (obs == ObsType::LIDAR || obs == ObsType::WHEEL_SPEED_AND_LIDAR) {
            lidar_obs_func_(x_, custom_obs_model_);
        } else if (obs == ObsType::WHEEL_SPEED) {
            wheelspeed_obs_func_(x_, custom_obs_model_);
        } else if (obs == ObsType::ACC_AS_GRAVITY) {
            acc_as_gravity_obs_func_(x_, custom_obs_model_);
        } else if (obs == ObsType::GPS) {
            gps_obs_func_(x_, custom_obs_model_);
        } else if (obs == ObsType::BIAS) {
            bias_obs_func_(x_, custom_obs_model_);
        }
        // Anderson加速收敛判断：当残差增大时回退到上一步状态
        // 如果使用Anderson加速且不是第一次迭代(i > -1)，并且当前是Lidar相关观测，
        // 且当前残差比上一次残差大1%以上，说明迭代可能发散，需要回退
        if (use_aa_ && i > -1 && (obs == ObsType::LIDAR || obs == ObsType::WHEEL_SPEED_AND_LIDAR) &&
            custom_obs_model_.lidar_residual_mean_ >= last_lidar_res * 1.01) {
            x_ = last_x;  // 回退到上一步的状态
            break;        // 跳出当前迭代循环
        }
        iterated_num += 1;

        if (!custom_obs_model_.valid_) {
            continue;
        }

        if (i == -1) {
            // 第一次迭代：记录初始残差作为收敛判断的基准
            init_res = custom_obs_model_.lidar_residual_mean_;
            if (init_res < 1e-9) {
                init_res = 1e-9;  // 防止除零，设置最小残差阈值
            }
        }

        iterations_ = i + 2;                                             // 迭代次数统计，i从-1开始，所以+2
        final_res_ = custom_obs_model_.lidar_residual_mean_ / init_res;  // 相对残差

        int dof_measurement = custom_obs_model_.h_x_.rows();  // 观测维度
        StateVecType dx = x_.boxminus(start_x);               // 计算当前状态相对于起始状态的变化量
        dx_current = dx;

        P_ = P_propagated;  // 重置为预测阶段的协方差

        // TODO.custom_obs_model_.h_x_中应该还有姿态部分名义状态关于误差状态的导数
        /// 更新P 和 dx
        /// P = J*P*J^T
        /// dx = J * dx
        // 处理SO3流形状态的协方差变换
        for (auto it : x_.SO3_states_) {
            int idx = it.idx_;                       // SO3状态在向量中的索引
            Vec3d seg_SO3 = dx.block<3, 1>(idx, 0);  // 提取SO3状态的变化量
            // [gj-2025-12-4] 右雅可比的转置！= 逆；虽然接近单位阵可以近似。
            Mat3d res_temp_SO3 = Mat3d::Identity() - 0.5 * SO3::hat(seg_SO3);
            // Mat3d res_temp_SO3 = math::A_matrix(seg_SO3).transpose();  // 先验误差映射到迭代点的雅可比矩阵
            // 对dx进行流形变换,这里理论是有负号的，后面dx_current迭代更新抵消了，所以就没写了
            dx_current.block<3, 1>(idx, 0) = res_temp_SO3 * dx.block<3, 1>(idx, 0);

            /// 更新协方差矩阵P的行：P_row = J * P_row
            for (int j = 0; j < state_dim_; j++) {
                P_.block<3, 1>(idx, j) = res_temp_SO3 * (P_.block<3, 1>(idx, j));
            }
            /// 更新协方差矩阵P的列：P_col = P_col * J^T
            for (int j = 0; j < state_dim_; j++) {
                P_.block<1, 3>(j, idx) = (P_.block<1, 3>(j, idx)) * res_temp_SO3.transpose();
            }
        }

        // 处理S2球面流形状态的协方差变换（重力向量）TODO，这里不懂
        for (auto it : x_.S2_states_) {
            int idx = it.idx_;  // S2状态在向量中的索引

            Vec2d seg_S2 = dx.block<2, 1>(idx, 0);  // 提取S2状态的变化量

            // 计算S2流形的局部坐标化雅可比矩阵
            Eigen::Matrix<double, 2, 3> Nx = x_.grav_.S2_Nx_yy();          // 当前重力状态的切空间投影矩阵
            Eigen::Matrix<double, 3, 2> Mx = start_x.grav_.S2_Mx(seg_S2);  // 起始状态的局部坐标映射
            Mat2d res_temp_S2 = Nx * Mx;                                   // 复合雅可比矩阵

            dx_current.block<2, 1>(idx, 0) = res_temp_S2 * dx.block<2, 1>(idx, 0);  // 对dx进行流形变换

            /// 更新协方差矩阵P的行：P_row = J * P_row
            for (int j = 0; j < state_dim_; j++) {
                P_.block<2, 1>(idx, j) = res_temp_S2 * (P_.block<2, 1>(idx, j));
            }

            /// 更新协方差矩阵P的列：P_col = P_col * J^T
            for (int j = 0; j < state_dim_; j++) {
                P_.block<1, 2>(j, idx) = (P_.block<1, 2>(j, idx)) * res_temp_S2.transpose();
            }
        }

        /// 处理各类观测模型：稀疏观测情况（观测维度 < 状态维度）
        if (state_dim_ > dof_measurement) {
            // 构建完整的雅可比矩阵，将激光雷达观测扩展到完整状态空间
            Eigen::MatrixXd h_x_cur = Eigen::MatrixXd::Zero(dof_measurement, state_dim_);
            h_x_cur.topLeftCorner(dof_measurement, 12) = custom_obs_model_.h_x_;                     // 只填充前12列
            custom_obs_model_.R_ = R * Eigen::MatrixXd::Identity(dof_measurement, dof_measurement);  // 观测噪声矩阵

            // 计算卡尔曼增益：K = P * H^T * (H * P * H^T + R)^(-1)
            Eigen::MatrixXd K =
                P_ * h_x_cur.transpose() * (h_x_cur * P_ * h_x_cur.transpose() + custom_obs_model_.R_).inverse();
            K_r = K * custom_obs_model_.residual_;  // 卡尔曼增益 × 残差
            K_H = K * h_x_cur;                      // 卡尔曼增益 × 雅可比矩阵
        } else {
            /// 纯雷达观测：观测维度 >= 状态维度 时，使用信息形式避免在大维度上求 (H P H^T + R)^-1
            /// 注意：这里假设 R = σ² I，传进来的标量 R 就是 σ²

            // H 为 12×state_dim_，只作用于前 12 维状态（姿态+位置+速度+陀螺偏置）
            // 这里先算 H^T H，对应于扩展到 23 维后的 H_full^T H_full 的左上 12×12 块
            Eigen::Matrix<double, 12, 12> HTH = custom_obs_model_.h_x_.transpose() * custom_obs_model_.h_x_;

            // 构造“先验信息矩阵”的一个缩放版本：
            // P_temp = (P / R)^(-1) = (P / σ²)^(-1) = σ² P^{-1} = R * P^{-1}
            // 也就是说，把先验协方差 P 转成先验信息矩阵（乘了一个 R 的尺度因子）
            CovType P_temp = (P_ / R).inverse();

            // 把观测信息叠加到信息矩阵上：
            // 左上 12×12 加上 H^T H，相当于
            // P_temp = R * (P^{-1} + H_full^T R^{-1} H_full)
            // 其中 H_full = [H  0] 为 12×23 的完整观测雅可比
            P_temp.block<12, 12>(0, 0) += HTH;

            // 取逆得到 Q_inv：
            // Q_inv = P_temp^{-1}
            //       = (1 / R) * (P^{-1} + H_full^T R^{-1} H_full)^{-1}
            // 也就是“真实后验协方差”的 1/R 缩放版本，这个尺度因子在下面计算 K_r、K_H 时会被抵消掉
            CovType Q_inv = P_temp.inverse();

            // 计算状态增量中的 K r 项：
            // K_r = Q_inv * H^T * r
            // 对应数学上 K r，其中
            //   K = P_post * H_full^T * R^{-1}
            //   P_post = (P^{-1} + H_full^T R^{-1} H_full)^{-1}
            // 由于 Q_inv = (1/R) * P_post，K 中的 R^{-1} 与 Q_inv 里的 1/R 正好抵消
            K_r = Q_inv.template block<23, 12>(0, 0) * custom_obs_model_.h_x_.transpose() * custom_obs_model_.residual_;

            // 计算 K H_full，用于后续 dx = K_r + (K_H - I) * dx_current：
            // K_H = Q_inv * H^T * H_full
            // 其中 H_full^T H_full 的左上 12×12 就是 HTH，其余列为 0
            K_H.setZero();
            K_H.template block<23, 12>(0, 0) = Q_inv.template block<23, 12>(0, 0) * HTH;
        }

        // 误差状态更新（dx_current前面省略了负号）
        dx_current = K_r + (K_H - Eigen::Matrix<double, 23, 23>::Identity()) * dx_current;

        // check nan
        for (int j = 0; j < 23; ++j) {
            if (std::isnan(dx_current(j, 0))) {
                return;
            }
        }

        if (!use_aa_) {
            x_ = x_.boxplus(dx_current);
        } else {
            // 转到起点的线性空间
            x_ = x_.boxplus(dx_current);

            if (i == -1) {
                aa_.init(dx_current);  // 初始化AA
            } else {
                // Anderson加速：基于历史迭代信息优化状态更新
                auto dx_all = x_.boxminus(start_x);     // 计算从起点到当前状态的总变化量
                auto new_dx_all = aa_.compute(dx_all);  // AA算法基于历史信息计算最优方向
                x_ = start_x.boxplus(new_dx_all);       // 应用AA优化后的状态更新
            }
        }

        last_x = x_;

        // update last res
        last_lidar_res = custom_obs_model_.lidar_residual_mean_;
        custom_obs_model_.converge_ = true;
        // 收敛性检查
        for (int j = 0; j < 23; j++) {
            if (std::fabs(dx_current[j]) > limit_[j]) {
                custom_obs_model_.converge_ = false;
                break;
            }
        }

        if (custom_obs_model_.converge_) {
            converged_times++;
        }

        // 兜底机制：如果从未收敛且到达倒数第二次迭代，强制标记为收敛
        // 防止优化无法结束，确保系统稳定性
        if (!converged_times && i == maximum_iter_ - 2) {
            custom_obs_model_.converge_ = true;
        }

        if (converged_times > 0 || i == maximum_iter_ - 1) {
            /// 结束条件：已经至少收敛过一次，或者已经到达最大迭代次数（兜底退出）
            /// 此时 P_、K_H 仍然是在“起始线性化点 start_x”的切空间下表达，
            /// 下面要把协方差 / 信息矩阵通过 dx_current 映射到最终状态 x_ 所在的流形切空间，再做最终的P更新（式(45)）

            L_ = P_;  ///< 临时拷贝一份当前协方差，用于做流形上的相似变换
            Mat3d res_temp_SO3;
            Vec3d seg_SO3;
            // 1) 处理所有 SO3 类型的状态块
            for (auto it : x_.SO3_states_) {
                int idx = it.idx_;  // 该SO3块在误差状态向量中的起始索引

                // seg_SO3 = 最终一次迭代得到的该 SO3 块的误差李代数（δθ）
                for (int j = 0; j < 3; j++) {
                    seg_SO3(j) = dx_current(j + idx);
                }

                // [gj-2025-12-4] 右雅可比的转置！= 逆；虽然接近单位阵可以近似。
                // res_temp_SO3 = A(δθ)^T，右雅可比的转置，对应“从 start_x 切空间 → 当前 x_ 切空间”的雅可比
                // res_temp_SO3 = math::A_matrix(seg_SO3).transpose();
                res_temp_SO3 = Mat3d::Identity() - 0.5 * SO3::hat(seg_SO3);

                // 先更新协方差的“行”：P_row = J * P_row
                for (int j = 0; j < 23; j++) {
                    L_.block<3, 1>(idx, j) = res_temp_SO3 * (P_.block<3, 1>(idx, j));
                }

                // 同样方式更新 K_H 中与该 SO3 状态对应的行（K_H_row = J * K_H_row）
                // 注意：这里只对参与观测的前 15 维做处理
                for (int j = 0; j < 15; j++) {
                    K_H.block<3, 1>(idx, j) = res_temp_SO3 * (K_H.block<3, 1>(idx, j));
                }

                // 再更新协方差的“列”：P_col = P_col * J^T
                // L_=  J * P_ * J^T
                // （这里既更新 L_ 也顺带把 P_ 本身做相同变换，以便后面使用）
                for (int j = 0; j < 23; j++) {
                    L_.block<1, 3>(j, idx) = (L_.block<1, 3>(j, idx)) * res_temp_SO3.transpose();
                    P_.block<1, 3>(j, idx) = (P_.block<1, 3>(j, idx)) * res_temp_SO3.transpose();
                }
            }

            Mat2d res_temp_S2;
            Vec2d seg_S2;
            // 2) 处理所有 S2（重力方向）类型的状态块
            for (auto it : x_.S2_states_) {
                int idx = it.idx_;  // S2 块在误差状态向量中的起始索引

                // seg_S2 = 最终一次迭代得到的该 S2 块的局部 2 维误差坐标
                for (int j = 0; j < 2; j++) {
                    seg_S2(j) = dx_current(j + idx);
                }

                // Nx：当前重力方向在切空间的投影矩阵
                // Mx：从起始重力方向（start_x）局部坐标 seg_S2 反映射到3维的矩阵
                // res_temp_S2 = J_S2：同样是“起点切空间 → 当前切空间”的雅可比
                Eigen::Matrix<double, 2, 3> Nx = x_.grav_.S2_Nx_yy();
                Eigen::Matrix<double, 3, 2> Mx = start_x.grav_.S2_Mx(seg_S2);
                res_temp_S2 = Nx * Mx;

                // 行变换：P_row = J * P_row
                for (auto j = 0; j < 23; j++) {
                    L_.block<2, 1>(idx, j) = res_temp_S2 * (P_.block<2, 1>(idx, j));
                }

                // 行变换：K_H_row = J * K_H_row（同上，只对前 15 维有意义）
                for (auto j = 0; j < 15; j++) {
                    K_H.block<2, 1>(idx, j) = res_temp_S2 * (K_H.block<2, 1>(idx, j));
                }

                // 列变换：P_col = P_col * J^T
                for (int j = 0; j < 23; j++) {
                    L_.block<1, 2>(j, idx) = (L_.block<1, 2>(j, idx)) * res_temp_S2.transpose();
                    P_.block<1, 2>(j, idx) = (P_.block<1, 2>(j, idx)) * res_temp_S2.transpose();
                }
            }

            // 3) 最终协方差更新
            // 这里的公式对应文中式 的一种实现形式：
            //   P_new = L - (K_H_full) * P_propagated_obs
            // 其中：
            //   - L 为经过流形相似变换后的协方差
            //   - K_H.block<23,15>(0,0) 相当于 K * H_full 的前 15 列（参与观测的状态块）
            //   - P_.block<15,23>(0,0) 为参与观测的先验协方差行块
            P_ = L_ - K_H.block<23, 15>(0, 0) * P_.template block<15, 23>(0, 0);

            break;  // 结束整个Update()的迭代
        }
    }
}

}  // namespace lightning