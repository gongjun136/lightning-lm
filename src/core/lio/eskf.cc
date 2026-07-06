//
// Created by xiang on 2022/2/15.
//

#include "core/lio/eskf.hpp"
#include "core/lightning_math.hpp"

#include <Eigen/Eigenvalues>
#include <algorithm>

namespace {

using CovType = lightning::ESKF::CovType;

void SymmetrizeAndFloorCovariance(CovType& P, double min_cov_diag) {
    P = 0.5 * (P + P.transpose()).eval();

    for (int i = 0; i < P.rows(); ++i) {
        if (P(i, i) < min_cov_diag) {
            P(i, i) = min_cov_diag;
        } else if (P(i, i) > 100.0) {
            P(i, i) = 100.0;
        }

        for (int j = 0; j < P.cols(); ++j) {
            if (std::isnan(P(i, j)) || std::isinf(P(i, j))) {
                LOG(WARNING) << "find nan or inf in P: " << P(i, j);
                P(i, j) = 1.0;
            }
        }
    }
}

}  // namespace

namespace lightning {

/**
 * @brief IMU预测：用当前IMU输入传播名义状态和误差状态协方差。
 *
 * gyro/acce一般由上层ImuProcess传入，可能是相邻两帧IMU量测的均值。这里虽然可以使用
 * 均值输入 u_avg，但状态导数和误差雅可比仍在当前状态 x_k 处计算，因此整体形式是
 * x_{k+1} = x_k ⊕ f(x_k, u_avg) * dt 的一阶欧拉传播，而不是严格中值积分。
 *
 * @param dt IMU积分时间间隔，单位秒。
 * @param Q 过程噪声协方差，噪声顺序需要和NavState::df_dw()保持一致。
 * @param gyro 当前积分区间使用的角速度输入。
 * @param acce 当前积分区间使用的加速度计输入。
 */
void ESKF::Predict(const double& dt, const ESKF::ProcessNoiseType& Q, const Vec3d& gyro, const Vec3d& acce) {
    // f_是名义状态连续时间导数 f(x,u)，例如位置导数、姿态角速度、速度导数。
    // f_x_是连续动力学对误差状态的雅可比 F_c，后面会通过 I + F_c * dt 做一阶离散化。
    Eigen::Matrix<double, NavState::full_dim, 1> f_ = x_.get_f(gyro, acce);
    if (!options_.propagate_velocity_) {
        f_.template segment<NavState::kBlockDim>(NavState::kVelIdx).setZero();
    }
    Eigen::Matrix<double, NavState::full_dim, state_dim_> f_x_ = x_.df_dx(acce);

    // f_w_是连续动力学对过程噪声的雅可比 G_c。
    // NavState为了兼容SO3/S2等流形，先在full_dim空间里组织导数；P_实际是state_dim_维。
    Eigen::Matrix<double, NavState::full_dim, process_noise_dim_> f_w_ = x_.df_dw();
    Eigen::Matrix<double, state_dim_, process_noise_dim_> f_w_final =
        Eigen::Matrix<double, state_dim_, process_noise_dim_>::Zero();

    NavState x_before = x_;  // 保存前一时刻状态（用于S2流形计算）
    x_.oplus(f_, dt);        // 名义状态积分: x_{k+1} = x_k ⊕ f(x_k,u)*dt

    // F_x1_最终表示离散误差状态转移矩阵 Phi，先放入单位阵，再叠加各状态块的一阶传播项。
    F_x1_ = CovType::Identity();

    // 普通向量状态块可以直接从full_dim导数空间拷贝到state_dim误差空间。
    // 这些块使用线性加法，不需要SO3那样额外的李群切空间映射。
    CovType f_x_final = CovType::Zero();
    for (auto st : x_.vect_states_) {
        int idx = st.idx_;
        int dim = st.dim_;
        int dof = st.dof_;

        for (int i = 0; i < state_dim_; i++) {
            for (int j = 0; j < dof; j++) {
                f_x_final(idx + j, i) = f_x_(dim + j, i);
            }
        }

        for (int i = 0; i < process_noise_dim_; i++) {
            for (int j = 0; j < dof; j++) {
                f_w_final(idx + j, i) = f_w_(dim + j, i);
            }
        }
    }

    // SO3状态块的误差定义在李代数切空间，不能像普通向量一样直接拷贝。
    // 这里需要显式构造旋转误差自身传播项，以及连续雅可比到切空间的映射。
    Mat3d res_temp_SO3;
    Vec3d seg_SO3;
    for (auto st : x_.SO3_states_) {
        int idx = st.idx_;
        int dim = st.dim_;
        for (int i = 0; i < 3; i++) {
            // 右扰动误差下，旋转误差自身传播项含有 -omega*dt。
            seg_SO3(i) = -1 * f_(dim + i) * dt;
        }

        // Phi_RR：姿态误差对上一时刻姿态误差的离散传播。
        F_x1_.block<3, 3>(idx, idx) = math::exp(seg_SO3, 0.5).matrix();

        // A_matrix(seg_SO3)用于把连续雅可比项映射到SO3误差切空间。
        res_temp_SO3 = math::A_matrix(seg_SO3);
        for (int i = 0; i < state_dim_; i++) {
            // [gj-2025-11-26] 为啥没有添加负号？？？
            // 答：这里不再额外加负号。seg_SO3 = -Omega，因此A_matrix(seg_SO3)=J_l(-Omega)=J_r(Omega)；
            // 姿态对gyro bias的负号已经在NavState::df_dx()的kRotIdx-kBgIdx块中，即f_x_里是-I。
            f_x_final.template block<3, 1>(idx, i) = res_temp_SO3 * (f_x_.block<3, 1>(dim, i));
        }

        for (int i = 0; i < process_noise_dim_; i++) {
            // [gj-2025-11-26] 为啥没有添加负号？？？
            // 答：同理，姿态对gyro noise的负号已经在NavState::df_dw()的kRotIdx噪声块中，即f_w_里是-I。
            // 此处只做SO3切空间映射；若再加负号，会把文档推导中的-J_r(Omega)dt变成错误的+J_r(Omega)dt。
            f_w_final.template block<3, 1>(idx, i) = res_temp_SO3 * (f_w_.block<3, 1>(dim, i));
        }
    }

    // 一阶离散化误差传播：
    //   Phi ~= I + F_c * dt
    //   P_{k+1} = Phi P_k Phi^T + (G_c dt) Q (G_c dt)^T
    // 这里没有构造中点状态处的F_c，因此属于起点状态线性化的欧拉式协方差传播。
    F_x1_ += f_x_final * dt;
    P_ = (F_x1_)*P_ * (F_x1_).transpose() + (dt * f_w_final) * Q * (dt * f_w_final).transpose();
    // 轻微膨胀协方差，给线性化误差、时间同步误差和未建模误差留余量。
    P_ *= options_.predict_cov_inflation_;
    // 数值保护：强制协方差对称，并限制对角线范围，避免后续更新阶段数值不稳定。
    SymmetrizeAndFloorCovariance(P_, options_.min_cov_diag_);
}

/**
 * @brief 根据指定观测模型迭代修正ESKF状态和协方差。
 *
 * Update()使用的是迭代误差状态更新。每轮迭代都会：
 * 1. 在当前名义状态x_处重新计算观测模型；
 * 2. 将观测函数返回的 H^T H 和 H^T r 与预测协方差组合；
 * 3. 求解当前误差增量dx_current，并通过boxplus()更新名义状态；
 * 4. 判断收敛或到达最大迭代次数后，更新协方差P_。
 *
 * 原版的迭代过程中，收敛次数大于1才会结果，所以需要两次收敛。
 * 在未收敛时，实际上不会计算最近邻，也就回避了一次ObsModel的计算。
 * 如果这边对每次迭代都计算最近邻的话，时间明显会变长一些，并不是非常合理。
 *
 * @param obs 观测类型，包括LIDAR、WHEEL_SPEED、GPS等，决定使用哪种观测模型
 * @param R   观测噪声方差，控制观测对状态估计的权重（越小权重越高，典型值1e-3）
 */
void ESKF::Update(ESKF::ObsType obs, const double& R) {
    // 每次Update开始前先把观测模型标志恢复为可用状态，具体质量检查交给观测函数设置。
    custom_obs_model_.valid_ = true;
    custom_obs_model_.converge_ = true;

    // P_propagated是IMU预测后的先验协方差。迭代过程中会多次临时变换P_，
    // 每一轮重新从这个先验协方差开始，避免把上一轮临时线性化结果重复叠加。
    CovType P_propagated = P_;  // 保存预测阶段的协方差矩阵

    Eigen::Matrix<double, state_dim_, 1> K_r;           // 卡尔曼增益与残差的乘积项 K*r
    Eigen::Matrix<double, state_dim_, state_dim_> K_H;  // 卡尔曼增益与雅可比的乘积项 K*H

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

        /// 计算observation function。这里不是直接返回H和r，而是返回已经累加好的 H^T H 和 H^T r。
        /// x_ 在每次迭代中都会被boxplus()更新，所以观测模型的线性化点也随之更新。
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

        if (custom_obs_model_.valid_ == false) {
            x_ = last_x;
            P_ = P_propagated;
            return;
        }

        if (use_aa_ && i > -1 && (obs == ObsType::LIDAR || obs == ObsType::WHEEL_SPEED_AND_LIDAR) &&
            custom_obs_model_.lidar_residual_mean_ >= last_lidar_res * 1.01) {
            x_ = last_x;  // 回退到上一步的状态
            break;        // 跳出当前迭代循环，停止继续迭代
        }
        iterated_num += 1;

        // 理论上上面valid=false已经return了，这里保留continue兼容早期逻辑。
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

        StateVecType dx = x_.boxminus(start_x);  // 当前x与起点之间的dx
        dx_current = dx;                         //

        P_ = P_propagated;  // 重置为预测阶段的协方差

        // TODO.custom_obs_model_.h_x_中应该还有姿态部分名义状态关于误差状态的导数
        /// 更新P 和 dx
        /// P = J*P*J^T
        /// dx = J * dx
        // 处理SO3流形状态的协方差变换。
        // 因为迭代起点start_x和当前线性化点x_不在同一个切空间，先验误差dx和协方差P_
        // 需要通过流形雅可比映射到当前x_的切空间，再和本轮观测模型组合。
        for (auto it : x_.SO3_states_) {
            int idx = it.idx_;                       // SO3状态在向量中的索引
            Vec3d seg_SO3 = dx.block<3, 1>(idx, 0);  // 提取SO3状态的变化量
            // 使用从start_x切空间到当前x_切空间的右雅可比J_r(δθ)。
            // math::A_matrix(δθ)=J_l(δθ)，由BCH.md中J_l(φ)^T=J_l(-φ)=J_r(φ)，所以转置后正好是J_r(δθ)。
            // I - 0.5*[δθ]x只是J_r(δθ)的小角度一阶近似；这里保留完整右雅可比。
            Mat3d res_temp_SO3 = math::A_matrix(seg_SO3).transpose();
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

        // 观测函数返回的是6维位姿约束的信息形式：
        // HTH = H^T H，HTr = H^T r。先做对称化，避免并行累加或数值误差导致特征分解不稳定。
        Mat6d HTH = custom_obs_model_.HTH_;
        Vec6d HTr = custom_obs_model_.HTr_;
        Mat6d HTH_sym = 0.5 * (HTH + HTH.transpose());

        // 对观测信息矩阵做特征分解，用特征值大小判断哪些方向可观、哪些方向退化。
        // 例如长直走廊中，某些平移/旋转方向可能缺乏几何约束，直接更新会过度相信不可靠残差。
        Eigen::SelfAdjointEigenSolver<Mat6d> eigen_solver(HTH_sym);
        if (eigen_solver.info() != Eigen::Success) {
            LOG(WARNING) << "Failed to decompose ESKF observation information matrix.";
            continue;
        }

        const Vec6d eigen_values = eigen_solver.eigenvalues();
        const Mat6d eigen_vectors = eigen_solver.eigenvectors();
        const double max_eigen_value = std::max(1e-12, eigen_values.maxCoeff());
        const double degeneracy_threshold = max_eigen_value * options_.degeneracy_threshold_ratio_;

        // LOG(INFO) << "eigen values of HTH: " << eigen_values.transpose();

        Vec6d observable_mask = Vec6d::Zero();
        int nullity = 0;
        // 特征值相对最大特征值太小的方向被认为是退化方向，后面会被投影掉。
        for (int k = 0; k < observable_mask.size(); ++k) {
            if (eigen_values(k) > degeneracy_threshold) {
                observable_mask(k) = 1.0;
            } else {
                nullity++;
            }
        }

        // 投影矩阵只保留可观方向。HTH_eff/HTr_eff是退化处理后的有效观测信息。
        const Mat6d observable_projector = eigen_vectors * observable_mask.asDiagonal() * eigen_vectors.transpose();
        const Mat6d HTH_eff = observable_projector * HTH_sym * observable_projector;
        const Vec6d HTr_eff = observable_projector * HTr;

        // 信息形式更新。P_ / R 等价于把观测噪声缩放合并进先验权重；
        // 再取逆得到先验信息矩阵 P^{-1} * R。
        CovType P_temp = (P_ / R).inverse();  // P阵上面已经更新

        /// 现在问题是这个权重太大，导致整体过于依赖先验 ...
        // P_temp.setIdentity();

        // 当前观测只约束前6维位姿，因此只把HTH_eff加到信息矩阵左上角位姿块。
        P_temp.block<pose_obs_dim_, pose_obs_dim_>(0, 0) += HTH_eff;
        CovType Q_inv = P_temp.inverse();  // Q inv

        // Q*H^T * R^-1 * r = K * r
        // <-- K ----->
        K_r = Q_inv.template block<state_dim_, pose_obs_dim_>(0, 0) * HTr_eff;

        // K_H = Q^-1 H^T R^-1 H
        //       <--  K     ->
        K_H.setZero();
        K_H.template block<state_dim_, pose_obs_dim_>(0, 0) =
            Q_inv.template block<state_dim_, pose_obs_dim_>(0, 0) * HTH_eff;

        // dx = Kr + (KH-I) dx
        // LOG(INFO) << "K_r: " << K_r.transpose()
        //           << ", prior: " << ((K_H - Eigen::Matrix<double, state_dim_, state_dim_>::Identity()) *
        //           dx_current).transpose();

        // 迭代误差状态更新公式。
        // K_r来自当前观测残差，(K_H-I)dx_current用于把先验误差项带入当前线性化点。
        dx_current = K_r + (K_H - Eigen::Matrix<double, state_dim_, state_dim_>::Identity()) * dx_current;

        // check nan
        for (int j = 0; j < state_dim_; ++j) {
            if (std::isnan(dx_current(j, 0))) {
                return;
            }
        }

        // Vec3d dv = dx_current.middleRows(NavState::kVelIdx, NavState::kBlockDim);
        // if (dv.norm() > options_.vel_clip_norm_) {
        //     dv = dv / dv.norm() * options_.vel_clip_norm_;
        // }

        // dv = dv * options_.dv_ratio_;
        // dx_current.middleRows(NavState::kVelIdx, NavState::kBlockDim) = dv;

        // dx_current.middleRows(18, 5).setZero();

        // LOG(INFO) << "iter " << iterations_ << ", dx: " << dx_current.transpose();
        const double dx_translation = dx_current.head<3>().norm();
        const double dx_rotation_deg = dx_current.segment<3>(3).norm() * 180.0 / M_PI;
        // 单次迭代修正过大通常意味着匹配错误、时间同步异常或初值偏差太大。
        // 这里直接拒绝本次观测更新，回退到预测状态，避免把错误观测注入滤波器。
        if (dx_translation > options_.max_update_translation_step_ ||
            dx_rotation_deg > options_.max_update_rotation_step_deg_) {
            LOG(ERROR) << "Reject ESKF iter update, dtrans: " << dx_translation << ", drot_deg: " << dx_rotation_deg
                       << ", dvel: " << dx_current.segment<NavState::kBlockDim>(NavState::kVelIdx).norm();
            x_ = start_x;
            P_ = P_propagated;
            return;
        }

        if (!use_aa_) {
            // 普通迭代：直接把误差状态增量施加到当前名义状态。
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

        // 记录本轮残差，下一轮Anderson加速时若残差明显变大，会回退到last_x。
        last_lidar_res = custom_obs_model_.lidar_residual_mean_;
        custom_obs_model_.converge_ = true;
        // 收敛性检查
        for (int j = 0; j < state_dim_; j++) {
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

                // 这里应使用从start_x切空间到当前x_切空间的右雅可比J_r(δθ)。
                // math::A_matrix(δθ)=J_l(δθ)，由BCH.md中J_l(φ)^T=J_l(-φ)=J_r(φ)，所以转置后正好是J_r(δθ)。
                // I - 0.5*[δθ]x只是J_r(δθ)的小角度一阶近似；这里保留完整右雅可比，避免较大迭代步下误差变大。
                res_temp_SO3 = math::A_matrix(seg_SO3).transpose();

                // 先更新协方差的“行”：P_row = J * P_row
                for (int j = 0; j < state_dim_; j++) {
                    L_.block<3, 1>(idx, j) = res_temp_SO3 * (P_.block<3, 1>(idx, j));
                }

                // 同样方式更新 K_H 中与该 SO3 状态对应的行（K_H_row = J * K_H_row）
                // 注意：这里只对参与观测的前 15 维做处理
                for (int j = 0; j < pose_obs_dim_; j++) {
                    K_H.block<3, 1>(idx, j) = res_temp_SO3 * (K_H.block<3, 1>(idx, j));
                }

                // 再更新协方差的“列”：P_col = P_col * J^T
                // L_=  J * P_ * J^T
                // （这里既更新 L_ 也顺带把 P_ 本身做相同变换，以便后面使用）
                for (int j = 0; j < state_dim_; j++) {
                    L_.block<1, 3>(j, idx) = (L_.block<1, 3>(j, idx)) * res_temp_SO3.transpose();
                    P_.block<1, 3>(j, idx) = (P_.block<1, 3>(j, idx)) * res_temp_SO3.transpose();
                }
            }

            // 最终协方差更新，对应信息形式下的 P = (I-KH)P，并结合上面的流形切空间映射。
            P_ = L_ - K_H.block<state_dim_, pose_obs_dim_>(0, 0) * P_.template block<pose_obs_dim_, state_dim_>(0, 0);

            if (nullity > 0) {
                // LOG_EVERY_N(INFO, 50) << "ESKF observation degeneracy rank " << (pose_obs_dim_ - nullity) << "/"
                //                      << pose_obs_dim_;
                // 如果观测存在退化方向，适当膨胀位姿协方差，避免滤波器对这些方向过度自信。
                P_.block<pose_obs_dim_, pose_obs_dim_>(0, 0) *= options_.degeneracy_cov_inflation_;
            }

            break;
        }
    }

    // 最后统一做协方差数值保护，确保P_对称并且对角线处于合理范围。
    SymmetrizeAndFloorCovariance(P_, options_.min_cov_diag_);
}

}  // namespace lightning
