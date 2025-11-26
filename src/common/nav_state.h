//
// Created by xiang on 2022/2/15.
//

#pragma once

#include "common/eigen_types.h"
#include "common/s2.hpp"

#include <glog/logging.h>
#include <iomanip>

namespace lightning {
/**
 * ESKF导航状态结构体
 *
 * 这是激光雷达SLAM系统中的核心状态表示，采用误差状态卡尔曼滤波(ESKF)架构。
 * 结构体支持李群流形优化，确保旋转和球面约束的数值稳定性。
 *
 * 状态向量维度分配 (共23维):
 * ┌─────────────────────────────────────────────────────────────┐
 * │ 状态维度    │ 物理含义           │ 数学表示          │ 用途      │
 * ├─────────────────────────────────────────────────────────────┤
 * │ [0-2]      │ 3D位置           │ (x, y, z)         │ 导航定位  │
 * │ [3-5]      │ 姿态角轴         │ 李代数so(3)       │ 旋转表示  │
 * │ [6-8]      │ 外参旋转         │ R_lidar←imu      │ 传感器标定│
 * │ [9-11]     │ 外参平移         │ t_lidar←imu      │ 传感器标定│
 * │ [12-14]    │ 世界系速度       │ (vx, vy, vz)      │ 运动估计  │
 * │ [15-17]    │ 陀螺仪偏置       │ bg (rad/s)        │ 偏置补偿  │
 * │ [18-20]    │ 加速度计偏置     │ ba (m/s²)        │ 偏置补偿  │
 * │ [21-22]    │ 重力向量(S2)     │ S2流形, |g|=9.81 │ 重力估计  │
 * └─────────────────────────────────────────────────────────────┘
 *
 * 数学特性：
 * - 支持SO(3)李群：避免万向锁，保持数值稳定性
 * - 支持S2球面流形：重力向量单位长度约束
 * - 误差状态分离：名义状态 + 误差状态的ESKF架构
 * - 外参在线标定：支持LiDAR-IMU外参数估计
 *
 * 总计: 3×7 + 2 = 23维名义状态
 *       full_dim = 24维误差状态(S2矢量化为3维用于误差计算)
 */
struct NavState {
    constexpr static int dim = 23;       //  状态变量维度
    constexpr static int full_dim = 24;  // 误差状态变量维度

    using VectState = Eigen::Matrix<double, dim, 1>;           // 矢量形式
    using FullVectState = Eigen::Matrix<double, full_dim, 1>;  // 全状态矢量形式

    NavState() = default;

    bool operator<(const NavState& other) { return timestamp_ < other.timestamp_; }

    FullVectState ToState() {
        FullVectState ret;
        ret.block<3, 1>(0, 0) = pos_;
        ret.block<3, 1>(3, 0) = rot_.log();
        ret.block<3, 1>(6, 0) = offset_R_lidar_.log();
        ret.block<3, 1>(9, 0) = offset_t_lidar_;
        ret.block<3, 1>(12, 0) = vel_;
        ret.block<3, 1>(15, 0) = bg_;
        ret.block<3, 1>(18, 0) = ba_;
        ret.block<3, 1>(21, 0) = grav_.vec_;
        return ret;
    }

    void FromVectState(const FullVectState& state) {
        pos_ = state.block<3, 1>(0, 0);
        rot_ = SO3::exp(state.block<3, 1>(3, 0));
        offset_R_lidar_ = SO3::exp(state.block<3, 1>(6, 0));
        offset_t_lidar_ = state.block<3, 1>(9, 0);
        vel_ = state.block<3, 1>(12, 0);
        bg_ = state.block<3, 1>(15, 0);
        ba_ = state.block<3, 1>(18, 0);
        grav_.vec_ = state.block<3, 1>(21, 0);
    }

    // 计算状态向量的时间导数（IMU不考虑测量噪声的离散运动学方程的雅可比矩阵）
    inline FullVectState get_f(const Vec3d& gyro, const Vec3d& acce) const {
        FullVectState res = FullVectState::Zero();
        // 减零偏
        Vec3d omega = gyro - bg_;
        Vec3d a_inertial = rot_ * (acce - ba_);  // 加计读数-ba 并转到 世界系下

        for (int i = 0; i < 3; i++) {
            res(i) = vel_[i];                        // 位置导数: dx/dt = velocity
            res(i + 3) = omega[i];                   // 旋转角轴导数: dθ/dt = angular_velocity (李代数空间)
            res(i + 12) = a_inertial[i] + grav_[i];  // 速度导数: dv/dt = acceleration + gravity
        }
        return res;
    }

    /// 运动方程对误差状态的雅可比矩阵: ∂f/∂(δx)
    /// 返回部分雅可比矩阵F，包含所有与dt相关的部分
    inline Eigen::Matrix<double, full_dim, dim> df_dx(const Vec3d& acce) const {
        Eigen::Matrix<double, full_dim, dim> cov = Eigen::Matrix<double, full_dim, dim>::Zero();
        cov.block<3, 3>(0, 12) = Mat3d::Identity();  // ∂(δẋ)/∂(δv) = I: 位置误差对速度误差的导数（缺少dt）
        Vec3d acc = acce - ba_;                      // 去偏置的加速度（用于科里奥利力计算）
        // Vec3d omega = gyro - bg_;
        cov.block<3, 3>(12, 3) = -rot_.matrix() * SO3::hat(acc);  // ∂(δv̇)/∂(δθ): 姿态误差对速度误差的影响（缺少dt）
        cov.block<3, 3>(12, 18) = -rot_.matrix();  // ∂(δv̇)/∂(δba): 加速度计偏置误差对速度误差的影响（缺少dt）

        Vec2d vec = Vec2d::Zero();
        Eigen::Matrix<double, 3, 2> grav_matrix = grav_.S2_Mx(vec);  // S2重力向量误差映射矩阵

        cov.block<3, 2>(12, 21) = grav_matrix;  // ∂(δv̇)/∂(δg): 重力向量误差对速度误差的影响（缺少dt）
        cov.block<3, 3>(3, 15) =
            -Eigen::Matrix3d::Identity();  // ∂(δθ̇)/∂(δbg): 陀螺仪偏置误差对姿态误差的影响（缺少dt，值也不对）
        return cov;
    }

    /// 运动方程对过程噪声的雅可比矩阵 ∂f/∂w (24×12)
    /// 返回部分雅可比矩阵F，包含所有与dt相关的部分
    inline Eigen::Matrix<double, 24, 12> df_dw() const {
        Eigen::Matrix<double, 24, 12> cov = Eigen::Matrix<double, 24, 12>::Zero();

        /// 噪声向量 w 的维度分配 (12维):
        /// w[0:2]   - 陀螺仪噪声 ω_noise (影响姿态变化率)
        /// w[3:5]   - 加速度计噪声 a_noise (影响速度变化率)
        /// w[6:8]   - 陀螺仪零偏噪声 bg_noise (零偏随机游走)
        /// w[9:11]  - 加速度计零偏噪声 ba_noise (零偏随机游走)

        // 速度 (12:14) 受加速度噪声影响: ẋ = R*(a - ba) + g, ∂ẋ/∂a_noise = -R
        cov.block<3, 3>(12, 3) = -rot_.matrix();
        // 姿态 (3:5) 受陀螺噪声影响: Ṙ = R*⌊ω-bg⌋×, ∂Ṙ/∂ω_noise = -I
        cov.block<3, 3>(3, 0) = -Eigen::Matrix3d::Identity();  // 值不对，外部赋值
        // 陀螺零偏 (15:17) 受零偏噪声影响: ḃg = bg_noise, ∂ḃg/∂bg_noise = I
        cov.block<3, 3>(15, 6) = Eigen::Matrix3d::Identity();
        // 加速度计零偏 (18:20) 受零偏噪声影响: ḃa = ba_noise, ∂ḃa/∂ba_noise = I
        cov.block<3, 3>(18, 9) = Eigen::Matrix3d::Identity();
        return cov;
    }

    /// 状态积分：基于状态导数进行欧拉积分更新
    void oplus(const FullVectState& vec, double dt) {
        timestamp_ += dt;                                   // 时间戳前进
        pos_ += vec.middleRows(0, 3) * dt;                  // 位置积分：p = p + ẋ*dt （ẋ为速度）
        rot_ = rot_ * SO3::exp(vec.middleRows(3, 3) * dt);  // 姿态积分：R = R * exp(ω*dt) （ω为角速度，SO3流形）
        offset_R_lidar_ =
            offset_R_lidar_ * SO3::exp(vec.middleRows(6, 3) * dt);      // 外参旋转积分：R_ext = R_ext * exp(ω_ext*dt)
        offset_t_lidar_ = offset_t_lidar_ + vec.middleRows(9, 3) * dt;  // 外参平移积分：t_ext = t_ext + ṫ_ext*dt
        vel_ += vec.middleRows(12, 3) * dt;                             // 速度积分：v = v + ẍ*dt （ẍ为加速度）
        bg_ += vec.middleRows(15, 3) * dt;        // 陀螺零偏积分：bg = bg + ω_bias*dt （一阶马尔可夫）
        ba_ += vec.middleRows(18, 3) * dt;        // 加计零偏积分：ba = ba + a_bias*dt （一阶马尔可夫）
        grav_.oplus(vec.middleRows(21, 3) * dt);  // 重力向量更新：g = g ⊕ ġ*dt （S2流形上的指数更新）
    }

    /**
     * 广义减法, this - other
     * @param result 减法结果
     * @param other 另一个状态变量
     */
    VectState boxminus(const NavState& other) {
        VectState result;
        result.block<3, 1>(0, 0) = pos_ - other.pos_;
        result.block<3, 1>(3, 0) = (other.rot_.inverse() * rot_).log();
        result.block<3, 1>(6, 0) = (other.offset_R_lidar_.inverse() * offset_R_lidar_).log();
        result.block<3, 1>(9, 0) = offset_t_lidar_ - other.offset_t_lidar_;
        result.block<3, 1>(12, 0) = vel_ - other.vel_;
        result.block<3, 1>(15, 0) = bg_ - other.bg_;
        result.block<3, 1>(18, 0) = ba_ - other.ba_;

        Vec2d dg = grav_.boxminus(other.grav_);
        result.block<2, 1>(21, 0) = dg;
        return result;
    }

    /**
     * 广义加法 this = this+dx
     * @param dx 增量
     */
    NavState boxplus(const VectState& dx) {
        NavState ret;
        ret.timestamp_ = timestamp_;
        ret.pos_ = pos_ + dx.middleRows(0, 3);
        ret.rot_ = rot_ * SO3::exp(dx.middleRows(3, 3));
        ret.offset_R_lidar_ = offset_R_lidar_ * SO3::exp(dx.middleRows(6, 3));
        ret.offset_t_lidar_ = offset_t_lidar_ + dx.middleRows(9, 3);
        ret.vel_ = vel_ + dx.middleRows(12, 3);
        ret.bg_ = bg_ + dx.middleRows(15, 3);
        ret.ba_ = ba_ + dx.middleRows(18, 3);
        ret.grav_ = grav_;
        ret.grav_.boxplus(dx.middleRows(21, 2));

        return ret;
    }

    /// 各个子变量所在维度信息，用于ESKF中状态空间的维度映射
    struct MetaInfo {
        MetaInfo(int idx, int vdim, int dof) : idx_(idx), dim_(vdim), dof_(dof) {}
        int idx_ = 0;  // 目标索引：在23维误差状态向量中的起始位置（用于f_x_final, f_w_final等结果矩阵）
        int dim_ = 0;  // 源维度索引：在24维全状态向量中的起始位置（用于访问f_x_, f_w_等原始雅可比矩阵）
        int dof_ = 0;  // 自由度：从24维空间读取多少维数据（向量状态为3，SO3状态为3，S2状态为2）
    };

    static const std::vector<MetaInfo> vect_states_;  // 矢量变量的维度
    static const std::vector<MetaInfo> SO3_states_;   // SO3 变量的维度
    static const std::vector<MetaInfo> S2_states_;    // S2 变量维度

    friend inline std::ostream& operator<<(std::ostream& os, const NavState& s) {
        os << std::setprecision(18) << s.pos_.transpose() << " " << s.rot_.unit_quaternion().coeffs().transpose() << " "
           << s.offset_R_lidar_.unit_quaternion().coeffs().transpose() << " " << s.offset_t_lidar_.transpose() << " "
           << s.vel_.transpose() << " " << s.bg_.transpose() << " " << s.ba_.transpose() << " "
           << s.grav_.vec_.transpose();
        return os;
    }

    inline SE3 GetPose() const { return SE3(rot_, pos_); }
    inline SO3 GetRot() const { return rot_; }
    inline void SetPose(const SE3& pose) {
        rot_ = pose.so3();
        pos_ = pose.translation();
    }

    inline Vec3d Getba() const { return ba_; }
    inline Vec3d Getbg() const { return bg_; }
    inline Vec3d GetVel() const { return vel_; }
    void SetVel(const Vec3d& v) { vel_ = v; }

    double timestamp_ = 0.0;           // 时间戳
    double confidence_ = 0.0;          // 定位置信度
    bool pose_is_ok_ = true;           // 定位是否有效
    bool lidar_odom_reliable_ = true;  // lio是否有效
    bool is_parking_ = false;          // 是否在停车
    // 名义状态：
    Vec3d pos_ = Vec3d::Zero();             // 位置
    SO3 rot_;                               // 旋转
    SO3 offset_R_lidar_;                    // 外参R
    Vec3d offset_t_lidar_ = Vec3d::Zero();  // 外参t
    Vec3d vel_ = Vec3d::Zero();             // 速度
    Vec3d bg_ = Vec3d::Zero();              // 陀螺零偏
    Vec3d ba_ = Vec3d::Zero();              // 加计零偏
    S2 grav_;                               // 重力
};

}  // namespace lightning
