//
// Created by xiang on 2022/2/15.
//

#pragma once

#include "common/eigen_types.h"

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
    constexpr static int kBlockDim = 3;
    constexpr static int kPosIdx = 0;
    constexpr static int kRotIdx = 3;
    constexpr static int kVelIdx = 6;
    constexpr static int kBgIdx = 9;

    using VectState = Eigen::Matrix<double, dim, 1>;           // 矢量形式
    using FullVectState = Eigen::Matrix<double, full_dim, 1>;  // 全状态矢量形式

    NavState() = default;

    bool operator<(const NavState& other) { return timestamp_ < other.timestamp_; }

    FullVectState ToState() {
        FullVectState ret;
        ret.block<kBlockDim, 1>(kPosIdx, 0) = pos_;
        ret.block<kBlockDim, 1>(kRotIdx, 0) = rot_.log();
        ret.block<kBlockDim, 1>(kVelIdx, 0) = vel_;
        ret.block<kBlockDim, 1>(kBgIdx, 0) = bg_;
        return ret;
    }

    void FromVectState(const FullVectState& state) {
        pos_ = state.block<kBlockDim, 1>(kPosIdx, 0);
        rot_ = SO3::exp(state.block<kBlockDim, 1>(kRotIdx, 0));
        vel_ = state.block<kBlockDim, 1>(kVelIdx, 0);
        bg_ = state.block<kBlockDim, 1>(kBgIdx, 0);
    }

    /**
     * @brief 计算名义状态连续时间导数 f(x,u)。
     *
     * gyro/acce通常来自当前IMU积分区间的输入，外部可能已经做均值、滤波和尺度修正。
     * 这里不显式引入测量噪声，只用陀螺零偏bg_修正角速度；加速度计零偏没有参与当前在线状态。
     *
     * 返回的FullVectState按full_dim布局保存各状态块导数：
     * - pos_dot = vel_
     * - rot_dot = gyro - bg_，后续在oplus()里通过SO3指数映射积分
     * - vel_dot = rot_ * acce + grav_
     */
    inline FullVectState get_f(const Vec3d& gyro, const Vec3d& acce) const {
        FullVectState res = FullVectState::Zero();
        // 陀螺仪量测先减去估计零偏，得到用于姿态积分的角速度。
        Vec3d omega = gyro - bg_;
        // 加速度计量测从IMU系转到世界系；当前实现没有在线估计加计零偏ba。
        Vec3d a_inertial = rot_ * acce;

        for (int i = 0; i < 3; i++) {
            res(i) = vel_[i];                         // p_dot
            res(i + kRotIdx) = omega[i];              // theta_dot
            res(i + kVelIdx) = a_inertial[i] + grav_[i];  // v_dot
        }
        return res;
    }

    /**
     * @brief 运动方程对误差状态的连续雅可比 F_c = ∂f/∂δx。
     *
     * 这里返回的是未乘dt的连续时间雅可比。ESKF::Predict()会将其映射到误差状态空间，
     * 再用 Phi ~= I + F_c * dt 做一阶离散化。
     *
     * 当前实现包含的主要项：
     * - δp_dot / δv = I
     * - δv_dot / δθ = -R * [a]_x
     * - δθ_dot / δbg = -I
     *
     * @param acce 当前积分区间使用的加速度计输入。
     */
    inline Eigen::Matrix<double, full_dim, dim> df_dx(const Vec3d& acce) const {
        Eigen::Matrix<double, full_dim, dim> cov = Eigen::Matrix<double, full_dim, dim>::Zero();
        // 位置误差导数由速度误差直接驱动：δp_dot = δv。
        cov.block<kBlockDim, kBlockDim>(kPosIdx, kVelIdx) = Mat3d::Identity();
        Vec3d acc = acce;
        // Vec3d omega = gyro - bg_;
        // 速度误差对姿态误差的敏感度：
        // R Exp(δθ) a ≈ R a - R [a]_x δθ，因此 δv_dot / δθ = -R [a]_x。
        cov.block<kBlockDim, kBlockDim>(kVelIdx, kRotIdx) = -rot_.matrix() * SO3::hat(acc);
        // 姿态误差由陀螺零偏误差驱动：δθ_dot ≈ -δbg。
        cov.block<kBlockDim, kBlockDim>(kRotIdx, kBgIdx) = -Eigen::Matrix3d::Identity();
        return cov;
    }

    /**
     * @brief 运动方程对过程噪声的连续雅可比 G_c = ∂f/∂w。
     *
     * 噪声向量按12维预留，当前代码实际使用的顺序为：
     * - w[0:3]：陀螺仪白噪声，进入姿态误差；
     * - w[3:6]：加速度计白噪声，进入速度误差；
     * - w[6:9]：陀螺零偏随机游走，进入bg；
     * - w[9:12]：加速度计零偏随机游走预留，当前实现未接入。
     */
    inline Eigen::Matrix<double, full_dim, 12> df_dw() const {
        Eigen::Matrix<double, full_dim, 12> cov = Eigen::Matrix<double, full_dim, 12>::Zero();
        // 加速度噪声在IMU系，转到世界系后影响速度误差，符号来自 a = a_meas - eta_a。
        cov.block<kBlockDim, kBlockDim>(kVelIdx, 3) = -rot_.matrix();
        // 陀螺噪声直接影响姿态误差，符号来自 omega = omega_meas - bg - eta_g。
        cov.block<kBlockDim, kBlockDim>(kRotIdx, 0) = -Eigen::Matrix3d::Identity();
        // 陀螺零偏按随机游走建模，bg_dot = eta_bg。
        cov.block<kBlockDim, kBlockDim>(kBgIdx, 6) = Eigen::Matrix3d::Identity();
        return cov;
    }

    /// 状态积分：基于状态导数进行欧拉积分更新
    void oplus(const FullVectState& vec, double dt) {
        timestamp_ += dt;
        pos_ += vec.middleRows(kPosIdx, kBlockDim) * dt;
        rot_ = rot_ * SO3::exp(vec.middleRows(kRotIdx, kBlockDim) * dt);
        // vel_ += vec.middleRows(kVelIdx, kBlockDim) * dt;
        bg_ += vec.middleRows(kBgIdx, kBlockDim) * dt;
    }

    /**
     * 广义减法, this - other
     * @param result 减法结果
     * @param other 另一个状态变量
     */
    VectState boxminus(const NavState& other) {
        VectState result;
        result.block<kBlockDim, 1>(kPosIdx, 0) = pos_ - other.pos_;
        result.block<kBlockDim, 1>(kRotIdx, 0) = (other.rot_.inverse() * rot_).log();
        result.block<kBlockDim, 1>(kVelIdx, 0) = vel_ - other.vel_;
        result.block<kBlockDim, 1>(kBgIdx, 0) = bg_ - other.bg_;

        return result;
    }

    /**
     * 广义加法 this = this+dx
     * @param dx 增量
     */
    NavState boxplus(const VectState& dx) {
        NavState ret;
        ret.timestamp_ = timestamp_;
        ret.pos_ = pos_ + dx.middleRows(kPosIdx, kBlockDim);
        ret.rot_ = rot_ * SO3::exp(dx.middleRows(kRotIdx, kBlockDim));
        ret.vel_ = vel_ + dx.middleRows(kVelIdx, kBlockDim);
        ret.bg_ = bg_ + dx.middleRows(kBgIdx, kBlockDim);
        ret.grav_ = grav_;

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

    friend inline std::ostream& operator<<(std::ostream& os, const NavState& s) {
        os << std::setprecision(18) << s.pos_.transpose() << " " << s.rot_.unit_quaternion().coeffs().transpose() << " "
           << s.vel_.transpose() << " " << s.bg_.transpose() << " " << s.grav_.transpose();
        return os;
    }

    inline SE3 GetPose() const { return SE3(rot_, pos_); }
    inline SO3 GetRot() const { return rot_; }
    inline void SetPose(const SE3& pose) {
        rot_ = pose.so3();
        pos_ = pose.translation();
    }

    inline Vec3d Getba() const { return Vec3d::Zero(); }
    inline Vec3d Getbg() const { return bg_; }
    inline Vec3d GetVel() const { return vel_; }
    void SetVel(const Vec3d& v) { vel_ = v; }

    double timestamp_ = 0.0;           // 时间戳
    double confidence_ = 0.0;          // 定位置信度
    bool pose_is_ok_ = true;           // 定位是否有效
    bool lidar_odom_reliable_ = true;  // lio是否有效
    bool is_parking_ = false;          // 是否在停车

    Vec3d pos_ = Vec3d::Zero();            // 位置（从IMU坐标系到世界坐标系，在世界坐标系下表示）
    SO3 rot_;                              // 旋转（从IMU坐标系到世界坐标系）
    Vec3d vel_ = Vec3d::Zero();            // 速度（世界坐标系）
    Vec3d bg_ = Vec3d::Zero();             // 陀螺零偏（在IMU原始测量中补偿，IMU坐标系）
    Vec3d grav_ = Vec3d(0.0, 0.0, -9.81);  // 重力向量（世界坐标系）
};

}  // namespace lightning
