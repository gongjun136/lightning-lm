#pragma once

#include "common/eigen_types.h"

#include <glog/logging.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <vector>

namespace lightning {

struct NavState {
    // Error-state layout: position, rotation, velocity, gyro bias, accel bias, gravity.
    constexpr static int dim = 18;
    constexpr static int full_dim = 18;
    constexpr static int kBlockDim = 3;
    constexpr static int kPosIdx = 0;
    constexpr static int kRotIdx = 3;
    constexpr static int kVelIdx = 6;
    constexpr static int kBgIdx = 9;
    constexpr static int kBaIdx = 12;
    constexpr static int kGravIdx = 15;
    constexpr static double kGravityNorm = 9.81;
    constexpr static double kMaxGravityDelta = 0.02;

    using VectState = Eigen::Matrix<double, dim, 1>;
    using FullVectState = Eigen::Matrix<double, full_dim, 1>;

    NavState() = default;

    bool operator<(const NavState& other) const { return timestamp_ < other.timestamp_; }

    static Vec3d NormalizeGravity(const Vec3d& grav) {
        const double norm = grav.norm();
        if (norm < 1e-9 || !std::isfinite(norm)) {
            return Vec3d(0.0, 0.0, -kGravityNorm);
        }
        return grav * (kGravityNorm / norm);
    }

    static Vec3d ApplyGravityDelta(const Vec3d& grav, const Vec3d& delta) {
        const Vec3d base = NormalizeGravity(grav);
        const Vec3d radial = base / kGravityNorm;
        Vec3d tangent_delta = delta - radial * radial.dot(delta);
        const double delta_norm = tangent_delta.norm();
        if (delta_norm > kMaxGravityDelta) {
            tangent_delta *= kMaxGravityDelta / delta_norm;
        }
        return NormalizeGravity(base + tangent_delta);
    }

    FullVectState ToState() const {
        FullVectState ret = FullVectState::Zero();
        ret.block<kBlockDim, 1>(kPosIdx, 0) = pos_;
        ret.block<kBlockDim, 1>(kRotIdx, 0) = rot_.log();
        ret.block<kBlockDim, 1>(kVelIdx, 0) = vel_;
        ret.block<kBlockDim, 1>(kBgIdx, 0) = bg_;
        ret.block<kBlockDim, 1>(kBaIdx, 0) = ba_;
        ret.block<kBlockDim, 1>(kGravIdx, 0) = grav_;
        return ret;
    }

    void FromVectState(const FullVectState& state) {
        pos_ = state.block<kBlockDim, 1>(kPosIdx, 0);
        rot_ = SO3::exp(state.block<kBlockDim, 1>(kRotIdx, 0));
        vel_ = state.block<kBlockDim, 1>(kVelIdx, 0);
        bg_ = state.block<kBlockDim, 1>(kBgIdx, 0);
        ba_ = state.block<kBlockDim, 1>(kBaIdx, 0);
        grav_ = NormalizeGravity(state.block<kBlockDim, 1>(kGravIdx, 0));
    }

    inline FullVectState get_f(const Vec3d& gyro, const Vec3d& acce) const {
        FullVectState res = FullVectState::Zero();
        const Vec3d omega = gyro - bg_;
        const Vec3d acc_unbiased = acce - ba_;
        const Vec3d a_inertial = rot_ * acc_unbiased;

        res.segment<kBlockDim>(kPosIdx) = vel_;
        res.segment<kBlockDim>(kRotIdx) = omega;
        res.segment<kBlockDim>(kVelIdx) = a_inertial + grav_;
        return res;
    }

    inline Eigen::Matrix<double, full_dim, dim> df_dx(const Vec3d& acce) const {
        Eigen::Matrix<double, full_dim, dim> cov = Eigen::Matrix<double, full_dim, dim>::Zero();
        const Vec3d acc_unbiased = acce - ba_;

        cov.block<kBlockDim, kBlockDim>(kPosIdx, kVelIdx) = Mat3d::Identity();
        cov.block<kBlockDim, kBlockDim>(kRotIdx, kBgIdx) = -Mat3d::Identity();
        cov.block<kBlockDim, kBlockDim>(kVelIdx, kRotIdx) = -rot_.matrix() * SO3::hat(acc_unbiased);
        cov.block<kBlockDim, kBlockDim>(kVelIdx, kBaIdx) = -rot_.matrix();
        cov.block<kBlockDim, kBlockDim>(kVelIdx, kGravIdx) = Mat3d::Identity();
        return cov;
    }

    inline Eigen::Matrix<double, full_dim, 12> df_dw() const {
        Eigen::Matrix<double, full_dim, 12> cov = Eigen::Matrix<double, full_dim, 12>::Zero();
        cov.block<kBlockDim, kBlockDim>(kRotIdx, 0) = -Mat3d::Identity();
        cov.block<kBlockDim, kBlockDim>(kVelIdx, 3) = -rot_.matrix();
        cov.block<kBlockDim, kBlockDim>(kBgIdx, 6) = Mat3d::Identity();
        cov.block<kBlockDim, kBlockDim>(kBaIdx, 9) = Mat3d::Identity();
        return cov;
    }

    void oplus(const FullVectState& vec, double dt) {
        timestamp_ += dt;
        pos_ += vec.middleRows(kPosIdx, kBlockDim) * dt;
        rot_ = rot_ * SO3::exp(vec.middleRows(kRotIdx, kBlockDim) * dt);
        vel_ += vec.middleRows(kVelIdx, kBlockDim) * dt;
        bg_ += vec.middleRows(kBgIdx, kBlockDim) * dt;
        ba_ += vec.middleRows(kBaIdx, kBlockDim) * dt;
        grav_ = ApplyGravityDelta(grav_, vec.middleRows(kGravIdx, kBlockDim) * dt);
    }

    VectState boxminus(const NavState& other) const {
        VectState result = VectState::Zero();
        result.block<kBlockDim, 1>(kPosIdx, 0) = pos_ - other.pos_;
        result.block<kBlockDim, 1>(kRotIdx, 0) = (other.rot_.inverse() * rot_).log();
        result.block<kBlockDim, 1>(kVelIdx, 0) = vel_ - other.vel_;
        result.block<kBlockDim, 1>(kBgIdx, 0) = bg_ - other.bg_;
        result.block<kBlockDim, 1>(kBaIdx, 0) = ba_ - other.ba_;
        const Vec3d other_grav = NormalizeGravity(other.grav_);
        const Vec3d radial = other_grav / kGravityNorm;
        const Vec3d grav_delta = NormalizeGravity(grav_) - other_grav;
        result.block<kBlockDim, 1>(kGravIdx, 0) = grav_delta - radial * radial.dot(grav_delta);
        return result;
    }

    NavState boxplus(const VectState& dx) const {
        NavState ret = *this;
        ret.pos_ = pos_ + dx.middleRows(kPosIdx, kBlockDim);
        ret.rot_ = rot_ * SO3::exp(dx.middleRows(kRotIdx, kBlockDim));
        ret.vel_ = vel_ + dx.middleRows(kVelIdx, kBlockDim);
        ret.bg_ = bg_ + dx.middleRows(kBgIdx, kBlockDim);
        ret.ba_ = ba_ + dx.middleRows(kBaIdx, kBlockDim);
        ret.grav_ = ApplyGravityDelta(grav_, dx.middleRows(kGravIdx, kBlockDim));
        return ret;
    }

    struct MetaInfo {
        MetaInfo(int idx, int vdim, int dof) : idx_(idx), dim_(vdim), dof_(dof) {}
        int idx_ = 0;
        int dim_ = 0;
        int dof_ = 0;
    };

    static const std::vector<MetaInfo> vect_states_;
    static const std::vector<MetaInfo> SO3_states_;

    friend inline std::ostream& operator<<(std::ostream& os, const NavState& s) {
        os << std::setprecision(18) << s.pos_.transpose() << " " << s.rot_.unit_quaternion().coeffs().transpose()
           << " " << s.vel_.transpose() << " " << s.bg_.transpose() << " " << s.ba_.transpose() << " "
           << s.grav_.transpose();
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

    double timestamp_ = 0.0;
    double confidence_ = 0.0;
    bool pose_is_ok_ = true;
    bool lidar_odom_reliable_ = true;
    bool is_parking_ = false;

    Vec3d pos_ = Vec3d::Zero();
    SO3 rot_;
    Vec3d vel_ = Vec3d::Zero();
    Vec3d bg_ = Vec3d::Zero();
    Vec3d ba_ = Vec3d::Zero();
    Vec3d grav_ = Vec3d(0.0, 0.0, -kGravityNorm);
};

}  // namespace lightning
