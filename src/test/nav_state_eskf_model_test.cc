#include "common/nav_state.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

using namespace lightning;

namespace {

template <typename DerivedA, typename DerivedB>
void ExpectNear(const Eigen::MatrixBase<DerivedA>& actual, const Eigen::MatrixBase<DerivedB>& expected, double tol,
                const std::string& name) {
    const double err = (actual - expected).norm();
    if (err > tol) {
        std::cerr << name << " mismatch, error=" << err << "\nactual:\n"
                  << actual << "\nexpected:\n"
                  << expected << std::endl;
        std::exit(1);
    }
}

}  // namespace

int main() {
    static_assert(NavState::dim == 18, "NavState error-state dimension should be p/R/v/bg/ba/g");
    static_assert(NavState::full_dim == 18, "NavState full derivative dimension should match the active state");
    static_assert(NavState::kBaIdx == 12, "Accelerometer bias block should follow gyro bias");
    static_assert(NavState::kGravIdx == 15, "Gravity block should follow accelerometer bias");

    NavState s;
    s.rot_ = SO3::exp(Vec3d(0.2, -0.1, 0.05));
    s.vel_ = Vec3d(0.4, -0.2, 0.1);
    s.bg_ = Vec3d(0.01, -0.02, 0.03);
    s.ba_ = Vec3d(0.12, -0.08, 0.04);
    s.grav_ = Vec3d(0.0, 0.0, -NavState::kGravityNorm);

    const Vec3d gyro(0.11, 0.03, -0.02);
    const Vec3d acc(0.5, -0.4, 9.7);
    const Vec3d acc_unbiased = acc - s.ba_;
    const Mat3d R = s.rot_.matrix();

    const auto f = s.get_f(gyro, acc);
    ExpectNear(f.segment<3>(NavState::kPosIdx), s.vel_, 1e-12, "pos derivative");
    ExpectNear(f.segment<3>(NavState::kRotIdx), gyro - s.bg_, 1e-12, "rot derivative");
    ExpectNear(f.segment<3>(NavState::kVelIdx), s.rot_ * acc_unbiased + s.grav_, 1e-12, "vel derivative");
    ExpectNear(f.segment<3>(NavState::kBgIdx), Vec3d::Zero(), 1e-12, "bg derivative");
    ExpectNear(f.segment<3>(NavState::kBaIdx), Vec3d::Zero(), 1e-12, "ba derivative");
    ExpectNear(f.segment<3>(NavState::kGravIdx), Vec3d::Zero(), 1e-12, "gravity derivative");

    const auto F = s.df_dx(acc);
    ExpectNear(F.block<3, 3>(NavState::kPosIdx, NavState::kVelIdx), Mat3d::Identity(), 1e-12, "F_p_v");
    ExpectNear(F.block<3, 3>(NavState::kRotIdx, NavState::kBgIdx), -Mat3d::Identity(), 1e-12, "F_rot_bg");
    ExpectNear(F.block<3, 3>(NavState::kVelIdx, NavState::kRotIdx), -R * SO3::hat(acc_unbiased), 1e-12,
               "F_v_rot");
    ExpectNear(F.block<3, 3>(NavState::kVelIdx, NavState::kBaIdx), -R, 1e-12, "F_v_ba");
    ExpectNear(F.block<3, 3>(NavState::kVelIdx, NavState::kGravIdx), Mat3d::Identity(), 1e-12, "F_v_g");

    const auto G = s.df_dw();
    ExpectNear(G.block<3, 3>(NavState::kRotIdx, 0), -Mat3d::Identity(), 1e-12, "G_rot_ng");
    ExpectNear(G.block<3, 3>(NavState::kVelIdx, 3), -R, 1e-12, "G_vel_na");
    ExpectNear(G.block<3, 3>(NavState::kBgIdx, 6), Mat3d::Identity(), 1e-12, "G_bg_nbg");
    ExpectNear(G.block<3, 3>(NavState::kBaIdx, 9), Mat3d::Identity(), 1e-12, "G_ba_nba");

    NavState::VectState dx = NavState::VectState::Zero();
    dx.segment<3>(NavState::kPosIdx) = Vec3d(0.1, -0.2, 0.3);
    dx.segment<3>(NavState::kRotIdx) = Vec3d(0.01, -0.02, 0.03);
    dx.segment<3>(NavState::kVelIdx) = Vec3d(-0.1, 0.2, -0.3);
    dx.segment<3>(NavState::kBgIdx) = Vec3d(0.001, -0.002, 0.003);
    dx.segment<3>(NavState::kBaIdx) = Vec3d(-0.004, 0.005, -0.006);
    dx.segment<3>(NavState::kGravIdx) = Vec3d(0.01, 0.015, 0.0);

    const NavState plus = s.boxplus(dx);
    ExpectNear(plus.ba_, s.ba_ + dx.segment<3>(NavState::kBaIdx), 1e-12, "boxplus ba");
    ExpectNear(plus.grav_, NavState::NormalizeGravity(s.grav_ + dx.segment<3>(NavState::kGravIdx)), 1e-12,
               "boxplus gravity");
    ExpectNear(plus.boxminus(s).head<NavState::kGravIdx>(), dx.head<NavState::kGravIdx>(), 1e-12, "boxminus linear");
    ExpectNear(plus.boxminus(s).segment<3>(NavState::kGravIdx), dx.segment<3>(NavState::kGravIdx), 1e-5,
               "boxminus gravity");

    NavState::VectState large_dx = NavState::VectState::Zero();
    large_dx.segment<3>(NavState::kGravIdx) = Vec3d(10.0, 20.0, 30.0);
    const NavState limited = s.boxplus(large_dx);
    if (std::abs(limited.grav_.norm() - NavState::kGravityNorm) > 1e-12) {
        std::cerr << "gravity norm should stay fixed, got " << limited.grav_.norm() << std::endl;
        return 1;
    }
    if (limited.boxminus(s).segment<3>(NavState::kGravIdx).norm() > NavState::kMaxGravityDelta + 1e-8) {
        std::cerr << "gravity update should be limited" << std::endl;
        return 1;
    }

    NavState::FullVectState full = s.ToState();
    NavState restored;
    restored.FromVectState(full);
    ExpectNear(restored.ba_, s.ba_, 1e-12, "FromVectState ba");
    ExpectNear(restored.grav_, s.grav_, 1e-12, "FromVectState gravity");

    std::cout << "nav_state_eskf_model_test passed" << std::endl;
    return 0;
}
