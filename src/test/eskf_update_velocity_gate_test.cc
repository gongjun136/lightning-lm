#include "core/lio/eskf.hpp"
#include "core/lio/imu_processing.hpp"

#include <cmath>
#include <iostream>
#include <memory>

namespace {
bool NearZeroVelocity(const lightning::NavState& state) {
    return state.vel_.norm() < 1e-9;
}

bool NearZeroVector(const lightning::Vec3d& value) {
    return value.norm() < 1e-9;
}

bool CovarianceIsPositiveSemidefinite(const lightning::ESKF::CovType& P) {
    lightning::ESKF::CovType sym = 0.5 * (P + P.transpose()).eval();
    Eigen::SelfAdjointEigenSolver<lightning::ESKF::CovType> solver(sym, Eigen::EigenvaluesOnly);
    return solver.info() == Eigen::Success && solver.eigenvalues().minCoeff() > -1e-9;
}

lightning::ESKF::Options MakeCrossCoupledLidarOptions() {
    using namespace lightning;

    ESKF::Options options;
    options.max_iterations_ = 1;
    options.epsi_ = ESKF::StateVecType::Constant(1e-12);
    options.max_update_translation_step_ = 0.5;
    options.max_update_rotation_step_deg_ = 5.0;
    options.lidar_obs_func_ = [](NavState&, ESKF::CustomObservationModel& obs) {
        obs.valid_ = true;
        obs.converge_ = true;
        obs.HTH_.setZero();
        obs.HTr_.setZero();
        obs.HTH_(0, 0) = 1.0;
        obs.HTr_(0) = 30.0;
        obs.lidar_residual_mean_ = 1.0;
    };
    return options;
}

void SetCrossCoupledCovariance(lightning::ESKF& eskf) {
    using namespace lightning;

    ESKF::CovType P = ESKF::CovType::Identity() * 1e-6;
    P(NavState::kPosIdx, NavState::kPosIdx) = 1e-4;
    P(NavState::kVelIdx, NavState::kVelIdx) = 100.0;
    P(NavState::kVelIdx, NavState::kPosIdx) = 0.09;
    P(NavState::kPosIdx, NavState::kVelIdx) = 0.09;
    eskf.ChangeP(P);
}

void SetBiasAndGravityCrossCoupledCovariance(lightning::ESKF& eskf) {
    using namespace lightning;

    ESKF::CovType P = ESKF::CovType::Identity() * 1e-6;
    P(NavState::kPosIdx, NavState::kPosIdx) = 1e-4;
    P(NavState::kBgIdx, NavState::kBgIdx) = 10.0;
    P(NavState::kBaIdx, NavState::kBaIdx) = 10.0;
    P(NavState::kGravIdx, NavState::kGravIdx) = 10.0;
    P(NavState::kBgIdx, NavState::kPosIdx) = 0.01;
    P(NavState::kPosIdx, NavState::kBgIdx) = 0.01;
    P(NavState::kBaIdx, NavState::kPosIdx) = 0.01;
    P(NavState::kPosIdx, NavState::kBaIdx) = 0.01;
    P(NavState::kGravIdx, NavState::kPosIdx) = 0.01;
    P(NavState::kPosIdx, NavState::kGravIdx) = 0.01;
    eskf.ChangeP(P);
}

bool FullLidarUpdateUsesCrossCovarianceForVelocity() {
    using namespace lightning;

    ESKF eskf;
    ESKF::Options options = MakeCrossCoupledLidarOptions();
    options.max_update_velocity_step_ = 100.0;
    eskf.Init(options);

    NavState state;
    state.pose_is_ok_ = true;
    eskf.ChangeX(state);
    SetCrossCoupledCovariance(eskf);

    eskf.Update(ESKF::ObsType::LIDAR, 1.0);

    const auto& updated = eskf.GetX();
    if (!eskf.LastUpdateAccepted()) {
        std::cerr << "Valid lidar update was not marked accepted." << std::endl;
        return false;
    }
    if (NearZeroVelocity(updated)) {
        std::cerr << "Full lidar update did not propagate pose information to velocity. velocity="
                  << updated.vel_.transpose() << std::endl;
        return false;
    }

    if (updated.pos_.norm() < 1e-6) {
        std::cerr << "Full lidar update did not update pose." << std::endl;
        return false;
    }

    return true;
}

bool DefaultLidarUpdateFallsBackToPoseOnlyForLargeVelocityStep() {
    using namespace lightning;

    ESKF eskf;
    ESKF::Options options = MakeCrossCoupledLidarOptions();
    eskf.Init(options);

    NavState state;
    state.pose_is_ok_ = true;
    eskf.ChangeX(state);
    SetCrossCoupledCovariance(eskf);

    eskf.Update(ESKF::ObsType::LIDAR, 1.0);

    const auto& updated = eskf.GetX();
    if (!NearZeroVelocity(updated)) {
        std::cerr << "Default lidar update accepted an unphysical velocity step. velocity="
                  << updated.vel_.transpose() << std::endl;
        return false;
    }

    if (updated.pos_.norm() < 1e-6) {
        std::cerr << "Default lidar update rejected the pose correction instead of falling back to pose-only."
                  << std::endl;
        return false;
    }

    return true;
}

bool LidarUpdateCanPreserveBiasAndGravityWhileUpdatingVelocity() {
    using namespace lightning;

    ESKF eskf;
    ESKF::Options options = MakeCrossCoupledLidarOptions();
    options.max_update_velocity_step_ = 100.0;
    options.lidar_update_inertial_states_ = false;
    eskf.Init(options);

    NavState state;
    state.pose_is_ok_ = true;
    eskf.ChangeX(state);

    ESKF::CovType P = ESKF::CovType::Identity() * 1e-6;
    P(NavState::kPosIdx, NavState::kPosIdx) = 1e-4;
    P(NavState::kVelIdx, NavState::kVelIdx) = 100.0;
    P(NavState::kVelIdx, NavState::kPosIdx) = 0.09;
    P(NavState::kPosIdx, NavState::kVelIdx) = 0.09;
    P(NavState::kBgIdx, NavState::kBgIdx) = 10.0;
    P(NavState::kBaIdx, NavState::kBaIdx) = 10.0;
    P(NavState::kGravIdx, NavState::kGravIdx) = 10.0;
    P(NavState::kBgIdx, NavState::kPosIdx) = 0.01;
    P(NavState::kPosIdx, NavState::kBgIdx) = 0.01;
    P(NavState::kBaIdx, NavState::kPosIdx) = 0.01;
    P(NavState::kPosIdx, NavState::kBaIdx) = 0.01;
    P(NavState::kGravIdx, NavState::kPosIdx) = 0.01;
    P(NavState::kPosIdx, NavState::kGravIdx) = 0.01;
    eskf.ChangeP(P);

    eskf.Update(ESKF::ObsType::LIDAR, 1.0);

    const auto& updated = eskf.GetX();
    if (NearZeroVelocity(updated)) {
        std::cerr << "Restricted lidar update did not update velocity through cross covariance." << std::endl;
        return false;
    }

    if (!NearZeroVector(updated.bg_)) {
        std::cerr << "Restricted lidar update changed gyro bias. bg=" << updated.bg_.transpose() << std::endl;
        return false;
    }

    if (!NearZeroVector(updated.ba_)) {
        std::cerr << "Restricted lidar update changed accel bias. ba=" << updated.ba_.transpose() << std::endl;
        return false;
    }

    if (!NearZeroVector(updated.grav_ - Vec3d(0.0, 0.0, -NavState::kGravityNorm))) {
        std::cerr << "Restricted lidar update changed gravity. grav=" << updated.grav_.transpose() << std::endl;
        return false;
    }

    if (!CovarianceIsPositiveSemidefinite(eskf.GetP())) {
        std::cerr << "Restricted lidar update produced a non-PSD covariance." << std::endl;
        return false;
    }

    return true;
}

bool FullLidarUpdateUsesCrossCovarianceForBiasAndGravity() {
    using namespace lightning;

    ESKF eskf;
    ESKF::Options options = MakeCrossCoupledLidarOptions();
    options.max_update_gyro_bias_step_ = 100.0;
    options.max_update_acc_bias_step_ = 100.0;
    options.max_update_gravity_step_ = 100.0;
    eskf.Init(options);

    NavState state;
    state.pose_is_ok_ = true;
    eskf.ChangeX(state);
    SetBiasAndGravityCrossCoupledCovariance(eskf);

    eskf.Update(ESKF::ObsType::LIDAR, 1.0);

    const auto& updated = eskf.GetX();
    if (NearZeroVector(updated.bg_)) {
        std::cerr << "Full lidar update did not update gyro bias through cross covariance." << std::endl;
        return false;
    }

    if (NearZeroVector(updated.ba_)) {
        std::cerr << "Full lidar update did not update accel bias through cross covariance." << std::endl;
        return false;
    }

    if (NearZeroVector(updated.grav_ - Vec3d(0.0, 0.0, -NavState::kGravityNorm))) {
        std::cerr << "Full lidar update did not update gravity direction through cross covariance." << std::endl;
        return false;
    }

    return true;
}

bool PoseOnlyLidarUpdatePreservesVelocity() {
    using namespace lightning;

    ESKF eskf;
    ESKF::Options options = MakeCrossCoupledLidarOptions();
    options.max_update_velocity_step_ = 100.0;
    options.lidar_update_pose_only_ = true;
    eskf.Init(options);

    NavState state;
    state.pose_is_ok_ = true;
    eskf.ChangeX(state);
    SetCrossCoupledCovariance(eskf);

    eskf.Update(ESKF::ObsType::LIDAR, 1.0);

    const auto& updated = eskf.GetX();
    if (!NearZeroVelocity(updated)) {
        std::cerr << "Pose-only lidar update changed velocity. velocity=" << updated.vel_.transpose() << std::endl;
        return false;
    }

    if (updated.pos_.norm() < 1e-6) {
        std::cerr << "Pose-only lidar update did not update pose." << std::endl;
        return false;
    }

    return true;
}

bool LidarUpdateFallsBackWhenInertialStepIsUnphysical() {
    using namespace lightning;

    ESKF eskf;
    ESKF::Options options = MakeCrossCoupledLidarOptions();
    options.max_update_velocity_step_ = 100.0;
    eskf.Init(options);

    NavState state;
    state.pose_is_ok_ = true;
    eskf.ChangeX(state);

    ESKF::CovType P = ESKF::CovType::Identity() * 1e-6;
    P(NavState::kPosIdx, NavState::kPosIdx) = 1e-4;
    P(NavState::kBgIdx, NavState::kBgIdx) = 1e4;
    P(NavState::kBaIdx, NavState::kBaIdx) = 1e4;
    P(NavState::kGravIdx, NavState::kGravIdx) = 1e4;
    P(NavState::kBgIdx, NavState::kPosIdx) = 0.9;
    P(NavState::kPosIdx, NavState::kBgIdx) = 0.9;
    P(NavState::kBaIdx, NavState::kPosIdx) = 0.9;
    P(NavState::kPosIdx, NavState::kBaIdx) = 0.9;
    P(NavState::kGravIdx, NavState::kPosIdx) = 0.9;
    P(NavState::kPosIdx, NavState::kGravIdx) = 0.9;
    eskf.ChangeP(P);

    eskf.Update(ESKF::ObsType::LIDAR, 1.0);

    const auto& updated = eskf.GetX();
    if (!NearZeroVector(updated.bg_)) {
        std::cerr << "Unphysical lidar inertial update changed gyro bias. bg=" << updated.bg_.transpose()
                  << std::endl;
        return false;
    }

    if (!NearZeroVector(updated.ba_)) {
        std::cerr << "Unphysical lidar inertial update changed accel bias. ba=" << updated.ba_.transpose()
                  << std::endl;
        return false;
    }

    if (!NearZeroVector(updated.grav_ - Vec3d(0.0, 0.0, -NavState::kGravityNorm))) {
        std::cerr << "Unphysical lidar inertial update changed gravity. grav=" << updated.grav_.transpose()
                  << std::endl;
        return false;
    }

    if (updated.pos_.norm() < 1e-6) {
        std::cerr << "Unphysical inertial fallback rejected the pose correction." << std::endl;
        return false;
    }

    return true;
}

bool InformationFormMatchesDenseKalmanUpdate() {
    using namespace lightning;

    ESKF::CovType a = ESKF::CovType::Zero();
    for (int r = 0; r < ESKF::state_dim_; ++r) {
        for (int c = 0; c < ESKF::state_dim_; ++c) {
            a(r, c) = 0.01 * static_cast<double>((r + 1) * (c + 3) % 7);
        }
    }

    ESKF::CovType P = a * a.transpose();
    P.diagonal().array() += 0.1;

    Mat6d hth = Mat6d::Zero();
    hth.diagonal() << 30.0, 20.0, 15.0, 10.0, 8.0, 5.0;
    hth(0, 3) = hth(3, 0) = 1.5;
    hth(1, 4) = hth(4, 1) = -0.7;

    Vec6d htr;
    htr << 0.4, -0.2, 0.1, 0.05, -0.03, 0.02;

    const double R = 0.25;

    ESKF::CovType info = (P / R).inverse();
    info.block<ESKF::pose_obs_dim_, ESKF::pose_obs_dim_>(0, 0) += hth;
    ESKF::CovType q_inv = info.inverse();
    ESKF::StateVecType kr_info = q_inv.block<ESKF::state_dim_, ESKF::pose_obs_dim_>(0, 0) * htr;
    ESKF::CovType kh_info = ESKF::CovType::Zero();
    kh_info.block<ESKF::state_dim_, ESKF::pose_obs_dim_>(0, 0) =
        q_inv.block<ESKF::state_dim_, ESKF::pose_obs_dim_>(0, 0) * hth;
    ESKF::CovType p_info =
        P - kh_info.block<ESKF::state_dim_, ESKF::pose_obs_dim_>(0, 0) *
                P.block<ESKF::pose_obs_dim_, ESKF::state_dim_>(0, 0);

    Eigen::Matrix<double, ESKF::pose_obs_dim_, ESKF::state_dim_> H =
        Eigen::Matrix<double, ESKF::pose_obs_dim_, ESKF::state_dim_>::Zero();
    H.block<ESKF::pose_obs_dim_, ESKF::pose_obs_dim_>(0, 0).setIdentity();
    ESKF::CovType p_kalman =
        (P.inverse() + H.transpose() * (hth / R) * H).inverse();
    ESKF::StateVecType kr_kalman =
        p_kalman * H.transpose() * (htr / R);

    if ((kr_info - kr_kalman).norm() > 1e-9) {
        std::cerr << "Information-form dx differs from dense Kalman dx: "
                  << (kr_info - kr_kalman).transpose() << std::endl;
        return false;
    }

    if ((p_info - p_kalman).norm() > 1e-9) {
        std::cerr << "Information-form covariance differs from dense Kalman covariance. diff norm="
                  << (p_info - p_kalman).norm() << std::endl;
        return false;
    }

    return true;
}

bool ZeroDurationPredictionDoesNotInflateCovariance() {
    using namespace lightning;

    ESKF eskf;
    ESKF::Options options;
    eskf.Init(options);

    const ESKF::CovType initial_covariance = ESKF::CovType::Identity() * 0.25;
    eskf.ChangeP(initial_covariance);
    const ESKF::ProcessNoiseType zero_process_noise = ESKF::ProcessNoiseType::Zero();
    for (int i = 0; i < 1000; ++i) {
        eskf.Predict(0.0, zero_process_noise, Vec3d::Zero(), Vec3d::Zero());
    }

    const double covariance_change = (eskf.GetP() - initial_covariance).norm();
    if (covariance_change > 1e-12) {
        std::cerr << "Zero-duration prediction inflated covariance. diff norm=" << covariance_change << std::endl;
        return false;
    }
    return true;
}

bool ForwardSpeedUpdateChangesOnlyVelocity() {
    using namespace lightning;

    NavState state;
    state.pose_is_ok_ = true;
    state.pos_ = Vec3d(3.0, -2.0, 1.0);
    state.bg_ = Vec3d(0.01, 0.02, 0.03);
    state.ba_ = Vec3d(0.1, 0.2, 0.3);
    ESKF eskf(state, ESKF::CovType::Identity());
    ESKF::Options options;
    eskf.Init(options);

    const NavState before = eskf.GetX();
    const auto result = eskf.UpdateBodyForwardSpeed(1.0, 0.25, 2.0, 9.0, 1.0);
    const NavState& after = eskf.GetX();
    if (!result.accepted || std::abs(after.vel_.x() - 0.8) > 1e-12) {
        std::cerr << "Forward speed update did not produce the expected velocity. velocity="
                  << after.vel_.transpose() << std::endl;
        return false;
    }
    if ((after.pos_ - before.pos_).norm() > 1e-12 ||
        (after.rot_.inverse() * before.rot_).log().norm() > 1e-12 ||
        (after.bg_ - before.bg_).norm() > 1e-12 ||
        (after.ba_ - before.ba_).norm() > 1e-12 ||
        (after.grav_ - before.grav_).norm() > 1e-12) {
        std::cerr << "Forward speed update changed a non-velocity state." << std::endl;
        return false;
    }
    return CovarianceIsPositiveSemidefinite(eskf.GetP());
}

bool ForwardSpeedUpdateUsesBodyHeading() {
    using namespace lightning;

    NavState state;
    state.rot_ = SO3::exp(Vec3d(0.0, 0.0, M_PI_2));
    ESKF eskf(state, ESKF::CovType::Identity());
    ESKF::Options options;
    eskf.Init(options);

    const auto result = eskf.UpdateBodyForwardSpeed(1.0, 0.25, 2.0, 9.0, 1.0);
    if (!result.accepted || std::abs(eskf.GetX().vel_.y() - 0.8) > 1e-12 ||
        std::abs(eskf.GetX().vel_.x()) > 1e-12) {
        std::cerr << "Forward speed update ignored body heading. velocity="
                  << eskf.GetX().vel_.transpose() << std::endl;
        return false;
    }
    return true;
}

bool ForwardSpeedUpdateRejectsAbsoluteInnovationOutlier() {
    using namespace lightning;

    ESKF eskf;
    ESKF::Options options;
    eskf.Init(options);
    const NavState before = eskf.GetX();
    const ESKF::CovType covariance_before = eskf.GetP();
    const auto result = eskf.UpdateBodyForwardSpeed(3.0, 0.25, 1.0, 100.0, 1.0);
    if (result.accepted || (eskf.GetX().boxminus(before)).norm() > 1e-12 ||
        (eskf.GetP() - covariance_before).norm() > 1e-12) {
        std::cerr << "Absolute innovation outlier changed the filter." << std::endl;
        return false;
    }
    return true;
}

bool ForwardSpeedUpdateRejectsNisOutlier() {
    using namespace lightning;

    ESKF eskf;
    ESKF::Options options;
    eskf.Init(options);
    ESKF::CovType covariance = ESKF::CovType::Identity() * 1e-6;
    eskf.ChangeP(covariance);
    const NavState before = eskf.GetX();
    const auto result = eskf.UpdateBodyForwardSpeed(0.5, 0.01, 2.0, 9.0, 1.0);
    if (result.accepted || result.normalized_innovation_squared <= 9.0 ||
        (eskf.GetX().boxminus(before)).norm() > 1e-12) {
        std::cerr << "NIS outlier was not rejected. nis="
                  << result.normalized_innovation_squared << std::endl;
        return false;
    }
    return true;
}

bool ForwardSpeedUpdateRejectsLargeVelocityStep() {
    using namespace lightning;

    ESKF eskf;
    ESKF::Options options;
    eskf.Init(options);
    const NavState before = eskf.GetX();
    const auto result = eskf.UpdateBodyForwardSpeed(1.0, 0.01, 2.0, 100.0, 0.2);
    if (result.accepted || (eskf.GetX().boxminus(before)).norm() > 1e-12) {
        std::cerr << "Large velocity step was not rejected." << std::endl;
        return false;
    }
    return true;
}

bool InvalidLidarUpdateIsNotMarkedAccepted() {
    using namespace lightning;

    ESKF eskf;
    ESKF::Options options = MakeCrossCoupledLidarOptions();
    options.lidar_obs_func_ = [](NavState&, ESKF::CustomObservationModel& obs) { obs.valid_ = false; };
    eskf.Init(options);
    eskf.Update(ESKF::ObsType::LIDAR, 1.0);
    if (eskf.LastUpdateAccepted()) {
        std::cerr << "Invalid lidar update was marked accepted." << std::endl;
        return false;
    }
    return true;
}

bool ImuPredictionPathAppliesConfiguredAccelerationScale() {
    using namespace lightning;

    ImuProcess imu;
    imu.SetAccCov(Vec3d::Ones());
    imu.SetGyrCov(Vec3d::Ones());

    ESKF eskf;
    MeasureGroup meas;
    for (int i = 0; i < 21; ++i) {
        auto sample = std::make_shared<IMU>();
        sample->timestamp = 0.01 * static_cast<double>(i);
        sample->linear_acceleration = Vec3d(0.0, 0.0, 1.0);
        sample->angular_velocity = Vec3d::Zero();
        meas.imu_.push_back(sample);
    }
    CloudPtr scan;
    imu.Process(meas, eskf, scan);

    const Vec3d raw_acc_in_g(0.0, 0.0, 1.0);
    const Vec3d scaled_acc = imu.ScaleAccelerationForPrediction(raw_acc_in_g);
    const Vec3d expected(0.0, 0.0, 9.81);
    if ((scaled_acc - expected).norm() > 1e-12) {
        std::cerr << "IMU prediction acceleration scale mismatch. got=" << scaled_acc.transpose()
                  << ", expected=" << expected.transpose() << std::endl;
        return false;
    }

    return true;
}

bool ImuInitializationRejectsInvalidMeanAccelerationNorm() {
    using namespace lightning;

    ImuProcess imu;
    ImuProcess::InitializationOptions options;
    options.min_duration = 0.1;
    options.min_samples = 20;
    options.min_mean_acc_norm = 0.5;
    options.max_mean_acc_norm = 1.5;
    imu.SetInitializationOptions(options);

    ESKF eskf;
    MeasureGroup meas;
    for (int i = 0; i < 21; ++i) {
        auto sample = std::make_shared<IMU>();
        sample->timestamp = 0.01 * static_cast<double>(i);
        sample->linear_acceleration = Vec3d(0.0, 0.0, 3.0);
        sample->angular_velocity = Vec3d::Zero();
        meas.imu_.push_back(sample);
    }
    CloudPtr scan;
    imu.Process(meas, eskf, scan);
    if (imu.IsIMUInited()) {
        std::cerr << "IMU initialization accepted an invalid mean acceleration norm." << std::endl;
        return false;
    }
    return true;
}

bool ImuInitializationAppliesConfiguredHeadingOffset() {
    using namespace lightning;

    ImuProcess imu;
    imu.SetAccCov(Vec3d::Ones());
    imu.SetGyrCov(Vec3d::Ones());
    ImuProcess::InitializationOptions options;
    options.min_duration = 0.1;
    options.min_samples = 20;
    options.min_mean_acc_norm = 0.5;
    options.max_mean_acc_norm = 1.5;
    options.initial_yaw_deg = 180.0;
    imu.SetInitializationOptions(options);

    constexpr double pitch = 0.6;
    const Vec3d stationary_acc(std::sin(pitch), 0.0, std::cos(pitch));
    ESKF eskf;
    MeasureGroup meas;
    for (int i = 0; i < 21; ++i) {
        auto sample = std::make_shared<IMU>();
        sample->timestamp = 0.01 * static_cast<double>(i);
        sample->linear_acceleration = stationary_acc;
        sample->angular_velocity = Vec3d::Zero();
        meas.imu_.push_back(sample);
    }
    CloudPtr scan;
    imu.Process(meas, eskf, scan);
    if (!imu.IsIMUInited()) {
        std::cerr << "IMU initialization did not finish for valid stationary samples." << std::endl;
        return false;
    }

    const SO3 gravity_alignment(Quatd::FromTwoVectors(stationary_acc, Vec3d::UnitZ()).normalized());
    const SO3 expected = SO3::exp(Vec3d(0.0, 0.0, M_PI)) * gravity_alignment;
    const SO3 actual = eskf.GetX().rot_;
    if ((expected.inverse() * actual).log().norm() > 1e-12) {
        std::cerr << "IMU initialization heading offset mismatch." << std::endl;
        return false;
    }
    if ((actual * stationary_acc - Vec3d::UnitZ()).norm() > 1e-12) {
        std::cerr << "Heading offset changed gravity alignment." << std::endl;
        return false;
    }
    return true;
}
}  // namespace

int main() {
    if (!InvalidLidarUpdateIsNotMarkedAccepted()) {
        return 1;
    }

    if (!InformationFormMatchesDenseKalmanUpdate()) {
        return 1;
    }

    if (!ZeroDurationPredictionDoesNotInflateCovariance()) {
        return 1;
    }

    if (!ForwardSpeedUpdateChangesOnlyVelocity()) {
        return 1;
    }

    if (!ForwardSpeedUpdateUsesBodyHeading()) {
        return 1;
    }

    if (!ForwardSpeedUpdateRejectsAbsoluteInnovationOutlier()) {
        return 1;
    }

    if (!ForwardSpeedUpdateRejectsNisOutlier()) {
        return 1;
    }

    if (!ForwardSpeedUpdateRejectsLargeVelocityStep()) {
        return 1;
    }

    if (!FullLidarUpdateUsesCrossCovarianceForVelocity()) {
        return 1;
    }

    if (!DefaultLidarUpdateFallsBackToPoseOnlyForLargeVelocityStep()) {
        return 1;
    }

    if (!LidarUpdateCanPreserveBiasAndGravityWhileUpdatingVelocity()) {
        return 1;
    }

    if (!FullLidarUpdateUsesCrossCovarianceForBiasAndGravity()) {
        return 1;
    }

    if (!PoseOnlyLidarUpdatePreservesVelocity()) {
        return 1;
    }

    if (!LidarUpdateFallsBackWhenInertialStepIsUnphysical()) {
        return 1;
    }

    if (!ImuPredictionPathAppliesConfiguredAccelerationScale()) {
        return 1;
    }

    if (!ImuInitializationRejectsInvalidMeanAccelerationNorm()) {
        return 1;
    }

    if (!ImuInitializationAppliesConfiguredHeadingOffset()) {
        return 1;
    }

    return 0;
}
