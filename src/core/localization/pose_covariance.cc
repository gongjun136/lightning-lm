#include "core/localization/pose_covariance.h"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <algorithm>

namespace lightning::loc {
bool MarginalPoseCovariance(const Eigen::MatrixXd& information, Mat6d& covariance) {
    const int n = information.rows();
    if (n < 6 || n != information.cols() || !information.allFinite()) return false;
    const Eigen::MatrixXd h = 0.5 * (information + information.transpose());
    Mat6d schur = h.bottomRightCorner<6, 6>();
    if (n > 6) {
        const Eigen::MatrixXd removed = h.topLeftCorner(n - 6, n - 6);
        Eigen::LDLT<Eigen::MatrixXd> solve(removed);
        if (solve.info() != Eigen::Success || !solve.isPositive() ||
            solve.vectorD().minCoeff() <= 1e-12) return false;
        schur -= h.bottomLeftCorner(6, n - 6) * solve.solve(h.topRightCorner(n - 6, 6));
    }
    schur = (0.5 * (schur + schur.transpose())).eval();
    Eigen::SelfAdjointEigenSolver<Mat6d> eigen(schur);
    if (eigen.info() != Eigen::Success || eigen.eigenvalues().minCoeff() <= 1e-12) return false;
    Eigen::LDLT<Mat6d> solve(schur);
    covariance = solve.solve(Mat6d::Identity());
    covariance = (0.5 * (covariance + covariance.transpose())).eval();
    return covariance.allFinite() && solve.info() == Eigen::Success;
}
Mat6d PropagatePoseCovariance(const Mat6d& covariance, const SE3& pose,
                            const SE3& increment, double dt,
                            double position_noise, double angle_noise) {
    Mat6d f = Mat6d::Identity();
    f.topRightCorner<3, 3>() = -pose.rotationMatrix() * SO3::hat(increment.translation());
    f.bottomRightCorner<3, 3>() = increment.rotationMatrix().transpose();
    Mat6d q = Mat6d::Zero();
    q.diagonal().head<3>().setConstant(position_noise * position_noise * std::max(0.0, dt));
    q.diagonal().tail<3>().setConstant(angle_noise * angle_noise * std::max(0.0, dt));
    return f * covariance * f.transpose() + q;
}
Mat6d TransformPoseCovariance(const Mat6d& covariance, const SE3& pose,
                            const SE3& right_transform, const SO3& left_rotation) {
    Mat6d j = Mat6d::Zero();
    j.topLeftCorner<3, 3>() = left_rotation.matrix();
    j.topRightCorner<3, 3>() = -left_rotation.matrix() * pose.rotationMatrix() * SO3::hat(right_transform.translation());
    j.bottomRightCorner<3, 3>() = right_transform.rotationMatrix().transpose();
    return j * covariance * j.transpose();
}
}
