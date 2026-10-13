#pragma once
#include "common/eigen_types.h"

namespace lightning::loc {
// VertexSE3 tangent: world translation, right body rotation. No LM damping.
bool MarginalPoseCovariance(const Eigen::MatrixXd& information, Mat6d& covariance);
Mat6d PropagatePoseCovariance(const Mat6d& covariance, const SE3& pose,
                            const SE3& increment, double dt,
                            double position_noise, double angle_noise);
Mat6d TransformPoseCovariance(const Mat6d& covariance, const SE3& pose,
                            const SE3& right_transform, const SO3& left_rotation);
}
