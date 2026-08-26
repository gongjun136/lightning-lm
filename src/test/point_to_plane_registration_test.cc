#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>

#include "core/localization/point_to_plane_registration.h"

namespace {

using lightning::CloudPtr;
using lightning::PointCloudType;
using lightning::PointType;
using lightning::SE3;
using lightning::SO3;
using lightning::Vec3d;
using lightning::Vec6d;
using lightning::loc::PointToPlaneRegistration;

void AddPoint(PointCloudType& cloud, const Vec3d& position) {
    PointType point;
    point.x = static_cast<float>(position.x());
    point.y = static_cast<float>(position.y());
    point.z = static_cast<float>(position.z());
    cloud.push_back(point);
}

CloudPtr MakeThreePlaneTarget() {
    CloudPtr target(new PointCloudType);
    for (double first = -4.0; first <= 4.0; first += 0.25) {
        for (double second = -4.0; second <= 4.0; second += 0.25) {
            AddPoint(*target, Vec3d(5.0, first, second));
            AddPoint(*target, Vec3d(first, 6.0, second));
            AddPoint(*target, Vec3d(first, second, 3.0));
        }
    }
    return target;
}

CloudPtr MakeSinglePlaneTarget() {
    CloudPtr target(new PointCloudType);
    for (double y = -4.0; y <= 4.0; y += 0.25) {
        for (double z = -4.0; z <= 4.0; z += 0.25) {
            AddPoint(*target, Vec3d(5.0, y, z));
        }
    }
    return target;
}

CloudPtr TransformToSource(const CloudPtr& target, const SE3& target_from_source) {
    CloudPtr source(new PointCloudType);
    source->reserve(target->size());
    const SE3 source_from_target = target_from_source.inverse();
    for (const auto& target_point : target->points) {
        AddPoint(*source, source_from_target *
                              Vec3d(target_point.x, target_point.y,
                                    target_point.z));
    }
    return source;
}

CloudPtr MakeNoisyThreePlaneSource(const CloudPtr& target,
                                   const SE3& target_from_source) {
    CloudPtr source(new PointCloudType);
    source->reserve(target->size());
    const SE3 source_from_target = target_from_source.inverse();
    for (std::size_t index = 0; index < target->size(); ++index) {
        const auto& point = target->points[index];
        Vec3d noisy_point(point.x, point.y, point.z);
        const double noise = (index / 3) % 2 == 0 ? 0.03 : -0.03;
        noisy_point[static_cast<int>(index % 3)] += noise;
        AddPoint(*source, source_from_target * noisy_point);
    }
    return source;
}

PointToPlaneRegistration::Options TestOptions() {
    PointToPlaneRegistration::Options options;
    options.target_voxel_size = 0.15;
    options.source_voxel_size = 0.15;
    options.plane_fit_threshold = 0.04;
    options.coarse_max_correspondence_distance = 1.0;
    options.fine_max_correspondence_distance = 0.4;
    options.huber_delta = 0.15;
    options.coarse_max_iterations = 20;
    options.fine_max_iterations = 20;
    options.min_matches = 300;
    options.min_inlier_ratio = 0.25;
    options.translation_convergence = 1e-5;
    options.rotation_convergence_deg = 0.001;
    options.min_normalized_hessian_eigenvalue = 1e-5;
    options.max_hessian_condition_number = 1e6;
    options.max_translation_correction = 1.0;
    options.max_rotation_correction_deg = 10.0;
    return options;
}

SE3 TruePose() {
    const Vec3d axis = Vec3d(0.2, -0.1, 0.3).normalized();
    return SE3(SO3::exp(axis * (5.0 * M_PI / 180.0)),
               Vec3d(1.0, -0.5, 0.3));
}

SE3 InitialPose(const SE3& true_pose) {
    Vec6d error = Vec6d::Zero();
    error.head<3>() = Vec3d(0.20, -0.15, 0.10);
    error.tail<3>() =
        Vec3d(-0.3, 0.4, 0.2).normalized() * (2.0 * M_PI / 180.0);
    return SE3::exp(error) * true_pose;
}

bool TestConvergence() {
    const CloudPtr target = MakeThreePlaneTarget();
    const SE3 true_pose = TruePose();
    const CloudPtr source = TransformToSource(target, true_pose);
    PointToPlaneRegistration registration(TestOptions());
    if (!registration.SetTarget(target)) {
        std::cerr << "failed to build three-plane target\n";
        return false;
    }

    SE3 estimate = InitialPose(true_pose);
    const auto result = registration.Refine(source, estimate);
    const SE3 error = true_pose.inverse() * estimate;
    if (!result.success || !result.converged || !result.attempted_pose_valid ||
        !result.attempted_pose.matrix().allFinite() ||
        (result.attempted_pose.matrix() - estimate.matrix())
                .cwiseAbs()
                .maxCoeff() != 0.0 ||
        result.iterations <= 0 ||
        result.matches < TestOptions().min_matches ||
        result.inlier_ratio < TestOptions().min_inlier_ratio ||
        !std::isfinite(result.rmse) || result.rmse > 0.01 ||
        error.translation().norm() > 0.01 ||
        error.so3().log().norm() > 0.1 * M_PI / 180.0) {
        std::cerr << "three-plane convergence failed: success=" << result.success
                  << " converged=" << result.converged
                  << " iterations=" << result.iterations
                  << " matches=" << result.matches
                  << " ratio=" << result.inlier_ratio
                  << " rmse=" << result.rmse
                  << " translation_error=" << error.translation().norm()
                  << " rotation_error_deg="
                  << error.so3().log().norm() * 180.0 / M_PI << '\n';
        return false;
    }
    return true;
}

bool TestCorrectionBoundDoesNotMutatePose() {
    const CloudPtr target = MakeThreePlaneTarget();
    const SE3 true_pose = TruePose();
    const CloudPtr source = TransformToSource(target, true_pose);
    auto options = TestOptions();
    options.max_translation_correction = 0.001;
    PointToPlaneRegistration registration(options);
    if (!registration.SetTarget(target)) return false;

    SE3 estimate = InitialPose(true_pose);
    const SE3 original = estimate;
    const auto result = registration.Refine(source, estimate);
    if (result.success || !result.attempted_pose_valid ||
        !result.attempted_pose.matrix().allFinite() ||
        (result.attempted_pose.matrix() - original.matrix())
                .cwiseAbs()
                .maxCoeff() == 0.0 ||
        (estimate.matrix() - original.matrix()).cwiseAbs().maxCoeff() != 0.0) {
        std::cerr << "correction-bound rejection pose diagnostics are invalid\n";
        return false;
    }
    return true;
}

bool TestRmseBoundRejectsWithoutMutatingPose() {
    const CloudPtr target = MakeThreePlaneTarget();
    const SE3 true_pose = TruePose();
    const CloudPtr source = MakeNoisyThreePlaneSource(target, true_pose);
    auto options = TestOptions();
    options.max_rmse = 0.01;
    PointToPlaneRegistration registration(options);
    if (!registration.SetTarget(target)) return false;

    SE3 estimate = InitialPose(true_pose);
    const SE3 original = estimate;
    const auto result = registration.Refine(source, estimate);
    if (result.success || !result.converged ||
        !std::isfinite(result.rmse) || result.rmse <= options.max_rmse ||
        !result.attempted_pose_valid ||
        !result.attempted_pose.matrix().allFinite() ||
        (result.attempted_pose.matrix() - original.matrix())
                .cwiseAbs()
                .maxCoeff() == 0.0 ||
        std::memcmp(estimate.data(), original.data(),
                    SE3::num_parameters * sizeof(double)) != 0) {
        std::cerr << "RMSE-bound rejection diagnostics are invalid: success="
                  << result.success << " converged=" << result.converged
                  << " rmse=" << result.rmse
                  << " threshold=" << options.max_rmse << '\n';
        return false;
    }
    return true;
}

bool TestMaxRmseValidation() {
    auto options = TestOptions();
    options.max_rmse = std::numeric_limits<double>::infinity();
    if (!PointToPlaneRegistration::ValidateOptions(options)) return false;

    const double invalid_values[] = {
        std::numeric_limits<double>::quiet_NaN(), 0.0, -1.0};
    for (double value : invalid_values) {
        options.max_rmse = value;
        if (PointToPlaneRegistration::ValidateOptions(options)) {
            std::cerr << "accepted invalid max_rmse=" << value << '\n';
            return false;
        }
    }
    return true;
}

bool TestDegeneratePlaneRejected() {
    const CloudPtr target = MakeSinglePlaneTarget();
    const SE3 true_pose = TruePose();
    const CloudPtr source = TransformToSource(target, true_pose);
    auto options = TestOptions();
    options.min_matches = 100;
    PointToPlaneRegistration registration(options);
    if (!registration.SetTarget(target)) return false;

    SE3 estimate = InitialPose(true_pose);
    const SE3 original = estimate;
    const auto result = registration.Refine(source, estimate);
    if (result.success ||
        (result.attempted_pose_valid &&
         !result.attempted_pose.matrix().allFinite()) ||
        (estimate.matrix() - original.matrix()).cwiseAbs().maxCoeff() != 0.0) {
        std::cerr << "degenerate single-plane rejection diagnostics are invalid\n";
        return false;
    }
    return true;
}

bool TestEmptySourceDoesNotAttemptOrMutatePose() {
    const CloudPtr target = MakeThreePlaneTarget();
    PointToPlaneRegistration registration(TestOptions());
    if (!registration.SetTarget(target)) return false;

    const CloudPtr source(new PointCloudType);
    SE3 estimate = InitialPose(TruePose());
    const SE3 original = estimate;
    const auto result = registration.Refine(source, estimate);
    if (result.success || result.attempted_pose_valid ||
        (estimate.matrix() - original.matrix()).cwiseAbs().maxCoeff() != 0.0) {
        std::cerr << "empty-source rejection diagnostics are invalid\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!TestConvergence()) return 1;
    if (!TestCorrectionBoundDoesNotMutatePose()) return 1;
    if (!TestRmseBoundRejectsWithoutMutatingPose()) return 1;
    if (!TestMaxRmseValidation()) return 1;
    if (!TestDegeneratePlaneRejected()) return 1;
    if (!TestEmptySourceDoesNotAttemptOrMutatePose()) return 1;
    return 0;
}
