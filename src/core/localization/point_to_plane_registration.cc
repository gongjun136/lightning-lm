#include "core/localization/point_to_plane_registration.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <pcl/filters/voxel_grid.h>

#include "core/lightning_math.hpp"

namespace lightning::loc {
namespace {

constexpr int kNeighborCount = 5;
constexpr double kDegreesPerRadian = 180.0 / M_PI;
constexpr int kMaximumTotalIterations = 100;

bool IsFinitePoint(const PointType& point) {
    return std::isfinite(point.x) && std::isfinite(point.y) &&
           std::isfinite(point.z);
}

void SetError(std::string* error, const std::string& message) {
    if (error) *error = message;
}

}  // namespace

PointToPlaneRegistration::PointToPlaneRegistration()
    : PointToPlaneRegistration(Options{}) {}

PointToPlaneRegistration::PointToPlaneRegistration(Options options)
    : options_(std::move(options)), target_(new PointCloudType) {
    std::string error;
    if (!ValidateOptions(options_, &error)) {
        throw std::invalid_argument(error);
    }
}

bool PointToPlaneRegistration::ValidateOptions(const Options& options,
                                               std::string* error) {
    const std::array<double, 13> finite_values = {
        options.target_voxel_size,
        options.source_voxel_size,
        options.plane_fit_threshold,
        options.coarse_max_correspondence_distance,
        options.fine_max_correspondence_distance,
        options.huber_delta,
        options.min_inlier_ratio,
        options.translation_convergence,
        options.rotation_convergence_deg,
        options.min_normalized_hessian_eigenvalue,
        options.max_hessian_condition_number,
        options.max_translation_correction,
        options.max_rotation_correction_deg};
    if (!std::all_of(finite_values.begin(), finite_values.end(),
                     [](double value) { return std::isfinite(value); })) {
        SetError(error, "relocalization.plane_icp values must be finite");
        return false;
    }
    if (std::isnan(options.max_rmse) || options.max_rmse <= 0.0) {
        SetError(error,
                 "relocalization.plane_icp max_rmse must be positive");
        return false;
    }
    if (options.target_voxel_size <= 0.0 ||
        options.source_voxel_size <= 0.0 ||
        options.plane_fit_threshold <= 0.0 ||
        options.coarse_max_correspondence_distance <= 0.0 ||
        options.fine_max_correspondence_distance <= 0.0 ||
        options.fine_max_correspondence_distance >
            options.coarse_max_correspondence_distance ||
        options.huber_delta <= 0.0) {
        SetError(error, "relocalization.plane_icp distances must be positive and fine_max_correspondence_distance must not exceed coarse_max_correspondence_distance");
        return false;
    }
    if (options.coarse_max_iterations <= 0 ||
        options.fine_max_iterations <= 0 ||
        options.coarse_max_iterations + options.fine_max_iterations >
            kMaximumTotalIterations) {
        SetError(error, "relocalization.plane_icp iteration counts must be positive and total at most 100");
        return false;
    }
    if (options.min_matches < 6 || options.min_inlier_ratio <= 0.0 ||
        options.min_inlier_ratio > 1.0) {
        SetError(error, "relocalization.plane_icp min_matches must be at least 6 and min_inlier_ratio must be in (0, 1]");
        return false;
    }
    if (options.translation_convergence <= 0.0 ||
        options.rotation_convergence_deg <= 0.0 ||
        options.min_normalized_hessian_eigenvalue <= 0.0 ||
        options.min_normalized_hessian_eigenvalue >= 1.0 ||
        options.max_hessian_condition_number <= 1.0 ||
        options.max_translation_correction <= 0.0 ||
        options.max_rotation_correction_deg <= 0.0 ||
        options.max_rotation_correction_deg > 180.0) {
        SetError(error, "relocalization.plane_icp convergence, observability, conditioning, and correction bounds are invalid");
        return false;
    }
    return true;
}

CloudPtr PointToPlaneRegistration::DownsampleFinite(const CloudPtr& cloud,
                                                    double voxel_size) const {
    CloudPtr finite(new PointCloudType);
    if (!cloud) return finite;
    finite->reserve(cloud->size());
    for (const auto& point : cloud->points) {
        if (IsFinitePoint(point)) finite->push_back(point);
    }
    if (finite->empty()) return finite;

    pcl::VoxelGrid<PointType> voxel;
    voxel.setLeafSize(static_cast<float>(voxel_size),
                      static_cast<float>(voxel_size),
                      static_cast<float>(voxel_size));
    voxel.setInputCloud(finite);
    CloudPtr downsampled(new PointCloudType);
    voxel.filter(*downsampled);
    return downsampled;
}

bool PointToPlaneRegistration::SetTarget(const CloudPtr& target) {
    CloudPtr downsampled = DownsampleFinite(target, options_.target_voxel_size);
    if (downsampled->size() < kNeighborCount) {
        target_.reset(new PointCloudType);
        target_kdtree_.reset();
        return false;
    }

    auto kdtree = std::make_unique<pcl::KdTreeFLANN<PointType>>();
    kdtree->setInputCloud(downsampled);
    target_ = std::move(downsampled);
    target_kdtree_ = std::move(kdtree);
    return true;
}

bool PointToPlaneRegistration::HasTarget() const {
    return target_kdtree_ && target_ && target_->size() >= kNeighborCount;
}

std::size_t PointToPlaneRegistration::TargetSize() const {
    return target_ ? target_->size() : 0;
}

PointToPlaneRegistration::Linearization PointToPlaneRegistration::Linearize(
    const CloudPtr& source, const SE3& pose,
    double max_correspondence_distance) const {
    Linearization linearization;
    if (!source || source->empty() || !HasTarget() ||
        !pose.matrix().allFinite()) {
        return linearization;
    }

    const double maximum_distance_squared =
        max_correspondence_distance * max_correspondence_distance;
    std::vector<int> neighbor_indices(kNeighborCount);
    std::vector<float> neighbor_distances_squared(kNeighborCount);
    PointVector neighbors;
    neighbors.reserve(kNeighborCount);

    for (const auto& source_point : source->points) {
        if (!IsFinitePoint(source_point)) continue;
        const Vec3d transformed = pose * ToVec3d(source_point);
        if (!transformed.allFinite()) continue;

        PointType query = source_point;
        query.x = static_cast<float>(transformed.x());
        query.y = static_cast<float>(transformed.y());
        query.z = static_cast<float>(transformed.z());
        if (!IsFinitePoint(query) ||
            target_kdtree_->nearestKSearch(
                query, kNeighborCount, neighbor_indices,
                neighbor_distances_squared) != kNeighborCount ||
            neighbor_distances_squared.back() > maximum_distance_squared) {
            continue;
        }

        neighbors.clear();
        for (int neighbor_index : neighbor_indices) {
            neighbors.push_back(target_->points[static_cast<std::size_t>(neighbor_index)]);
        }
        Vec4f plane;
        if (!math::esti_plane<float>(
                plane, neighbors,
                static_cast<float>(options_.plane_fit_threshold)) ||
            !plane.allFinite()) {
            continue;
        }

        const Vec3d normal = plane.head<3>().cast<double>();
        const double residual = normal.dot(transformed) + plane.w();
        if (!normal.allFinite() || !std::isfinite(residual) ||
            std::abs(residual) > max_correspondence_distance) {
            continue;
        }

        H6d jacobian;
        jacobian.head<3>() = normal.transpose();
        jacobian.tail<3>() = transformed.cross(normal).transpose();
        if (!jacobian.allFinite()) continue;

        const double absolute_residual = std::abs(residual);
        const double weight = absolute_residual <= options_.huber_delta
                                  ? 1.0
                                  : options_.huber_delta / absolute_residual;
        linearization.hessian.noalias() +=
            weight * jacobian.transpose() * jacobian;
        linearization.gradient.noalias() +=
            weight * jacobian.transpose() * residual;
        linearization.squared_error += residual * residual;
        ++linearization.matches;
    }

    linearization.inlier_ratio =
        static_cast<double>(linearization.matches) /
        static_cast<double>(source->size());
    if (linearization.matches < options_.min_matches ||
        linearization.inlier_ratio < options_.min_inlier_ratio ||
        !linearization.hessian.allFinite() ||
        !linearization.gradient.allFinite()) {
        return linearization;
    }

    Vec6d inverse_scale;
    for (int index = 0; index < 6; ++index) {
        const double diagonal = linearization.hessian(index, index);
        if (!std::isfinite(diagonal) || diagonal <= 0.0) return linearization;
        inverse_scale[index] = 1.0 / std::sqrt(diagonal);
    }
    const Mat6d normalized_hessian =
        inverse_scale.asDiagonal() * linearization.hessian *
        inverse_scale.asDiagonal();
    Eigen::SelfAdjointEigenSolver<Mat6d> eigen_solver(normalized_hessian);
    if (eigen_solver.info() != Eigen::Success ||
        !eigen_solver.eigenvalues().allFinite()) {
        return linearization;
    }
    const double minimum_eigenvalue = eigen_solver.eigenvalues().minCoeff();
    const double maximum_eigenvalue = eigen_solver.eigenvalues().maxCoeff();
    if (minimum_eigenvalue <= 0.0 || maximum_eigenvalue <= 0.0) {
        return linearization;
    }
    linearization.condition_number = maximum_eigenvalue / minimum_eigenvalue;
    linearization.observable =
        minimum_eigenvalue >= options_.min_normalized_hessian_eigenvalue &&
        linearization.condition_number <= options_.max_hessian_condition_number;
    return linearization;
}

bool PointToPlaneRegistration::WithinCorrectionBounds(
    const SE3& initial_pose, const SE3& pose, Result& result) const {
    if (!pose.matrix().allFinite()) return false;
    result.translation_correction =
        (pose.translation() - initial_pose.translation()).norm();
    result.rotation_correction_deg =
        (initial_pose.so3().inverse() * pose.so3()).log().norm() *
        kDegreesPerRadian;
    return std::isfinite(result.translation_correction) &&
           std::isfinite(result.rotation_correction_deg) &&
           result.translation_correction <=
               options_.max_translation_correction &&
           result.rotation_correction_deg <=
               options_.max_rotation_correction_deg;
}

bool PointToPlaneRegistration::RunStage(
    const CloudPtr& source, double max_correspondence_distance,
    int max_iterations, const SE3& initial_pose, SE3& pose,
    Result& result) const {
    const double rotation_convergence =
        options_.rotation_convergence_deg / kDegreesPerRadian;
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        const Linearization linearization =
            Linearize(source, pose, max_correspondence_distance);
        result.matches = linearization.matches;
        result.inlier_ratio = linearization.inlier_ratio;
        result.hessian_condition_number = linearization.condition_number;
        result.rmse = linearization.matches > 0
                          ? std::sqrt(linearization.squared_error /
                                      linearization.matches)
                          : std::numeric_limits<double>::infinity();
        if (!linearization.observable || !std::isfinite(result.rmse)) {
            return false;
        }

        const Eigen::LDLT<Mat6d> solver(linearization.hessian);
        if (solver.info() != Eigen::Success) return false;
        const Vec6d increment = solver.solve(-linearization.gradient);
        if (solver.info() != Eigen::Success || !increment.allFinite()) {
            return false;
        }

        const SE3 proposed_pose = SE3::exp(increment) * pose;
        ++result.iterations;
        if (proposed_pose.matrix().allFinite()) {
            result.attempted_pose = proposed_pose;
            result.attempted_pose_valid = true;
        }
        if (!WithinCorrectionBounds(initial_pose, proposed_pose, result)) {
            return false;
        }
        pose = proposed_pose;
        if (increment.head<3>().norm() <= options_.translation_convergence &&
            increment.tail<3>().norm() <= rotation_convergence) {
            return true;
        }
    }
    return false;
}

PointToPlaneRegistration::Result PointToPlaneRegistration::Refine(
    const CloudPtr& source, SE3& pose) const {
    Result result;
    if (!HasTarget() || !source || source->empty() ||
        !pose.matrix().allFinite()) {
        return result;
    }
    const CloudPtr downsampled_source =
        DownsampleFinite(source, options_.source_voxel_size);
    if (downsampled_source->size() <
        static_cast<std::size_t>(options_.min_matches)) {
        return result;
    }

    const SE3 initial_pose = pose;
    SE3 refined_pose = initial_pose;
    if (!RunStage(downsampled_source,
                  options_.coarse_max_correspondence_distance,
                  options_.coarse_max_iterations, initial_pose, refined_pose,
                  result) ||
        !RunStage(downsampled_source,
                  options_.fine_max_correspondence_distance,
                  options_.fine_max_iterations, initial_pose, refined_pose,
                  result)) {
        return result;
    }

    const Linearization final_linearization =
        Linearize(downsampled_source, refined_pose,
                  options_.fine_max_correspondence_distance);
    result.matches = final_linearization.matches;
    result.inlier_ratio = final_linearization.inlier_ratio;
    result.hessian_condition_number = final_linearization.condition_number;
    result.rmse = final_linearization.matches > 0
                      ? std::sqrt(final_linearization.squared_error /
                                  final_linearization.matches)
                      : std::numeric_limits<double>::infinity();
    result.converged = final_linearization.observable &&
                       std::isfinite(result.rmse) &&
                       WithinCorrectionBounds(initial_pose, refined_pose, result);
    if (!result.converged) return result;

    result.attempted_pose = refined_pose;
    result.attempted_pose_valid = true;
    result.success = result.rmse <= options_.max_rmse;
    if (!result.success) return result;

    pose = refined_pose;
    return result;
}

}  // namespace lightning::loc
