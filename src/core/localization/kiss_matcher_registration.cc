#include "core/localization/kiss_matcher_registration.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <set>
#include <utility>

#include <Eigen/SVD>
#include <pcl/features/fpfh_omp.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/search/kdtree.h>

namespace lightning::loc {
namespace {

using Correspondence = KissMatcherRegistration::Correspondence;
using Correspondences =
    std::vector<Correspondence, Eigen::aligned_allocator<Correspondence>>;

void SetError(std::string* error, const std::string& message) {
    if (error) *error = message;
}

bool DescriptorFinite(const pcl::FPFHSignature33& descriptor) {
    for (const float value : descriptor.histogram) {
        if (!std::isfinite(value)) return false;
    }
    return true;
}

bool EstimateTransform(const Correspondences& correspondences,
                       const std::vector<std::size_t>& indices,
                       const std::vector<double>& weights, SE3& transform) {
    if (indices.size() < 3 || weights.size() != indices.size()) return false;
    double weight_sum = 0.0;
    std::size_t active_weights = 0;
    Vec3d source_mean = Vec3d::Zero();
    Vec3d target_mean = Vec3d::Zero();
    for (std::size_t i = 0; i < indices.size(); ++i) {
        const double weight = weights[i];
        if (!std::isfinite(weight) || weight <= 0.0) continue;
        const auto& correspondence = correspondences[indices[i]];
        ++active_weights;
        weight_sum += weight;
        source_mean += weight * correspondence.source;
        target_mean += weight * correspondence.target;
    }
    // GNC weights are fractional.  Requiring their sum to be at least three
    // incorrectly rejects a valid three-or-more-point weighted estimate when
    // all current weights are small.
    if (active_weights < 3 || !std::isfinite(weight_sum) ||
        weight_sum <= std::numeric_limits<double>::epsilon()) {
        return false;
    }
    source_mean /= weight_sum;
    target_mean /= weight_sum;

    Mat3d covariance = Mat3d::Zero();
    for (std::size_t i = 0; i < indices.size(); ++i) {
        const double weight = weights[i];
        if (!std::isfinite(weight) || weight <= 0.0) continue;
        const auto& correspondence = correspondences[indices[i]];
        covariance += weight * (correspondence.source - source_mean) *
                      (correspondence.target - target_mean).transpose();
    }
    const Eigen::JacobiSVD<Mat3d> svd(
        covariance, Eigen::ComputeFullU | Eigen::ComputeFullV);
    if (svd.info() != Eigen::Success || !svd.singularValues().allFinite() ||
        svd.singularValues()[1] < 1e-9) {
        return false;
    }
    Mat3d correction = Mat3d::Identity();
    correction(2, 2) =
        (svd.matrixV() * svd.matrixU().transpose()).determinant() < 0.0
            ? -1.0
            : 1.0;
    const Mat3d rotation =
        svd.matrixV() * correction * svd.matrixU().transpose();
    const Vec3d translation = target_mean - rotation * source_mean;
    if (!rotation.allFinite() || !translation.allFinite()) return false;
    transform = SE3(Quatd(rotation).normalized(), translation);
    return true;
}

double SquaredResidual(const Correspondence& correspondence,
                       const SE3& transform) {
    return (transform * correspondence.source - correspondence.target)
        .squaredNorm();
}

std::vector<std::size_t> MaximumCore(
    const std::vector<std::vector<int>>& adjacency, int& maximum_core_degree) {
    const std::size_t count = adjacency.size();
    std::vector<int> degrees(count, 0);
    std::vector<int> core_numbers(count, 0);
    std::vector<bool> removed(count, false);
    using QueueEntry = std::pair<int, int>;
    std::priority_queue<QueueEntry, std::vector<QueueEntry>,
                        std::greater<QueueEntry>> queue;
    for (std::size_t i = 0; i < count; ++i) {
        degrees[i] = static_cast<int>(adjacency[i].size());
        queue.emplace(degrees[i], static_cast<int>(i));
    }

    maximum_core_degree = 0;
    while (!queue.empty()) {
        const auto [degree, vertex] = queue.top();
        queue.pop();
        if (removed[vertex] || degrees[vertex] != degree) continue;
        removed[vertex] = true;
        core_numbers[vertex] = degree;
        maximum_core_degree = std::max(maximum_core_degree, degree);
        for (const int neighbor : adjacency[vertex]) {
            if (removed[neighbor] || degrees[neighbor] <= degree) continue;
            --degrees[neighbor];
            queue.emplace(degrees[neighbor], neighbor);
        }
    }

    std::vector<std::size_t> selected;
    for (std::size_t i = 0; i < count; ++i) {
        if (core_numbers[i] >= maximum_core_degree) selected.push_back(i);
    }
    return selected;
}

}  // namespace

KissMatcherRegistration::KissMatcherRegistration() = default;

KissMatcherRegistration::KissMatcherRegistration(Options options)
    : options_(std::move(options)) {}

bool KissMatcherRegistration::ValidateOptions(const Options& options,
                                              std::string* error) {
    if (!std::isfinite(options.voxel_size) || options.voxel_size <= 0.0 ||
        !std::isfinite(options.normal_radius) ||
        options.normal_radius <= options.voxel_size ||
        !std::isfinite(options.feature_radius) ||
        options.feature_radius <= options.normal_radius) {
        SetError(error,
                 "kiss_matcher radii must satisfy 0 < voxel_size < normal_radius < feature_radius");
        return false;
    }
    if (!std::isfinite(options.descriptor_ratio) ||
        options.descriptor_ratio <= 0.0 || options.descriptor_ratio >= 1.0 ||
        options.max_feature_points < 100 || options.max_correspondences < 3) {
        SetError(error, "invalid kiss_matcher feature matching options");
        return false;
    }
    if (!std::isfinite(options.compatibility_tolerance) ||
        options.compatibility_tolerance <= 0.0 ||
        !std::isfinite(options.minimum_pair_distance) ||
        options.minimum_pair_distance < 0.0 || options.minimum_core_degree < 1 ||
        options.minimum_inliers < 3) {
        SetError(error, "invalid kiss_matcher compatibility graph options");
        return false;
    }
    if (!std::isfinite(options.inlier_threshold) ||
        options.inlier_threshold <= 0.0 ||
        !std::isfinite(options.maximum_rmse) || options.maximum_rmse <= 0.0 ||
        options.gnc_max_iterations <= 0 || !std::isfinite(options.gnc_factor) ||
        options.gnc_factor <= 1.0) {
        SetError(error, "invalid kiss_matcher robust solver options");
        return false;
    }
    if (!std::isfinite(options.spatial_cell_size) ||
        options.spatial_cell_size <= 0.0 || options.minimum_spatial_cells <= 0 ||
        options.feature_threads <= 0) {
        SetError(error, "invalid kiss_matcher coverage or threading options");
        return false;
    }
    return true;
}

bool KissMatcherRegistration::BuildFeatures(const CloudPtr& input,
                                            PreparedCloud& output,
                                            std::string* error) const {
    output.points->clear();
    output.descriptors->clear();
    if (!input || input->empty()) {
        SetError(error, "input cloud is empty");
        return false;
    }

    CloudPtr voxelized(new PointCloudType);
    pcl::VoxelGrid<PointType> voxel;
    const float leaf = static_cast<float>(options_.voxel_size);
    voxel.setLeafSize(leaf, leaf, leaf);
    voxel.setInputCloud(input);
    voxel.filter(*voxelized);
    if (voxelized->size() < 20) {
        SetError(error, "too few points remain after voxel filtering");
        return false;
    }

    const std::size_t stride = std::max<std::size_t>(
        1, (voxelized->size() +
            static_cast<std::size_t>(options_.max_feature_points) - 1) /
               static_cast<std::size_t>(options_.max_feature_points));
    pcl::PointCloud<pcl::PointXYZ>::Ptr points(
        new pcl::PointCloud<pcl::PointXYZ>);
    points->reserve(std::min<std::size_t>(
        voxelized->size(), static_cast<std::size_t>(options_.max_feature_points)));
    for (std::size_t i = 0; i < voxelized->size(); i += stride) {
        const auto& source = voxelized->points[i];
        if (!std::isfinite(source.x) || !std::isfinite(source.y) ||
            !std::isfinite(source.z)) {
            continue;
        }
        points->emplace_back(source.x, source.y, source.z);
    }
    if (points->size() < 20) {
        SetError(error, "too few finite feature points");
        return false;
    }

    pcl::PointCloud<pcl::Normal>::Ptr normals(
        new pcl::PointCloud<pcl::Normal>);
    pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> normal_estimator;
    normal_estimator.setNumberOfThreads(options_.feature_threads);
    normal_estimator.setInputCloud(points);
    normal_estimator.setSearchMethod(
        pcl::search::KdTree<pcl::PointXYZ>::Ptr(
            new pcl::search::KdTree<pcl::PointXYZ>));
    normal_estimator.setRadiusSearch(options_.normal_radius);
    normal_estimator.compute(*normals);

    pcl::PointCloud<pcl::FPFHSignature33>::Ptr descriptors(
        new pcl::PointCloud<pcl::FPFHSignature33>);
    pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal,
                           pcl::FPFHSignature33>
        feature_estimator;
    feature_estimator.setNumberOfThreads(options_.feature_threads);
    feature_estimator.setInputCloud(points);
    feature_estimator.setInputNormals(normals);
    feature_estimator.setSearchMethod(
        pcl::search::KdTree<pcl::PointXYZ>::Ptr(
            new pcl::search::KdTree<pcl::PointXYZ>));
    feature_estimator.setRadiusSearch(options_.feature_radius);
    feature_estimator.compute(*descriptors);

    output.points->reserve(points->size());
    output.descriptors->reserve(descriptors->size());
    for (std::size_t i = 0; i < points->size() && i < descriptors->size() &&
                            i < normals->size();
         ++i) {
        const auto& normal = normals->points[i];
        if (!pcl::isFinite(normal) || !DescriptorFinite(descriptors->points[i])) {
            continue;
        }
        output.points->push_back(points->points[i]);
        output.descriptors->push_back(descriptors->points[i]);
    }
    output.points->width = static_cast<std::uint32_t>(output.points->size());
    output.points->height = 1;
    output.descriptors->width =
        static_cast<std::uint32_t>(output.descriptors->size());
    output.descriptors->height = 1;
    if (output.points->size() < 20) {
        SetError(error, "too few valid FPFH descriptors");
        return false;
    }
    return true;
}

bool KissMatcherRegistration::SetTarget(const CloudPtr& target,
                                        std::string* error) {
    target_points_.reset();
    target_features_.reset();
    target_feature_tree_.reset();
    if (!ValidateOptions(options_, error)) return false;

    PreparedCloud features;
    if (!BuildFeatures(target, features, error)) return false;
    target_points_ = std::move(features.points);
    target_features_ = std::move(features.descriptors);
    target_feature_tree_.reset(new pcl::KdTreeFLANN<pcl::FPFHSignature33>);
    target_feature_tree_->setInputCloud(target_features_);
    return true;
}

bool KissMatcherRegistration::PrepareSource(
    const CloudPtr& source, PreparedCloud& prepared, std::string* error) const {
    if (!ValidateOptions(options_, error)) return false;
    return BuildFeatures(source, prepared, error);
}

std::vector<Correspondence, Eigen::aligned_allocator<Correspondence>>
KissMatcherRegistration::MatchFeatures(const PreparedCloud& source) const {
    Correspondences correspondences;
    if (!IsReady() || !source.descriptors || source.descriptors->empty()) {
        return correspondences;
    }

    pcl::KdTreeFLANN<pcl::FPFHSignature33> source_tree;
    source_tree.setInputCloud(source.descriptors);
    const double maximum_ratio_squared =
        options_.descriptor_ratio * options_.descriptor_ratio;
    std::vector<int> target_indices(2);
    std::vector<float> target_distances(2);
    std::vector<int> source_indices(1);
    std::vector<float> source_distances(1);
    correspondences.reserve(source.descriptors->size());
    for (std::size_t source_index = 0;
         source_index < source.descriptors->size(); ++source_index) {
        if (target_feature_tree_->nearestKSearch(
                source.descriptors->points[source_index], 2, target_indices,
                target_distances) != 2 ||
            target_distances[1] <= std::numeric_limits<float>::epsilon() ||
            target_distances[0] >
                maximum_ratio_squared * target_distances[1]) {
            continue;
        }
        const int target_index = target_indices[0];
        if (source_tree.nearestKSearch(target_features_->points[target_index], 1,
                                       source_indices,
                                       source_distances) != 1 ||
            source_indices[0] != static_cast<int>(source_index)) {
            continue;
        }
        const auto& source_point = source.points->points[source_index];
        const auto& target_point = target_points_->points[target_index];
        Correspondence correspondence;
        correspondence.source =
            Vec3d(source_point.x, source_point.y, source_point.z);
        correspondence.target =
            Vec3d(target_point.x, target_point.y, target_point.z);
        correspondence.descriptor_distance =
            std::sqrt(std::max(0.0F, target_distances[0]));
        correspondences.push_back(std::move(correspondence));
    }

    std::stable_sort(correspondences.begin(), correspondences.end(),
                     [](const Correspondence& left,
                        const Correspondence& right) {
                         return left.descriptor_distance <
                                right.descriptor_distance;
                     });
    if (correspondences.size() >
        static_cast<std::size_t>(options_.max_correspondences)) {
        correspondences.resize(
            static_cast<std::size_t>(options_.max_correspondences));
    }
    return correspondences;
}

KissMatcherRegistration::Result KissMatcherRegistration::SolveCorrespondences(
    const Correspondences& correspondences, const Options& options) {
    Result result;
    result.rough_correspondences = correspondences.size();
    std::string validation_error;
    if (!ValidateOptions(options, &validation_error)) {
        result.reason = validation_error;
        return result;
    }
    if (correspondences.size() <
        static_cast<std::size_t>(options.minimum_inliers)) {
        result.reason = "too_few_correspondences";
        return result;
    }

    std::vector<std::vector<int>> adjacency(correspondences.size());
    for (std::size_t left = 0; left < correspondences.size(); ++left) {
        for (std::size_t right = left + 1; right < correspondences.size();
             ++right) {
            const double source_distance =
                (correspondences[left].source - correspondences[right].source)
                    .norm();
            const double target_distance =
                (correspondences[left].target - correspondences[right].target)
                    .norm();
            if (std::min(source_distance, target_distance) <
                    options.minimum_pair_distance ||
                std::abs(source_distance - target_distance) >
                    options.compatibility_tolerance) {
                continue;
            }
            adjacency[left].push_back(static_cast<int>(right));
            adjacency[right].push_back(static_cast<int>(left));
        }
    }

    int maximum_core_degree = 0;
    const std::vector<std::size_t> core_indices =
        MaximumCore(adjacency, maximum_core_degree);
    result.core_correspondences = core_indices.size();
    if (maximum_core_degree < options.minimum_core_degree ||
        core_indices.size() < static_cast<std::size_t>(options.minimum_inliers)) {
        result.reason = "compatibility_core_too_small";
        return result;
    }

    std::vector<double> weights(core_indices.size(), 1.0);
    SE3 transform;
    if (!EstimateTransform(correspondences, core_indices, weights, transform)) {
        result.reason = "degenerate_correspondences";
        return result;
    }

    const double threshold_squared =
        options.inlier_threshold * options.inlier_threshold;
    double maximum_residual_squared = threshold_squared;
    for (const std::size_t index : core_indices) {
        maximum_residual_squared = std::max(
            maximum_residual_squared,
            SquaredResidual(correspondences[index], transform));
    }
    double mu = 1.0;
    if (2.0 * maximum_residual_squared > threshold_squared) {
        mu = threshold_squared /
             (2.0 * maximum_residual_squared - threshold_squared);
    }
    mu = std::max(mu, 1e-6);

    for (int iteration = 0; iteration < options.gnc_max_iterations;
         ++iteration) {
        result.iterations = iteration + 1;
        if (!EstimateTransform(correspondences, core_indices, weights,
                               transform)) {
            result.reason = "robust_pose_estimation_failed";
            return result;
        }
        const double lower =
            (mu / (mu + 1.0)) * threshold_squared;
        const double upper =
            ((mu + 1.0) / mu) * threshold_squared;
        double maximum_weight_change = 0.0;
        std::vector<double> updated_weights(core_indices.size(), 0.0);
        std::size_t active_weights = 0;
        for (std::size_t i = 0; i < core_indices.size(); ++i) {
            const double residual_squared = SquaredResidual(
                correspondences[core_indices[i]], transform);
            double updated_weight = 0.0;
            if (residual_squared <= lower) {
                updated_weight = 1.0;
            } else if (residual_squared < upper) {
                updated_weight =
                    std::sqrt(threshold_squared * mu * (mu + 1.0) /
                              residual_squared) -
                    mu;
                updated_weight = std::clamp(updated_weight, 0.0, 1.0);
            }
            updated_weights[i] = updated_weight;
            if (updated_weight > 0.0) ++active_weights;
            maximum_weight_change =
                std::max(maximum_weight_change,
                         std::abs(updated_weight - weights[i]));
        }
        mu *= options.gnc_factor;
        // At the beginning of GNC, a poor least-squares pose can place every
        // residual outside the current truncated basin.  Increasing mu
        // expands the basin; failing immediately here makes the solver least
        // reliable precisely when the initial pose is hardest.
        if (active_weights < 3) continue;
        weights.swap(updated_weights);
        if (maximum_weight_change < 1e-4 && mu > 1.0) {
            result.converged = true;
            break;
        }
    }

    std::vector<std::size_t> inlier_indices;
    for (const std::size_t index : core_indices) {
        if (SquaredResidual(correspondences[index], transform) <=
            threshold_squared) {
            inlier_indices.push_back(index);
        }
    }
    result.inliers = inlier_indices.size();
    result.inlier_ratio = core_indices.empty()
                              ? 0.0
                              : static_cast<double>(inlier_indices.size()) /
                                    static_cast<double>(core_indices.size());
    if (inlier_indices.size() <
        static_cast<std::size_t>(options.minimum_inliers)) {
        result.reason = "too_few_robust_inliers";
        return result;
    }

    std::vector<double> inlier_weights(inlier_indices.size(), 1.0);
    if (!EstimateTransform(correspondences, inlier_indices, inlier_weights,
                           transform)) {
        result.reason = "inlier_refinement_failed";
        return result;
    }
    double squared_error_sum = 0.0;
    std::set<std::pair<int, int>> inlier_cells;
    std::set<std::pair<int, int>> core_cells;
    for (const std::size_t index : core_indices) {
        const Vec3d& point = correspondences[index].target;
        core_cells.emplace(
            static_cast<int>(std::floor(point.x() / options.spatial_cell_size)),
            static_cast<int>(std::floor(point.y() / options.spatial_cell_size)));
    }
    for (const std::size_t index : inlier_indices) {
        squared_error_sum +=
            SquaredResidual(correspondences[index], transform);
        const Vec3d& point = correspondences[index].target;
        inlier_cells.emplace(
            static_cast<int>(std::floor(point.x() / options.spatial_cell_size)),
            static_cast<int>(std::floor(point.y() / options.spatial_cell_size)));
    }
    result.rmse =
        std::sqrt(squared_error_sum / static_cast<double>(inlier_indices.size()));
    result.spatial_cells = inlier_cells.size();
    result.spatial_coverage =
        core_cells.empty()
            ? 0.0
            : static_cast<double>(inlier_cells.size()) /
                  static_cast<double>(core_cells.size());
    result.T_target_source = transform;
    result.converged = result.converged || result.iterations > 0;
    result.score = result.inlier_ratio *
                   std::exp(-result.rmse / options.inlier_threshold) *
                   std::min(1.0, static_cast<double>(result.spatial_cells) /
                                     static_cast<double>(
                                         options.minimum_spatial_cells));
    if (!std::isfinite(result.rmse) || result.rmse > options.maximum_rmse) {
        result.reason = "rmse_above_threshold";
        return result;
    }
    if (result.spatial_cells <
        static_cast<std::size_t>(options.minimum_spatial_cells)) {
        result.reason = "insufficient_spatial_coverage";
        return result;
    }
    result.success = true;
    result.reason = "accepted";
    return result;
}

KissMatcherRegistration::Result KissMatcherRegistration::Align(
    const CloudPtr& source) const {
    Result result;
    if (!IsReady()) {
        result.reason = "target_not_ready";
        return result;
    }
    PreparedCloud features;
    std::string feature_error;
    if (!PrepareSource(source, features, &feature_error)) {
        result.reason = feature_error;
        return result;
    }
    return Align(features);
}

KissMatcherRegistration::Result KissMatcherRegistration::Align(
    const PreparedCloud& source) const {
    Result result;
    if (!IsReady()) {
        result.reason = "target_not_ready";
        return result;
    }
    if (!source.points || !source.descriptors || source.points->empty() ||
        source.points->size() != source.descriptors->size()) {
        result.reason = "prepared_source_invalid";
        return result;
    }
    result.source_points = source.points->size();
    result.source_descriptors = source.descriptors->size();
    const Correspondences correspondences = MatchFeatures(source);
    Result solved = SolveCorrespondences(correspondences, options_);
    solved.source_points = result.source_points;
    solved.source_descriptors = result.source_descriptors;
    return solved;
}

}  // namespace lightning::loc
