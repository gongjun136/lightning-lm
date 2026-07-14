#include "core/backend/voxel_bundle_adjustment.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <utility>

#include <Eigen/Eigenvalues>

namespace lightning::backend {
namespace {

constexpr double kRadToDeg = 180.0 / M_PI;

Eigen::Matrix3d Hat(const Eigen::Vector3d& value) {
    Eigen::Matrix3d result;
    result << 0.0, -value.z(), value.y(), value.z(), 0.0, -value.x(), -value.y(), value.x(), 0.0;
    return result;
}

Eigen::Matrix3d ExpSO3(const Eigen::Vector3d& delta) {
    const double angle = delta.norm();
    if (angle < 1e-12) return Eigen::Matrix3d::Identity() + Hat(delta);
    const Eigen::Matrix3d axis_hat = Hat(delta / angle);
    return Eigen::Matrix3d::Identity() + std::sin(angle) * axis_hat +
           (1.0 - std::cos(angle)) * axis_hat * axis_hat;
}

Eigen::Vector3d LogSO3(const Eigen::Matrix3d& rotation) {
    Eigen::AngleAxisd angle_axis(rotation);
    if (!std::isfinite(angle_axis.angle()) || angle_axis.angle() < 1e-12) return Eigen::Vector3d::Zero();
    return angle_axis.axis() * angle_axis.angle();
}

struct PoseState {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();
    Eigen::Vector3d translation = Eigen::Vector3d::Zero();
};

struct PointCluster {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Eigen::Matrix3d second_moment = Eigen::Matrix3d::Zero();
    Eigen::Vector3d sum = Eigen::Vector3d::Zero();
    int count = 0;

    void Add(const Eigen::Vector3d& point) {
        second_moment.noalias() += point * point.transpose();
        sum += point;
        ++count;
    }

    void Add(const PointCluster& other) {
        second_moment += other.second_moment;
        sum += other.sum;
        count += other.count;
    }

    PointCluster Transformed(const PoseState& pose) const {
        PointCluster result;
        result.count = count;
        result.sum = pose.rotation * sum + static_cast<double>(count) * pose.translation;
        const Eigen::Matrix3d cross = pose.rotation * sum * pose.translation.transpose();
        result.second_moment = pose.rotation * second_moment * pose.rotation.transpose() + cross +
                               cross.transpose() + static_cast<double>(count) * pose.translation *
                                                       pose.translation.transpose();
        return result;
    }

    Eigen::Matrix3d Covariance() const {
        if (count <= 0) return Eigen::Matrix3d::Zero();
        const Eigen::Vector3d center = sum / static_cast<double>(count);
        return second_moment / static_cast<double>(count) - center * center.transpose();
    }
};

struct PointRecord {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    int frame = 0;
    Eigen::Vector3d local = Eigen::Vector3d::Zero();
};

struct VoxelKey {
    int64_t x = 0;
    int64_t y = 0;
    int64_t z = 0;

    bool operator==(const VoxelKey& other) const { return x == other.x && y == other.y && z == other.z; }
};

struct VoxelKeyHash {
    std::size_t operator()(const VoxelKey& key) const {
        std::size_t seed = std::hash<int64_t>{}(key.x);
        seed ^= std::hash<int64_t>{}(key.y) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<int64_t>{}(key.z) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

VoxelKey ToVoxel(const Eigen::Vector3d& point, double size) {
    return {static_cast<int64_t>(std::floor(point.x() / size)),
            static_cast<int64_t>(std::floor(point.y() / size)),
            static_cast<int64_t>(std::floor(point.z() / size))};
}

struct PlaneFactor {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    std::vector<PointCluster, Eigen::aligned_allocator<PointCluster>> local_clusters;
};

using PoseStates = std::vector<PoseState, Eigen::aligned_allocator<PoseState>>;
using PointRecords = std::vector<PointRecord, Eigen::aligned_allocator<PointRecord>>;

struct FactorBuildResult {
    std::vector<PlaneFactor, Eigen::aligned_allocator<PlaneFactor>> factors;
    std::size_t sampled_points = 0;
};

void ProcessBuckets(const PointRecords& records, const std::vector<std::size_t>& indices, const PoseStates& poses,
                    int frame_count, double voxel_size, int depth, const BundleAdjustmentOptions& options,
                    std::vector<PlaneFactor, Eigen::aligned_allocator<PlaneFactor>>& factors) {
    std::unordered_map<VoxelKey, std::vector<std::size_t>, VoxelKeyHash> buckets;
    buckets.reserve(indices.size());
    for (const std::size_t index : indices) {
        const auto& record = records[index];
        const Eigen::Vector3d world = poses[record.frame].rotation * record.local + poses[record.frame].translation;
        buckets[ToVoxel(world, voxel_size)].push_back(index);
    }

    for (auto& [key, bucket] : buckets) {
        (void)key;
        if (bucket.size() < static_cast<std::size_t>(options.min_points_per_voxel)) continue;

        PlaneFactor factor;
        factor.local_clusters.resize(frame_count);
        PointCluster world_cluster;
        int contributing_frames = 0;
        for (const std::size_t index : bucket) {
            const auto& record = records[index];
            factor.local_clusters[record.frame].Add(record.local);
            world_cluster.Add(poses[record.frame].rotation * record.local + poses[record.frame].translation);
        }
        for (const auto& cluster : factor.local_clusters) contributing_frames += cluster.count > 0 ? 1 : 0;

        bool planar = false;
        if (contributing_frames >= options.min_frames_per_voxel) {
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(world_cluster.Covariance());
            if (eigen_solver.info() == Eigen::Success) {
                const Eigen::Vector3d values = eigen_solver.eigenvalues();
                const double denominator = std::max(values.y(), options.eigenvalue_floor);
                planar = values.x() >= -options.eigenvalue_floor && values.x() / denominator < options.planarity_ratio;
            }
        }

        if (planar) {
            factors.push_back(std::move(factor));
        } else if (depth + 1 < options.max_voxel_depth && voxel_size > 1e-3) {
            ProcessBuckets(records, bucket, poses, frame_count, voxel_size * 0.5, depth + 1, options, factors);
        }
    }
}

FactorBuildResult BuildFactors(const std::vector<BundleFrame>& frames, const PoseStates& poses,
                               const BundleAdjustmentOptions& options) {
    PointRecords records;
    for (std::size_t frame_index = 0; frame_index < frames.size(); ++frame_index) {
        const auto& cloud = frames[frame_index].cloud;
        if (!cloud || cloud->empty()) continue;
        const std::size_t maximum = std::max(1, options.max_points_per_frame);
        const std::size_t stride = std::max<std::size_t>(1, (cloud->size() + maximum - 1) / maximum);
        records.reserve(records.size() + std::min(cloud->size(), maximum));
        for (std::size_t point_index = 0; point_index < cloud->size(); point_index += stride) {
            const auto& point = cloud->points[point_index];
            if (!std::isfinite(point.x) || !std::isfinite(point.y) || !std::isfinite(point.z)) continue;
            PointRecord record;
            record.frame = static_cast<int>(frame_index);
            record.local = Eigen::Vector3d(point.x, point.y, point.z);
            records.push_back(record);
        }
    }

    FactorBuildResult result;
    result.sampled_points = records.size();
    if (records.empty()) return result;
    std::vector<std::size_t> indices(records.size());
    for (std::size_t index = 0; index < indices.size(); ++index) indices[index] = index;
    ProcessBuckets(records, indices, poses, static_cast<int>(frames.size()), options.voxel_size, 0, options,
                   result.factors);
    return result;
}

double EvaluateFactors(const std::vector<PlaneFactor, Eigen::aligned_allocator<PlaneFactor>>& factors,
                       const PoseStates& poses, double eigenvalue_floor, Eigen::MatrixXd* hessian,
                       Eigen::VectorXd* gradient) {
    const int frame_count = static_cast<int>(poses.size());
    const int dimension = frame_count * 6;
    if (hessian) hessian->setZero(dimension, dimension);
    if (gradient) gradient->setZero(dimension);
    double cost = 0.0;

    for (const auto& factor : factors) {
        PointCluster total;
        for (int frame = 0; frame < frame_count; ++frame) {
            if (factor.local_clusters[frame].count > 0) total.Add(factor.local_clusters[frame].Transformed(poses[frame]));
        }
        if (total.count <= 0) continue;

        const Eigen::Vector3d center = total.sum / static_cast<double>(total.count);
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(total.Covariance());
        if (eigen_solver.info() != Eigen::Success) continue;
        const Eigen::Vector3d eigenvalues = eigen_solver.eigenvalues();
        const Eigen::Matrix3d eigenvectors = eigen_solver.eigenvectors();
        if (!eigenvalues.allFinite() || !eigenvectors.allFinite()) continue;
        cost += std::max(0.0, eigenvalues.x());
        if (!hessian || !gradient) continue;

        const Eigen::Vector3d normal = eigenvectors.col(0);
        const Eigen::Matrix3d normal_outer = normal * normal.transpose();
        Eigen::Matrix3d eigen_coupling = Eigen::Matrix3d::Zero();
        for (int eigen_index = 1; eigen_index < 3; ++eigen_index) {
            double denominator = eigenvalues.x() - eigenvalues[eigen_index];
            if (std::abs(denominator) < eigenvalue_floor) denominator = -eigenvalue_floor;
            eigen_coupling += 2.0 / denominator * eigenvectors.col(eigen_index) *
                              eigenvectors.col(eigen_index).transpose();
        }

        std::vector<Eigen::Matrix<double, 3, 6>,
                    Eigen::aligned_allocator<Eigen::Matrix<double, 3, 6>>>
            pose_jacobians(frame_count, Eigen::Matrix<double, 3, 6>::Zero());
        std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> skew_sum_rotated(frame_count);
        std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>>
            skew_sum_rotated_normal(frame_count);

        for (int frame = 0; frame < frame_count; ++frame) {
            const PointCluster& local = factor.local_clusters[frame];
            if (local.count <= 0) continue;
            const Eigen::Matrix3d& rotation = poses[frame].rotation;
            const Eigen::Matrix3d sum_hat = Hat(local.sum);
            const Eigen::Vector3d rotated_normal = rotation.transpose() * normal;
            const Eigen::Matrix3d rotated_normal_hat = Hat(rotated_normal);
            const Eigen::Vector3d second_rotated_normal = local.second_moment * rotated_normal;
            skew_sum_rotated[frame] = sum_hat * rotated_normal;
            skew_sum_rotated_normal[frame] = skew_sum_rotated[frame] * normal.transpose();

            const Eigen::Vector3d translation_from_center = poses[frame].translation - center;
            const double normal_translation = normal.dot(translation_from_center);
            const Eigen::Matrix3d combo_rotation = Hat(second_rotated_normal) + sum_hat * normal_translation;
            const Eigen::Vector3d combo_translation =
                rotation * local.sum + static_cast<double>(local.count) * translation_from_center;

            auto& pose_jacobian = pose_jacobians[frame];
            pose_jacobian.block<3, 3>(0, 0) =
                (rotation * local.second_moment + translation_from_center * local.sum.transpose()) *
                    rotated_normal_hat -
                rotation * combo_rotation;
            pose_jacobian.block<3, 3>(0, 3) =
                combo_translation * normal.transpose() + combo_translation.dot(normal) * Eigen::Matrix3d::Identity();
            pose_jacobian /= static_cast<double>(total.count);

            const Eigen::Matrix<double, 6, 1> local_gradient = pose_jacobian.transpose() * normal;
            gradient->segment<6>(6 * frame) += local_gradient;

            const Eigen::Matrix3d rotation_translation =
                2.0 / static_cast<double>(total.count) *
                (1.0 - static_cast<double>(local.count) / static_cast<double>(total.count)) *
                skew_sum_rotated_normal[frame];
            Eigen::Matrix<double, 6, 6> local_hessian =
                pose_jacobian.transpose() * eigen_coupling * pose_jacobian;
            local_hessian.block<3, 3>(0, 0) +=
                2.0 / static_cast<double>(total.count) *
                    (combo_rotation - rotated_normal_hat * local.second_moment) * rotated_normal_hat -
                2.0 / (static_cast<double>(total.count) * static_cast<double>(total.count)) *
                    skew_sum_rotated[frame] * skew_sum_rotated[frame].transpose() -
                0.5 * Hat(local_gradient.head<3>());
            local_hessian.block<3, 3>(0, 3) += rotation_translation;
            local_hessian.block<3, 3>(3, 0) += rotation_translation.transpose();
            local_hessian.block<3, 3>(3, 3) +=
                2.0 / static_cast<double>(total.count) *
                (static_cast<double>(local.count) -
                 static_cast<double>(local.count * local.count) / static_cast<double>(total.count)) *
                normal_outer;
            hessian->block<6, 6>(6 * frame, 6 * frame) +=
                0.5 * (local_hessian + local_hessian.transpose());
        }

        for (int first = 0; first + 1 < frame_count; ++first) {
            const auto& first_cluster = factor.local_clusters[first];
            if (first_cluster.count <= 0) continue;
            for (int second = first + 1; second < frame_count; ++second) {
                const auto& second_cluster = factor.local_clusters[second];
                if (second_cluster.count <= 0) continue;
                Eigen::Matrix<double, 6, 6> cross =
                    pose_jacobians[first].transpose() * eigen_coupling * pose_jacobians[second];
                const double count = static_cast<double>(total.count);
                cross.block<3, 3>(0, 0) +=
                    -2.0 / (count * count) *
                    skew_sum_rotated[first] * skew_sum_rotated[second].transpose();
                cross.block<3, 3>(0, 3) +=
                    -2.0 * static_cast<double>(second_cluster.count) / (count * count) *
                    skew_sum_rotated_normal[first];
                cross.block<3, 3>(3, 0) +=
                    -2.0 * static_cast<double>(first_cluster.count) / (count * count) *
                    skew_sum_rotated_normal[second].transpose();
                cross.block<3, 3>(3, 3) +=
                    -2.0 * static_cast<double>(first_cluster.count * second_cluster.count) / (count * count) *
                    normal_outer;
                hessian->block<6, 6>(6 * first, 6 * second) += cross;
            }
        }
    }

    if (hessian) *hessian = hessian->selfadjointView<Eigen::Upper>();
    return cost;
}

Eigen::Matrix<double, 6, 1> Correction(const PoseState& initial, const PoseState& current) {
    Eigen::Matrix<double, 6, 1> correction;
    correction.head<3>() = LogSO3(initial.rotation.transpose() * current.rotation);
    correction.tail<3>() = current.translation - initial.translation;
    return correction;
}

double AddSmoothness(const PoseStates& initial, const PoseStates& current, const BundleAdjustmentOptions& options,
                     Eigen::MatrixXd* hessian, Eigen::VectorXd* gradient) {
    double cost = 0.0;
    Eigen::Matrix<double, 6, 6> prior_information = Eigen::Matrix<double, 6, 6>::Zero();
    prior_information.diagonal().head<3>().setConstant(options.correction_prior_rotation);
    prior_information.diagonal().tail<3>().setConstant(options.correction_prior_translation);
    for (std::size_t index = 1; index < current.size(); ++index) {
        const Eigen::Matrix<double, 6, 1> residual = Correction(initial[index], current[index]);
        cost += 0.5 * residual.dot(prior_information * residual);
        if (!hessian || !gradient) continue;
        const int offset = static_cast<int>(index * 6);
        hessian->block<6, 6>(offset, offset) += prior_information;
        gradient->segment<6>(offset) += prior_information * residual;
    }

    Eigen::Matrix<double, 6, 6> information = Eigen::Matrix<double, 6, 6>::Zero();
    information.diagonal().head<3>().setConstant(options.correction_smoothness_rotation);
    information.diagonal().tail<3>().setConstant(options.correction_smoothness_translation);
    for (std::size_t index = 1; index < current.size(); ++index) {
        const Eigen::Matrix<double, 6, 1> residual =
            Correction(initial[index], current[index]) - Correction(initial[index - 1], current[index - 1]);
        cost += 0.5 * residual.dot(information * residual);
        if (!hessian || !gradient) continue;
        const int previous_offset = static_cast<int>((index - 1) * 6);
        const int current_offset = static_cast<int>(index * 6);
        hessian->block<6, 6>(previous_offset, previous_offset) += information;
        hessian->block<6, 6>(current_offset, current_offset) += information;
        hessian->block<6, 6>(previous_offset, current_offset) -= information;
        hessian->block<6, 6>(current_offset, previous_offset) -= information;
        gradient->segment<6>(previous_offset) -= information * residual;
        gradient->segment<6>(current_offset) += information * residual;
    }
    return cost;
}

double Evaluate(const std::vector<PlaneFactor, Eigen::aligned_allocator<PlaneFactor>>& factors,
                const PoseStates& initial, const PoseStates& current, const BundleAdjustmentOptions& options,
                Eigen::MatrixXd* hessian, Eigen::VectorXd* gradient) {
    double cost = EvaluateFactors(factors, current, options.eigenvalue_floor, hessian, gradient);
    cost += AddSmoothness(initial, current, options, hessian, gradient);
    return cost;
}

bool WithinTotalCorrectionLimit(const PoseStates& initial, const PoseStates& candidate,
                                const BundleAdjustmentOptions& options) {
    for (std::size_t index = 0; index < candidate.size(); ++index) {
        const auto correction = Correction(initial[index], candidate[index]);
        if (correction.head<3>().norm() * kRadToDeg > options.max_total_rotation_deg ||
            correction.tail<3>().norm() > options.max_total_translation) {
            return false;
        }
    }
    return true;
}

}  // namespace

VoxelBundleAdjuster::VoxelBundleAdjuster(BundleAdjustmentOptions options) : options_(std::move(options)) {}

BundleAdjustmentSummary VoxelBundleAdjuster::Optimize(std::vector<BundleFrame>& frames) const {
    BundleAdjustmentSummary summary;
    summary.frame_count = frames.size();
    if (!options_.enabled) {
        summary.reason = "disabled";
        return summary;
    }
    if (frames.size() < 2) {
        summary.reason = "too_few_frames";
        return summary;
    }
    summary.attempted = true;

    PoseStates initial(frames.size());
    PoseStates prior(frames.size());
    for (std::size_t index = 0; index < frames.size(); ++index) {
        initial[index].rotation = frames[index].pose.rotationMatrix();
        initial[index].translation = frames[index].pose.translation();
        const SE3& prior_pose = frames[index].has_prior_pose ? frames[index].prior_pose : frames[index].pose;
        prior[index].rotation = prior_pose.rotationMatrix();
        prior[index].translation = prior_pose.translation();
    }
    PoseStates current = initial;
    const FactorBuildResult build = BuildFactors(frames, current, options_);
    summary.factor_count = build.factors.size();
    summary.sampled_point_count = build.sampled_points;
    if (build.factors.empty()) {
        summary.reason = "no_shared_planar_voxels";
        return summary;
    }

    Eigen::MatrixXd hessian;
    Eigen::VectorXd gradient;
    double cost = Evaluate(build.factors, prior, current, options_, &hessian, &gradient);
    summary.initial_cost = cost;
    if (!std::isfinite(cost)) {
        summary.reason = "non_finite_initial_cost";
        return summary;
    }

    double damping = std::max(1e-9, options_.initial_damping);
    bool accepted_any = false;
    for (int iteration = 0; iteration < options_.max_iterations; ++iteration) {
        const int dimension = static_cast<int>(frames.size() * 6);
        hessian.topRows(6).setZero();
        hessian.leftCols(6).setZero();
        hessian.block<6, 6>(0, 0).setIdentity();
        gradient.head<6>().setZero();

        Eigen::MatrixXd diagonal = Eigen::MatrixXd::Zero(dimension, dimension);
        diagonal.diagonal() = hessian.diagonal().cwiseAbs().cwiseMax(1e-9);
        Eigen::LDLT<Eigen::MatrixXd> solver(hessian + damping * diagonal);
        if (solver.info() != Eigen::Success) {
            damping *= 10.0;
            continue;
        }
        Eigen::VectorXd step = solver.solve(-gradient);
        if (solver.info() != Eigen::Success || !step.allFinite()) {
            damping *= 10.0;
            continue;
        }

        const double max_rotation_step = options_.max_rotation_step_deg / kRadToDeg;
        for (std::size_t frame = 1; frame < frames.size(); ++frame) {
            auto rotation_step = step.segment<3>(static_cast<int>(6 * frame));
            auto translation_step = step.segment<3>(static_cast<int>(6 * frame + 3));
            const double rotation_norm = rotation_step.norm();
            const double translation_norm = translation_step.norm();
            if (rotation_norm > max_rotation_step) rotation_step *= max_rotation_step / rotation_norm;
            if (translation_norm > options_.max_translation_step)
                translation_step *= options_.max_translation_step / translation_norm;
        }

        PoseStates candidate = current;
        for (std::size_t frame = 1; frame < frames.size(); ++frame) {
            candidate[frame].rotation =
                current[frame].rotation * ExpSO3(step.segment<3>(static_cast<int>(6 * frame)));
            candidate[frame].translation =
                current[frame].translation + step.segment<3>(static_cast<int>(6 * frame + 3));
        }
        if (!WithinTotalCorrectionLimit(prior, candidate, options_)) {
            damping *= 4.0;
            continue;
        }

        const double candidate_cost = Evaluate(build.factors, prior, candidate, options_, nullptr, nullptr);
        const double predicted_reduction = 0.5 * step.dot(damping * diagonal * step - gradient);
        if (std::isfinite(candidate_cost) && candidate_cost < cost && predicted_reduction > 0.0) {
            const double relative_improvement = (cost - candidate_cost) / std::max(std::abs(cost), 1e-12);
            current = std::move(candidate);
            cost = candidate_cost;
            accepted_any = true;
            ++summary.iterations;
            damping = std::max(1e-9, damping * 0.5);
            if (relative_improvement < 1e-6 || step.norm() < 1e-6) break;
            Evaluate(build.factors, prior, current, options_, &hessian, &gradient);
        } else {
            damping *= 4.0;
        }
    }

    summary.final_cost = cost;
    if (!accepted_any) {
        summary.reason = "no_cost_decreasing_step";
        return summary;
    }

    for (std::size_t index = 0; index < frames.size(); ++index) {
        const auto correction = Correction(initial[index], current[index]);
        summary.max_rotation_update_deg =
            std::max(summary.max_rotation_update_deg, correction.head<3>().norm() * kRadToDeg);
        summary.max_translation_update = std::max(summary.max_translation_update, correction.tail<3>().norm());
        frames[index].pose = SE3(Quatd(current[index].rotation).normalized(), current[index].translation);
    }
    summary.accepted = true;
    summary.reason = "accepted";
    return summary;
}

BundleAdjustmentSummary VoxelBundleAdjuster::OptimizeKeyframes(const std::vector<Keyframe::Ptr>& keyframes,
                                                               const SE3& T_imu_lidar,
                                                               bool apply_results,
                                                               bool use_lio_prior) const {
    std::vector<BundleFrame> frames;
    std::vector<Keyframe::Ptr> selected_keyframes;
    frames.reserve(keyframes.size());
    selected_keyframes.reserve(keyframes.size());
    for (const auto& keyframe : keyframes) {
        if (!keyframe) continue;
        BundleFrame frame;
        frame.id = keyframe->GetID();
        frame.timestamp = keyframe->GetState().timestamp_;
        frame.cloud = keyframe->GetCloud();
        frame.pose = keyframe->GetOptPose() * T_imu_lidar;
        if (use_lio_prior) {
            frame.prior_pose = keyframe->GetLIOPose() * T_imu_lidar;
            frame.has_prior_pose = true;
        }
        frames.push_back(std::move(frame));
        selected_keyframes.push_back(keyframe);
    }
    BundleAdjustmentSummary summary = Optimize(frames);
    if (!summary.accepted || !apply_results) return summary;
    const SE3 T_lidar_imu = T_imu_lidar.inverse();
    for (std::size_t index = 0; index < frames.size(); ++index) {
        selected_keyframes[index]->SetOptPose(frames[index].pose * T_lidar_imu);
    }
    return summary;
}

}  // namespace lightning::backend
