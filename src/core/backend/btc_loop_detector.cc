#include "core/backend/btc_loop_detector.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <utility>

#include <Eigen/Eigenvalues>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>

namespace lightning::backend {
namespace {

constexpr double kRadToDeg = 180.0 / M_PI;

Eigen::Matrix3d Hat(const Eigen::Vector3d& value) {
    Eigen::Matrix3d result;
    result << 0.0, -value.z(), value.y(), value.z(), 0.0, -value.x(), -value.y(), value.x(), 0.0;
    return result;
}

double ElapsedMilliseconds(const std::chrono::steady_clock::time_point& begin) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - begin).count();
}

}  // namespace

BtcLoopDetector::BtcLoopDetector(BtcLoopDetectorOptions options)
    : options_(std::move(options)), manager_(options_.descriptor) {}

void BtcLoopDetector::Reset() {
    manager_ = STDescManager(options_.descriptor);
    pending_keyframes_.clear();
    entries_.clear();
    last_keyframe_.reset();
    journey_ = 0.0;
    last_confirmation_current_ = -1;
    last_confirmation_history_ = -1;
    confirmation_count_ = 0;
    last_accepted_descriptor_ = -1000000;
}

std::optional<BtcLoopResult> BtcLoopDetector::AddKeyframe(const Keyframe::Ptr& keyframe,
                                                          const SE3& T_imu_lidar) {
    if (!options_.enabled || !keyframe) return std::nullopt;
    if (last_keyframe_) {
        journey_ += (keyframe->GetLIOPose().translation() - last_keyframe_->GetLIOPose().translation()).norm();
    }
    last_keyframe_ = keyframe;
    pending_keyframes_.push_back(keyframe);
    if (pending_keyframes_.size() < static_cast<std::size_t>(std::max(1, options_.descriptor_submap_size))) {
        return std::nullopt;
    }
    std::vector<Keyframe::Ptr> keyframes;
    keyframes.swap(pending_keyframes_);
    return ProcessSubmap(keyframes, T_imu_lidar);
}

pcl::PointCloud<pcl::PointXYZI>::Ptr BtcLoopDetector::BuildSubmap(
    const std::vector<Keyframe::Ptr>& keyframes, const SE3& T_imu_lidar) const {
    pcl::PointCloud<pcl::PointXYZI>::Ptr combined(new pcl::PointCloud<pcl::PointXYZI>);
    if (keyframes.empty() || !keyframes.back()) return combined;
    const SE3 T_world_current_lidar = keyframes.back()->GetOptPose() * T_imu_lidar;
    const SE3 T_current_lidar_world = T_world_current_lidar.inverse();

    std::size_t available_points = 0;
    for (const auto& keyframe : keyframes) {
        if (keyframe && keyframe->GetCloud()) available_points += keyframe->GetCloud()->size();
    }
    const std::size_t maximum = static_cast<std::size_t>(std::max(1, options_.max_points_per_submap));
    const std::size_t stride = std::max<std::size_t>(1, (available_points + maximum - 1) / maximum);
    combined->reserve(std::min(available_points, maximum));

    std::size_t source_index = 0;
    for (const auto& keyframe : keyframes) {
        if (!keyframe || !keyframe->GetCloud()) continue;
        const SE3 T_current_lidar_frame_lidar =
            T_current_lidar_world * keyframe->GetOptPose() * T_imu_lidar;
        for (const auto& source : keyframe->GetCloud()->points) {
            if (source_index++ % stride != 0) continue;
            if (!std::isfinite(source.x) || !std::isfinite(source.y) || !std::isfinite(source.z)) continue;
            const Eigen::Vector3d transformed =
                T_current_lidar_frame_lidar * Eigen::Vector3d(source.x, source.y, source.z);
            pcl::PointXYZI point;
            point.x = static_cast<float>(transformed.x());
            point.y = static_cast<float>(transformed.y());
            point.z = static_cast<float>(transformed.z());
            point.intensity = source.intensity;
            combined->push_back(point);
        }
    }

    if (options_.downsample_leaf_size <= 0.0 || combined->empty()) return combined;
    pcl::PointCloud<pcl::PointXYZI>::Ptr filtered(new pcl::PointCloud<pcl::PointXYZI>);
    pcl::VoxelGrid<pcl::PointXYZI> voxel_grid;
    voxel_grid.setLeafSize(options_.downsample_leaf_size, options_.downsample_leaf_size,
                           options_.downsample_leaf_size);
    voxel_grid.setInputCloud(combined);
    voxel_grid.filter(*filtered);
    return filtered;
}

BtcLoopDetector::RefineSummary BtcLoopDetector::RefinePlaneTransform(
    const pcl::PointCloud<pcl::PointXYZINormal>::ConstPtr& current,
    const pcl::PointCloud<pcl::PointXYZINormal>::ConstPtr& history,
    Eigen::Matrix3d& rotation, Eigen::Vector3d& translation) const {
    RefineSummary summary;
    if (!current || !history || current->empty() || history->empty()) return summary;

    pcl::PointCloud<pcl::PointXYZ>::Ptr history_xyz(new pcl::PointCloud<pcl::PointXYZ>);
    history_xyz->reserve(history->size());
    for (const auto& point : history->points) history_xyz->push_back(pcl::PointXYZ(point.x, point.y, point.z));
    pcl::KdTreeFLANN<pcl::PointXYZ> kd_tree;
    kd_tree.setInputCloud(history_xyz);

    bool fine_stage = false;
    Eigen::Matrix3d normal_information = Eigen::Matrix3d::Zero();
    for (int iteration = 0; iteration < std::max(1, options_.plane_icp_iterations); ++iteration) {
        Eigen::Matrix<double, 6, 6> hessian = Eigen::Matrix<double, 6, 6>::Zero();
        Eigen::Matrix<double, 6, 1> gradient = Eigen::Matrix<double, 6, 1>::Zero();
        normal_information.setZero();
        int matches = 0;
        const double normal_threshold = fine_stage ? 0.10 : 0.20;
        const double plane_threshold = fine_stage ? 0.10 : 0.50;
        const double distance_threshold = fine_stage ? 1.0 : 3.0;

        for (const auto& source : current->points) {
            const Eigen::Vector3d local(source.x, source.y, source.z);
            const Eigen::Vector3d point = rotation * local + translation;
            pcl::PointXYZ query(static_cast<float>(point.x()), static_cast<float>(point.y()),
                                static_cast<float>(point.z()));
            std::vector<int> indices(1);
            std::vector<float> distances(1);
            if (kd_tree.nearestKSearch(query, 1, indices, distances) <= 0) continue;
            const auto& target = history->points[indices[0]];
            const Eigen::Vector3d target_point(target.x, target.y, target.z);
            const Eigen::Vector3d source_normal =
                rotation * Eigen::Vector3d(source.normal_x, source.normal_y, source.normal_z);
            const Eigen::Vector3d target_normal(target.normal_x, target.normal_y, target.normal_z);
            if (std::min((source_normal - target_normal).norm(), (source_normal + target_normal).norm()) >=
                normal_threshold) {
                continue;
            }
            const Eigen::Vector3d difference = point - target_point;
            const double residual = target_normal.dot(difference);
            if (std::abs(residual) >= plane_threshold || difference.norm() >= distance_threshold) continue;

            Eigen::Matrix<double, 6, 1> jacobian;
            jacobian.head<3>() = Hat(local) * rotation.transpose() * target_normal;
            jacobian.tail<3>() = target_normal;
            hessian.noalias() += jacobian * jacobian.transpose();
            gradient.noalias() += jacobian * residual;
            normal_information.noalias() += target_normal * target_normal.transpose();
            ++matches;
        }

        summary.matches = matches;
        if (matches < std::max(6, options_.plane_icp_min_matches)) return summary;
        Eigen::LDLT<Eigen::Matrix<double, 6, 6>> solver(hessian);
        if (solver.info() != Eigen::Success) return summary;
        const Eigen::Matrix<double, 6, 1> step = solver.solve(-gradient);
        if (solver.info() != Eigen::Success || !step.allFinite()) return summary;
        rotation = rotation * Sophus::SO3d::exp(step.head<3>()).matrix();
        translation += step.tail<3>();

        if (step.head<3>().norm() < 1e-3 && step.tail<3>().norm() < 1e-3) {
            if (fine_stage) {
                summary.converged = true;
                break;
            }
            fine_stage = true;
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(normal_information);
    if (eigen_solver.info() != Eigen::Success) return summary;
    summary.observability = eigen_solver.eigenvalues().x();
    summary.converged = summary.converged || fine_stage;
    summary.accepted = summary.converged && summary.matches >= std::max(6, options_.plane_icp_min_matches) &&
                       summary.observability >= options_.plane_icp_min_observability;
    return summary;
}

BtcLoopResult BtcLoopDetector::ProcessSubmap(const std::vector<Keyframe::Ptr>& keyframes,
                                             const SE3& T_imu_lidar) {
    BtcLoopResult result;
    result.current_descriptor_id = static_cast<int>(entries_.size());
    if (!options_.enabled) {
        result.rejection_reason = "disabled";
        return result;
    }
    if (keyframes.empty() || !keyframes.front() || !keyframes.back()) {
        result.rejection_reason = "invalid_keyframes";
        return result;
    }

    const Keyframe::Ptr& current_endpoint = keyframes.back();
    result.current_keyframe_id = current_endpoint->GetID();
    result.current_timestamp = current_endpoint->GetState().timestamp_;
    auto cloud = BuildSubmap(keyframes, T_imu_lidar);
    result.point_count = cloud->size();
    if (cloud->size() < static_cast<std::size_t>(std::max(1, options_.min_points_per_submap))) {
        result.rejection_reason = "too_few_submap_points";
        return result;
    }

    const auto generation_begin = std::chrono::steady_clock::now();
    std::vector<STD> descriptors;
    manager_.GenerateSTDescs(cloud, descriptors, static_cast<int>(current_endpoint->GetID()));
    result.generation_time_ms = ElapsedMilliseconds(generation_begin);
    result.descriptor_generated = true;
    result.descriptor_count = descriptors.size();

    std::pair<int, double> search_result(-1, 0.0);
    std::pair<Eigen::Vector3d, Eigen::Matrix3d> transform;
    std::vector<std::pair<STD, STD>> matches;
    const auto search_begin = std::chrono::steady_clock::now();
    if (!entries_.empty() && !manager_.plane_cloud_vec_.empty()) {
        manager_.SearchLoop(descriptors, search_result, transform, matches, manager_.plane_cloud_vec_.back());
    }
    result.search_time_ms = ElapsedMilliseconds(search_begin);

    BtcDescriptorEntry current_entry;
    current_entry.descriptor_id = result.current_descriptor_id;
    current_entry.first_keyframe_id = keyframes.front()->GetID();
    current_entry.last_keyframe_id = current_endpoint->GetID();
    current_entry.timestamp = result.current_timestamp;
    current_entry.journey = journey_;
    current_entry.endpoint = current_endpoint;

    if (search_result.first >= 0 && search_result.first < static_cast<int>(entries_.size())) {
        result.candidate_found = true;
        result.history_descriptor_id = search_result.first;
        result.score = search_result.second;
        const BtcDescriptorEntry& history_entry = entries_[search_result.first];
        result.history_keyframe_id = history_entry.last_keyframe_id;
        result.history_timestamp = history_entry.timestamp;
        if (history_entry.endpoint) {
            result.odom_revisit_distance =
                (current_endpoint->GetLIOPose().translation() -
                 history_entry.endpoint->GetLIOPose().translation())
                    .norm();
        }

        Eigen::Vector3d translation = transform.first;
        Eigen::Matrix3d rotation = transform.second;
        if (result.score < options_.min_loop_score) {
            result.rejection_reason = "score_below_threshold";
        } else if (options_.max_odom_revisit_distance > 0.0 &&
                   result.odom_revisit_distance > options_.max_odom_revisit_distance) {
            result.rejection_reason = "odom_revisit_distance_exceeded";
        } else {
            RefineSummary refine;
            if (options_.refine_with_plane_icp) {
                refine = RefinePlaneTransform(manager_.plane_cloud_vec_.back(),
                                              manager_.plane_cloud_vec_[search_result.first], rotation,
                                              translation);
                result.plane_icp_observability = refine.observability;
                result.plane_icp_matches = refine.matches;
                result.plane_icp_converged = refine.converged;
            } else {
                refine.accepted = true;
                refine.converged = true;
            }

            const bool consistent_with_pending =
                last_confirmation_current_ >= 0 &&
                result.current_descriptor_id - last_confirmation_current_ <=
                    std::max(1, options_.confirmation_max_current_gap) &&
                std::abs(result.history_descriptor_id - last_confirmation_history_) <=
                    std::max(0, options_.confirmation_max_history_gap);
            const bool converged_degenerate_fallback =
                refine.converged && refine.matches >= options_.degenerate_min_matches;
            const bool temporally_confirmed_fallback =
                consistent_with_pending && refine.matches >= options_.degenerate_min_matches &&
                refine.observability >= options_.plane_icp_min_observability;
            const bool degenerate_fallback =
                options_.refine_with_plane_icp && options_.allow_degenerate_plane_icp &&
                result.score >= options_.degenerate_min_loop_score &&
                (converged_degenerate_fallback || temporally_confirmed_fallback);
            result.used_degenerate_plane_fallback = !refine.accepted && degenerate_fallback;
            if (!refine.accepted && !degenerate_fallback) {
                result.rejection_reason = "plane_icp_rejected";
            } else if (!history_entry.endpoint) {
                result.rejection_reason = "missing_history_endpoint";
            } else {
                result.T_history_lidar_current_lidar =
                    SE3(Quatd(rotation).normalized(), translation);
                const SE3 T_world_history_lidar = history_entry.endpoint->GetOptPose() * T_imu_lidar;
                const SE3 T_world_current_lidar = current_endpoint->GetOptPose() * T_imu_lidar;
                const SE3 T_world_current_from_loop =
                    T_world_history_lidar * result.T_history_lidar_current_lidar;
                const SE3 correction = T_world_current_from_loop.inverse() * T_world_current_lidar;
                result.drift_translation = correction.translation().norm();
                result.drift_rotation_deg = correction.so3().log().norm() * kRadToDeg;
                result.journey_span = std::max(0.0, journey_ - history_entry.journey);
                result.drift_ratio = result.drift_translation / std::max(1e-6, result.journey_span);

                if (result.journey_span <= 1e-6) {
                    result.rejection_reason = "invalid_journey_span";
                } else if (result.drift_ratio >= options_.max_drift_ratio) {
                    result.rejection_reason = "drift_ratio_exceeded";
                } else if (result.drift_rotation_deg >= options_.max_rotation_correction_deg) {
                    result.rejection_reason = "rotation_correction_exceeded";
                } else if (result.current_descriptor_id - last_accepted_descriptor_ <=
                           options_.loop_cooldown_descriptors) {
                    result.rejection_reason = "loop_cooldown";
                } else {
                    confirmation_count_ =
                        consistent_with_pending ? confirmation_count_ + 1 : 1;
                    last_confirmation_current_ = result.current_descriptor_id;
                    last_confirmation_history_ = result.history_descriptor_id;
                    result.confirmation_count = confirmation_count_;
                    if (confirmation_count_ < std::max(1, options_.confirmation_count)) {
                        result.rejection_reason = "awaiting_confirmation";
                    } else {
                        result.accepted = true;
                        result.optimization_warranted =
                            result.drift_translation >= options_.min_optimization_translation ||
                            result.drift_rotation_deg >= options_.min_optimization_rotation_deg;
                        result.rejection_reason = result.optimization_warranted
                                                      ? "accepted"
                                                      : "accepted_no_optimization_needed";
                        last_accepted_descriptor_ = result.current_descriptor_id;
                        confirmation_count_ = 0;
                        last_confirmation_current_ = -1;
                        last_confirmation_history_ = -1;
                    }
                }
            }
        }
    } else {
        result.rejection_reason = descriptors.empty() ? "no_descriptors" : "no_candidate";
    }

    manager_.AddSTDescs(descriptors);
    entries_.push_back(std::move(current_entry));
    return result;
}

}  // namespace lightning::backend
