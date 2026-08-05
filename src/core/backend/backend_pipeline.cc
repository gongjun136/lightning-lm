#include "core/backend/backend_pipeline.h"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <utility>

#include <glog/logging.h>
#include <yaml-cpp/yaml.h>

#include "core/graph/optimizer.h"
#include "core/lightning_math.hpp"
#include "core/opti_algo/algo_select.h"
#include "core/robust_kernel/cauchy.h"
#include "core/types/edge_se3.h"
#include "core/types/vertex_se3.h"

namespace lightning::backend {
namespace {

constexpr double kDegToRad = M_PI / 180.0;

template <typename T>
T ValueOr(const YAML::Node& node, const char* key, const T& fallback) {
    if (!node || !node[key]) return fallback;
    try {
        return node[key].as<T>();
    } catch (const YAML::Exception&) {
        return fallback;
    }
}

BundleAdjustmentOptions ReadBaOptions(const YAML::Node& node, BundleAdjustmentOptions options) {
    options.enabled = ValueOr(node, "enabled", options.enabled);
    options.window_size = ValueOr(node, "window_size", options.window_size);
    options.window_stride = ValueOr(node, "window_stride", options.window_stride);
    options.max_iterations = ValueOr(node, "max_iterations", options.max_iterations);
    options.max_points_per_frame = ValueOr(node, "max_points_per_frame", options.max_points_per_frame);
    options.min_points_per_voxel = ValueOr(node, "min_points_per_voxel", options.min_points_per_voxel);
    options.min_frames_per_voxel = ValueOr(node, "min_frames_per_voxel", options.min_frames_per_voxel);
    options.max_voxel_depth = ValueOr(node, "max_voxel_depth", options.max_voxel_depth);
    options.voxel_size = ValueOr(node, "voxel_size", options.voxel_size);
    options.planarity_ratio = ValueOr(node, "planarity_ratio", options.planarity_ratio);
    options.correction_prior_rotation =
        ValueOr(node, "correction_prior_rotation", options.correction_prior_rotation);
    options.correction_prior_translation =
        ValueOr(node, "correction_prior_translation", options.correction_prior_translation);
    options.correction_smoothness_rotation =
        ValueOr(node, "correction_smoothness_rotation", options.correction_smoothness_rotation);
    options.correction_smoothness_translation =
        ValueOr(node, "correction_smoothness_translation", options.correction_smoothness_translation);
    options.max_rotation_step_deg = ValueOr(node, "max_rotation_step_deg", options.max_rotation_step_deg);
    options.max_translation_step = ValueOr(node, "max_translation_step", options.max_translation_step);
    options.max_total_rotation_deg = ValueOr(node, "max_total_rotation_deg", options.max_total_rotation_deg);
    options.max_total_translation = ValueOr(node, "max_total_translation", options.max_total_translation);
    options.initial_damping = ValueOr(node, "initial_damping", options.initial_damping);
    return options;
}

ConfigSetting ReadBtcDescriptorOptions(const YAML::Node& node, ConfigSetting options) {
    options.useful_corner_num_ = ValueOr(node, "useful_corner_num", options.useful_corner_num_);
    options.plane_merge_normal_thre_ = ValueOr(node, "plane_merge_normal_threshold", options.plane_merge_normal_thre_);
    options.plane_merge_dis_thre_ = ValueOr(node, "plane_merge_distance_threshold", options.plane_merge_dis_thre_);
    options.plane_detection_thre_ = ValueOr(node, "plane_detection_threshold", options.plane_detection_thre_);
    options.voxel_size_ = ValueOr(node, "voxel_size", options.voxel_size_);
    options.voxel_init_num_ = ValueOr(node, "voxel_init_points", options.voxel_init_num_);
    options.proj_plane_num_ = ValueOr(node, "projection_plane_count", options.proj_plane_num_);
    options.proj_image_resolution_ = ValueOr(node, "projection_resolution", options.proj_image_resolution_);
    options.proj_image_high_inc_ = ValueOr(node, "projection_height_increment", options.proj_image_high_inc_);
    options.proj_dis_min_ = ValueOr(node, "projection_min_distance", options.proj_dis_min_);
    options.proj_dis_max_ = ValueOr(node, "projection_max_distance", options.proj_dis_max_);
    options.summary_min_thre_ = ValueOr(node, "summary_min_threshold", options.summary_min_thre_);
    options.line_filter_enable_ = ValueOr(node, "line_filter", options.line_filter_enable_);
    options.touch_filter_enable_ = ValueOr(node, "touch_filter", options.touch_filter_enable_);
    options.descriptor_near_num_ = ValueOr(node, "descriptor_near_count", options.descriptor_near_num_);
    options.descriptor_min_len_ = ValueOr(node, "descriptor_min_length", options.descriptor_min_len_);
    options.descriptor_max_len_ = ValueOr(node, "descriptor_max_length", options.descriptor_max_len_);
    options.non_max_suppression_radius_ =
        ValueOr(node, "non_max_suppression_radius", options.non_max_suppression_radius_);
    options.std_side_resolution_ = ValueOr(node, "triangle_side_resolution", options.std_side_resolution_);
    options.skip_near_num_ = ValueOr(node, "skip_near_descriptors", options.skip_near_num_);
    options.candidate_num_ = ValueOr(node, "candidate_count", options.candidate_num_);
    options.candidate_min_votes_ =
        ValueOr(node, "candidate_min_votes", options.candidate_min_votes_);
    options.verification_threads_ =
        ValueOr(node, "verification_threads", options.verification_threads_);
    options.rough_dis_threshold_ = ValueOr(node, "rough_distance_threshold", options.rough_dis_threshold_);
    options.similarity_threshold_ = ValueOr(node, "similarity_threshold", options.similarity_threshold_);
    options.icp_threshold_ = ValueOr(node, "internal_icp_threshold", options.icp_threshold_);
    options.normal_threshold_ = ValueOr(node, "normal_threshold", options.normal_threshold_);
    options.dis_threshold_ = ValueOr(node, "plane_distance_threshold", options.dis_threshold_);
    return options;
}

BackendPipelineOptions ReadOptions(const YAML::Node& root, bool online_mode) {
    BackendPipelineOptions options;
    options.online_mode = online_mode;
    const YAML::Node backend = root["backend"];
    options.verbose = ValueOr(backend, "verbose", options.verbose);
    const YAML::Node graph = backend["pose_graph"];
    options.pose_graph_iterations = ValueOr(graph, "iterations", options.pose_graph_iterations);
    options.pose_graph_outlier_iterations =
        ValueOr(graph, "outlier_iterations", options.pose_graph_outlier_iterations);
    options.motion_translation_noise =
        ValueOr(graph, "motion_translation_noise", options.motion_translation_noise);
    options.motion_rotation_noise_deg =
        ValueOr(graph, "motion_rotation_noise_deg", options.motion_rotation_noise_deg);
    options.loop_translation_noise = ValueOr(graph, "loop_translation_noise", options.loop_translation_noise);
    options.loop_rotation_noise_deg = ValueOr(graph, "loop_rotation_noise_deg", options.loop_rotation_noise_deg);
    options.robust_kernel_delta = ValueOr(graph, "robust_kernel_delta", options.robust_kernel_delta);
    options.loop_outlier_chi2 = ValueOr(graph, "loop_outlier_chi2", options.loop_outlier_chi2);

    options.local_ba = ReadBaOptions(backend["local_ba"], options.local_ba);
    const YAML::Node btc = backend["btc"];
    options.btc.enabled = ValueOr(btc, "enabled", options.btc.enabled);
    options.btc.descriptor_submap_size =
        ValueOr(btc, "descriptor_submap_size", options.btc.descriptor_submap_size);
    options.btc.descriptor_submap_stride =
        ValueOr(btc, "descriptor_submap_stride", options.btc.descriptor_submap_stride);
    options.btc.max_points_per_submap = ValueOr(btc, "max_points_per_submap", options.btc.max_points_per_submap);
    options.btc.min_points_per_submap = ValueOr(btc, "min_points_per_submap", options.btc.min_points_per_submap);
    options.btc.downsample_leaf_size = ValueOr(btc, "downsample_leaf_size", options.btc.downsample_leaf_size);
    options.btc.min_loop_score = ValueOr(btc, "min_loop_score", options.btc.min_loop_score);
    options.btc.max_drift_ratio = ValueOr(btc, "max_drift_ratio", options.btc.max_drift_ratio);
    options.btc.max_rotation_correction_deg =
        ValueOr(btc, "max_rotation_correction_deg", options.btc.max_rotation_correction_deg);
    options.btc.refine_with_plane_icp =
        ValueOr(btc, "refine_with_plane_icp", options.btc.refine_with_plane_icp);
    options.btc.plane_icp_iterations = ValueOr(btc, "plane_icp_iterations", options.btc.plane_icp_iterations);
    options.btc.plane_icp_min_matches = ValueOr(btc, "plane_icp_min_matches", options.btc.plane_icp_min_matches);
    options.btc.plane_icp_min_observability =
        ValueOr(btc, "plane_icp_min_observability", options.btc.plane_icp_min_observability);
    options.btc.allow_degenerate_plane_icp =
        ValueOr(btc, "allow_degenerate_plane_icp", options.btc.allow_degenerate_plane_icp);
    options.btc.degenerate_min_loop_score =
        ValueOr(btc, "degenerate_min_loop_score", options.btc.degenerate_min_loop_score);
    options.btc.degenerate_min_matches =
        ValueOr(btc, "degenerate_min_matches", options.btc.degenerate_min_matches);
    options.btc.confirmation_count = ValueOr(btc, "confirmation_count", options.btc.confirmation_count);
    options.btc.confirmation_max_current_gap =
        ValueOr(btc, "confirmation_max_current_gap", options.btc.confirmation_max_current_gap);
    options.btc.confirmation_max_history_gap =
        ValueOr(btc, "confirmation_max_history_gap", options.btc.confirmation_max_history_gap);
    options.btc.loop_cooldown_descriptors =
        ValueOr(btc, "loop_cooldown_descriptors", options.btc.loop_cooldown_descriptors);
    options.btc.enable_odom_revisit_fallback =
        ValueOr(btc, "enable_odom_revisit_fallback", options.btc.enable_odom_revisit_fallback);
    options.btc.odom_revisit_search_radius =
        ValueOr(btc, "odom_revisit_search_radius", options.btc.odom_revisit_search_radius);
    options.btc.odom_revisit_min_journey =
        ValueOr(btc, "odom_revisit_min_journey", options.btc.odom_revisit_min_journey);
    options.btc.odom_revisit_degenerate_min_matches = ValueOr(
        btc, "odom_revisit_degenerate_min_matches", options.btc.odom_revisit_degenerate_min_matches);
    options.btc.odom_revisit_max_drift_ratio =
        ValueOr(btc, "odom_revisit_max_drift_ratio", options.btc.odom_revisit_max_drift_ratio);
    options.btc.max_odom_revisit_distance =
        ValueOr(btc, "max_odom_revisit_distance", options.btc.max_odom_revisit_distance);
    options.btc.min_optimization_translation =
        ValueOr(btc, "min_optimization_translation", options.btc.min_optimization_translation);
    options.btc.min_optimization_rotation_deg =
        ValueOr(btc, "min_optimization_rotation_deg", options.btc.min_optimization_rotation_deg);
    options.btc.descriptor = ReadBtcDescriptorOptions(btc["descriptor"], options.btc.descriptor);

    const YAML::Node hba = backend["hba"];
    options.hba.enabled = ValueOr(hba, "enabled", options.hba.enabled);
    options.hba.leaf_submap_size = ValueOr(hba, "leaf_submap_size", options.hba.leaf_submap_size);
    options.hba.branching_factor = ValueOr(hba, "branching_factor", options.hba.branching_factor);
    options.hba.max_levels = ValueOr(hba, "max_levels", options.hba.max_levels);
    options.hba.max_points_per_submap = ValueOr(hba, "max_points_per_submap", options.hba.max_points_per_submap);
    options.hba.level_voxel_scale = ValueOr(hba, "level_voxel_scale", options.hba.level_voxel_scale);
    options.hba.global_top_level = ValueOr(hba, "global_top_level", options.hba.global_top_level);
    options.hba.final_local_refinement =
        ValueOr(hba, "final_local_refinement", options.hba.final_local_refinement);
    options.hba.require_applied_loop_for_commit = ValueOr(
        hba, "require_applied_loop_for_commit", options.hba.require_applied_loop_for_commit);
    options.hba.optimizer = ReadBaOptions(hba["optimizer"], options.hba.optimizer);
    return options;
}

double ElapsedMs(const std::chrono::steady_clock::time_point& begin) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - begin).count();
}

}  // namespace

BackendMode ReadBackendMode(const std::string& yaml_path) {
    try {
        const YAML::Node root = YAML::LoadFile(yaml_path);
        const std::string mode = ValueOr<std::string>(root["backend"], "mode", "ba_btc_hba");
        if (mode == "legacy") return BackendMode::kLegacy;
        if (mode == "disabled" || mode == "none") return BackendMode::kDisabled;
        return BackendMode::kBaBtcHba;
    } catch (const YAML::Exception&) {
        return BackendMode::kBaBtcHba;
    }
}

const char* BackendModeName(BackendMode mode) {
    if (mode == BackendMode::kLegacy) return "legacy";
    if (mode == BackendMode::kDisabled) return "disabled";
    return "ba_btc_hba";
}

BackendPipeline::~BackendPipeline() { Shutdown(); }

bool BackendPipeline::Init(const std::string& yaml_path, bool online_mode) {
    if (initialized_) return true;
    try {
        const YAML::Node root = YAML::LoadFile(yaml_path);
        options_ = ReadOptions(root, online_mode);
        const std::vector<double> translation = root["fasterlio"]["extrinsic_T"].as<std::vector<double>>();
        const std::vector<double> rotation = root["fasterlio"]["extrinsic_R"].as<std::vector<double>>();
        T_imu_lidar_ = SE3(Quatd(math::MatFromArray<double>(rotation)).normalized(),
                           math::VecFromArray<double>(translation));
    } catch (const YAML::Exception& error) {
        LOG(ERROR) << "failed to read ba_btc_hba backend configuration: " << error.what();
        return false;
    }

    local_ba_ = std::make_unique<VoxelBundleAdjuster>(options_.local_ba);
    btc_ = std::make_unique<BtcLoopDetector>(options_.btc);
    hba_ = std::make_unique<HierarchicalBundleAdjuster>(options_.hba);
    worker_stop_ = false;
    hba_stop_ = false;
    if (options_.online_mode) worker_thread_ = std::thread(&BackendPipeline::WorkerLoop, this);
    if (options_.hba.enabled) hba_thread_ = std::thread(&BackendPipeline::HbaLoop, this);
    initialized_ = true;
    LOG(INFO) << "backend initialized: mode=ba_btc_hba, online=" << options_.online_mode
              << ", local_ba=" << options_.local_ba.enabled << ", btc=" << options_.btc.enabled
              << ", hba=" << options_.hba.enabled;
    return true;
}

void BackendPipeline::AddKeyframe(const Keyframe::Ptr& keyframe) {
    if (!initialized_ || !keyframe) return;
    if (!options_.online_mode) {
        HandleKeyframe(keyframe);
        return;
    }
    {
        std::lock_guard<std::mutex> lock(queue_mutex_);
        queue_.push_back(keyframe);
    }
    queue_cv_.notify_one();
}

void BackendPipeline::HandleKeyframe(const Keyframe::Ptr& keyframe) {
    std::vector<Keyframe::Ptr> local_window;
    {
        std::lock_guard<std::mutex> lock(data_mutex_);
        if (!keyframes_.empty() && keyframes_.back() == keyframe) return;
        keyframes_.push_back(keyframe);
        runtime_.keyframes = keyframes_.size();
        const std::size_t window = static_cast<std::size_t>(std::max(2, options_.local_ba.window_size));
        if (keyframes_.size() >= window &&
            keyframes_.size() - last_local_ba_keyframe_count_ >=
                static_cast<std::size_t>(std::max(1, options_.local_ba.window_stride))) {
            local_window.assign(keyframes_.end() - window, keyframes_.end());
            last_local_ba_keyframe_count_ = keyframes_.size();
        }
    }

    if (options_.local_ba.enabled && !local_window.empty()) {
        const auto begin = std::chrono::steady_clock::now();
        BundleAdjustmentSummary summary;
        {
            std::lock_guard<std::mutex> lock(optimization_mutex_);
            summary = local_ba_->OptimizeKeyframes(local_window, T_imu_lidar_, true);
        }
        std::lock_guard<std::mutex> lock(data_mutex_);
        ++runtime_.local_ba_attempts;
        runtime_.local_ba_accepted += summary.accepted ? 1 : 0;
        runtime_.local_ba_time_ms += ElapsedMs(begin);
    }

    const auto btc_begin = std::chrono::steady_clock::now();
    const std::optional<BtcLoopResult> detection = btc_->AddKeyframe(keyframe, T_imu_lidar_);
    if (detection) {
        const bool accepted = detection->accepted;
        const bool apply_constraint = accepted && detection->optimization_warranted;
        {
            std::lock_guard<std::mutex> lock(data_mutex_);
            loop_results_.push_back(*detection);
            ++runtime_.btc_descriptors;
            runtime_.btc_candidates += detection->candidate_found ? 1 : 0;
            runtime_.loops_accepted += accepted ? 1 : 0;
            runtime_.loops_applied += apply_constraint ? 1 : 0;
            runtime_.btc_time_ms += ElapsedMs(btc_begin);
            if (apply_constraint) {
                LoopConstraint constraint;
                constraint.detection = *detection;
                constraint.T_history_imu_current_imu =
                    T_imu_lidar_ * detection->T_history_lidar_current_lidar * T_imu_lidar_.inverse();
                loop_constraints_.push_back(std::move(constraint));
            }
        }
        if (apply_constraint) {
            {
                std::lock_guard<std::mutex> lock(optimization_mutex_);
                OptimizePoseGraph();
                UpdateMapToOdomAndNotify();
            }
        }
        // HBA is useful only after a loop constraint has actually entered the
        // graph. A verified revisit whose correction is below the graph-update
        // gate must not launch an expensive point-level optimization that is
        // guaranteed to be rolled back by the transactional commit guard.
        if (apply_constraint) {
            RequestHba(false, "accepted_btc_loop");
        }
    }
    UpdateMapToOdomAndNotify();
}

bool BackendPipeline::OptimizePoseGraph() {
    std::vector<Keyframe::Ptr> keyframes;
    std::vector<LoopConstraint> constraints;
    {
        std::lock_guard<std::mutex> lock(data_mutex_);
        keyframes = keyframes_;
        constraints = loop_constraints_;
    }
    if (keyframes.size() < 2 || constraints.empty()) return false;
    const auto begin = std::chrono::steady_clock::now();

    miao::OptimizerConfig config(miao::AlgorithmType::LEVENBERG_MARQUARDT,
                                 miao::LinearSolverType::LINEAR_SOLVER_SPARSE_EIGEN, false);
    auto optimizer = miao::SetupOptimizer<6, 3>(config);
    std::unordered_map<unsigned long, Keyframe::Ptr> keyframes_by_id;
    for (std::size_t index = 0; index < keyframes.size(); ++index) {
        if (!keyframes[index]) continue;
        auto vertex = std::make_shared<miao::VertexSE3>();
        vertex->SetId(static_cast<int>(keyframes[index]->GetID()));
        vertex->SetEstimate(keyframes[index]->GetOptPose());
        if (index == 0) vertex->SetFixed(true);
        optimizer->AddVertex(vertex);
        keyframes_by_id[keyframes[index]->GetID()] = keyframes[index];
    }

    Mat6d motion_information = Mat6d::Zero();
    motion_information.diagonal().head<3>().setConstant(
        1.0 / std::pow(options_.motion_translation_noise, 2));
    motion_information.diagonal().tail<3>().setConstant(
        1.0 / std::pow(options_.motion_rotation_noise_deg * kDegToRad, 2));
    for (std::size_t index = 1; index < keyframes.size(); ++index) {
        if (!keyframes[index - 1] || !keyframes[index]) continue;
        auto edge = std::make_shared<miao::EdgeSE3>();
        edge->SetVertex(0, optimizer->GetVertex(keyframes[index - 1]->GetID()));
        edge->SetVertex(1, optimizer->GetVertex(keyframes[index]->GetID()));
        edge->SetMeasurement(keyframes[index - 1]->GetLIOPose().inverse() * keyframes[index]->GetLIOPose());
        edge->SetInformation(motion_information);
        optimizer->AddEdge(edge);
    }

    Mat6d loop_information = Mat6d::Zero();
    loop_information.diagonal().head<3>().setConstant(1.0 / std::pow(options_.loop_translation_noise, 2));
    loop_information.diagonal().tail<3>().setConstant(
        1.0 / std::pow(options_.loop_rotation_noise_deg * kDegToRad, 2));
    std::vector<std::shared_ptr<miao::EdgeSE3>> loop_edges;
    std::vector<std::size_t> constraint_indices;
    for (std::size_t index = 0; index < constraints.size(); ++index) {
        const auto& constraint = constraints[index];
        auto history = optimizer->GetVertex(constraint.detection.history_keyframe_id);
        auto current = optimizer->GetVertex(constraint.detection.current_keyframe_id);
        if (!history || !current) continue;
        auto edge = std::make_shared<miao::EdgeSE3>();
        edge->SetVertex(0, history);
        edge->SetVertex(1, current);
        edge->SetMeasurement(constraint.T_history_imu_current_imu);
        edge->SetInformation(loop_information);
        auto robust = std::make_shared<miao::RobustKernelCauchy>();
        robust->SetDelta(options_.robust_kernel_delta);
        edge->SetRobustKernel(robust);
        optimizer->AddEdge(edge);
        loop_edges.push_back(edge);
        constraint_indices.push_back(index);
    }
    if (loop_edges.empty() || !optimizer->InitializeOptimization()) return false;
    optimizer->SetVerbose(false);
    optimizer->Optimize(std::max(1, options_.pose_graph_iterations));
    optimizer->ComputeActiveErrors();

    std::size_t inliers = 0;
    for (std::size_t index = 0; index < loop_edges.size(); ++index) {
        const bool inlier = std::isfinite(loop_edges[index]->Chi2()) &&
                            loop_edges[index]->Chi2() <= options_.loop_outlier_chi2;
        constraints[constraint_indices[index]].graph_inlier = inlier;
        if (inlier) {
            loop_edges[index]->SetRobustKernel(nullptr);
            ++inliers;
        } else {
            loop_edges[index]->SetLevel(1);
        }
    }
    if (inliers > 0 && optimizer->InitializeOptimization(0)) {
        optimizer->Optimize(std::max(1, options_.pose_graph_outlier_iterations));
    }
    for (const auto& [id, keyframe] : keyframes_by_id) {
        auto vertex = std::dynamic_pointer_cast<miao::VertexSE3>(optimizer->GetVertex(static_cast<int>(id)));
        if (vertex && keyframe) keyframe->SetOptPose(vertex->Estimate());
    }

    {
        std::lock_guard<std::mutex> lock(data_mutex_);
        for (std::size_t index = 0; index < constraints.size() && index < loop_constraints_.size(); ++index) {
            loop_constraints_[index].graph_inlier = constraints[index].graph_inlier;
        }
        runtime_.loops_graph_inliers = inliers;
        runtime_.pose_graph_time_ms += ElapsedMs(begin);
    }
    if (options_.verbose) {
        LOG(INFO) << "BTC pose graph optimized: keyframes=" << keyframes.size()
                  << ", loop_inliers=" << inliers << "/" << loop_edges.size();
    }
    return inliers > 0;
}

void BackendPipeline::RequestGlobalOptimization(const std::string& reason) {
    RequestHba(true, reason);
}

void BackendPipeline::RequestHba(bool rerun_pose_graph, const std::string& reason) {
    if (!initialized_) return;
    if (!options_.hba.enabled) {
        if (rerun_pose_graph) {
            std::lock_guard<std::mutex> lock(optimization_mutex_);
            OptimizePoseGraph();
            UpdateMapToOdomAndNotify();
        }
        return;
    }
    {
        std::lock_guard<std::mutex> lock(hba_mutex_);
        hba_requested_ = true;
        hba_rerun_pose_graph_ = hba_rerun_pose_graph_ || rerun_pose_graph;
        hba_reason_ = reason;
    }
    hba_cv_.notify_one();
}

void BackendPipeline::WorkerLoop() {
    while (true) {
        Keyframe::Ptr keyframe;
        {
            std::unique_lock<std::mutex> lock(queue_mutex_);
            queue_cv_.wait(lock, [this]() { return worker_stop_ || !queue_.empty(); });
            if (worker_stop_ && queue_.empty()) break;
            keyframe = queue_.front();
            queue_.pop_front();
            worker_processing_ = true;
        }
        HandleKeyframe(keyframe);
        {
            std::lock_guard<std::mutex> lock(queue_mutex_);
            worker_processing_ = false;
            if (queue_.empty()) queue_idle_cv_.notify_all();
        }
    }
}

void BackendPipeline::HbaLoop() {
    while (true) {
        bool rerun_pose_graph = false;
        std::string reason;
        {
            std::unique_lock<std::mutex> lock(hba_mutex_);
            hba_cv_.wait(lock, [this]() { return hba_stop_ || hba_requested_; });
            if (hba_stop_ && !hba_requested_) break;
            hba_requested_ = false;
            rerun_pose_graph = hba_rerun_pose_graph_;
            hba_rerun_pose_graph_ = false;
            hba_running_ = true;
            reason = hba_reason_;
        }

        std::vector<Keyframe::Ptr> keyframes;
        bool has_applied_loop = false;
        {
            std::lock_guard<std::mutex> lock(data_mutex_);
            keyframes = keyframes_;
            has_applied_loop = runtime_.loops_applied > 0;
        }
        const auto begin = std::chrono::steady_clock::now();
        HierarchicalBundleAdjustmentSummary summary;
        {
            std::lock_guard<std::mutex> lock(optimization_mutex_);
            if (rerun_pose_graph) OptimizePoseGraph();
            std::vector<SE3> poses_before_hba;
            poses_before_hba.reserve(keyframes.size());
            for (const auto& keyframe : keyframes) {
                poses_before_hba.push_back(keyframe ? keyframe->GetOptPose() : SE3());
            }
            summary = hba_->Optimize(keyframes, T_imu_lidar_);
            if (summary.accepted && options_.hba.require_applied_loop_for_commit &&
                !has_applied_loop) {
                for (std::size_t index = 0; index < keyframes.size(); ++index) {
                    if (keyframes[index]) keyframes[index]->SetOptPose(poses_before_hba[index]);
                }
                summary.accepted = false;
                summary.reason = "rolled_back_without_applied_loop";
            }
            UpdateMapToOdomAndNotify();
        }
        {
            std::lock_guard<std::mutex> lock(data_mutex_);
            ++runtime_.hba_runs;
            runtime_.hba_accepted += summary.accepted ? 1 : 0;
            runtime_.hba_time_ms += ElapsedMs(begin);
        }
        if (options_.verbose) {
            LOG(INFO) << "HBA finished: reason=" << reason << ", keyframes=" << keyframes.size()
                      << ", levels=" << summary.levels_accepted << "/" << summary.levels_attempted
                      << ", global_top=" << summary.global_top_level_accepted << "/"
                      << summary.global_top_level_attempted
                      << ", status=" << summary.reason << ", elapsed_ms=" << summary.elapsed_ms;
        }
        {
            std::lock_guard<std::mutex> lock(hba_mutex_);
            hba_running_ = false;
            if (!hba_requested_) hba_idle_cv_.notify_all();
        }
    }
}

void BackendPipeline::UpdateMapToOdomAndNotify() {
    OptimizedCallback callback;
    {
        std::lock_guard<std::mutex> lock(data_mutex_);
        if (!keyframes_.empty() && keyframes_.back()) {
            T_map_odom_ = keyframes_.back()->GetOptPose() * keyframes_.back()->GetLIOPose().inverse();
        }
        callback = optimized_callback_;
    }
    if (callback) callback();
}

void BackendPipeline::WaitUntilIdle(bool force_global_optimization) {
    if (!initialized_) return;
    if (options_.online_mode) {
        std::unique_lock<std::mutex> lock(queue_mutex_);
        queue_idle_cv_.wait(lock, [this]() { return queue_.empty() && !worker_processing_; });
    }
    if (force_global_optimization) RequestGlobalOptimization("data_end");
    if (options_.hba.enabled) {
        std::unique_lock<std::mutex> lock(hba_mutex_);
        hba_idle_cv_.wait(lock, [this]() { return !hba_requested_ && !hba_running_; });
    }
}

void BackendPipeline::Shutdown() {
    if (!initialized_) return;
    WaitUntilIdle(false);
    if (worker_thread_.joinable()) {
        {
            std::lock_guard<std::mutex> lock(queue_mutex_);
            worker_stop_ = true;
        }
        queue_cv_.notify_all();
        worker_thread_.join();
    }
    if (hba_thread_.joinable()) {
        {
            std::lock_guard<std::mutex> lock(hba_mutex_);
            hba_stop_ = true;
        }
        hba_cv_.notify_all();
        hba_thread_.join();
    }
    initialized_ = false;
}

void BackendPipeline::SetOptimizedCallback(OptimizedCallback callback) {
    std::lock_guard<std::mutex> lock(data_mutex_);
    optimized_callback_ = std::move(callback);
}

SE3 BackendPipeline::GetMapToOdom() const {
    std::lock_guard<std::mutex> lock(data_mutex_);
    return T_map_odom_;
}

std::vector<BtcLoopResult> BackendPipeline::GetLoopResults() const {
    std::lock_guard<std::mutex> lock(data_mutex_);
    return loop_results_;
}

BackendRuntimeSummary BackendPipeline::GetRuntimeSummary() const {
    std::lock_guard<std::mutex> lock(data_mutex_);
    return runtime_;
}

bool BackendPipeline::SaveDiagnostics(
    const std::string& directory,
    const map_frame::Metadata* map_metadata) const {
    std::vector<BtcLoopResult> loops;
    std::vector<LoopConstraint> constraints;
    BackendRuntimeSummary runtime;
    {
        std::lock_guard<std::mutex> lock(data_mutex_);
        loops = loop_results_;
        constraints = loop_constraints_;
        runtime = runtime_;
    }
    std::error_code error;
    std::filesystem::create_directories(directory, error);
    if (error) return false;

    std::ofstream loop_csv(std::filesystem::path(directory) / "btc_loop_candidates.csv");
    if (!loop_csv.is_open()) return false;
    loop_csv << "current_descriptor,history_descriptor,current_keyframe,history_keyframe,current_timestamp,"
                "history_timestamp,candidate,candidate_source,accepted,optimization_warranted,score,drift_translation,"
                "drift_rotation_deg,journey_span,drift_ratio,odom_revisit_distance,observability,"
                "plane_icp_matches,plane_icp_converged,degenerate_fallback,"
                "confirmation_count,points,descriptors,"
                "generation_ms,search_ms,reason\n";
    for (const auto& loop : loops) {
        loop_csv << loop.current_descriptor_id << ',' << loop.history_descriptor_id << ','
                 << loop.current_keyframe_id << ',' << loop.history_keyframe_id << ',' << std::setprecision(12)
                 << loop.current_timestamp << ',' << loop.history_timestamp << ',' << loop.candidate_found << ','
                 << loop.candidate_source << ',' << loop.accepted << ',' << loop.optimization_warranted << ','
                 << loop.score << ','
                 << loop.drift_translation << ',' << loop.drift_rotation_deg << ',' << loop.journey_span << ','
                 << loop.drift_ratio << ',' << loop.odom_revisit_distance << ','
                 << loop.plane_icp_observability << ',' << loop.plane_icp_matches << ','
                 << loop.plane_icp_converged << ',' << loop.used_degenerate_plane_fallback << ','
                 << loop.confirmation_count << ',' << loop.point_count << ',' << loop.descriptor_count << ','
                 << loop.generation_time_ms << ',' << loop.search_time_ms << ',' << loop.rejection_reason << '\n';
    }

    std::ofstream summary(std::filesystem::path(directory) / "backend_summary.yaml");
    if (!summary.is_open()) return false;
    summary << "mode: ba_btc_hba\n"
            << "keyframes: " << runtime.keyframes << '\n'
            << "local_ba_attempts: " << runtime.local_ba_attempts << '\n'
            << "local_ba_accepted: " << runtime.local_ba_accepted << '\n'
            << "btc_descriptors: " << runtime.btc_descriptors << '\n'
            << "btc_candidates: " << runtime.btc_candidates << '\n'
            << "loops_accepted: " << runtime.loops_accepted << '\n'
            << "loops_applied: " << runtime.loops_applied << '\n'
            << "loops_graph_inliers: " << runtime.loops_graph_inliers << '\n'
            << "hba_runs: " << runtime.hba_runs << '\n'
            << "hba_accepted: " << runtime.hba_accepted << '\n'
            << "local_ba_time_ms: " << runtime.local_ba_time_ms << '\n'
            << "btc_time_ms: " << runtime.btc_time_ms << '\n'
            << "pose_graph_time_ms: " << runtime.pose_graph_time_ms << '\n'
            << "hba_time_ms: " << runtime.hba_time_ms << '\n';
    if (map_metadata && map_metadata->normalized) {
        summary << "map_frame_transform_id: " << map_metadata->transform_id << '\n'
                << "map_frame_z_offset: "
                << map_metadata->T_export_slam.translation().z() << '\n';
    }
    return true;
}

bool BackendPipeline::SaveRelocalizationDatabase(
    const std::string& directory,
    const map_frame::Metadata* map_metadata) const {
    std::lock_guard<std::mutex> lock(data_mutex_);
    if (!btc_) {
        LOG(ERROR) << "cannot save BTC relocalization database: BTC is not initialized";
        return false;
    }
    return btc_->SaveRelocalizationDatabase(directory, T_imu_lidar_, map_metadata);
}

}  // namespace lightning::backend
