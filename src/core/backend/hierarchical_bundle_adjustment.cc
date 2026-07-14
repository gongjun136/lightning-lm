#include "core/backend/hierarchical_bundle_adjustment.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <utility>

namespace lightning::backend {
namespace {

struct HierarchyNode {
    BundleFrame frame;
    std::vector<Keyframe::Ptr> members;
};

HierarchyNode BuildNode(const std::vector<Keyframe::Ptr>& members, const SE3& T_imu_lidar,
                        int max_points) {
    HierarchyNode node;
    node.members = members;
    if (members.empty() || !members.back()) return node;
    const Keyframe::Ptr& anchor = members.back();
    const SE3 T_world_anchor_lidar = anchor->GetOptPose() * T_imu_lidar;
    const SE3 T_anchor_lidar_world = T_world_anchor_lidar.inverse();
    node.frame.id = anchor->GetID();
    node.frame.timestamp = anchor->GetState().timestamp_;
    node.frame.pose = T_world_anchor_lidar;
    node.frame.cloud = std::make_shared<PointCloudType>();

    std::size_t available = 0;
    for (const auto& keyframe : members) {
        if (keyframe && keyframe->GetCloud()) available += keyframe->GetCloud()->size();
    }
    const std::size_t maximum = static_cast<std::size_t>(std::max(1, max_points));
    const std::size_t stride = std::max<std::size_t>(1, (available + maximum - 1) / maximum);
    node.frame.cloud->reserve(std::min(available, maximum));
    std::size_t source_index = 0;

    for (const auto& keyframe : members) {
        if (!keyframe || !keyframe->GetCloud()) continue;
        const SE3 T_anchor_lidar_frame_lidar =
            T_anchor_lidar_world * keyframe->GetOptPose() * T_imu_lidar;
        for (const auto& source : keyframe->GetCloud()->points) {
            if (source_index++ % stride != 0) continue;
            if (!std::isfinite(source.x) || !std::isfinite(source.y) || !std::isfinite(source.z)) continue;
            const Eigen::Vector3d point =
                T_anchor_lidar_frame_lidar * Eigen::Vector3d(source.x, source.y, source.z);
            PointType target = source;
            target.x = static_cast<float>(point.x());
            target.y = static_cast<float>(point.y());
            target.z = static_cast<float>(point.z());
            node.frame.cloud->push_back(target);
        }
    }
    return node;
}

std::vector<HierarchyNode> BuildLeafNodes(const std::vector<Keyframe::Ptr>& keyframes,
                                          const SE3& T_imu_lidar, int leaf_size,
                                          int max_points) {
    std::vector<HierarchyNode> nodes;
    const std::size_t size = static_cast<std::size_t>(std::max(2, leaf_size));
    for (std::size_t begin = 0; begin < keyframes.size(); begin += size) {
        const std::size_t end = std::min(keyframes.size(), begin + size);
        std::vector<Keyframe::Ptr> members(keyframes.begin() + begin, keyframes.begin() + end);
        if (members.size() >= 2) nodes.push_back(BuildNode(members, T_imu_lidar, max_points));
    }
    return nodes;
}

std::vector<HierarchyNode> GroupNodes(const std::vector<HierarchyNode>& nodes,
                                      const SE3& T_imu_lidar, int branching,
                                      int max_points) {
    std::vector<HierarchyNode> grouped;
    const std::size_t group_size = static_cast<std::size_t>(std::max(2, branching));
    for (std::size_t begin = 0; begin < nodes.size(); begin += group_size) {
        const std::size_t end = std::min(nodes.size(), begin + group_size);
        std::vector<Keyframe::Ptr> members;
        for (std::size_t index = begin; index < end; ++index) {
            members.insert(members.end(), nodes[index].members.begin(), nodes[index].members.end());
        }
        if (members.size() >= 2) grouped.push_back(BuildNode(members, T_imu_lidar, max_points));
    }
    return grouped;
}

void ApplyCorrections(const std::vector<BundleFrame>& frames, std::vector<HierarchyNode>& nodes,
                      std::size_t begin, const SE3& T_imu_lidar) {
    const SE3 T_lidar_imu = T_imu_lidar.inverse();
    for (std::size_t offset = 0; offset < frames.size(); ++offset) {
        const std::size_t index = begin + offset;
        const SE3 correction = frames[offset].pose * nodes[index].frame.pose.inverse();
        nodes[index].frame.pose = frames[offset].pose;
        for (const auto& keyframe : nodes[index].members) {
            if (!keyframe) continue;
            const SE3 updated_lidar_pose = correction * keyframe->GetOptPose() * T_imu_lidar;
            keyframe->SetOptPose(updated_lidar_pose * T_lidar_imu);
        }
    }
}

}  // namespace

HierarchicalBundleAdjuster::HierarchicalBundleAdjuster(HierarchicalBundleAdjustmentOptions options)
    : options_(std::move(options)) {}

HierarchicalBundleAdjustmentSummary HierarchicalBundleAdjuster::Optimize(
    const std::vector<Keyframe::Ptr>& keyframes, const SE3& T_imu_lidar) const {
    HierarchicalBundleAdjustmentSummary summary;
    summary.keyframe_count = keyframes.size();
    if (!options_.enabled) {
        summary.reason = "disabled";
        return summary;
    }
    if (keyframes.size() < static_cast<std::size_t>(std::max(2, options_.leaf_submap_size))) {
        summary.reason = "too_few_keyframes";
        return summary;
    }
    summary.attempted = true;
    const auto begin_time = std::chrono::steady_clock::now();
    std::vector<HierarchyNode> nodes = BuildLeafNodes(
        keyframes, T_imu_lidar, options_.leaf_submap_size, options_.max_points_per_submap);

    for (int level = 0; level < std::max(1, options_.max_levels) && nodes.size() >= 2; ++level) {
        ++summary.levels_attempted;
        BundleAdjustmentOptions optimizer_options = options_.optimizer;
        optimizer_options.voxel_size *= std::pow(std::max(1.0, options_.level_voxel_scale), level);
        optimizer_options.max_points_per_frame =
            std::min(options_.max_points_per_submap, optimizer_options.max_points_per_frame);
        optimizer_options.min_frames_per_voxel = 2;
        VoxelBundleAdjuster optimizer(optimizer_options);
        bool level_accepted = false;
        const std::size_t group_size = static_cast<std::size_t>(std::max(2, options_.branching_factor));
        for (std::size_t begin = 0; begin < nodes.size(); begin += group_size) {
            const std::size_t end = std::min(nodes.size(), begin + group_size);
            if (end - begin < 2) continue;
            std::vector<BundleFrame> frames;
            frames.reserve(end - begin);
            for (std::size_t index = begin; index < end; ++index) frames.push_back(nodes[index].frame);
            const BundleAdjustmentSummary level_summary = optimizer.Optimize(frames);
            summary.max_translation_update =
                std::max(summary.max_translation_update, level_summary.max_translation_update);
            summary.max_rotation_update_deg =
                std::max(summary.max_rotation_update_deg, level_summary.max_rotation_update_deg);
            if (!level_summary.accepted) continue;
            level_accepted = true;
            summary.accepted = true;
            ApplyCorrections(frames, nodes, begin, T_imu_lidar);
        }
        if (level_accepted) {
            ++summary.levels_accepted;
        }
        nodes = GroupNodes(nodes, T_imu_lidar, options_.branching_factor,
                           options_.max_points_per_submap);
    }

    // The lower levels deliberately use bounded temporal groups.  The remaining
    // coarse nodes must be optimized together so spatially overlapping parts of
    // a loop can share plane voxels even when they are far apart in time.  This
    // is the point-cloud counterpart of Voxel-SLAM HBA's global top layer.
    if (options_.global_top_level && nodes.size() >= 2) {
        summary.global_top_level_attempted = true;
        BundleAdjustmentOptions optimizer_options = options_.optimizer;
        optimizer_options.voxel_size *=
            std::pow(std::max(1.0, options_.level_voxel_scale),
                     std::max(0, summary.levels_attempted - 1));
        optimizer_options.max_points_per_frame =
            std::min(options_.max_points_per_submap, optimizer_options.max_points_per_frame);
        optimizer_options.min_frames_per_voxel = 2;
        VoxelBundleAdjuster optimizer(optimizer_options);
        std::vector<BundleFrame> frames;
        frames.reserve(nodes.size());
        for (const auto& node : nodes) frames.push_back(node.frame);
        const BundleAdjustmentSummary global_summary = optimizer.Optimize(frames);
        summary.max_translation_update =
            std::max(summary.max_translation_update, global_summary.max_translation_update);
        summary.max_rotation_update_deg =
            std::max(summary.max_rotation_update_deg, global_summary.max_rotation_update_deg);
        if (global_summary.accepted) {
            summary.global_top_level_accepted = true;
            summary.accepted = true;
            ApplyCorrections(frames, nodes, 0, T_imu_lidar);
        }
    }

    if (options_.final_local_refinement) {
        const std::size_t window = static_cast<std::size_t>(std::max(2, options_.leaf_submap_size));
        VoxelBundleAdjuster local_optimizer(options_.optimizer);
        for (std::size_t begin = 0; begin + 1 < keyframes.size(); begin += window) {
            const std::size_t end = std::min(keyframes.size(), begin + window);
            std::vector<Keyframe::Ptr> local(keyframes.begin() + begin, keyframes.begin() + end);
            const auto local_summary = local_optimizer.OptimizeKeyframes(local, T_imu_lidar, true, false);
            if (local_summary.accepted) {
                ++summary.final_local_windows_accepted;
                summary.accepted = true;
            }
        }
    }

    summary.elapsed_ms = std::chrono::duration<double, std::milli>(
                             std::chrono::steady_clock::now() - begin_time)
                             .count();
    summary.reason = summary.accepted ? "accepted" : "no_cost_decreasing_level";
    return summary;
}

}  // namespace lightning::backend
