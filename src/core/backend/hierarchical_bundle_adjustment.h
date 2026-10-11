#ifndef LIGHTNING_BACKEND_HIERARCHICAL_BUNDLE_ADJUSTMENT_H
#define LIGHTNING_BACKEND_HIERARCHICAL_BUNDLE_ADJUSTMENT_H

#include <cstddef>
#include <string>
#include <vector>

#include "core/backend/voxel_bundle_adjustment.h"

namespace lightning::backend {

struct HierarchicalBundleAdjustmentOptions {
    bool enabled = true;  ///< 是否启用分层 BA。
    int leaf_submap_size = 10;  ///< 叶子子地图包含的关键帧数。
    int branching_factor = 4;  ///< 每个上层子地图聚合的下层子地图数。
    int max_levels = 3;  ///< 分层优化层数上限。
    int max_points_per_submap = 50000;  ///< 每个子地图进入优化的点数上限。
    double level_voxel_scale = 2.0;  ///< 随层级增大的体素尺度倍数。
    bool global_top_level = true;  ///< 是否在顶层执行全局优化。
    bool final_local_refinement = true;  ///< 全局修正后是否再执行局部窗口精化。
    bool require_applied_loop_for_commit = true;  ///< 提交全局修正前是否要求已有实际应用的回环。
    BundleAdjustmentOptions optimizer;  ///< 各层复用的体素 BA 优化配置。
};

struct HierarchicalBundleAdjustmentSummary {
    bool attempted = false;  ///< 是否实际尝试分层优化。
    bool accepted = false;  ///< 分层优化结果是否被接受。
    std::size_t keyframe_count = 0;  ///< 本次参与分层优化的关键帧数。
    int levels_attempted = 0;  ///< 实际尝试的优化层数。
    int levels_accepted = 0;  ///< 接受优化结果的层数。
    bool global_top_level_attempted = false;  ///< 是否尝试顶层全局优化。
    bool global_top_level_accepted = false;  ///< 是否接受顶层全局优化结果。
    int final_local_windows_accepted = 0;  ///< 最终精化中接受结果的局部窗口数。
    double max_translation_update = 0.0;  ///< 各关键帧中最大的平移修正（m）。
    double max_rotation_update_deg = 0.0;  ///< 各关键帧中最大的旋转修正（deg）。
    double elapsed_ms = 0.0;  ///< 本次分层优化总耗时（ms）。
    std::string reason;  ///< 本次优化的接受或拒绝原因。
};

// Multi-resolution plane-voxel BA. Leaf scan groups are optimized first,
// grouped recursively into coarser submaps, and their corrections are applied
// top-down to the member keyframes before an optional local refinement pass.
class HierarchicalBundleAdjuster {
   public:
    explicit HierarchicalBundleAdjuster(HierarchicalBundleAdjustmentOptions options = {});

    HierarchicalBundleAdjustmentSummary Optimize(const std::vector<Keyframe::Ptr>& keyframes,
                                                  const SE3& T_imu_lidar) const;
    const HierarchicalBundleAdjustmentOptions& GetOptions() const { return options_; }

   private:
    HierarchicalBundleAdjustmentOptions options_;  ///< 分层子地图构建、优化及提交配置。
};

}  // namespace lightning::backend

#endif  // LIGHTNING_BACKEND_HIERARCHICAL_BUNDLE_ADJUSTMENT_H
