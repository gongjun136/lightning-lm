#ifndef LIGHTNING_BACKEND_HIERARCHICAL_BUNDLE_ADJUSTMENT_H
#define LIGHTNING_BACKEND_HIERARCHICAL_BUNDLE_ADJUSTMENT_H

#include <cstddef>
#include <string>
#include <vector>

#include "core/backend/voxel_bundle_adjustment.h"

namespace lightning::backend {

struct HierarchicalBundleAdjustmentOptions {
    bool enabled = true;
    int leaf_submap_size = 10;
    int branching_factor = 4;
    int max_levels = 3;
    int max_points_per_submap = 50000;
    double level_voxel_scale = 2.0;
    bool global_top_level = true;
    bool final_local_refinement = true;
    bool require_applied_loop_for_commit = true;
    BundleAdjustmentOptions optimizer;
};

struct HierarchicalBundleAdjustmentSummary {
    bool attempted = false;
    bool accepted = false;
    std::size_t keyframe_count = 0;
    int levels_attempted = 0;
    int levels_accepted = 0;
    bool global_top_level_attempted = false;
    bool global_top_level_accepted = false;
    int final_local_windows_accepted = 0;
    double max_translation_update = 0.0;
    double max_rotation_update_deg = 0.0;
    double elapsed_ms = 0.0;
    std::string reason;
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
    HierarchicalBundleAdjustmentOptions options_;
};

}  // namespace lightning::backend

#endif  // LIGHTNING_BACKEND_HIERARCHICAL_BUNDLE_ADJUSTMENT_H
