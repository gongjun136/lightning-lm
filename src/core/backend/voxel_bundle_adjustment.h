#ifndef LIGHTNING_BACKEND_VOXEL_BUNDLE_ADJUSTMENT_H
#define LIGHTNING_BACKEND_VOXEL_BUNDLE_ADJUSTMENT_H

#include <cstddef>
#include <string>
#include <vector>

#include "common/keyframe.h"

namespace lightning::backend {

struct BundleAdjustmentOptions {
    bool enabled = true;
    int window_size = 10;
    int window_stride = 1;
    int max_iterations = 3;
    int max_points_per_frame = 4000;
    int min_points_per_voxel = 20;
    int min_frames_per_voxel = 2;
    int max_voxel_depth = 2;
    int num_threads = 4;
    double voxel_size = 1.0;
    double planarity_ratio = 0.10;
    double eigenvalue_floor = 1e-9;
    // LIO already fuses the built-in IMU. These information weights preserve
    // that motion estimate in directions that planar LiDAR geometry cannot
    // observe, playing the same stabilizing role as the IMU term in LI-BA.
    double correction_prior_rotation = 400.0;
    double correction_prior_translation = 100.0;
    double correction_smoothness_rotation = 400.0;
    double correction_smoothness_translation = 100.0;
    double max_rotation_step_deg = 2.0;
    double max_translation_step = 0.30;
    double max_total_rotation_deg = 10.0;
    double max_total_translation = 2.0;
    double initial_damping = 1e-2;
};

struct BundleFrame {
    unsigned long id = 0;
    double timestamp = 0.0;
    CloudPtr cloud;
    SE3 pose;
    // Optional immutable LIO/IMU reference. Sliding-window BA must not use the
    // previous window's optimized pose as a new prior, or corrections random-walk.
    SE3 prior_pose;
    bool has_prior_pose = false;
};

struct BundleAdjustmentSummary {
    bool attempted = false;
    bool accepted = false;
    std::size_t frame_count = 0;
    std::size_t factor_count = 0;
    std::size_t sampled_point_count = 0;
    int iterations = 0;
    double initial_cost = 0.0;
    double final_cost = 0.0;
    double max_translation_update = 0.0;
    double max_rotation_update_deg = 0.0;
    std::string reason;
};

// Plane-voxel LiDAR bundle adjustment adapted from Voxel-SLAM's LidarFactor
// and Lidar_BA_Optimizer. The first pose is fixed; all other poses are refined
// against the smallest-eigenvalue plane residual over shared voxels.
class VoxelBundleAdjuster {
   public:
    explicit VoxelBundleAdjuster(BundleAdjustmentOptions options = {});

    void SetOptions(const BundleAdjustmentOptions& options) { options_ = options; }
    const BundleAdjustmentOptions& GetOptions() const { return options_; }

    BundleAdjustmentSummary Optimize(std::vector<BundleFrame>& frames) const;

    BundleAdjustmentSummary OptimizeKeyframes(const std::vector<Keyframe::Ptr>& keyframes,
                                               const SE3& T_imu_lidar, bool apply_results,
                                               bool use_lio_prior = true) const;

   private:
    BundleAdjustmentOptions options_;
};

}  // namespace lightning::backend

#endif  // LIGHTNING_BACKEND_VOXEL_BUNDLE_ADJUSTMENT_H
