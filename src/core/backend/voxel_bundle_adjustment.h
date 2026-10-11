#ifndef LIGHTNING_BACKEND_VOXEL_BUNDLE_ADJUSTMENT_H
#define LIGHTNING_BACKEND_VOXEL_BUNDLE_ADJUSTMENT_H

#include <cstddef>
#include <string>
#include <vector>

#include "common/keyframe.h"

namespace lightning::backend {

struct BundleAdjustmentOptions {
    bool enabled = true;  ///< 是否启用体素 BA。
    int window_size = 10;  ///< 局部优化窗口包含的关键帧数。
    int window_stride = 1;  ///< 相邻优化窗口的关键帧滑动步长。
    int max_iterations = 3;  ///< 优化迭代次数上限。
    int max_points_per_frame = 4000;  ///< 每帧进入体素优化的采样点数上限。
    int min_points_per_voxel = 20;  ///< 形成体素平面因子所需的最少点数。
    int min_frames_per_voxel = 2;  ///< 形成体素平面因子所需的最少不同帧数。
    int max_voxel_depth = 2;  ///< 体素递归细分深度上限。
    int num_threads = 4;  ///< 体素因子构建与计算的并行线程数。
    double voxel_size = 1.0;  ///< 基础体素边长（m）。
    double planarity_ratio = 0.10;  ///< 判定体素点集为平面的特征值比阈值。
    double eigenvalue_floor = 1e-9;  ///< 特征值数值稳定下限。
    // LIO already fuses the built-in IMU. These information weights preserve
    // that motion estimate in directions that planar LiDAR geometry cannot
    // observe, playing the same stabilizing role as the IMU term in LI-BA.
    double correction_prior_rotation = 400.0;  ///< 保持 LIO 旋转参考的修正先验信息权重。
    double correction_prior_translation = 100.0;  ///< 保持 LIO 平移参考的修正先验信息权重。
    double correction_smoothness_rotation = 400.0;  ///< 相邻帧旋转修正的平滑信息权重。
    double correction_smoothness_translation = 100.0;  ///< 相邻帧平移修正的平滑信息权重。
    double max_rotation_step_deg = 2.0;  ///< 单次迭代旋转修正上限（deg）。
    double max_translation_step = 0.30;  ///< 单次迭代平移修正上限（m）。
    double max_total_rotation_deg = 10.0;  ///< 整体优化旋转修正上限（deg）。
    double max_total_translation = 2.0;  ///< 整体优化平移修正上限（m）。
    double initial_damping = 1e-2;  ///< 优化器初始阻尼系数。
};

struct BundleFrame {
    unsigned long id = 0;  ///< 待优化帧的关键帧 ID。
    double timestamp = 0.0;  ///< 待优化帧的传感器时间戳（s）。
    CloudPtr cloud;  ///< 雷达坐标系下的该帧点云。
    SE3 pose;  ///< 待优化和回写的帧位姿。
    // Optional immutable LIO/IMU reference. Sliding-window BA must not use the
    // previous window's optimized pose as a new prior, or corrections random-walk.
    SE3 prior_pose;  ///< 不可变的 LIO 参考位姿，避免把上一窗口优化结果当作新先验。
    bool has_prior_pose = false;  ///< 该帧是否提供独立 LIO 参考位姿。
};

struct BundleAdjustmentSummary {
    bool attempted = false;  ///< 是否实际尝试本次优化。
    bool accepted = false;  ///< 优化结果是否通过改善与修正幅度门控。
    std::size_t frame_count = 0;  ///< 参与优化的帧数。
    std::size_t factor_count = 0;  ///< 构建的体素平面因子数。
    std::size_t sampled_point_count = 0;  ///< 进入优化的采样点总数。
    int iterations = 0;  ///< 实际执行的迭代次数。
    double initial_cost = 0.0;  ///< 优化前的目标函数值。
    double final_cost = 0.0;  ///< 优化后的目标函数值。
    double max_translation_update = 0.0;  ///< 各帧中最大的平移修正（m）。
    double max_rotation_update_deg = 0.0;  ///< 各帧中最大的旋转修正（deg）。
    std::string reason;  ///< 本次优化的接受或拒绝原因。
};

// Plane-voxel LiDAR bundle adjustment adapted from Voxel-SLAM's LidarFactor
// and Lidar_BA_Optimizer. The first pose is fixed; all other poses are refined
// against the smallest-eigenvalue plane residual over shared voxels.
/** @see @ref backend_optimization */
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
    BundleAdjustmentOptions options_;  ///< 体素 BA 优化及修正门控参数。
};

}  // namespace lightning::backend

#endif  // LIGHTNING_BACKEND_VOXEL_BUNDLE_ADJUSTMENT_H
