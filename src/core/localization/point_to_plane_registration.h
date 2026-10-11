#pragma once

#include <cstddef>
#include <limits>
#include <memory>
#include <string>

#include <pcl/kdtree/kdtree_flann.h>

#include "common/eigen_types.h"
#include "common/point_def.h"

namespace lightning::loc {

/** @see @ref lidar_residuals */
class PointToPlaneRegistration {
   public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    struct Options {
        double target_voxel_size = 0.5;  ///< 目标地图体素降采样分辨率（m）。
        double source_voxel_size = 0.5;  ///< 查询点云体素降采样分辨率（m）。
        double plane_fit_threshold = 0.15;  ///< 局部平面拟合允许的距离误差（m）。
        double coarse_max_correspondence_distance = 2.0;  ///< 粗精配阶段的最大近邻距离（m）。
        double fine_max_correspondence_distance = 0.75;  ///< 细精配阶段的最大近邻距离（m）。
        double huber_delta = 0.2;  ///< 点面残差鲁棒核的尺度（m）。
        int coarse_max_iterations = 10;  ///< 粗精配阶段的迭代次数上限。
        int fine_max_iterations = 15;  ///< 细精配阶段的迭代次数上限。
        int min_matches = 500;  ///< 接受精配所需的最少点面匹配数。
        double min_inlier_ratio = 0.30;  ///< 接受精配所需的最小内点比例。
        double translation_convergence = 0.005;  ///< 判断迭代收敛的平移变化阈值（m）。
        double rotation_convergence_deg = 0.05;  ///< 判断迭代收敛的旋转变化阈值（deg）。
        double min_normalized_hessian_eigenvalue = 1e-4;  ///< 归一化 Hessian 最小特征值的可观测性门槛。
        double max_hessian_condition_number = 1e5;  ///< 接受精配所允许的 Hessian 条件数上限。
        double max_translation_correction = 2.0;  ///< 相对初始猜测允许的最大平移修正（m）。
        double max_rotation_correction_deg = 10.0;  ///< 相对初始猜测允许的最大旋转修正（deg）。
        double max_rmse = std::numeric_limits<double>::infinity();  ///< 接受精配所允许的点面均方根误差（m）。
    };

    struct Result {
        SE3 attempted_pose;  ///< 本次尝试得到的精配位姿，失败时可用于诊断。
        int iterations = 0;  ///< 实际执行的精配迭代次数。
        int matches = 0;  ///< 进入线性化或最终精配的点面匹配数。
        double inlier_ratio = 0.0;  ///< 有效点面匹配占查询点的比例。
        double rmse = std::numeric_limits<double>::infinity();  ///< 最终点面残差的均方根（m）。
        double translation_correction = 0.0;  ///< 相对初始猜测的平移修正幅度（m）。
        double rotation_correction_deg = 0.0;  ///< 相对初始猜测的旋转修正幅度（deg）。
        double hessian_condition_number = std::numeric_limits<double>::infinity();  ///< 最终观测 Hessian 的条件数。
        bool attempted_pose_valid = false;  ///< 尝试位姿是否有限且可用于诊断。
        bool converged = false;  ///< 精配迭代是否达到收敛条件。
        bool success = false;  ///< 精配结果是否通过全部质量门控。
    };

    PointToPlaneRegistration();
    explicit PointToPlaneRegistration(Options options);

    static bool ValidateOptions(const Options& options, std::string* error = nullptr);

    bool SetTarget(const CloudPtr& target);
    bool HasTarget() const;
    std::size_t TargetSize() const;

    // On failure, pose is left exactly as supplied by the caller.
    Result Refine(const CloudPtr& source, SE3& pose) const;

   private:
    struct Linearization {
        Mat6d hessian = Mat6d::Zero();  ///< 当前迭代累计的六维位姿观测 Hessian。
        Vec6d gradient = Vec6d::Zero();  ///< 当前迭代累计的六维位姿梯度。
        int matches = 0;  ///< 进入线性化或最终精配的点面匹配数。
        double squared_error = 0.0;  ///< 有效点面残差平方和。
        double inlier_ratio = 0.0;  ///< 有效点面匹配占查询点的比例。
        double condition_number = std::numeric_limits<double>::infinity();  ///< 当前线性化 Hessian 的条件数。
        bool observable = false;  ///< 当前点面几何是否通过可观测性门控。
    };

    CloudPtr DownsampleFinite(const CloudPtr& cloud, double voxel_size) const;
    Linearization Linearize(const CloudPtr& source, const SE3& pose,
                            double max_correspondence_distance) const;
    bool RunStage(const CloudPtr& source, double max_correspondence_distance,
                  int max_iterations, const SE3& initial_pose, SE3& pose,
                  Result& result) const;
    bool WithinCorrectionBounds(const SE3& initial_pose, const SE3& pose,
                                Result& result) const;

    Options options_;  ///< 点面精配的采样、迭代、可观测性与修正门控参数。
    CloudPtr target_;  ///< 已经过滤和降采样的目标地图点云。
    std::unique_ptr<pcl::KdTreeFLANN<PointType>> target_kdtree_;  ///< 目标地图的最近邻索引。
};

}  // namespace lightning::loc
