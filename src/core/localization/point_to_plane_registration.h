#pragma once

#include <cstddef>
#include <limits>
#include <memory>
#include <string>

#include <pcl/kdtree/kdtree_flann.h>

#include "common/eigen_types.h"
#include "common/point_def.h"

namespace lightning::loc {

class PointToPlaneRegistration {
   public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    struct Options {
        double target_voxel_size = 0.5;
        double source_voxel_size = 0.5;
        double plane_fit_threshold = 0.15;
        double coarse_max_correspondence_distance = 2.0;
        double fine_max_correspondence_distance = 0.75;
        double huber_delta = 0.2;
        int coarse_max_iterations = 10;
        int fine_max_iterations = 15;
        int min_matches = 500;
        double min_inlier_ratio = 0.30;
        double translation_convergence = 0.005;
        double rotation_convergence_deg = 0.05;
        double min_normalized_hessian_eigenvalue = 1e-4;
        double max_hessian_condition_number = 1e5;
        double max_translation_correction = 2.0;
        double max_rotation_correction_deg = 10.0;
        double max_rmse = std::numeric_limits<double>::infinity();
    };

    struct Result {
        SE3 attempted_pose;
        int iterations = 0;
        int matches = 0;
        double inlier_ratio = 0.0;
        double rmse = std::numeric_limits<double>::infinity();
        double translation_correction = 0.0;
        double rotation_correction_deg = 0.0;
        double hessian_condition_number = std::numeric_limits<double>::infinity();
        bool attempted_pose_valid = false;
        bool converged = false;
        bool success = false;
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
        Mat6d hessian = Mat6d::Zero();
        Vec6d gradient = Vec6d::Zero();
        int matches = 0;
        double squared_error = 0.0;
        double inlier_ratio = 0.0;
        double condition_number = std::numeric_limits<double>::infinity();
        bool observable = false;
    };

    CloudPtr DownsampleFinite(const CloudPtr& cloud, double voxel_size) const;
    Linearization Linearize(const CloudPtr& source, const SE3& pose,
                            double max_correspondence_distance) const;
    bool RunStage(const CloudPtr& source, double max_correspondence_distance,
                  int max_iterations, const SE3& initial_pose, SE3& pose,
                  Result& result) const;
    bool WithinCorrectionBounds(const SE3& initial_pose, const SE3& pose,
                                Result& result) const;

    Options options_;
    CloudPtr target_;
    std::unique_ptr<pcl::KdTreeFLANN<PointType>> target_kdtree_;
};

}  // namespace lightning::loc
