#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include "common/eigen_types.h"
#include "common/point_def.h"

namespace lightning::loc {

// Initial-guess-free feature registration inspired by KISS-Matcher.  The
// implementation intentionally uses the PCL/Eigen stack already shipped with
// Lightning-LM: FPFH matching, maximum k-core pruning, and GNC-TLS pose
// estimation.  Fine alignment remains the responsibility of the existing
// relocalization refinement backend.
class KissMatcherRegistration {
   public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    struct Options {
        double voxel_size = 0.50;
        double normal_radius = 1.75;
        double feature_radius = 2.50;
        double descriptor_ratio = 0.90;
        int max_feature_points = 60000;
        int max_correspondences = 800;
        double compatibility_tolerance = 1.00;
        double minimum_pair_distance = 1.00;
        int minimum_core_degree = 4;
        int minimum_inliers = 12;
        double inlier_threshold = 0.75;
        double maximum_rmse = 0.50;
        int gnc_max_iterations = 50;
        double gnc_factor = 1.40;
        double spatial_cell_size = 5.0;
        int minimum_spatial_cells = 4;
        int feature_threads = 8;
    };

    struct Correspondence {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Vec3d source = Vec3d::Zero();
        Vec3d target = Vec3d::Zero();
        double descriptor_distance = 0.0;
    };

    struct Result {
        bool success = false;
        bool converged = false;
        SE3 T_target_source;
        std::size_t source_points = 0;
        std::size_t source_descriptors = 0;
        std::size_t rough_correspondences = 0;
        std::size_t core_correspondences = 0;
        std::size_t inliers = 0;
        std::size_t spatial_cells = 0;
        double inlier_ratio = 0.0;
        double spatial_coverage = 0.0;
        double rmse = 0.0;
        double score = 0.0;
        int iterations = 0;
        std::string reason;
    };

    // A query cloud's normals and FPFH descriptors do not depend on the
    // candidate target.  Preparing them once avoids repeating the dominant
    // feature extraction work for every SOLiD shortlist entry.
    struct PreparedCloud {
        pcl::PointCloud<pcl::PointXYZ>::Ptr points{
            new pcl::PointCloud<pcl::PointXYZ>};
        pcl::PointCloud<pcl::FPFHSignature33>::Ptr descriptors{
            new pcl::PointCloud<pcl::FPFHSignature33>};
    };

    KissMatcherRegistration();
    explicit KissMatcherRegistration(Options options);

    static bool ValidateOptions(const Options& options, std::string* error = nullptr);

    bool SetTarget(const CloudPtr& target, std::string* error = nullptr);
    bool PrepareSource(const CloudPtr& source, PreparedCloud& prepared,
                       std::string* error = nullptr) const;
    bool IsReady() const { return target_features_ && !target_features_->empty(); }
    std::size_t TargetPointCount() const {
        return target_points_ ? target_points_->size() : 0;
    }

    Result Align(const CloudPtr& source) const;
    Result Align(const PreparedCloud& source) const;

    // Public to make the robust estimator independently testable and useful
    // for future descriptor front-ends such as DualReg-style anchors.
    static Result SolveCorrespondences(
        const std::vector<Correspondence,
                          Eigen::aligned_allocator<Correspondence>>& correspondences,
        const Options& options);

   private:
    bool BuildFeatures(const CloudPtr& input, PreparedCloud& output,
                       std::string* error) const;
    std::vector<Correspondence, Eigen::aligned_allocator<Correspondence>>
    MatchFeatures(const PreparedCloud& source) const;

    Options options_;
    pcl::PointCloud<pcl::PointXYZ>::Ptr target_points_;
    pcl::PointCloud<pcl::FPFHSignature33>::Ptr target_features_;
    pcl::KdTreeFLANN<pcl::FPFHSignature33>::Ptr target_feature_tree_;
};

}  // namespace lightning::loc
