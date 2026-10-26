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
        double voxel_size = 0.50;  ///< 特征提取前的体素降采样分辨率（m）。
        double normal_radius = 1.75;  ///< 估计局部法向的搜索半径（m）。
        double feature_radius = 2.50;  ///< 计算 FPFH 特征的搜索半径（m）。
        double descriptor_ratio = 0.90;  ///< 特征最近邻匹配的比值筛选阈值。
        int max_feature_points = 60000;  ///< 参与特征计算的点数上限。
        int max_correspondences = 800;  ///< 进入一致性图与鲁棒求解的对应数上限。
        double compatibility_tolerance = 1.00;  ///< 判断点对距离一致性的容差（m）。
        double minimum_pair_distance = 1.00;  ///< 进入几何一致性判断的最小点对间距（m）。
        int minimum_core_degree = 4;  ///< 保留一致性图顶点所需的最小核心度。
        int minimum_inliers = 12;  ///< 接受配准所需的最少内点数。
        double inlier_threshold = 0.75;  ///< 配准内点的距离阈值（m）。
        double maximum_rmse = 0.50;  ///< 接受配准所允许的均方根误差（m）。
        int gnc_max_iterations = 50;  ///< GNC 鲁棒估计的迭代次数上限。
        double gnc_factor = 1.40;  ///< GNC 逐步调整鲁棒性参数的倍率。
        double spatial_cell_size = 5.0;  ///< 评估内点空间覆盖的网格边长（m）。
        int minimum_spatial_cells = 4;  ///< 接受配准所需的最少覆盖网格数。
        int feature_threads = 8;  ///< 特征计算使用的并行线程数。
    };

    struct Correspondence {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Vec3d source = Vec3d::Zero();  ///< 源点云中参与对应的点坐标（m）。
        Vec3d target = Vec3d::Zero();  ///< 目标点云中参与对应的点坐标（m）。
        double descriptor_distance = 0.0;  ///< 该对应的特征距离。
    };

    struct Result {
        bool success = false;  ///< 粗配准是否通过质量门控。
        bool converged = false;  ///< 达到迭代收敛条件，或完成末端统计且 iterations>0；接受以 success 为准。
        SE3 T_target_source;  ///< 将源点云坐标变换到目标点云坐标的配准结果。
        std::size_t source_points = 0;  ///< 查询点云预处理后的点数。
        std::size_t source_descriptors = 0;  ///< 查询点云生成的特征数量。
        std::size_t rough_correspondences = 0;  ///< 特征检索建立的初始对应数。
        std::size_t core_correspondences = 0;  ///< 一致性图剪枝后保留的对应数。
        std::size_t inliers = 0;  ///< 末次内点重拟合前选出的内点数，重拟合后不重新分类。
        std::size_t spatial_cells = 0;  ///< 所选内点在目标坐标系中覆盖的 XY 网格数。
        double inlier_ratio = 0.0;  ///< 所选内点数 / 兼容图保留的 core 对应数。
        double spatial_coverage = 0.0;  ///< 内点占据的 target XY 网格数 / core 占据的网格数。
        double rmse = 0.0;  ///< 最终内点配准均方根误差（m）。
        double score = 0.0;  ///< 由配准质量计算的候选分值。
        int iterations = 0;  ///< 实际执行的鲁棒求解迭代次数。
        std::string reason;  ///< 配准成功或失败的原因。
    };

    // A query cloud's normals and FPFH descriptors do not depend on the
    // candidate target.  Preparing them once avoids repeating the dominant
    // feature extraction work for every SOLiD shortlist entry.
    struct PreparedCloud {
        pcl::PointCloud<pcl::PointXYZ>::Ptr points{
            new pcl::PointCloud<pcl::PointXYZ>};  ///< 预处理点云，可供多个候选复用。
        pcl::PointCloud<pcl::FPFHSignature33>::Ptr descriptors{
            new pcl::PointCloud<pcl::FPFHSignature33>};  ///< 与预处理点云对应的 FPFH 特征。
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

    Options options_;  ///< FPFH 对应、图剪枝与 GNC 配准参数。
    pcl::PointCloud<pcl::PointXYZ>::Ptr target_points_;  ///< 预处理后的目标点云。
    pcl::PointCloud<pcl::FPFHSignature33>::Ptr target_features_;  ///< 与目标点对应的 FPFH 描述子。
    pcl::KdTreeFLANN<pcl::FPFHSignature33>::Ptr target_feature_tree_;  ///< 目标 FPFH 特征的最近邻检索树。
};

}  // namespace lightning::loc
