#pragma once

#include <deque>
#include <memory>
#include <string>
#include <vector>

#include "core/localization/global_relocalizer.h"
#include "core/localization/kiss_matcher_registration.h"
#include "core/localization/solid_descriptor.h"

namespace lightning::loc {

// Project map<-odom to a gravity-aligned yaw rotation while preserving the
// matched current lidar position exactly. This keeps the result independent
// of the absolute odometry translation after a long localization outage.
SE3 MakePlanarMapOdom(const SE3& T_map_current_lidar,
                      const SE3& T_odom_current_lidar);

class SolidRelocalizer : public GlobalRelocalizer {
   public:
    SolidRelocalizer();

    bool Init(const std::string& config_path, const std::string& map_path,
              const SE3& T_imu_lidar) override;
    void ResetQuery() override;
    bool IsReady() const override { return ready_; }
    std::size_t DatabaseSize() const override { return entries_.size(); }
    const char* Name() const override {
        return options_.coarse_registration_backend == "kiss_matcher"
                   ? "solid_kiss"
                   : "solid";
    }
    std::optional<RelocalizationResult> AddFrame(
        const CloudPtr& cloud, const SE3& T_odom_imu, double timestamp) override;

    // Build a SOLiD database from the optimized BTC submap clouds already
    // exported with a Lightning-LM map package.
    static bool BuildDatabase(const std::string& config_path,
                              const std::string& map_path);

   private:
    struct Options {
        std::vector<int> query_submap_sizes = {1, 3, 5};  ///< 用于查询的多尺度点云帧数列表。
        int query_stride = 1;  ///< 相邻查询的点云帧间隔。
        int top_k = 10;  ///< 精配阶段保留的候选数。
        int retrieval_pool_size = 200;  ///< 初步检索候选池的数量上限。
        double min_similarity = 0.50;  ///< 候选描述子相似度下限。
        double candidate_dedup_radius = 2.0;  ///< 候选位置去重半径（m）。
        double candidate_dedup_yaw_deg = 15.0;  ///< 候选航向去重角度（deg）。
        int max_points_per_submap = 50000;  ///< 查询子地图点数上限。
        int min_points_per_submap = 100;  ///< 有效查询需要的最少点数。
        double downsample_leaf_size = 0.20;  ///< 子地图体素降采样分辨率（m）。
        bool refine_with_icp = true;  ///< 是否启用候选粗配准细化流程。
        std::string coarse_registration_backend = "icp";  ///< 候选粗配准方法名称，如 ICP 或 KISS。
        KissMatcherRegistration::Options kiss_registration;  ///< KISS 粗配准参数或历史子地图的预处理配准器。
        int icp_max_iterations = 30;  ///< 候选 ICP 粗配准迭代次数上限。
        double icp_max_correspondence_distance = 5.0;  ///< 候选 ICP 允许的最大对应距离（m）。
        double icp_max_fitness_score = 2.0;  ///< 接受候选 ICP 结果的分值上限。
        double icp_max_translation_correction = 10.0;  ///< 候选粗配准允许的最大平移修正（m）。
        int icp_batch_size = 8;  ///< 一批执行粗配准的候选数。
        int icp_workers = 8;  ///< 候选粗配准工作线程数。
        int icp_worker_nice = 0;  ///< 候选工作线程的 nice 优先级配置。
        std::vector<int> icp_cpu_affinity;  ///< 候选工作线程允许使用的 CPU 编号。
        std::vector<double> icp_yaw_hypothesis_offsets_deg = {0.0};  ///< 相对描述子航向估计追加的航向假设（deg）。
        std::string database_subdirectory = "solid_relocalization";  ///< SOLiD 数据库在地图包中的相对目录。
        std::string source_database_subdirectory = "btc_relocalization";  ///< 构建 SOLiD 数据库使用的源子地图相对目录。
        std::string compute_backend = "cpu";  ///< 描述子计算后端名称。
        std::vector<int> debug_candidate_ids;  ///< 指定输出调试信息的候选 ID 列表。
    };

    struct DatabaseEntry {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        int descriptor_id = -1;  ///< 历史子地图的描述子数据库 ID。
        SE3 T_world_lidar;  ///< 历史子地图雷达在地图世界坐标系中的位姿。
        SolidDescriptor descriptor;  ///< 历史子地图的 SOLiD 描述子。
        std::string cloud_path;  ///< 历史子地图点云文件路径。
        pcl::PointCloud<pcl::PointXYZI>::Ptr cloud;  ///< 雷达坐标系下的查询或历史子地图点云。
        std::unique_ptr<KissMatcherRegistration> kiss_registration;  ///< KISS 粗配准参数或历史子地图的预处理配准器。
    };

    struct QueryFrame {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        CloudPtr cloud;  ///< 雷达坐标系下的查询或历史子地图点云。
        SE3 T_odom_lidar;  ///< 查询帧雷达在里程计坐标系中的位姿。
        double timestamp = 0.0;  ///< 查询帧的传感器时间戳（s）。
    };

    pcl::PointCloud<pcl::PointXYZI>::Ptr BuildQuerySubmap(
        std::size_t frame_count) const;
    bool RefineCandidateWithIcp(
        const pcl::PointCloud<pcl::PointXYZI>::ConstPtr& query,
        DatabaseEntry& entry, RelocalizationCandidate& candidate,
        double& fitness_score);
    bool RefineCandidateWithKissMatcher(
        const KissMatcherRegistration::PreparedCloud& query,
        DatabaseEntry& entry, RelocalizationCandidate& candidate,
        double& fitness_score);

    Options options_;  ///< SOLiD 检索、粗配准与算力配置。
    SolidDescriptorOptions descriptor_options_;  ///< SOLiD 距离、角度与高度分箱参数。
    SolidDescriptorEngine descriptor_engine_;  ///< 生成查询描述子并估计候选航向的引擎。
    std::vector<DatabaseEntry, Eigen::aligned_allocator<DatabaseEntry>> entries_;  ///< 历史描述子、位姿及候选点云索引。
    std::deque<QueryFrame, Eigen::aligned_allocator<QueryFrame>> query_frames_;  ///< 近期查询点云及里程计位姿历史。
    std::size_t frames_since_query_ = 0;  ///< 距离上次查询已积累的帧数。
    std::size_t icp_batch_cursor_ = 0;  ///< 候选粗配准批次的处理游标。
    SE3 T_imu_lidar_;  ///< 将雷达坐标变换到 IMU 坐标的固定外参。
    bool ready_ = false;  ///< 地图描述子数据库是否已就绪。
};

}  // namespace lightning::loc
