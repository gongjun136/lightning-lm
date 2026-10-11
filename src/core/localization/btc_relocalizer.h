#pragma once

#include <deque>
#include <optional>
#include <string>
#include <vector>

#include "common/eigen_types.h"
#include "common/point_def.h"
#include "core/backend/third_party/voxel_slam_btc/BTC.h"
#include "core/localization/global_relocalizer.h"

namespace lightning::loc {

using BtcRelocalizationCandidate = RelocalizationCandidate;
using BtcRelocalizationResult = RelocalizationResult;

// Global place recognition against BTC submaps saved together with a map.
// Incoming scans are accumulated in their current lidar frame using lidar
// odometry, then BTC supplies a coarse SE(3) pose for NDT verification.
class BtcRelocalizer : public GlobalRelocalizer {
   public:
    BtcRelocalizer();

    bool Init(const std::string& config_path, const std::string& map_path,
              const SE3& T_imu_lidar) override;
    void ResetQuery() override;
    bool IsReady() const override { return ready_; }
    std::size_t DatabaseSize() const override { return entries_.size(); }
    const char* Name() const override { return "btc"; }

    std::optional<RelocalizationResult> AddFrame(
        const CloudPtr& cloud, const SE3& T_odom_imu, double timestamp) override;

   private:
    struct Options {
        bool enabled = false;  ///< 是否启用 BTC 全局重定位。
        int query_submap_size = 10;  ///< 默认查询子地图的点云帧数。
        std::vector<int> query_submap_sizes = {5, 10, 20};  ///< 用于多尺度检索的查询帧数列表。
        int query_stride = 1;  ///< 相邻全局查询的点云帧间隔。
        int top_k = 10;  ///< 最多返回的候选数量。
        double candidate_dedup_radius = 2.0;  ///< 候选位置去重半径（m）。
        double candidate_dedup_yaw_deg = 15.0;  ///< 候选航向去重角度（deg）。
        int max_points_per_submap = 50000;  ///< 查询子地图点数上限。
        int min_points_per_submap = 100;  ///< 触发有效查询所需的最少点数。
        double downsample_leaf_size = 0.20;  ///< 查询子地图体素分辨率（m）。
        double min_btc_score = 0.45;  ///< 返回候选所需的最低 BTC 分值。
        std::string database_subdirectory = "btc_relocalization";  ///< 地图包中 BTC 描述子数据库的相对目录。
    };

    struct DatabaseEntry {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        int descriptor_id = -1;  ///< 地图数据库中描述子的 ID。
        SE3 T_world_lidar;  ///< 历史子地图雷达在地图世界坐标系中的位姿。
    };

    struct QueryFrame {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        CloudPtr cloud;  ///< 查询帧的雷达坐标系点云。
        SE3 T_odom_lidar;  ///< 查询帧雷达在里程计坐标系中的位姿。
        double timestamp = 0.0;  ///< 查询帧的传感器时间戳（s）。
    };

    pcl::PointCloud<pcl::PointXYZI>::Ptr BuildQuerySubmap(std::size_t frame_count) const;

    Options options_;  ///< BTC 全局重定位查询配置。
    ConfigSetting descriptor_config_;  ///< BTC 描述子生成和检索参数。
    STDescManager manager_;  ///< 加载地图描述子并执行 BTC 检索的管理器。
    std::vector<DatabaseEntry, Eigen::aligned_allocator<DatabaseEntry>> entries_;  ///< 地图描述子 ID 与雷达世界位姿的索引。
    std::deque<QueryFrame, Eigen::aligned_allocator<QueryFrame>> query_frames_;  ///< 等待构成查询子地图的近期点云与里程计位姿。
    std::size_t frames_since_query_ = 0;  ///< 距离上次查询已积累的帧数。
    SE3 T_imu_lidar_;  ///< 将雷达坐标变换到 IMU 坐标的固定外参。
    bool ready_ = false;  ///< 地图数据库是否已加载并可查询。
};

}  // namespace lightning::loc
