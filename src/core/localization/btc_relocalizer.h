#pragma once

#include <deque>
#include <optional>
#include <string>
#include <vector>

#include "common/eigen_types.h"
#include "common/point_def.h"
#include "core/backend/third_party/voxel_slam_btc/BTC.h"

namespace lightning::loc {

struct BtcRelocalizationResult {
    bool attempted = false;
    bool candidate_found = false;
    bool accepted = false;
    int candidate_id = -1;
    double score = 0.0;
    double timestamp = 0.0;
    std::size_t point_count = 0;
    std::size_t descriptor_count = 0;
    SE3 T_world_imu;
    std::string reason;
};

// Global place recognition against BTC submaps saved together with a map.
// Incoming scans are accumulated in their current lidar frame using lidar
// odometry, then BTC supplies a coarse SE(3) pose for NDT verification.
class BtcRelocalizer {
   public:
    BtcRelocalizer();

    bool Init(const std::string& config_path, const std::string& map_path,
              const SE3& T_imu_lidar);
    void ResetQuery();
    bool IsReady() const { return ready_; }
    std::size_t DatabaseSize() const { return entries_.size(); }

    std::optional<BtcRelocalizationResult> AddFrame(
        const CloudPtr& cloud, const SE3& T_odom_imu, double timestamp);

   private:
    struct Options {
        bool enabled = false;
        int query_submap_size = 10;
        int max_points_per_submap = 50000;
        int min_points_per_submap = 100;
        double downsample_leaf_size = 0.20;
        double min_btc_score = 0.45;
        std::string database_subdirectory = "btc_relocalization";
    };

    struct DatabaseEntry {
        int descriptor_id = -1;
        SE3 T_world_lidar;
    };

    struct QueryFrame {
        CloudPtr cloud;
        SE3 T_odom_lidar;
        double timestamp = 0.0;
    };

    pcl::PointCloud<pcl::PointXYZI>::Ptr BuildQuerySubmap() const;

    Options options_;
    ConfigSetting descriptor_config_;
    STDescManager manager_;
    std::vector<DatabaseEntry> entries_;
    std::deque<QueryFrame> query_frames_;
    SE3 T_imu_lidar_;
    bool ready_ = false;
};

}  // namespace lightning::loc
