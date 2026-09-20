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
        bool enabled = false;
        int query_submap_size = 10;
        std::vector<int> query_submap_sizes = {5, 10, 20};
        int query_stride = 1;
        int top_k = 10;
        double candidate_dedup_radius = 2.0;
        double candidate_dedup_yaw_deg = 15.0;
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

    pcl::PointCloud<pcl::PointXYZI>::Ptr BuildQuerySubmap(std::size_t frame_count) const;

    Options options_;
    ConfigSetting descriptor_config_;
    STDescManager manager_;
    std::vector<DatabaseEntry> entries_;
    std::deque<QueryFrame> query_frames_;
    std::size_t frames_since_query_ = 0;
    SE3 T_imu_lidar_;
    bool ready_ = false;
};

}  // namespace lightning::loc
