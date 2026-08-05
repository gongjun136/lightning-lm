#pragma once

#include <deque>
#include <string>
#include <vector>

#include "core/localization/global_relocalizer.h"
#include "core/localization/solid_descriptor.h"

namespace lightning::loc {

class SolidRelocalizer : public GlobalRelocalizer {
   public:
    SolidRelocalizer();

    bool Init(const std::string& config_path, const std::string& map_path,
              const SE3& T_imu_lidar) override;
    void ResetQuery() override;
    bool IsReady() const override { return ready_; }
    std::size_t DatabaseSize() const override { return entries_.size(); }
    const char* Name() const override { return "solid"; }
    std::optional<RelocalizationResult> AddFrame(
        const CloudPtr& cloud, const SE3& T_odom_imu, double timestamp) override;

    // Build a SOLiD database from the optimized BTC submap clouds already
    // exported with a Lightning-LM map package.
    static bool BuildDatabase(const std::string& config_path,
                              const std::string& map_path);

   private:
    struct Options {
        std::vector<int> query_submap_sizes = {1, 3, 5};
        int query_stride = 1;
        int top_k = 10;
        int retrieval_pool_size = 200;
        double min_similarity = 0.50;
        double candidate_dedup_radius = 2.0;
        double candidate_dedup_yaw_deg = 15.0;
        int max_points_per_submap = 50000;
        int min_points_per_submap = 100;
        double downsample_leaf_size = 0.20;
        bool refine_with_icp = true;
        int icp_max_iterations = 30;
        double icp_max_correspondence_distance = 5.0;
        double icp_max_fitness_score = 2.0;
        double icp_max_translation_correction = 10.0;
        int icp_batch_size = 8;
        int icp_workers = 8;
        std::vector<double> icp_yaw_hypothesis_offsets_deg = {0.0};
        std::string database_subdirectory = "solid_relocalization";
        std::string source_database_subdirectory = "btc_relocalization";
        std::string compute_backend = "cpu";
        std::vector<int> debug_candidate_ids;
    };

    struct DatabaseEntry {
        int descriptor_id = -1;
        SE3 T_world_lidar;
        SolidDescriptor descriptor;
        std::string cloud_path;
        pcl::PointCloud<pcl::PointXYZI>::Ptr cloud;
    };

    struct QueryFrame {
        CloudPtr cloud;
        SE3 T_odom_lidar;
        double timestamp = 0.0;
    };

    pcl::PointCloud<pcl::PointXYZI>::Ptr BuildQuerySubmap(
        std::size_t frame_count) const;
    bool RefineCandidateWithIcp(
        const pcl::PointCloud<pcl::PointXYZI>::ConstPtr& query,
        DatabaseEntry& entry, RelocalizationCandidate& candidate,
        double& fitness_score);

    Options options_;
    SolidDescriptorOptions descriptor_options_;
    SolidDescriptorEngine descriptor_engine_;
    std::vector<DatabaseEntry> entries_;
    std::deque<QueryFrame> query_frames_;
    std::size_t frames_since_query_ = 0;
    std::size_t icp_batch_cursor_ = 0;
    SE3 T_imu_lidar_;
    bool ready_ = false;
};

}  // namespace lightning::loc
