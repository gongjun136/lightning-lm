#ifndef LIGHTNING_BACKEND_BTC_LOOP_DETECTOR_H
#define LIGHTNING_BACKEND_BTC_LOOP_DETECTOR_H

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "common/keyframe.h"
#include "core/backend/third_party/voxel_slam_btc/BTC.h"

namespace lightning::backend {

struct BtcLoopDetectorOptions {
    bool enabled = true;
    int descriptor_submap_size = 10;
    int max_points_per_submap = 100000;
    int min_points_per_submap = 100;
    double downsample_leaf_size = 0.20;
    double min_loop_score = 0.45;
    double max_drift_ratio = 0.05;
    double max_rotation_correction_deg = 45.0;
    bool refine_with_plane_icp = true;
    int plane_icp_iterations = 20;
    int plane_icp_min_matches = 20;
    double plane_icp_min_observability = 14.0;
    bool allow_degenerate_plane_icp = true;
    double degenerate_min_loop_score = 0.40;
    int degenerate_min_matches = 100;
    int confirmation_count = 2;
    int confirmation_max_current_gap = 2;
    int confirmation_max_history_gap = 2;
    int loop_cooldown_descriptors = 10;
    // Recover revisits that are close in LIO odometry but whose BTC
    // descriptor is absent or below the global retrieval threshold. This is
    // candidate recall only; point-to-plane verification remains mandatory.
    bool enable_odom_revisit_fallback = false;
    double odom_revisit_search_radius = 5.0;
    double odom_revisit_min_journey = 50.0;
    int odom_revisit_degenerate_min_matches = 60;
    double odom_revisit_max_drift_ratio = 0.05;
    double max_odom_revisit_distance = 0.0;
    double min_optimization_translation = 0.20;
    double min_optimization_rotation_deg = 2.0;

    ConfigSetting descriptor;
};

struct BtcDescriptorEntry {
    int descriptor_id = -1;
    unsigned long first_keyframe_id = 0;
    unsigned long last_keyframe_id = 0;
    double timestamp = 0.0;
    double journey = 0.0;
    Keyframe::Ptr endpoint;
};

struct BtcLoopResult {
    bool descriptor_generated = false;
    bool candidate_found = false;
    bool accepted = false;
    bool optimization_warranted = false;
    int current_descriptor_id = -1;
    int history_descriptor_id = -1;
    unsigned long current_keyframe_id = 0;
    unsigned long history_keyframe_id = 0;
    double current_timestamp = 0.0;
    double history_timestamp = 0.0;
    double score = 0.0;
    double drift_translation = 0.0;
    double drift_rotation_deg = 0.0;
    double journey_span = 0.0;
    double drift_ratio = 0.0;
    double odom_revisit_distance = 0.0;
    double plane_icp_observability = 0.0;
    int plane_icp_matches = 0;
    bool plane_icp_converged = false;
    bool used_degenerate_plane_fallback = false;
    int confirmation_count = 0;
    std::size_t point_count = 0;
    std::size_t descriptor_count = 0;
    double generation_time_ms = 0.0;
    double search_time_ms = 0.0;
    std::string candidate_source = "none";
    std::string rejection_reason;
    SE3 T_history_lidar_current_lidar;
};

// BTC place recognition adapter. As in Voxel-SLAM, one descriptor is built
// from a non-overlapping aggregate of locally optimized scans instead of from
// a single scan. BTC remains the global retrieval path. An optional odometry
// proximity path can recover a start/end revisit, but it never bypasses point-
// cloud verification or the existing safety gates.
class BtcLoopDetector {
   public:
    explicit BtcLoopDetector(BtcLoopDetectorOptions options = {});

    void Reset();
    std::optional<BtcLoopResult> AddKeyframe(const Keyframe::Ptr& keyframe,
                                             const SE3& T_imu_lidar);
    BtcLoopResult ProcessSubmap(const std::vector<Keyframe::Ptr>& keyframes,
                                const SE3& T_imu_lidar);
    bool SaveRelocalizationDatabase(const std::string& directory,
                                    const SE3& T_imu_lidar) const;

    const BtcLoopDetectorOptions& GetOptions() const { return options_; }
    const std::vector<BtcDescriptorEntry>& Entries() const { return entries_; }
    std::size_t PendingKeyframes() const { return pending_keyframes_.size(); }

   private:
    struct RefineSummary {
        bool accepted = false;
        bool converged = false;
        double observability = 0.0;
        int matches = 0;
    };

    pcl::PointCloud<pcl::PointXYZI>::Ptr BuildSubmap(
        const std::vector<Keyframe::Ptr>& keyframes, const SE3& T_imu_lidar) const;
    RefineSummary RefinePlaneTransform(
        const pcl::PointCloud<pcl::PointXYZINormal>::ConstPtr& current,
        const pcl::PointCloud<pcl::PointXYZINormal>::ConstPtr& history,
        Eigen::Matrix3d& rotation, Eigen::Vector3d& translation) const;
    int FindOdomRevisitCandidate(const BtcDescriptorEntry& current) const;

    BtcLoopDetectorOptions options_;
    STDescManager manager_;
    std::vector<Keyframe::Ptr> pending_keyframes_;
    std::vector<BtcDescriptorEntry> entries_;
    std::vector<std::vector<Keyframe::Ptr>> descriptor_keyframes_;
    Keyframe::Ptr last_keyframe_;
    double journey_ = 0.0;
    int last_confirmation_current_ = -1;
    int last_confirmation_history_ = -1;
    int confirmation_count_ = 0;
    int last_accepted_descriptor_ = -1000000;
};

}  // namespace lightning::backend

#endif  // LIGHTNING_BACKEND_BTC_LOOP_DETECTOR_H
