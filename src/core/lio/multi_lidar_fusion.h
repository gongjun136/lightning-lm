#pragma once

#include <yaml-cpp/yaml.h>

#include <atomic>
#include <cstdint>
#include <deque>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "common/eigen_types.h"
#include "common/point_def.h"

namespace lightning {

struct MultiLidarSensorConfig {
    int id = 0;
    std::string lidar_topic;
    std::string imu_topic;
    Mat3d R_lidar_to_primary = Mat3d::Identity();
    Vec3d t_lidar_to_primary = Vec3d::Zero();
};

struct AdaptiveLidarLoadConfig {
    bool enabled = false;
    int tracking_min_lidars = 1;
    int relocalization_min_lidars = 2;
    int cloud_publish_min_lidars = 3;
    bool cloud_publish_require_primary = true;
    double target_latency_sec = 0.2;
    double hard_latency_sec = 0.3;
    double degrade_processing_ratio = 0.8;
    double recover_processing_ratio = 0.6;
    int degrade_consecutive_frames = 3;
    int recover_consecutive_frames = 20;
    std::vector<int> point_strides{1, 2, 3};
    // Empty keeps legacy stride-only behavior. Caps apply after voxel filtering,
    // and only while both map localization and LIO tracking are healthy.
    std::vector<int> lio_point_budgets;
    int tracking_lidar_count = 0;  // 0: legacy adaptive source count
    bool rotate_secondary_lidars = false;
    bool predictive = false;
};

struct MultiLidarConfig {
    bool enabled = false;
    int primary_lidar_id = 0;
    double frame_period = 0.1;
    double match_tolerance = 0.002;
    double reorder_window = 0.5;
    int min_lidars = 1;
    bool online_extrinsic_estimation = false;
    AdaptiveLidarLoadConfig adaptive_load;
    std::vector<MultiLidarSensorConfig> lidars;

    const MultiLidarSensorConfig* FindLidar(int id) const;
    const MultiLidarSensorConfig* PrimaryLidar() const { return FindLidar(primary_lidar_id); }
};

struct SelfPointFilterConfig {
    bool enabled = false;
    Vec3d min_body = Vec3d::Zero();
    Vec3d max_body = Vec3d::Zero();
};

struct MultiLidarFrameStats {
    double begin_time = 0.0;
    double end_time = 0.0;
    std::vector<int> present_lidar_ids;
    std::vector<int> missing_lidar_ids;
    std::map<int, std::size_t> points_by_lidar;
    std::size_t merged_points = 0;
    bool partial = false;
};

struct FusedLidarFrame {
    CloudPtr cloud;
    MultiLidarFrameStats stats;
};

bool LoadMultiLidarConfig(const YAML::Node& root, MultiLidarConfig& config, std::string* error = nullptr);

/// Load a body-aligned ego box whose origin is the primary LiDAR.
bool LoadSelfPointFilterConfig(const YAML::Node& root, SelfPointFilterConfig& config,
                               std::string* error = nullptr);

/// Remove points inside the ego box after expressing them in vehicle forward-left-up axes.
std::size_t FilterSelfPoints(PointCloudType& cloud, const SelfPointFilterConfig& config,
                             const Mat3d& R_primary_to_body);

/// Voxel downsampling that keeps a real input point, so lidar_id is never averaged.
CloudPtr DownsamplePreservingSource(const CloudPtr& cloud, double leaf_size);

/// Deterministic hard cap with equal source/azimuth/elevation/range-bin quotas.
/// Keeps real points and their original order, timestamps and source IDs.
/// A zero cap disables sampling. Never mutates the input cloud.
/// preserve_density allocates proportional quotas (NDT score preservation).
CloudPtr SampleSpatiallyBalanced(const CloudPtr& cloud, std::size_t max_points,
                                bool preserve_density = false);

struct AdaptiveLidarSelection {
    std::vector<int> lidar_ids;
    int point_stride = 1;
    int degradation_step = 0;
    std::size_t max_points = 0;
};

class AdaptiveLidarLoadController {
   public:
    void Reset(const MultiLidarConfig& config);
    void SetLocalizationGood(bool good) { localization_good_ = good; }
    void Observe(double processing_sec, double latency_sec, bool tracking_healthy);

    AdaptiveLidarSelection Select(const MultiLidarFrameStats& stats, double latency_sec = 0.0);
    bool CanPublishCloud(const MultiLidarFrameStats& stats) const;
    bool IsHardStale(double latency_sec) const;

    int DegradationStep() const { return degradation_step_; }
    int TargetLidarCount() const;
    int PointStride() const;

   private:
    int MaximumDegradationStep() const;

    MultiLidarConfig config_;
    std::atomic_bool localization_good_{false};
    int degradation_step_ = 0;
    int overload_frames_ = 0;
    int recovery_frames_ = 0;
    bool tracking_healthy_ = false;
    double predicted_processing_sec_ = 0.0;
    int last_secondary_id_ = -1;
};

CloudPtr SelectLidarPoints(const CloudPtr& cloud, const AdaptiveLidarSelection& selection);

class MultiLidarFrameAssembler {
   public:
    MultiLidarFrameAssembler() = default;
    explicit MultiLidarFrameAssembler(MultiLidarConfig config) { Reset(std::move(config)); }

    void Reset(MultiLidarConfig config);
    bool AddCloud(int lidar_id, double timestamp, CloudPtr cloud);
    bool PopReady(FusedLidarFrame& frame);
    void Flush();

    std::size_t PendingFrameCount() const { return frames_.size(); }
    std::size_t ReadyFrameCount() const { return ready_frames_.size(); }
    std::size_t LateDropCount() const { return late_drop_count_; }
    std::size_t DuplicateDropCount() const { return duplicate_drop_count_; }
    std::size_t ToleranceDropCount() const { return tolerance_drop_count_; }
    std::size_t InvalidDropCount() const { return invalid_drop_count_; }
    std::size_t EmittedFrameCount() const { return emitted_frame_count_; }
    std::size_t InsufficientLidarDropCount() const { return insufficient_lidar_drop_count_; }

   private:
    struct PartialFrame {
        std::map<int, CloudPtr> clouds;
        std::map<int, double> timestamps;
    };

    long long BucketFor(double timestamp) const;
    double BucketTime(long long bucket) const;
    bool IsComplete(const PartialFrame& frame) const;
    bool IsExpired(long long bucket) const;
    void PromoteReady(bool force);
    FusedLidarFrame Assemble(const PartialFrame& frame) const;

    MultiLidarConfig config_;
    std::set<int> expected_lidar_ids_;
    std::map<long long, PartialFrame> frames_;
    std::deque<FusedLidarFrame> ready_frames_;
    bool anchor_initialized_ = false;
    double anchor_time_ = 0.0;
    double max_seen_time_ = -std::numeric_limits<double>::infinity();
    long long last_emitted_bucket_ = std::numeric_limits<long long>::min();
    std::size_t late_drop_count_ = 0;
    std::size_t duplicate_drop_count_ = 0;
    std::size_t tolerance_drop_count_ = 0;
    std::size_t invalid_drop_count_ = 0;
    std::size_t emitted_frame_count_ = 0;
    std::size_t insufficient_lidar_drop_count_ = 0;
};

}  // namespace lightning
