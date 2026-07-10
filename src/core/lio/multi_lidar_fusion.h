#pragma once

#include <yaml-cpp/yaml.h>

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

struct MultiLidarConfig {
    bool enabled = false;
    int primary_lidar_id = 0;
    double frame_period = 0.1;
    double match_tolerance = 0.002;
    double reorder_window = 0.5;
    int min_lidars = 1;
    bool online_extrinsic_estimation = false;
    std::vector<MultiLidarSensorConfig> lidars;

    const MultiLidarSensorConfig* FindLidar(int id) const;
    const MultiLidarSensorConfig* PrimaryLidar() const { return FindLidar(primary_lidar_id); }
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

/// Voxel downsampling that keeps a real input point, so lidar_id is never averaged.
CloudPtr DownsamplePreservingSource(const CloudPtr& cloud, double leaf_size);

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
