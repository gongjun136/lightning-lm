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
    int id = 0;  ///< 在配置、点来源标记与统计中使用的雷达 ID。
    std::string lidar_topic;  ///< 该雷达的点云输入话题。
    std::string imu_topic;  ///< 与该雷达配置关联的 IMU 输入话题。
    Mat3d R_lidar_to_primary = Mat3d::Identity();  ///< 将该雷达点坐标旋转到主雷达坐标系的固定外参。
    Vec3d t_lidar_to_primary = Vec3d::Zero();  ///< 将该雷达点坐标平移到主雷达坐标系的固定外参（m）。
};

struct AdaptiveLidarLoadConfig {
    bool enabled = false;  ///< 是否启用对应的组帧、降载或自反射过滤功能。
    int tracking_min_lidars = 1;  ///< 健康跟踪时允许选用的最少雷达数。
    int relocalization_min_lidars = 2;  ///< 重定位阶段选用的最少雷达数。
    int cloud_publish_min_lidars = 3;  ///< 允许发布完整点云的最少雷达数。
    bool cloud_publish_require_primary = true;  ///< 点云发布是否必须包含主雷达。
    double target_latency_sec = 0.2;  ///< 期望维持的点云处理滞后目标（s）。
    double hard_latency_sec = 0.3;  ///< 拒绝过期点云的硬滞后阈值（s）。
    double degrade_processing_ratio = 0.8;  ///< 触发降载的处理耗时与 target_latency_sec 的比值阈值。
    double recover_processing_ratio = 0.6;  ///< 恢复时处理耗时和滞后各自相对 target_latency_sec 的上限比例。
    int degrade_consecutive_frames = 3;  ///< 触发降载所需的连续过载帧数。
    int recover_consecutive_frames = 20;  ///< 恢复一个级别所需的连续健康帧数。
    std::vector<int> point_strides{1, 2, 3};  ///< 各降载级别的点采样步长。
    // Empty keeps legacy stride-only behavior. Caps apply after voxel filtering,
    // and only while both map localization and LIO tracking are healthy.
    std::vector<int> lio_point_budgets;  ///< 各降载级别的体素过滤后点数上限，空值保留仅步长降载。
    int tracking_lidar_count = 0;  ///< 0: legacy adaptive source count
    bool rotate_secondary_lidars = false;  ///< 是否轮换参与跟踪的辅助雷达。
    bool predictive = false;  ///< 是否结合预测处理耗时提前调整预算。
};

struct MultiLidarConfig {
    // [config-multi-default-disabled]
    bool enabled = false;  ///< 是否启用多雷达组帧；缺少整个配置块时保持关闭。
    // [config-multi-default-disabled]
    int primary_lidar_id = 0;  ///< 组帧与输出坐标系使用的主雷达 ID。
    double frame_period = 0.1;  ///< 组帧时间桶的周期（s）。
    double match_tolerance = 0.002;  ///< 同一时间桶内各雷达扫描开始时间戳的最大允许跨度（s）。
    double reorder_window = 0.5;  ///< 等待乱序或晚到雷达输入的时间窗口（s）。
    int min_lidars = 1;  ///< 排出融合帧所需的最少雷达数。
    bool online_extrinsic_estimation = false;  ///< 是否配置在线外参估计。
    AdaptiveLidarLoadConfig adaptive_load;  ///< 多雷达选择与点数预算配置。
    std::vector<MultiLidarSensorConfig> lidars;  ///< 各雷达的话题、ID 与相对主雷达外参。

    const MultiLidarSensorConfig* FindLidar(int id) const;
    const MultiLidarSensorConfig* PrimaryLidar() const { return FindLidar(primary_lidar_id); }
};

struct SelfPointFilterConfig {
    bool enabled = false;  ///< 是否启用对应的组帧、降载或自反射过滤功能。
    Vec3d min_body = Vec3d::Zero();  ///< 自反射过滤盒的下界（m），车体轴对齐且原点在主雷达。
    Vec3d max_body = Vec3d::Zero();  ///< 自反射过滤盒的上界（m），车体轴对齐且原点在主雷达。
};

struct MultiLidarFrameStats {
    double begin_time = 0.0;  ///< 融合帧最早的扫描起始时间戳（s）。
    double end_time = 0.0;  ///< 融合帧最晚的扫描结束时间戳（s）。
    std::vector<int> present_lidar_ids;  ///< 融合帧实际包含的雷达 ID。
    std::vector<int> missing_lidar_ids;  ///< 融合帧缺失的已配置雷达 ID。
    std::map<int, std::size_t> points_by_lidar;  ///< 融合帧中各雷达贡献的点数。
    std::size_t merged_points = 0;  ///< 融合点云的总点数。
    bool partial = false;  ///< 该帧是否缺少部分已配置雷达。
};

struct FusedLidarFrame {
    CloudPtr cloud;  ///< 已经统一到主雷达坐标系的融合点云。
    MultiLidarFrameStats stats;  ///< 融合帧的时间范围、来源与完整性统计。
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
    std::vector<int> lidar_ids;  ///< 当前选择参与定位计算的雷达 ID。
    int point_stride = 1;  ///< 当前定位点云的采样步长。
    int degradation_step = 0;  ///< 该次选择对应的降载级别。
    std::size_t max_points = 0;  ///< 体素过滤后的点数上限，0 表示不限制。
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

    MultiLidarConfig config_;  ///< 雷达组帧及自适应降载配置。
    std::atomic_bool localization_good_{false};  ///< 地图定位是否健康，可否使用跟踪阶段的降载策略。
    int degradation_step_ = 0;  ///< 当前雷达与点数预算的降载级别。
    int overload_frames_ = 0;  ///< 连续满足降载条件的帧数。
    int recovery_frames_ = 0;  ///< 连续满足恢复条件的帧数。
    bool tracking_healthy_ = false;  ///< 最近前端跟踪健康度。
    double predicted_processing_sec_ = 0.0;  ///< 预测的下一帧处理耗时（s）。
    int last_secondary_id_ = -1;  ///< 最近选中的辅助雷达 ID，用于轮换选择。
};

CloudPtr SelectLidarPoints(const CloudPtr& cloud, const AdaptiveLidarSelection& selection);

/** @see @ref sensor_pipeline */
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
        std::map<int, CloudPtr> clouds;  ///< 当前时间桶内按雷达 ID 保存的点云。
        std::map<int, double> timestamps;  ///< 当前时间桶内各雷达的传感器时间戳（s）。
    };

    long long BucketFor(double timestamp) const;
    double BucketTime(long long bucket) const;
    bool IsComplete(const PartialFrame& frame) const;
    bool IsExpired(long long bucket) const;
    void PromoteReady(bool force);
    FusedLidarFrame Assemble(const PartialFrame& frame) const;

    MultiLidarConfig config_;  ///< 雷达组帧及自适应降载配置。
    std::set<int> expected_lidar_ids_;  ///< 组帧器期望接收的雷达 ID 集合。
    std::map<long long, PartialFrame> frames_;  ///< 按时间桶保存尚未排出的部分组帧。
    std::deque<FusedLidarFrame> ready_frames_;  ///< 已完成组帧、等待调用者取出的点云队列。
    bool anchor_initialized_ = false;  ///< 是否已建立组帧时间桶的参考原点。
    double anchor_time_ = 0.0;  ///< 组帧时间桶的参考时间戳（s）。
    double max_seen_time_ = -std::numeric_limits<double>::infinity();  ///< 已观测到的最新雷达时间戳（s）。
    long long last_emitted_bucket_ = std::numeric_limits<long long>::min();  ///< 最近已排出的时间桶编号，用于拒绝晚到输入。
    std::size_t late_drop_count_ = 0;  ///< 晚于对应时间桶排出后才到达的点云数。
    std::size_t duplicate_drop_count_ = 0;  ///< 同一雷达在同一时间桶重复输入的丢弃数。
    std::size_t tolerance_drop_count_ = 0;  ///< 时间戳偏离桶时刻超过容差的丢弃数。
    std::size_t invalid_drop_count_ = 0;  ///< 输入 ID、时间戳或点云无效的丢弃数。
    std::size_t emitted_frame_count_ = 0;  ///< 已经排出的融合点云帧数。
    std::size_t insufficient_lidar_drop_count_ = 0;  ///< 因雷达数量不足而丢弃的时间桶数。
};

}  // namespace lightning
