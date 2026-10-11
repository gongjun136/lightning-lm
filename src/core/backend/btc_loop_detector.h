#ifndef LIGHTNING_BACKEND_BTC_LOOP_DETECTOR_H
#define LIGHTNING_BACKEND_BTC_LOOP_DETECTOR_H

#include <cstddef>
#include <deque>
#include <optional>
#include <string>
#include <vector>

#include "common/keyframe.h"
#include "core/backend/third_party/voxel_slam_btc/BTC.h"
#include "core/maps/map_frame.h"

namespace lightning::backend {

struct BtcLoopDetectorOptions {
    bool enabled = true;  ///< 是否启用 BTC 回环检测。
    int descriptor_submap_size = 10;  ///< 组成一个描述子子地图的关键帧数。
    // Values <= 0 retain the legacy non-overlapping behavior.
    int descriptor_submap_stride = 0;  ///< 相邻子地图的关键帧滑动步长；非正值保留不重叠组帧。
    int max_points_per_submap = 100000;  ///< 子地图采样点数上限。
    int min_points_per_submap = 100;  ///< 生成描述子所需的最少点数。
    double downsample_leaf_size = 0.20;  ///< 子地图体素降采样分辨率（m）。
    double min_loop_score = 0.45;  ///< 接受全局检索候选所需的最低 BTC 分值。
    double max_drift_ratio = 0.05;  ///< 回环平移修正与行程之比的上限。
    double max_rotation_correction_deg = 45.0;  ///< 允许的最大回环旋转修正（deg）。
    bool refine_with_plane_icp = true;  ///< 是否用点面 ICP 验证和精配回环。
    int plane_icp_iterations = 20;  ///< 点面 ICP 迭代次数上限。
    int plane_icp_min_matches = 20;  ///< 点面 ICP 验证所需的最少匹配数。
    double plane_icp_min_observability = 14.0;  ///< 点面 ICP 的最低可观测性要求。
    bool allow_degenerate_plane_icp = true;  ///< 是否允许通过专门门控的退化点面回退。
    double degenerate_min_loop_score = 0.40;  ///< 退化点面回退所需的最低 BTC 分值。
    int degenerate_min_matches = 100;  ///< 退化点面回退所需的最少匹配数。
    int confirmation_count = 2;  ///< 接受回环所需或当前已达到的连续确认次数。
    int confirmation_max_current_gap = 2;  ///< 连续确认允许的当前描述子编号间隔。
    int confirmation_max_history_gap = 2;  ///< 连续确认允许的历史描述子编号间隔。
    int loop_cooldown_descriptors = 10;  ///< 接受回环后再次检索的描述子冷却间隔。
    // Recover revisits that are close in LIO odometry but whose BTC
    // descriptor is absent or below the global retrieval threshold. This is
    // candidate recall only; point-to-plane verification remains mandatory.
    bool enable_odom_revisit_fallback = false;  ///< 是否以里程计邻近关系补充候选，仍需几何验证。
    double odom_revisit_search_radius = 5.0;  ///< 里程计重访候选的搜索半径（m）。
    double odom_revisit_min_journey = 50.0;  ///< 允许里程计重访候选的最小行程间隔（m）。
    int odom_revisit_degenerate_min_matches = 60;  ///< 里程计重访退化回退所需的最少匹配数。
    double odom_revisit_max_drift_ratio = 0.05;  ///< 里程计重访允许的最大修正与行程比。
    double max_odom_revisit_distance = 0.0;  ///< 里程计重访距离门控配置（m）。
    double min_optimization_translation = 0.20;  ///< 判定回环值得触发优化的最小平移修正（m）。
    double min_optimization_rotation_deg = 2.0;  ///< 判定回环值得触发优化的最小旋转修正（deg）。

    ConfigSetting descriptor;  ///< BTC 特征提取与检索算法配置。
};

struct BtcDescriptorEntry {
    int descriptor_id = -1;  ///< 子地图描述子的数据库 ID。
    unsigned long first_keyframe_id = 0;  ///< 子地图包含的首个关键帧 ID。
    unsigned long last_keyframe_id = 0;  ///< 子地图包含的末个关键帧 ID。
    double timestamp = 0.0;  ///< 描述子端点对应的传感器时间戳（s）。
    double journey = 0.0;  ///< 到该子地图端点的累计行程（m）。
    Keyframe::Ptr endpoint;  ///< 该子地图的末端关键帧。
};

struct BtcLoopResult {
    bool descriptor_generated = false;  ///< 本次是否成功生成描述子。
    bool candidate_found = false;  ///< 本次是否找到回环候选。
    bool accepted = false;  ///< 候选是否通过相应的几何或回环门控。
    bool optimization_warranted = false;  ///< 该回环修正是否达到触发优化的幅度。
    int current_descriptor_id = -1;  ///< 当前查询子地图的描述子 ID。
    int history_descriptor_id = -1;  ///< 历史候选子地图的描述子 ID。
    unsigned long current_keyframe_id = 0;  ///< 当前子地图端点的关键帧 ID。
    unsigned long history_keyframe_id = 0;  ///< 历史子地图端点的关键帧 ID。
    double current_timestamp = 0.0;  ///< 当前子地图端点时间戳（s）。
    double history_timestamp = 0.0;  ///< 历史子地图端点时间戳（s）。
    double score = 0.0;  ///< 候选的 BTC 检索分值。
    double drift_translation = 0.0;  ///< 回环相对里程计的平移修正幅度（m）。
    double drift_rotation_deg = 0.0;  ///< 回环相对里程计的旋转修正幅度（deg）。
    double journey_span = 0.0;  ///< 当前与历史子地图之间的累计行程差（m）。
    double drift_ratio = 0.0;  ///< 平移修正与行程差的比值。
    double odom_revisit_distance = 0.0;  ///< 里程计坐标中当前与历史端点距离（m）。
    double plane_icp_observability = 0.0;  ///< 精配后的点面可观测性指标。
    int plane_icp_matches = 0;  ///< 精配后的点面匹配数。
    bool plane_icp_converged = false;  ///< 点面精配是否收敛。
    bool used_degenerate_plane_fallback = false;  ///< 是否使用退化点面回退路径。
    int confirmation_count = 0;  ///< 接受回环所需或当前已达到的连续确认次数。
    std::size_t point_count = 0;  ///< 当前查询子地图的有效点数。
    std::size_t descriptor_count = 0;  ///< 当前查询生成的描述子数量。
    double generation_time_ms = 0.0;  ///< 子地图与描述子生成耗时（ms）。
    double search_time_ms = 0.0;  ///< 回环检索与验证耗时（ms）。
    std::string candidate_source = "none";  ///< 候选由 BTC 检索或里程计重访产生的来源标识。
    std::string rejection_reason;  ///< 候选被拒绝的原因。
    SE3 T_history_lidar_current_lidar;  ///< 将当前雷达坐标变换到历史雷达坐标的回环相对位姿。
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
                                    const SE3& T_imu_lidar,
                                    const map_frame::Metadata* map_metadata = nullptr) const;

    const BtcLoopDetectorOptions& GetOptions() const { return options_; }
    const std::vector<BtcDescriptorEntry>& Entries() const { return entries_; }
    std::size_t PendingKeyframes() const { return pending_keyframes_.size(); }

   private:
    struct RefineSummary {
        bool accepted = false;  ///< 候选是否通过相应的几何或回环门控。
        bool converged = false;  ///< 点面 ICP 是否收敛。
        double observability = 0.0;  ///< 点面验证的可观测性指标。
        int matches = 0;  ///< 点面验证的匹配数量。
    };

    pcl::PointCloud<pcl::PointXYZI>::Ptr BuildSubmap(
        const std::vector<Keyframe::Ptr>& keyframes, const SE3& T_imu_lidar) const;
    RefineSummary RefinePlaneTransform(
        const pcl::PointCloud<pcl::PointXYZINormal>::ConstPtr& current,
        const pcl::PointCloud<pcl::PointXYZINormal>::ConstPtr& history,
        Eigen::Matrix3d& rotation, Eigen::Vector3d& translation) const;
    int FindOdomRevisitCandidate(const BtcDescriptorEntry& current) const;

    BtcLoopDetectorOptions options_;  ///< 回环检索、几何验证与连续确认参数。
    STDescManager manager_;  ///< 管理 BTC 描述子数据库与候选检索。
    std::deque<Keyframe::Ptr> pending_keyframes_;  ///< 等待组成描述子子地图的关键帧。
    std::vector<BtcDescriptorEntry> entries_;  ///< 描述子 ID 到子地图端点与行程的索引。
    std::vector<std::vector<Keyframe::Ptr>> descriptor_keyframes_;  ///< 每个描述子对应的关键帧组。
    Keyframe::Ptr last_keyframe_;  ///< 最近加入的关键帧，用于累计行程。
    double journey_ = 0.0;  ///< 累计里程计行程（m）。
    int last_confirmation_current_ = -1;  ///< 最近连续确认对应的当前描述子 ID。
    int last_confirmation_history_ = -1;  ///< 最近连续确认对应的历史描述子 ID。
    int confirmation_count_ = 0;  ///< 当前候选连续确认次数。
    int last_accepted_descriptor_ = -1000000;  ///< 最近接受回环的当前描述子 ID，用于冷却门控。
};

}  // namespace lightning::backend

#endif  // LIGHTNING_BACKEND_BTC_LOOP_DETECTOR_H
