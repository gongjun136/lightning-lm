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

/// @brief 单个雷达的话题及到主雷达的固定外参。
struct MultiLidarSensorConfig {
    int id = 0;  ///< 在配置、点来源标记与统计中使用的雷达 ID。
    std::string lidar_topic;  ///< 该雷达的点云输入话题。
    std::string imu_topic;  ///< 与该雷达配置关联的 IMU 输入话题。
    Mat3d R_lidar_to_primary = Mat3d::Identity();  ///< 将该雷达点坐标旋转到主雷达坐标系的固定外参。
    Vec3d t_lidar_to_primary = Vec3d::Zero();  ///< 将该雷达点坐标平移到主雷达坐标系的固定外参（m）。
};

/// @brief 依据处理负载选择雷达和点数预算的配置。
struct AdaptiveLidarLoadConfig {
    bool enabled = false;  ///< 是否启用自适应降载。
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
    std::vector<int> lio_point_budgets;  ///< 体素过滤后的空间点数上限；空配置关闭此上限但传统步长仍可生效，有配置时仅定位和跟踪均健康才启用。
    int tracking_lidar_count = 0;  ///< 0 表示按原有自适应策略决定雷达数量。
    bool rotate_secondary_lidars = false;  ///< 是否轮换参与跟踪的辅助雷达。
    bool predictive = false;  ///< 是否结合预测处理耗时提前调整预算。
};

/// @brief 多雷达组帧、外参与自适应负载配置。
struct MultiLidarConfig {
    // [config-multi-default-disabled]
    bool enabled = false;  ///< 是否启用多雷达组帧；缺少整个配置块时保持关闭。
    // [config-multi-default-disabled]
    int primary_lidar_id = 0;  ///< 组帧与输出坐标系使用的主雷达 ID。
    double frame_period = 0.1;  ///< 组帧时间桶的周期（s）。
    double match_tolerance = 0.002;  ///< 同一时间桶内各雷达扫描开始时间戳的最大允许跨度（s）。
    double reorder_window = 0.5;  ///< 等待乱序或晚到输入的传感器时间进度窗口（s），不是墙上时间。
    int min_lidars = 1;  ///< 排出融合帧所需的最少雷达数；自适应启用时组帧使用 tracking_min_lidars。
    bool online_extrinsic_estimation = false;  ///< 是否配置在线外参估计；启用时当前实现会报错。
    AdaptiveLidarLoadConfig adaptive_load;  ///< 多雷达选择与点数预算配置。
    std::vector<MultiLidarSensorConfig> lidars;  ///< 各雷达的话题、ID 与相对主雷达外参。

    /// 按 ID 查找配置；不存在时返回空指针。
    const MultiLidarSensorConfig* FindLidar(int id) const;
    /// 获取主雷达配置；不存在时返回空指针。
    const MultiLidarSensorConfig* PrimaryLidar() const { return FindLidar(primary_lidar_id); }
};

/// @brief 以主雷达为原点、车体轴对齐的自车过滤盒。
struct SelfPointFilterConfig {
    bool enabled = false;  ///< 是否启用自车点过滤。
    Vec3d min_body = Vec3d::Zero();  ///< 自反射过滤盒的下界（m），车体轴对齐且原点在主雷达。
    Vec3d max_body = Vec3d::Zero();  ///< 自反射过滤盒的上界（m），车体轴对齐且原点在主雷达。
};

/// @brief 融合帧的时间范围、来源和点数统计。
struct MultiLidarFrameStats {
    double begin_time = 0.0;  ///< 融合帧最早的扫描起始时间戳（s）。
    double end_time = 0.0;  ///< 融合帧最晚的扫描结束时间戳（s）。
    std::vector<int> present_lidar_ids;  ///< 融合帧实际包含的雷达 ID。
    std::vector<int> missing_lidar_ids;  ///< 融合帧缺失的已配置雷达 ID。
    std::map<int, std::size_t> points_by_lidar;  ///< 融合帧中各雷达贡献的点数。
    std::size_t merged_points = 0;  ///< 融合点云的总点数。
    bool partial = false;  ///< 该帧是否缺少部分已配置雷达。
};

/// @brief 主雷达坐标系中的融合点云及其统计。
struct FusedLidarFrame {
    CloudPtr cloud;  ///< 已经统一到主雷达坐标系的融合点云。
    MultiLidarFrameStats stats;  ///< 融合帧的时间范围、来源与完整性统计。
};

/**
 * @brief 从 YAML 根节点读取并校验多雷达配置。
 * @param root 包含可选 multi_lidar 配置块的根节点。
 * @param config 输出配置；调用时先重置为默认值，失败时可能已部分写入。
 * @param error 非空时接收错误说明。
 * @return 配置有效或配置块缺失时返回 true。
 * @note 缺少配置块或显式 enabled=false 时保持禁用；配置块存在且缺少 enabled 时默认启用。
 */
bool LoadMultiLidarConfig(const YAML::Node& root, MultiLidarConfig& config, std::string* error = nullptr);

/**
 * @brief 读取以主雷达为原点的车体轴对齐过滤盒。
 * @param root 包含可选 self_point_filter 配置块的根节点。
 * @param config 输出配置；调用时先重置为默认值。
 * @param error 非空时接收错误说明。
 * @return 配置有效或配置块缺失时返回 true。
 * @note 缺少配置块或显式 enabled=false 时保持禁用；配置块存在且缺少 enabled 时默认启用。
 */
bool LoadSelfPointFilterConfig(const YAML::Node& root, SelfPointFilterConfig& config,
                               std::string* error = nullptr);

/**
 * @brief 原地移除过滤盒内及边界上的点。
 * @param cloud 主雷达坐标系点云，过滤后更新点云元数据。
 * @param config 车体轴对齐的过滤盒，单位 m。
 * @param R_primary_to_body 主雷达坐标到车体前左上轴的旋转。
 * @return 移除的点数；过滤关闭或点云为空时为 0。
 */
std::size_t FilterSelfPoints(PointCloudType& cloud, const SelfPointFilterConfig& config,
                             const Mat3d& R_primary_to_body);

/**
 * @brief 每个体素保留离中心最近的原始点，保留来源 ID。
 * @param cloud 输入点云，可为空指针。
 * @param leaf_size 体素边长（m）；非正值时复制输入。
 * @return 新点云；空指针输入产生空点云。
 */
CloudPtr DownsamplePreservingSource(const CloudPtr& cloud, double leaf_size);

/**
 * @brief 按来源、方位、俯仰和距离分桶，确定性地限制点数。
 * @par 主要流程
 * 1. 未超上限时直接返回输入。
 * 2. 将有限坐标点按来源和空间特征分桶。
 * 3. 按配置分配各桶配额。
 * 4. 按输入顺序选取原始点并设置元数据。
 * @param cloud 输入点云；不会被修改。
 * @param max_points 点数上限，0 表示不限制。
 * @param preserve_density true 按各桶点数比例分配配额，false 尽量等额分配。
 * @return 保留原始顺序、时间和来源 ID 的点云；输入为空指针、上限为 0 或点数未超上限时直接返回原指针。
 */
CloudPtr SampleSpatiallyBalanced(const CloudPtr& cloud, std::size_t max_points,
                                bool preserve_density = false);

/// @brief 一帧 LIO 所选的雷达来源和点数预算。
struct AdaptiveLidarSelection {
    std::vector<int> lidar_ids;  ///< 当前选择参与定位计算的雷达 ID。
    int point_stride = 1;  ///< 当前定位点云的采样步长。
    int degradation_step = 0;  ///< 该次选择对应的降载级别。
    std::size_t max_points = 0;  ///< 体素过滤后的点数上限，0 表示不限制。
};

/** @brief 根据负载与健康状态选择雷达和点数预算。
 * @note 仅定位健康标志使用原子变量；其他状态的并发访问需由调用者同步。
 */
class AdaptiveLidarLoadController {
   public:
    /// 重置配置、健康状态和降载统计。
    void Reset(const MultiLidarConfig& config);
    /// 设置地图定位是否健康。
    void SetLocalizationGood(bool good) { localization_good_ = good; }
    // 根据处理耗时、滞后（均为秒）和跟踪健康度更新降载状态，详见定义。
    void Observe(double processing_sec, double latency_sec, bool tracking_healthy);

    // 根据帧来源及滞后（秒）选择预算；可能更新降载和轮换状态，详见定义。
    AdaptiveLidarSelection Select(const MultiLidarFrameStats& stats, double latency_sec = 0.0);
    /// 判断帧是否满足点云发布的来源要求。
    bool CanPublishCloud(const MultiLidarFrameStats& stats) const;
    /// 判断滞后（秒）是否达到或超过硬阈值。
    bool IsHardStale(double latency_sec) const;

    /// 获取当前降载级别。
    int DegradationStep() const { return degradation_step_; }
    /// 获取当前目标雷达数量。
    int TargetLidarCount() const;
    /// 获取当前点采样步长。
    int PointStride() const;

   private:
    /// 计算当前配置允许的最大降载级别。
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

/// 根据选择的来源和步长筛选点云；不实施 selection.max_points 上限。
CloudPtr SelectLidarPoints(const CloudPtr& cloud, const AdaptiveLidarSelection& selection);

/**
 * @brief 按传感器时间桶组装多雷达帧，输出统一到主雷达坐标系。
 * @note 本类没有内部锁，调用者须串行化所有访问，包括统计读取。
 * @see @ref sensor_pipeline
 */
class MultiLidarFrameAssembler {
   public:
    enum class Mode { kClouds, kInputStamps };
    /// 构造空组帧器；使用前须以已校验配置调用 Reset。
    MultiLidarFrameAssembler() = default;
    /// 使用已校验配置构造；本类不主动校验配置。
    explicit MultiLidarFrameAssembler(MultiLidarConfig config) { Reset(std::move(config)); }

    /// 用已校验配置重置缓存和统计；本方法不主动校验配置。
    void Reset(MultiLidarConfig config, Mode mode = Mode::kClouds);
    // 接纳一个雷达的点云；扫描开始时间以秒计，详细契约见定义。
    bool AddCloud(int lidar_id, double timestamp, CloudPtr cloud);
    // Input-stamp mode uses identical bucket/source/time rules without touching points.
    bool AddStamp(int lidar_id, double timestamp);
    /// 取出最早的就绪帧；队列为空时返回 false，frame 不变。
    bool PopReady(FusedLidarFrame& frame);
    /// 强制处理待排出的桶，仍按最低雷达数量丢弃不足的帧。
    void Flush();

    /// 获取待处理时间桶数。
    std::size_t PendingFrameCount() const { return frames_.size(); }
    /// 获取就绪帧数。
    std::size_t ReadyFrameCount() const { return ready_frames_.size(); }
    /// 获取晚到输入丢弃数。
    std::size_t LateDropCount() const { return late_drop_count_; }
    /// 获取重复输入丢弃数。
    std::size_t DuplicateDropCount() const { return duplicate_drop_count_; }
    /// 获取超出匹配容差的输入丢弃数。
    std::size_t ToleranceDropCount() const { return tolerance_drop_count_; }
    /// 获取无效输入丢弃数。
    std::size_t InvalidDropCount() const { return invalid_drop_count_; }
    /// 获取已排出的融合帧数。
    std::size_t EmittedFrameCount() const { return emitted_frame_count_; }
    /// 获取雷达数量不足而丢弃的桶数。
    std::size_t InsufficientLidarDropCount() const { return insufficient_lidar_drop_count_; }

   private:
    struct PartialFrame {
        std::map<int, CloudPtr> clouds;  ///< 当前时间桶内按雷达 ID 保存的点云。
        std::map<int, double> timestamps;  ///< 当前时间桶内各雷达的传感器时间戳（s）。
    };

    /// 将传感器时间映射到时间桶。
    long long BucketFor(double timestamp) const;
    /// 计算时间桶对应的参考时间（s）。
    double BucketTime(long long bucket) const;
    /// 判断桶是否包含全部期望雷达。
    bool IsComplete(const PartialFrame& frame) const;
    /// 按已见传感器时间判断桶是否过期。
    bool IsExpired(long long bucket) const;
    bool AddInput(int lidar_id, double timestamp, CloudPtr cloud);
    FusedLidarFrame DescribeInputFrame(const PartialFrame& frame) const;
    // 按时间顺序排出就绪或强制处理的桶。
    void PromoteReady(bool force);
    // 将桶内点云变换并合并到主雷达系。
    FusedLidarFrame Assemble(const PartialFrame& frame) const;

    MultiLidarConfig config_;  ///< 雷达组帧及自适应降载配置。
    Mode mode_ = Mode::kClouds;
    std::set<int> expected_lidar_ids_;  ///< 组帧器期望接收的雷达 ID 集合。
    std::map<long long, PartialFrame> frames_;  ///< 按时间桶保存尚未排出的部分组帧。
    std::deque<FusedLidarFrame> ready_frames_;  ///< 已完成组帧、等待调用者取出的点云队列。
    bool anchor_initialized_ = false;  ///< 是否已建立组帧时间桶的参考原点。
    double anchor_time_ = 0.0;  ///< 组帧时间桶的参考时间戳（s）。
    double max_seen_time_ = -std::numeric_limits<double>::infinity();  ///< 已观测到的最新雷达时间戳（s）。
    long long last_emitted_bucket_ = std::numeric_limits<long long>::min();  ///< 最近已排出的时间桶编号，用于拒绝晚到输入。
    std::size_t late_drop_count_ = 0;  ///< 晚于对应时间桶排出后才到达的点云数。
    std::size_t duplicate_drop_count_ = 0;  ///< 同一雷达在同一时间桶重复输入的丢弃数。
    std::size_t tolerance_drop_count_ = 0;  ///< 同桶输入的最大时间跨度超过匹配容差的丢弃数。
    std::size_t invalid_drop_count_ = 0;  ///< 输入 ID、时间戳或点云无效的丢弃数。
    std::size_t emitted_frame_count_ = 0;  ///< 已经排出的融合点云帧数。
    std::size_t insufficient_lidar_drop_count_ = 0;  ///< 因雷达数量不足而丢弃的时间桶数。
};

}  // namespace lightning
