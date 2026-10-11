#pragma once

#include <optional>
#include <string>
#include <vector>

#include <Eigen/StdVector>

#include "common/eigen_types.h"
#include "common/point_def.h"

namespace lightning::loc {

struct RelocalizationCandidate {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    int candidate_id = -1;  ///< 候选在地图描述子数据库中的 ID。
    double score = 0.0;  ///< 后端给出的候选检索分值。
    int query_submap_size = 0;  ///< 查询子地图使用的点云帧数。
    std::size_t point_count = 0;  ///< 查询子地图的有效点数。
    std::size_t descriptor_count = 0;  ///< 查询生成的描述子数量。
    std::size_t rough_match_count = 0;  ///< 粗配准建立的对应关系数。
    double spatial_coverage = 0.0;  ///< 粗配准内点的空间覆盖指标。
    SE3 T_world_imu;  ///< 候选 IMU 在地图世界坐标系中的位姿。
};

using RelocalizationCandidateVector =
    std::vector<RelocalizationCandidate,
                Eigen::aligned_allocator<RelocalizationCandidate>>;

struct RelocalizationResult {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    bool attempted = false;  ///< 是否执行了候选检索。
    bool candidate_found = false;  ///< 是否检索到有效候选。
    bool accepted = false;  ///< 当前后端是否接受候选。
    int candidate_id = -1;  ///< 候选在地图描述子数据库中的 ID。
    double score = 0.0;  ///< 后端给出的候选检索分值。
    double timestamp = 0.0;  ///< 查询末帧的传感器时间戳（s）。
    std::size_t point_count = 0;  ///< 查询子地图的有效点数。
    std::size_t descriptor_count = 0;  ///< 查询生成的描述子数量。
    double search_time_ms = 0.0;  ///< 本次全局候选检索耗时（ms）。
    SE3 T_world_imu;  ///< 候选 IMU 在地图世界坐标系中的位姿。
    RelocalizationCandidateVector candidates;  ///< 供后续精配与地图一致性验证的候选列表。
    std::string reason;  ///< 检索成功或未能返回候选的原因。
};

class GlobalRelocalizer {
   public:
    virtual ~GlobalRelocalizer() = default;
    virtual bool Init(const std::string& config_path, const std::string& map_path,
                      const SE3& T_imu_lidar) = 0;
    virtual void ResetQuery() = 0;
    virtual bool IsReady() const = 0;
    virtual std::size_t DatabaseSize() const = 0;
    virtual const char* Name() const = 0;
    virtual std::optional<RelocalizationResult> AddFrame(
        const CloudPtr& cloud, const SE3& T_odom_imu, double timestamp) = 0;
};

}  // namespace lightning::loc
