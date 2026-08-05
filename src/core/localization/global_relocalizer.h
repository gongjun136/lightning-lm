#pragma once

#include <optional>
#include <string>
#include <vector>

#include "common/eigen_types.h"
#include "common/point_def.h"

namespace lightning::loc {

struct RelocalizationCandidate {
    int candidate_id = -1;
    double score = 0.0;
    int query_submap_size = 0;
    std::size_t point_count = 0;
    std::size_t descriptor_count = 0;
    std::size_t rough_match_count = 0;
    double spatial_coverage = 0.0;
    SE3 T_world_imu;
};

struct RelocalizationResult {
    bool attempted = false;
    bool candidate_found = false;
    bool accepted = false;
    int candidate_id = -1;
    double score = 0.0;
    double timestamp = 0.0;
    std::size_t point_count = 0;
    std::size_t descriptor_count = 0;
    double search_time_ms = 0.0;
    SE3 T_world_imu;
    std::vector<RelocalizationCandidate> candidates;
    std::string reason;
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
