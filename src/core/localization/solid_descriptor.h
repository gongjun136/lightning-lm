#pragma once

#include <cstddef>
#include <optional>

#include "common/eigen_types.h"
#include "common/point_def.h"

namespace lightning::loc {

struct SolidDescriptorOptions {
    int range_bins = 40;  ///< 距离分箱数量。
    int angle_bins = 60;  ///< 水平角度分箱数量。
    int elevation_bins = 32;  ///< 垂直角度分箱数量。
    double min_distance = 3.0;  ///< 参与描述子计算的最小距离（m）。
    double max_distance = 80.0;  ///< 参与描述子计算的最大距离（m）。
    double elevation_min_deg = -90.0;  ///< 参与描述子计算的仰角下界（deg）。
    double elevation_max_deg = 90.0;  ///< 参与描述子计算的仰角上界（deg）。
};

struct SolidDescriptor {
    VecXd range;  ///< 用于地点检索的距离描述向量。
    VecXd angle;  ///< 用于候选航向估计的角度描述向量。
    std::size_t accepted_points = 0;  ///< 进入描述子计算的有效点数。
};

// Spatially Organized and Lightweight (SOLiD) global descriptor.  This is a
// bounds-safe implementation of the public SOLiD formulation.  Input clouds
// must already be expressed in one common sensor frame; for SANY this is the
// four-LiDAR fused cloud in the primary/front LiDAR frame.
/** @see @ref relocalization_theory */
class SolidDescriptorEngine {
   public:
    explicit SolidDescriptorEngine(SolidDescriptorOptions options = {});

    const SolidDescriptorOptions& GetOptions() const { return options_; }
    std::optional<SolidDescriptor> Compute(const PointCloudType& cloud) const;

    static double Similarity(const SolidDescriptor& query,
                             const SolidDescriptor& candidate);
    // Returns the yaw rotation that maps query-frame points into the
    // candidate frame.  The result is in [-pi, pi).
    static std::optional<double> EstimateCandidateFromQueryYaw(
        const SolidDescriptor& query, const SolidDescriptor& candidate);

   private:
    SolidDescriptorOptions options_;  ///< 描述子空间分箱与点云范围配置。
};

}  // namespace lightning::loc
