#pragma once

#include <cstddef>
#include <optional>

#include "common/eigen_types.h"
#include "common/point_def.h"

namespace lightning::loc {

struct SolidDescriptorOptions {
    int range_bins = 40;
    int angle_bins = 60;
    int elevation_bins = 32;
    double min_distance = 3.0;
    double max_distance = 80.0;
    double elevation_min_deg = -90.0;
    double elevation_max_deg = 90.0;
};

struct SolidDescriptor {
    VecXd range;
    VecXd angle;
    std::size_t accepted_points = 0;
};

// Spatially Organized and Lightweight (SOLiD) global descriptor.  This is a
// bounds-safe implementation of the public SOLiD formulation.  Input clouds
// must already be expressed in one common sensor frame; for SANY this is the
// four-LiDAR fused cloud in the primary/front LiDAR frame.
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
    SolidDescriptorOptions options_;
};

}  // namespace lightning::loc
