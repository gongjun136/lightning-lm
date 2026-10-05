#include "core/navigation/ins_cloud.h"
#include <algorithm>
#include <cmath>

namespace lightning::ins {
bool Deskew(const CloudPtr& input, double begin, double end, const SE3& rear_from_primary,
            const PoseBuffer& poses, CloudPtr& output, SE3& map_from_rear) {
    output.reset();
    if (!input || input->empty() || !poses.Covers(begin, end)) return false;
    Pose reference;
    if (!poses.At(end, reference)) return false;
    map_from_rear = reference.map_from_rear;
    const SE3 rear_end_from_map = map_from_rear.inverse();
    CloudPtr result(new PointCloudType);
    result->reserve(input->size());
    for (const auto& point : input->points) {
        if (!std::isfinite(point.time) || point.time < 0 ||
            !point.getVector3fMap().allFinite()) return false;
        const double stamp = begin + point.time * 1e-3;
        if (stamp < begin || stamp > end + 1e-7) return false;
        Pose at;
        if (!poses.At(std::min(stamp,end), at)) return false;
        const Vec3d p = rear_end_from_map * at.map_from_rear * rear_from_primary *
                        point.getVector3fMap().cast<double>();
        PointType transformed = point;
        transformed.getVector3fMap() = p.cast<float>();
        result->push_back(transformed);
    }
    result->header = input->header;
    result->is_dense = false;
    output = result;
    return true;
}
}
