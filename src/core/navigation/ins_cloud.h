#pragma once
#include "core/navigation/ins_navigation.h"
#include "common/point_def.h"

namespace lightning::ins {
// Input points are in primary-LiDAR axes at their individual acquisition
// times; output is rear-body FLU at scan end. Point times remain in ms.
bool Deskew(const CloudPtr& input, double begin, double end,
            const SE3& rear_from_primary, const PoseBuffer& poses,
            CloudPtr& output, SE3& map_from_rear);
}
