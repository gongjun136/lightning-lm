//
// Created by xiang on 25-7-7.
//

#ifndef LIGHTNING_TIMED_POSE_H
#define LIGHTNING_TIMED_POSE_H

#include "common/eigen_types.h"

namespace lightning {

/// 带时间的pose
struct TimedPose {
    TimedPose() {}
    TimedPose(double t, const SE3& pose) : timestamp_(t), pose_(pose) {}
    double timestamp_ = 0;  ///< 位姿对应的传感器时间戳（s）。
    SE3 pose_;  ///< 该时刻的刚体位姿。

    bool operator<(const TimedPose& timed_pose) const { return ((timestamp_ - timed_pose.timestamp_) < 0); };
};

}  // namespace lightning

#endif  // LIGHTNING_TIMED_POSE_H
