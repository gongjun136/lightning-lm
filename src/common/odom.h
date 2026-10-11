//
// Created by xiang on 25-3-12.
//

#ifndef LIGHTNING_ODOM_H
#define LIGHTNING_ODOM_H

#include "common/eigen_types.h"

namespace lightning {

/// 内部定义的里程计观测，如果需要ROS，则从odom的消息转换过来
struct Odom {
    double timestamp_ = 0;  ///< 里程计观测的时间戳（s）。
    SE3 pose;  ///< 该条里程计观测的位姿。

    Vec3d linear;  ///< 里程计线速度（m/s）。
    Vec3d angular;  ///< 角速度（若有）
};

using OdomPtr = std::shared_ptr<Odom>;

}  // namespace lightning

#endif  // LIGHTNING_ODOM_H
