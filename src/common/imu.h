//
// Created by xiang on 25-3-12.
//

#ifndef LIGHTNING_IMU_H
#define LIGHTNING_IMU_H

#include "common/eigen_types.h"

namespace lightning {

/// 内部定义的里程计观测，如果需要ROS，则从sensor_msg/Imu 的消息转换过来
struct IMU {
    double timestamp = 0;  ///< IMU 观测的传感器时间戳（s）。
    Vec3d angular_velocity = Vec3d::Zero();  ///< IMU 传感器坐标系的角速度（rad/s）。
    Vec3d linear_acceleration = Vec3d::Zero();  ///< IMU 传感器坐标系的加速度量测（m/s²），包含重力。
};

using IMUPtr = std::shared_ptr<IMU>;

}  // namespace lightning

#endif  // LIGHTNING_IMU_H
