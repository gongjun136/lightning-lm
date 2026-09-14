#pragma once
#include "common/eigen_types.h"
namespace lightning {
// Rotation only: velocity remains at the IMU reference point.
inline Vec3d ImuVelocityToBody(const SO3& imu_to_body, const Vec3d& velocity_in_imu) {
    return imu_to_body * velocity_in_imu;
}
inline Vec3d BodyForwardAxisInImu(const SO3& imu_to_body) {
    return imu_to_body.inverse() * Vec3d::UnitX();
}
// r is rear axle -> IMU, expressed in body coordinates.
inline double RearAxleSpeedOffset(const SO3& imu_to_body,
                                 const Vec3d& unbiased_gyro_imu,
                                 const Vec3d& rear_to_imu_body) {
    return (imu_to_body * unbiased_gyro_imu).cross(rear_to_imu_body).x();
}
}
