#pragma once

#include "common/nav_state.h"

namespace lightning {

/// Converts the estimator IMU pose to a rear-axle pose and rebases the first valid body pose to identity.
class RearAxlePoseTransformer {
   public:
    RearAxlePoseTransformer() = default;
    RearAxlePoseTransformer(const Mat3d& R_lidar_to_imu, const Vec3d& t_lidar_to_imu,
                            const Vec3d& primary_lidar_position_in_body)
        : T_imu_lidar_(SO3(R_lidar_to_imu), t_lidar_to_imu),
          primary_lidar_position_in_body_(primary_lidar_position_in_body) {}

    void Reset() { initialized_ = false; }
    bool IsInitialized() const { return initialized_; }

    SE3 PrimaryLidarPose(const NavState& state) const { return state.GetPose() * T_imu_lidar_; }

    SE3 Transform(const NavState& state) {
        const SE3 T_estimator_lidar = PrimaryLidarPose(state);
        if (!initialized_) {
            // Gravity initialization fixes roll/pitch and the user defines initial lidar yaw as vehicle yaw zero.
            const SE3 T_body_lidar(T_estimator_lidar.so3(), primary_lidar_position_in_body_);
            T_lidar_body_ = T_body_lidar.inverse();
            T_estimator_body_initial_ = T_estimator_lidar * T_lidar_body_;
            initialized_ = true;
        }
        return T_estimator_body_initial_.inverse() * T_estimator_lidar * T_lidar_body_;
    }

   private:
    SE3 T_imu_lidar_;
    Vec3d primary_lidar_position_in_body_ = Vec3d::Zero();
    SE3 T_lidar_body_;
    SE3 T_estimator_body_initial_;
    bool initialized_ = false;
};

}  // namespace lightning
