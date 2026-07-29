#pragma once

#include <cstddef>
#include <mutex>
#include <string>

#include <geometry_msgs/msg/pose_stamped.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>

#include "common/eigen_types.h"
#include "common/point_def.h"
#include "geosun_msgs/msg/pos_res.hpp"

namespace lightning::sany_output {

geosun_msgs::msg::PosRes MakePosResMessage(const SE3& map_rear_axle_pose, double vehicle_speed, double stamp,
                                           const std::string& frame_id);

geometry_msgs::msg::PoseStamped MakePoseMessage(const geosun_msgs::msg::PosRes& position);

sensor_msgs::msg::PointCloud2 MakeCloudMessage(const CloudPtr& cloud, double begin_time, double end_time,
                                               const SE3& T_output_lidar, const std::string& frame_id);

SE3 MakeLivoxLidarTransform(const SO3& initial_lidar_rotation);
SE3 MakeMapLivoxPose(const SE3& map_lidar_pose, const SO3& initial_lidar_rotation);
SO3 ApplyRearAxleYawCompensation(const SO3& initial_lidar_rotation, double compensation_deg);
SE3 MakeRearAxleLidarTransform(const SO3& initial_lidar_rotation,
                               const Vec3d& primary_lidar_position_in_body);
SE3 MakeMapRearAxlePose(const SE3& map_lidar_pose, const SO3& initial_lidar_rotation,
                        const Vec3d& primary_lidar_position_in_body);

class LocalizationPublicationGate {
   public:
    explicit LocalizationPublicationGate(std::size_t lost_frame_threshold = 5);

    void SetLostFrameThreshold(std::size_t lost_frame_threshold);
    void ObserveLidarMatch(bool valid);
    bool MapOutputsEnabled() const;
    std::size_t ConsecutiveLostFrames() const;

   private:
    mutable std::mutex mutex_;
    std::size_t lost_frame_threshold_ = 5;
    std::size_t consecutive_lost_frames_ = 0;
    bool has_valid_match_ = false;
};

class FrameDecimator {
   public:
    explicit FrameDecimator(std::size_t interval) : interval_(interval) {}

    bool Tick() { return interval_ > 0 && ++frame_count_ % interval_ == 0; }

   private:
    std::size_t interval_;
    std::size_t frame_count_ = 0;
};

}  // namespace lightning::sany_output
