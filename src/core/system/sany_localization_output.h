#pragma once

#include <cstddef>
#include <string>

#include <geometry_msgs/msg/pose_stamped.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>

#include "common/eigen_types.h"
#include "common/point_def.h"
#include "geosun_msgs/msg/pos_res.hpp"

namespace lightning::sany_output {

geosun_msgs::msg::PosRes MakePosResMessage(const SE3& rear_axle_pose, double vehicle_speed, double stamp,
                                           const std::string& frame_id);

geometry_msgs::msg::PoseStamped MakePoseMessage(const geosun_msgs::msg::PosRes& position);

sensor_msgs::msg::PointCloud2 MakeCloudMessage(const CloudPtr& cloud, double begin_time, double end_time,
                                               const SE3& T_output_lidar, const std::string& frame_id);

class FrameDecimator {
   public:
    explicit FrameDecimator(std::size_t interval) : interval_(interval) {}

    bool Tick() { return interval_ > 0 && ++frame_count_ % interval_ == 0; }

   private:
    std::size_t interval_;
    std::size_t frame_count_ = 0;
};

}  // namespace lightning::sany_output
