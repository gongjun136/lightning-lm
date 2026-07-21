#include "core/system/sany_localization_output.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

#include <sensor_msgs/point_cloud2_iterator.hpp>

namespace {

void Require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << std::endl;
        std::exit(1);
    }
}

bool Near(double left, double right, double tolerance = 1e-6) {
    return std::abs(left - right) <= tolerance;
}

}  // namespace

int main() {
    using namespace lightning;
    using namespace lightning::sany_output;

    const double roll = -0.4 * M_PI / 180.0;
    const double pitch = -1.8 * M_PI / 180.0;
    const double yaw = -7.0 * M_PI / 180.0;
    const Quatd quaternion = Eigen::AngleAxisd(yaw, Vec3d::UnitZ()) *
                             Eigen::AngleAxisd(pitch, Vec3d::UnitY()) *
                             Eigen::AngleAxisd(roll, Vec3d::UnitX());
    const SE3 rear_pose(quaternion, Vec3d(0.3, -18.2, 0.9));

    const auto position = MakePosResMessage(rear_pose, -1.25, 123.5, "map");
    Require(Near(position.f8enh[0], 0.3) && Near(position.f8enh[1], -18.2) && Near(position.f8enh[2], 0.9),
            "PosRes ENH");
    Require(Near(position.f8pry[0], -0.4) && Near(position.f8pry[1], -1.8) && Near(position.f8pry[2], 353.0),
            "PosRes roll/pitch/yaw degrees");
    Require(Near(position.f8vehiclespeed, -1.25), "PosRes signed vehicle speed");
    Require(Near(position.f8venh[0], 0.0) && Near(position.f8venh[1], 0.0) && Near(position.f8venh[2], 0.0),
            "PosRes only publishes f8vehiclespeed");

    const auto pose_message = MakePoseMessage(position);
    const Quatd reconstructed(pose_message.pose.orientation.w, pose_message.pose.orientation.x,
                              pose_message.pose.orientation.y, pose_message.pose.orientation.z);
    Require(std::abs(std::abs(reconstructed.dot(quaternion)) - 1.0) < 1e-9, "PosRes PRY to pose quaternion");

    CloudPtr cloud(new PointCloudType());
    PointType point;
    point.x = 1.0F;
    point.y = 2.0F;
    point.z = 3.0F;
    point.intensity = 4.0F;
    point.time = 50.0;
    point.lidar_id = 3;
    cloud->push_back(point);
    cloud->is_dense = true;

    const SE3 T_rear_lidar(SO3(), Vec3d(2.0, 0.0, 1.0));
    const auto rear_cloud = MakeCloudMessage(cloud, 10.0, 10.1, T_rear_lidar, "rear_axle");
    Require(rear_cloud.point_step == 26 && rear_cloud.fields.size() == 7, "reference PointCloud2 layout");
    sensor_msgs::PointCloud2ConstIterator<float> rear_x(rear_cloud, "x"), rear_y(rear_cloud, "y"),
        rear_z(rear_cloud, "z");
    sensor_msgs::PointCloud2ConstIterator<std::uint8_t> rear_tag(rear_cloud, "tag"), rear_line(rear_cloud, "line");
    sensor_msgs::PointCloud2ConstIterator<double> rear_timestamp(rear_cloud, "timestamp");
    Require(Near(*rear_x, 3.0) && Near(*rear_y, 2.0) && Near(*rear_z, 4.0), "lidar to rear axle cloud");
    Require(*rear_tag == 0 && *rear_line == 3, "tag and source line fields");
    Require(Near(*rear_timestamp, 10.05), "absolute point timestamp");

    const SE3 T_map_rear(SO3::exp(Vec3d(0.0, 0.0, M_PI_2)), Vec3d(5.0, 6.0, 7.0));
    const auto map_cloud = MakeCloudMessage(cloud, 10.0, 10.1, T_map_rear * T_rear_lidar, "map");
    sensor_msgs::PointCloud2ConstIterator<float> map_x(map_cloud, "x"), map_y(map_cloud, "y"), map_z(map_cloud, "z");
    Require(Near(*map_x, 3.0) && Near(*map_y, 9.0) && Near(*map_z, 11.0), "rear pose maps cloud to local map");

    FrameDecimator decimator(10);
    for (int frame = 1; frame < 10; ++frame) Require(!decimator.Tick(), "no early decimated frame");
    Require(decimator.Tick(), "publish every tenth frame");

    std::cout << "sany_localization_output_test passed" << std::endl;
    return 0;
}
