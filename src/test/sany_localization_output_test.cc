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
    const SE3 map_livox_pose(quaternion, Vec3d(0.3, -18.2, 0.9));

    const auto position = MakePosResMessage(map_livox_pose, -1.25, 123.5, "map");
    Require(position.header.frame_id == "map" && position.header.stamp.sec == 123 &&
                position.header.stamp.nanosec == 500000000,
            "PosRes uses the map frame and sensor timestamp");
    Require(Near(position.f8enh[0], 0.3) && Near(position.f8enh[1], -18.2) && Near(position.f8enh[2], 0.9),
            "PosRes ENH");
    Require(Near(position.f8pry[0], -0.4) && Near(position.f8pry[1], -1.8) && Near(position.f8pry[2], 353.0),
            "PosRes roll/pitch/yaw degrees");
    Require(Near(position.f8vehiclespeed, -1.25), "PosRes signed vehicle speed");
    Require(Near(position.f8venh[0], 0.0) && Near(position.f8venh[1], 0.0) && Near(position.f8venh[2], 0.0),
            "PosRes only publishes f8vehiclespeed");

    const auto pose_message = MakePoseMessage(position);
    Require(pose_message.header.frame_id == "map" && pose_message.header.stamp == position.header.stamp,
            "PoseStamped and PosRes share frame and timestamp");
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

    const SO3 initial_lidar_rotation = SO3::exp(Vec3d(0.0, 0.0, M_PI_2));
    const SE3 T_livox_lidar = MakeLivoxLidarTransform(initial_lidar_rotation);
    const auto livox_cloud = MakeCloudMessage(cloud, 10.0, 10.1, T_livox_lidar, "livox_frame");
    Require(livox_cloud.header.frame_id == "livox_frame" && livox_cloud.header.stamp.sec == 10 &&
                livox_cloud.header.stamp.nanosec == 100000000,
            "inverse cloud uses livox_frame and scan end timestamp");
    Require(livox_cloud.point_step == 26 && livox_cloud.fields.size() == 7, "reference PointCloud2 layout");
    sensor_msgs::PointCloud2ConstIterator<float> livox_x(livox_cloud, "x"), livox_y(livox_cloud, "y"),
        livox_z(livox_cloud, "z");
    sensor_msgs::PointCloud2ConstIterator<std::uint8_t> livox_tag(livox_cloud, "tag"),
        livox_line(livox_cloud, "line");
    sensor_msgs::PointCloud2ConstIterator<double> livox_timestamp(livox_cloud, "timestamp");
    Require(Near(*livox_x, -2.0) && Near(*livox_y, 1.0) && Near(*livox_z, 3.0),
            "fixed initial rotation maps lidar cloud to livox frame without translation");
    Require(*livox_tag == 0 && *livox_line == 3, "tag and source line fields");
    Require(Near(*livox_timestamp, 10.05), "absolute point timestamp");

    const SE3 T_map_lidar(SO3::exp(Vec3d(0.1, -0.2, 0.3)), Vec3d(5.0, 6.0, 7.0));
    const SE3 T_map_livox = MakeMapLivoxPose(T_map_lidar, initial_lidar_rotation);
    const Vec3d raw_point(1.0, 2.0, 3.0);
    Require((T_map_livox * (T_livox_lidar * raw_point) - T_map_lidar * raw_point).norm() < 1e-9,
            "map-livox pose composes to the original map-lidar registration");
    const auto map_cloud = MakeCloudMessage(cloud, 10.0, 10.1, T_map_lidar, "map");
    Require(map_cloud.header.frame_id == "map" && map_cloud.header.stamp == livox_cloud.header.stamp,
            "same lidar batch uses the same timestamp in livox and map frames");
    sensor_msgs::PointCloud2ConstIterator<float> map_x(map_cloud, "x"), map_y(map_cloud, "y"), map_z(map_cloud, "z");
    const Vec3d expected_map_point = T_map_lidar * raw_point;
    Require(Near(*map_x, expected_map_point.x(), 1e-5) && Near(*map_y, expected_map_point.y(), 1e-5) &&
                Near(*map_z, expected_map_point.z(), 1e-5),
            "map cloud preserves the lidar registration result");

    FrameDecimator decimator(10);
    for (int frame = 1; frame < 10; ++frame) Require(!decimator.Tick(), "no early decimated frame");
    Require(decimator.Tick(), "publish every tenth frame");

    LocalizationPublicationGate publication_gate(5);
    Require(!publication_gate.MapOutputsEnabled(), "map outputs wait for the first valid match");
    publication_gate.ObserveLidarMatch(false);
    Require(!publication_gate.MapOutputsEnabled(), "startup failures do not enable map outputs");
    publication_gate.ObserveLidarMatch(true);
    Require(publication_gate.MapOutputsEnabled(), "first valid match enables map outputs");
    for (int lost = 1; lost < 5; ++lost) {
        publication_gate.ObserveLidarMatch(false);
        Require(publication_gate.MapOutputsEnabled(), "grace frames keep map outputs enabled");
    }
    publication_gate.ObserveLidarMatch(false);
    Require(!publication_gate.MapOutputsEnabled(), "fifth consecutive failure disables map outputs");
    publication_gate.ObserveLidarMatch(true);
    Require(publication_gate.MapOutputsEnabled(), "one valid match immediately resumes map outputs");

    std::cout << "sany_localization_output_test passed" << std::endl;
    return 0;
}
