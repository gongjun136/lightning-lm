#include "core/system/sany_localization_output.h"
#include "common/timestamp_gate.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

#include <sensor_msgs/point_cloud2_iterator.hpp>

#include "core/lightning_math.hpp"

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

    MonotonicTimestampGate live_output_gate;
    Require(live_output_gate.Observe(100.0).accepted,
            "monotonic output gate accepts the first timestamp");
    Require(!live_output_gate.Observe(100.0).accepted,
            "monotonic output gate rejects a duplicate timestamp");
    Require(!live_output_gate.Observe(99.5).accepted,
            "monotonic output gate rejects a timestamp rollback");
    Require(live_output_gate.Observe(100.1).accepted,
            "monotonic output gate resumes on a newer timestamp");
    Require(live_output_gate.RejectedCount() == 2 &&
                Near(live_output_gate.WorstRollbackSec(), 0.5) &&
                Near(live_output_gate.LastAcceptedTimestamp(), 100.1),
            "monotonic output gate reports duplicate and rollback diagnostics");

    MaximumLagTimestampGate lidar_input_gate;
    lidar_input_gate.SetMaximumLag(0.3);
    Require(lidar_input_gate.Observe(200.0).accepted &&
                lidar_input_gate.Observe(199.75).accepted,
            "lidar timestamp gate permits bounded cross-topic reordering");
    const auto stale_lidar = lidar_input_gate.Observe(198.0);
    Require(!stale_lidar.accepted && Near(stale_lidar.lag_sec, 2.0),
            "lidar timestamp gate rejects a cross-topic stale frame");
    Require(lidar_input_gate.Observe(200.1).accepted,
            "lidar timestamp gate continues with current frames after a stale drop");

    YAML::Node transform_root;
    auto fixed = transform_root["output"]["fixed_map_transform"];
    fixed["enabled"] = true;
    fixed["convention"] = "target_from_localization";
    fixed["source_frame"] = "localization_map";
    fixed["target_frame"] = "map";
    fixed["translation_xyz"] = std::vector<double>{10.0, -2.0, 1.0};
    const lightning::Quatd fixed_q(
        Eigen::AngleAxisd(M_PI / 2.0, lightning::Vec3d::UnitZ()));
    fixed["quaternion_xyzw"] =
        std::vector<double>{fixed_q.x(), fixed_q.y(), fixed_q.z(), fixed_q.w()};
    FixedMapTransform fixed_transform;
    std::string fixed_error;
    Require(LoadFixedMapTransform(transform_root, fixed_transform, fixed_error),
            "load fixed output map transform");
    const lightning::SE3 localization_pose(lightning::Quatd::Identity(),
                                            lightning::Vec3d(2.0, 3.0, 4.0));
    const auto output_pose = TransformPoseForOutput(localization_pose, fixed_transform);
    Require(Near(output_pose.translation().x(), 7.0) &&
                Near(output_pose.translation().y(), 0.0) &&
                Near(output_pose.translation().z(), 5.0),
            "fixed map transform is left-multiplied at output");

    geometry_msgs::msg::TransformStamped localization_tf;
    localization_tf.header.frame_id = "map";
    localization_tf.child_frame_id = "base_link";
    localization_tf.transform.translation.x = 2.0;
    localization_tf.transform.translation.y = 3.0;
    localization_tf.transform.translation.z = 4.0;
    localization_tf.transform.rotation.w = 1.0;
    const auto output_tf = TransformTfForOutput(localization_tf, fixed_transform);
    Require(output_tf.header.frame_id == "map" && output_tf.child_frame_id == "base_link" &&
                Near(output_tf.transform.translation.x, 7.0) &&
                Near(output_tf.transform.translation.y, 0.0) &&
                Near(output_tf.transform.translation.z, 5.0),
            "TF uses the same fixed output transform and preserves frame names");

    const double roll = -0.4 * M_PI / 180.0;
    const double pitch = -1.8 * M_PI / 180.0;
    const double yaw = -7.0 * M_PI / 180.0;
    const Quatd quaternion = Eigen::AngleAxisd(yaw, Vec3d::UnitZ()) *
                             Eigen::AngleAxisd(pitch, Vec3d::UnitY()) *
                             Eigen::AngleAxisd(roll, Vec3d::UnitX());
    const SE3 map_rear_axle_pose(quaternion, Vec3d(0.3, -18.2, 0.9));

    const auto position = MakePosResMessage(map_rear_axle_pose, -1.25, 123.5, "map");
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
    const Vec3d primary_lidar_position_in_body(2.199, -0.25, 2.740);
    const SE3 T_rear_lidar =
        MakeRearAxleLidarTransform(initial_lidar_rotation, primary_lidar_position_in_body);
    const auto rear_cloud = MakeCloudMessage(cloud, 10.0, 10.1, T_rear_lidar, "rear_axle");
    Require(rear_cloud.header.frame_id == "rear_axle" && rear_cloud.header.stamp.sec == 10 &&
                rear_cloud.header.stamp.nanosec == 100000000,
            "inverse cloud uses rear_axle and scan end timestamp");
    Require(rear_cloud.point_step == 26 && rear_cloud.fields.size() == 7, "reference PointCloud2 layout");
    sensor_msgs::PointCloud2ConstIterator<float> rear_x(rear_cloud, "x"), rear_y(rear_cloud, "y"),
        rear_z(rear_cloud, "z");
    sensor_msgs::PointCloud2ConstIterator<std::uint8_t> rear_tag(rear_cloud, "tag"),
        rear_line(rear_cloud, "line");
    sensor_msgs::PointCloud2ConstIterator<double> rear_timestamp(rear_cloud, "timestamp");
    const Vec3d raw_point(1.0, 2.0, 3.0);
    const Vec3d expected_rear_point = T_livox_lidar * raw_point + primary_lidar_position_in_body;
    Require(Near(*rear_x, expected_rear_point.x()) && Near(*rear_y, expected_rear_point.y()) &&
                Near(*rear_z, expected_rear_point.z()),
            "rear cloud rotates first and adds the positive lidar lever arm");
    Require(*rear_tag == 0 && *rear_line == 3, "tag and source line fields");
    Require(Near(*rear_timestamp, 10.05), "absolute point timestamp");

    const SE3 T_map_lidar(SO3::exp(Vec3d(0.1, -0.2, 0.3)), Vec3d(5.0, 6.0, 7.0));
    const SE3 T_map_livox = MakeMapLivoxPose(T_map_lidar, initial_lidar_rotation);
    Require((T_map_livox * (T_livox_lidar * raw_point) - T_map_lidar * raw_point).norm() < 1e-9,
            "map-livox pose composes to the original map-lidar registration");
    const SE3 T_map_rear =
        MakeMapRearAxlePose(T_map_lidar, initial_lidar_rotation, primary_lidar_position_in_body);
    Require((T_map_rear * (T_rear_lidar * raw_point) - T_map_lidar * raw_point).norm() < 1e-9,
            "map-rear pose composes rear cloud to the original map-lidar registration");
    Require((T_map_rear.translation() -
             (T_map_livox.translation() - T_map_livox.so3() * primary_lidar_position_in_body))
                .norm() < 1e-9,
            "rear pose subtracts the lidar lever arm in the map orientation");
    Require((T_map_rear.so3().inverse() * T_map_livox.so3()).log().norm() < 1e-9,
            "parallel livox and rear frames keep the same map orientation");

    const SO3 initial_heading = SO3::exp(Vec3d(0.0, 0.0, M_PI));
    const SE3 T_rear_lidar_at_start =
        MakeRearAxleLidarTransform(initial_heading, primary_lidar_position_in_body);
    const SE3 published_start_pose =
        MakeMapRearAxlePose(T_rear_lidar_at_start, initial_heading, primary_lidar_position_in_body);
    Require(published_start_pose.translation().norm() < 1e-9 &&
                published_start_pose.so3().log().norm() < 1e-9,
            "configured 180 degree lidar heading publishes vehicle yaw zero at the map origin");

    const SE3 expected_vehicle_pose(
        SO3::exp(Vec3d(0.0, 0.0, 40.0 * M_PI / 180.0)), Vec3d(8.0, -3.0, 0.6));
    const SE3 raw_map_lidar_pose = expected_vehicle_pose * T_rear_lidar_at_start;
    const SE3 published_vehicle_pose =
        MakeMapRearAxlePose(raw_map_lidar_pose, initial_heading, primary_lidar_position_in_body);
    Require((published_vehicle_pose.inverse() * expected_vehicle_pose).log().norm() < 1e-9,
            "vehicle yaw is not shifted by 180 degrees in the published map pose");

    const auto map_cloud = MakeCloudMessage(cloud, 10.0, 10.1, T_map_lidar, "map");
    Require(map_cloud.header.frame_id == "map" && map_cloud.header.stamp == rear_cloud.header.stamp,
            "same lidar batch uses the same timestamp in rear and map frames");
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

    LocalizationPublicationGate freshness_gate(5);
    freshness_gate.SetMaxLidarMatchAge(0.6);
    freshness_gate.ObserveLidarMatch(true, 10.0);
    Require(freshness_gate.MapOutputsEnabled(10.6), "lidar match remains fresh at the age boundary");
    Require(!freshness_gate.MapOutputsEnabled(10.61), "stale lidar worker disables map outputs");
    Require(freshness_gate.LidarMatchStale(10.61), "stale lidar worker is diagnosed explicitly");
    Require(Near(freshness_gate.LidarMatchAgeSec(10.7), 0.7), "lidar match age uses sensor time");
    Require(Near(freshness_gate.LastLidarMatchStamp(), 10.0), "last completed lidar match is retained");
    freshness_gate.ObserveLidarMatch(false, 10.7);
    Require(!freshness_gate.MapOutputsEnabled(10.7), "an invalid frame cannot clear a stale-output latch");
    freshness_gate.ObserveLidarMatch(true, 10.8);
    Require(freshness_gate.MapOutputsEnabled(10.8), "a fresh valid match clears the stale-output latch");

    LocalizationTelemetryState telemetry(5, 500, 0.1);
    telemetry.Start();
    const auto telemetry_stamp = math::FromSec(20.0);
    Require(telemetry.MakeLocalizationStatus(telemetry_stamp).status ==
                lightning::msg::LocalizationStatus::STATUS_INITIALIZING,
            "telemetry starts in INITIALIZING");
    Require(telemetry.MakeFaultStatus(telemetry_stamp).level ==
                lightning::msg::FaultStatus::LEVEL_NO_FAULT,
            "initialization is not a fault");
    LocalizationTelemetryState startup_telemetry(5);
    startup_telemetry.Start();
    startup_telemetry.ObserveLocalization(loc::LocalizationStatus::FOLLOWING_DR, 5);
    Require(startup_telemetry.MakeFaultStatus(telemetry_stamp).level ==
                lightning::msg::FaultStatus::LEVEL_P1,
            "startup failures do not report localization lost before the first GOOD state");

    telemetry.ObserveLocalization(loc::LocalizationStatus::GOOD, 0);
    telemetry.ObserveLocalization(loc::LocalizationStatus::FOLLOWING_DR, 1);
    auto fault = telemetry.MakeFaultStatus(telemetry_stamp);
    Require(fault.level == lightning::msg::FaultStatus::LEVEL_P1 &&
                fault.fault_type == static_cast<std::int32_t>(
                    LocalizationFaultType::LOCALIZATION_DEGRADED),
            "FOLLOWING_DR immediately reports localization degradation");
    telemetry.ObserveLocalization(loc::LocalizationStatus::FOLLOWING_DR, 5);
    fault = telemetry.MakeFaultStatus(telemetry_stamp);
    Require(fault.level == lightning::msg::FaultStatus::LEVEL_P0 &&
                fault.fault_type == static_cast<std::int32_t>(
                    LocalizationFaultType::LOCALIZATION_LOST),
            "five consecutive failures latch localization lost");
    telemetry.ObserveLocalization(loc::LocalizationStatus::INITIALIZING, 0);
    Require(telemetry.MakeFaultStatus(telemetry_stamp).level ==
                lightning::msg::FaultStatus::LEVEL_P0,
            "localization lost remains latched during relocalization");
    telemetry.ObserveLocalization(loc::LocalizationStatus::GOOD, 0);
    fault = telemetry.MakeFaultStatus(telemetry_stamp);
    Require(fault.level == lightning::msg::FaultStatus::LEVEL_NO_FAULT &&
                fault.fault_type == 0 && fault.description.empty(),
            "GOOD clears localization faults");
    telemetry.ObserveLocalizationStale(true);
    fault = telemetry.MakeFaultStatus(telemetry_stamp);
    Require(fault.level == lightning::msg::FaultStatus::LEVEL_P0 &&
                telemetry.MakeLocalizationStatus(telemetry_stamp).status ==
                    lightning::msg::LocalizationStatus::STATUS_FAIL,
            "a stalled lidar worker changes GOOD to localization lost");
    telemetry.ObserveLocalization(loc::LocalizationStatus::GOOD, 0);
    Require(telemetry.MakeFaultStatus(telemetry_stamp).level ==
                lightning::msg::FaultStatus::LEVEL_NO_FAULT,
            "a new GOOD match recovers from lidar freshness timeout");

    geometry_msgs::msg::PoseStamped path_pose;
    path_pose.header.frame_id = "map";
    telemetry.ObserveLocalization(loc::LocalizationStatus::FOLLOWING_DR, 1);
    path_pose.header.stamp = math::FromSec(29.9);
    telemetry.ObservePose(path_pose);
    Require(telemetry.PathSize() == 0, "FOLLOWING_DR poses are excluded from the path");
    telemetry.ObserveLocalization(loc::LocalizationStatus::GOOD, 0);
    for (int index = 0; index < 501; ++index) {
        path_pose.header.stamp = math::FromSec(30.0 + 0.1 * index);
        path_pose.pose.position.x = static_cast<double>(index);
        telemetry.ObservePose(path_pose);
    }
    Require(telemetry.PathSize() == 500, "path keeps a 500-pose sliding window");
    const auto path = telemetry.MakePath(math::FromSec(80.0));
    Require(path.header.frame_id == "map" && path.poses.size() == 500 &&
                Near(path.poses.front().pose.position.x, 1.0) &&
                Near(path.poses.back().pose.position.x, 500.0),
            "path drops the oldest sampled pose");
    Require(telemetry.OfflineHealthPublishDue(100.0), "offline health publishes immediately");
    Require(!telemetry.OfflineHealthPublishDue(100.05), "offline health is capped at 10 Hz");
    Require(telemetry.OfflineHealthPublishDue(100.1), "offline health publishes after 0.1 seconds");
    Require(!telemetry.OfflinePathPublishDue(100.0), "offline path waits for its first interval");
    Require(!telemetry.OfflinePathPublishDue(101.9), "offline path waits two seconds");
    Require(telemetry.OfflinePathPublishDue(102.0), "offline path publishes every two seconds");

    std::cout << "sany_localization_output_test passed" << std::endl;
    return 0;
}
