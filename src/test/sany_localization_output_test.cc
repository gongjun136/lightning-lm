#include "core/system/sany_localization_output.h"
#include "core/system/speed_smoothing_shadow.h"
#include "common/timestamp_gate.h"
#include "core/localization/localization_result.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

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
    using Shadow = lightning::sany_output::SpeedSmoothingShadow;
    Require(Shadow::Enabled("diagnostic", "1", false), "shadow explicit diagnostic opt-in");
    Require(!Shadow::Enabled("production", "1", false), "production blocks shadow");
    Require(!Shadow::Enabled("diagnostic", "1", true), "reduced IO blocks shadow");
    Require(!Shadow::Enabled(nullptr, "1", false), "missing mode blocks shadow");
    Require(!Shadow::Enabled("diagnostic", nullptr, false), "shadow defaults off");
    Shadow shadow;
    Require(Near(shadow.Observe(1.0, 0.0, false), 0.0), "shadow initial raw");
    Require(Near(shadow.Observe(1.005, 1.0, false), -std::expm1(-.005/.03)), "shadow exact causal gain");
    Require(shadow.Observe(1.01, 1.0, true) == 0.0, "parking immediate zero");
    Require(Near(shadow.Observe(1.015, -1.0, false), std::expm1(-.005/.03)), "signed reverse release");
    Require(shadow.Observe(1.2, 2.0, false) == 2.0 && shadow.ResetLast(), "gap resets shadow");
    Require(shadow.Observe(1.2, 3.0, false) == 3.0 && shadow.ResetLast(), "duplicate resets shadow");
    Require(shadow.Observe(1.1, 4.0, false) == 4.0, "rollback resets shadow");
    shadow.Reset();
    Require(shadow.Observe(1.105, 5.0, false) == 5.0, "publication suppression reset");
    Require(std::isnan(shadow.Observe(1.11, NAN, false)), "invalid speed remains visible");
    Require(shadow.Observe(1.115, 6.0, false) == 6.0, "invalid sample clears history");
    lightning::loc::LocalizationResult speed_result;
    speed_result.timestamp_ = 1.0;
    speed_result.rear_axle_speed_timestamp_ = 1.0;
    speed_result.rear_axle_speed_offset_ = 0.3;
    if (std::abs(speed_result.RearAxleSpeedOffset() - 0.3) > 1e-12) return 1;
    speed_result.is_parking_ = true;
    if (speed_result.RearAxleSpeedOffset() != 0.0) return 1;
    speed_result.is_parking_ = false;
    speed_result.timestamp_ = 1.01;
    if (speed_result.RearAxleSpeedOffset() != 0.0) return 1;
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

    const auto position = MakeVehiclePoseMessage(map_rear_axle_pose, -1.25, 123.5, "map");
    Require(position.header.frame_id == "map" && position.header.stamp.sec == 123 &&
                position.header.stamp.nanosec == 500000000,
            "VehiclePose uses the map frame and sensor timestamp");
    Require(Near(position.x, 0.3) && Near(position.y, -18.2) && Near(position.z, 0.9),
            "VehiclePose ENH");
    Require(Near(position.roll, -0.4) && Near(position.pitch, -1.8) && Near(position.yaw, 353.0),
            "VehiclePose roll/pitch/yaw degrees");
    Require(Near(position.speed, -1.25), "VehiclePose signed vehicle speed");
    const auto pose_message = MakePoseMessage(position);
    Require(pose_message.header.frame_id == "map" && pose_message.header.stamp == position.header.stamp,
            "PoseStamped and VehiclePose share frame and timestamp");
    const Quatd reconstructed(pose_message.pose.orientation.w, pose_message.pose.orientation.x,
                              pose_message.pose.orientation.y, pose_message.pose.orientation.z);
    Require(std::abs(std::abs(reconstructed.dot(quaternion)) - 1.0) < 1e-9, "VehiclePose PRY to pose quaternion");

    std::ostringstream published_tum;
    double last_published_tum_timestamp = 0.0;
    Require(WriteTumPoseLine(published_tum, pose_message, last_published_tum_timestamp),
            "published PoseStamped is written to TUM");
    Require(Near(last_published_tum_timestamp, 123.5),
            "published TUM tracks the last timestamp");
    Require(!WriteTumPoseLine(published_tum, pose_message, last_published_tum_timestamp),
            "published TUM rejects duplicate timestamps");
    std::istringstream published_tum_input(published_tum.str());
    double tum_stamp = 0.0;
    double tum_x = 0.0, tum_y = 0.0, tum_z = 0.0;
    double tum_qx = 0.0, tum_qy = 0.0, tum_qz = 0.0, tum_qw = 0.0;
    published_tum_input >> tum_stamp >> tum_x >> tum_y >> tum_z >> tum_qx >> tum_qy >> tum_qz >> tum_qw;
    Require(Near(tum_stamp, 123.5) && Near(tum_x, position.x) &&
                Near(tum_y, position.y) && Near(tum_z, position.z),
            "published TUM preserves the downstream timestamp and rear-axle position");
    Require(Near(tum_qx, pose_message.pose.orientation.x) &&
                Near(tum_qy, pose_message.pose.orientation.y) &&
                Near(tum_qz, pose_message.pose.orientation.z) &&
                Near(tum_qw, pose_message.pose.orientation.w),
            "published TUM preserves the downstream quaternion");

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
    Require(!publication_gate.PoseOutputsEnabled(), "pose outputs wait for the first valid match");
    Require(!publication_gate.MapOutputsEnabled(), "map outputs wait for the first valid match");
    publication_gate.ObserveLidarMatch(false);
    Require(!publication_gate.MapOutputsEnabled(), "startup failures do not enable map outputs");
    Require(!publication_gate.TrackingState().initialized && !publication_gate.TrackingState().confirmed_lost,
            "startup without a valid match is pending, not confirmed runtime loss");
    publication_gate.ObserveLidarMatch(true);
    Require(publication_gate.PoseOutputsEnabled(), "first valid match enables pose outputs");
    Require(publication_gate.MapOutputsEnabled(), "first valid match enables map outputs");
    for (int lost = 1; lost < 5; ++lost) {
        publication_gate.ObserveLidarMatch(false);
        Require(!publication_gate.PoseOutputsEnabled(), "first rejected match stops pose outputs");
        Require(!publication_gate.MapOutputsEnabled(), "rejected match stops map outputs during loss debounce");
        Require(!publication_gate.TrackingState().confirmed_lost &&
                !publication_gate.TrackingState().latest_match_valid,
                "rejected matches below the confirmation threshold are degraded, not lost");
    }
    publication_gate.ObserveLidarMatch(false);
    Require(!publication_gate.PoseOutputsEnabled(), "fifth failure disables pose outputs");
    Require(!publication_gate.MapOutputsEnabled(), "fifth consecutive failure disables map outputs");
    Require(publication_gate.TrackingState().confirmed_lost,
            "the existing consecutive-rejection threshold confirms tracking loss");
    publication_gate.ObserveLidarMatch(true);
    Require(publication_gate.PoseOutputsEnabled(), "one valid match resumes pose outputs");
    Require(publication_gate.MapOutputsEnabled(), "one valid match immediately resumes map outputs");
    Require(!publication_gate.TrackingState().confirmed_lost, "accepted match clears confirmed tracking loss");

    const auto status_time = LocalizationPublicationGate::Clock::now();
    publication_gate.ObserveLidarMatch(false, 30.0, status_time, loc::LocalizationStatus::INITIALIZING);
    Require(publication_gate.TrackingState(30.0, status_time).confirmed_lost,
            "runtime transition back to relocalization confirms loss without a new eight-second grace");
    publication_gate.ObserveLidarMatch(true, 30.1, status_time, loc::LocalizationStatus::GOOD);
    publication_gate.ObserveLidarMatch(false, 30.2, status_time, loc::LocalizationStatus::FAIL);
    Require(publication_gate.TrackingState(30.2, status_time).confirmed_lost,
            "explicit algorithm failure confirms loss without waiting for five more frames");
    publication_gate.ObserveLidarMatch(true, 30.3, status_time, loc::LocalizationStatus::GOOD);
    publication_gate.ObserveLidarMatch(false, 30.3, status_time, loc::LocalizationStatus::FAIL);
    Require(publication_gate.TrackingState(30.3, status_time).confirmed_lost &&
            publication_gate.ConsecutiveLostFrames() == 0 && !publication_gate.PoseOutputsEnabled(30.3, status_time),
            "same-stamp explicit failure is retained without counting a duplicate rejection");
    publication_gate.ObserveLidarMatch(true, 30.3, status_time, loc::LocalizationStatus::GOOD);
    Require(publication_gate.TrackingState(30.3, status_time).confirmed_lost,
            "a duplicate good frame cannot clear the confirmed-loss latch");
    publication_gate.ObserveLidarMatch(true, 30.4, status_time, loc::LocalizationStatus::GOOD);
    publication_gate.ObserveLidarMatch(false, 30.3, status_time, loc::LocalizationStatus::INITIALIZING);
    Require(!publication_gate.TrackingState(30.4, status_time).confirmed_lost,
            "older cached relocalization status cannot roll back a newer accepted match");
    publication_gate.ObserveLidarMatch(false, 30.4, status_time, loc::LocalizationStatus::INITIALIZING);
    Require(publication_gate.TrackingState(30.4, status_time).confirmed_lost,
            "same-stamp runtime relocalization is explicit loss, independent of frame counting");

    LocalizationPublicationGate freshness_gate(5);
    freshness_gate.SetMaxLidarMatchAge(0.6);
    freshness_gate.ObserveLidarMatch(true, 10.0);
    Require(freshness_gate.MapOutputsEnabled(10.6), "lidar match remains fresh at the age boundary");
    Require(!freshness_gate.MapOutputsEnabled(10.61), "stale lidar worker disables map outputs");
    Require(!freshness_gate.PoseOutputsEnabled(),
            "stale lidar worker stops IMU/CAN DR business pose outputs");
    Require(freshness_gate.LidarMatchStale(10.61), "stale lidar worker is diagnosed explicitly");
    Require(!freshness_gate.TrackingState(10.61).confirmed_lost &&
            freshness_gate.TrackingState(10.61).valid_match_age_sec > .6,
            "map-match age alone is not confirmed tracking loss");
    Require(Near(freshness_gate.LidarMatchAgeSec(10.7), 0.7), "lidar match age uses sensor time");
    Require(Near(freshness_gate.LastLidarMatchStamp(), 10.0), "last completed lidar match is retained");
    freshness_gate.ObserveLidarMatch(false, 10.7);
    Require(!freshness_gate.MapOutputsEnabled(10.7), "an invalid frame cannot clear a stale-output latch");
    freshness_gate.ObserveLidarMatch(true, 10.8);
    Require(freshness_gate.MapOutputsEnabled(10.8), "a fresh valid match clears the stale-output latch");

    LocalizationPublicationGate silent_gate(5);
    silent_gate.SetMaxLidarMatchAge(0.5);
    const auto arrival = LocalizationPublicationGate::Clock::now();
    silent_gate.ObserveLidarMatch(true, 20.0, arrival);
    Require(silent_gate.PoseOutputsEnabled(20.5, arrival), "pose age boundary is inclusive");
    Require(!silent_gate.PoseOutputsEnabled(20.51, arrival),
            "pose callback enforces freshness without waiting for the health timer");
    silent_gate.ObserveLidarMatch(true, 20.0, arrival);
    Require(!silent_gate.PoseOutputsEnabled(20.0, arrival), "duplicate match cannot clear stale latch");
    silent_gate.ObserveLidarMatch(true, 20.6, arrival);
    const auto silent_now = arrival + std::chrono::milliseconds(501);
    Require(silent_gate.LidarMatchStale(20.6, silent_now),
            "steady clock detects loss when every sensor timestamp stops");
    Require(!silent_gate.PoseOutputsEnabled(20.6, silent_now), "all-input silence stops business poses");
    silent_gate.ObserveLidarMatch(false, 21.2, silent_now);
    Require(silent_gate.LidarMatchStale(21.2, silent_now),
            "invalid fresh matches cannot refresh valid-localization age");
    silent_gate.ObserveLidarMatch(true, 21.3, silent_now);
    Require(silent_gate.PoseOutputsEnabled(21.3, silent_now), "fresh valid match recovers publication");
    Require(silent_gate.LidarMatchStale(23.0, silent_now), "newer input exposes a processing backlog");
    silent_gate.ObserveLidarMatch(true, 21.4, silent_now);
    Require(!silent_gate.PoseOutputsEnabled(21.4, silent_now),
            "fresh completion of an old scan cannot reopen outputs behind the known sensor watermark");
    silent_gate.ObserveLidarMatch(true, 23.1, silent_now);
    Require(silent_gate.PoseOutputsEnabled(23.1, silent_now), "caught-up match can reopen outputs");

    const auto request_base = LocalizationPublicationGate::Clock::time_point{} + std::chrono::seconds(1);
    const auto steady_ns = [](LocalizationPublicationGate::Clock::time_point time) {
        return std::chrono::duration_cast<std::chrono::nanoseconds>(time.time_since_epoch()).count();
    };
    LocalizationPublicationGate request_gate(5);
    request_gate.SetMaxLidarMatchAge(.5);
    loc::dtc::Monitor request_monitor(loc::dtc::Config{}, request_base);
    Require(!request_gate.ObserveRelocalizationRequested(), "initial pose before first localization stays startup");
    const auto initial_request = request_monitor.Evaluate(request_gate.TrackingState(100, request_base), request_base);
    Require(initial_request.startup_pending && !initial_request.localization_lost &&
            loc::dtc::SafetyValue(initial_request, false, false, true) == 4,
            "initial external-pose request does not invent runtime loss");
    request_gate.CompleteRelocalizationRequest(request_base + std::chrono::milliseconds(10));
    Require(request_gate.ObserveLidarMatch(true, 100, request_base + std::chrono::milliseconds(20),
                loc::LocalizationStatus::GOOD, 100, steady_ns(request_base + std::chrono::milliseconds(15))),
            "fresh match after initial request can establish localization");
    Require(request_gate.ObserveRelocalizationRequested(), "runtime external-pose request remembers prior initialization");
    const auto requested = request_monitor.Evaluate(
        request_gate.TrackingState(100, request_base + std::chrono::milliseconds(30)),
        request_base + std::chrono::milliseconds(30));
    Require(requested.localization_lost && !requested.startup_pending &&
            loc::dtc::SafetyValue(requested, false, false, true) == -9,
            "runtime request reports confirmed loss immediately without waiting for another scan");
    Require(!request_gate.ObserveLidarMatch(true, 100.1, request_base + std::chrono::milliseconds(40),
                loc::LocalizationStatus::GOOD, 100.1, steady_ns(request_base + std::chrono::milliseconds(31))) &&
            request_gate.TrackingState(100.1, request_base + std::chrono::milliseconds(40)).confirmed_lost,
            "GOOD arriving while external reset is pending cannot clear confirmed loss");
    request_gate.CompleteRelocalizationRequest(request_base + std::chrono::milliseconds(50));
    Require(!request_gate.ObserveLidarMatch(true, 100.2, request_base + std::chrono::milliseconds(60),
                loc::LocalizationStatus::GOOD, 100.2, steady_ns(request_base + std::chrono::milliseconds(20))) &&
            request_gate.TrackingState(100.2, request_base + std::chrono::milliseconds(60)).confirmed_lost,
            "pre-request in-flight GOOD cannot clear loss even when its callback arrives after reset");
    Require(!request_gate.ObserveLidarMatch(true, 100.3, request_base + std::chrono::milliseconds(65),
                loc::LocalizationStatus::GOOD, 100.3, steady_ns(request_base + std::chrono::milliseconds(49))),
            "a match begun before reset completion is not recovery evidence");
    Require(request_gate.ObserveLidarMatch(true, 100.4, request_base + std::chrono::milliseconds(70),
                loc::LocalizationStatus::GOOD, 100.4, steady_ns(request_base + std::chrono::milliseconds(51))) &&
            request_gate.PoseOutputsEnabled(100.4, request_base + std::chrono::milliseconds(70)),
            "genuinely new accepted match after reset can recover the gate");

    LocalizationPublicationGate cached_gate(5);
    cached_gate.SetMaxLidarMatchAge(.5);
    cached_gate.ObserveLidarMatch(true, 200, request_base, loc::LocalizationStatus::GOOD, 200, steady_ns(request_base));
    for (int ms : {100, 200, 300, 400, 501}) {
        Require(!cached_gate.ObserveLidarMatch(true, 200 + ms * .001,
                    request_base + std::chrono::milliseconds(ms), loc::LocalizationStatus::GOOD,
                    200, steady_ns(request_base)),
                "advancing cached GOOD timestamps do not refresh real matching evidence");
    }
    const auto cached = cached_gate.TrackingState(200.501, request_base + std::chrono::milliseconds(501));
    Require(cached.match_stale && !cached.confirmed_lost && !cached.output_gate_open &&
            cached.valid_match_age_sec > .5 && Near(cached_gate.LastValidLidarMatchStamp(), 200) &&
            Near(cached_gate.LastLidarMatchStamp(), 200),
            "parked skipped matching becomes stale without fabricating a tracking failure or a fresh match");
    const auto cached_dtc = request_monitor.Evaluate(cached, request_base + std::chrono::milliseconds(501));
    Require(loc::dtc::SafetyValue(cached_dtc, false, true, false) == 4 &&
            !loc::dtc::HealthyLocalizationForDegradedRadar(cached_dtc),
            "cached parked GOOD cannot authorize partial-radar suppression");

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
    Require(fault.level == lightning::msg::FaultStatus::LEVEL_P1 &&
                fault.fault_type == static_cast<std::int32_t>(
                    LocalizationFaultType::LOCALIZATION_DEGRADED) &&
                telemetry.MakeLocalizationStatus(telemetry_stamp).status ==
                    lightning::msg::LocalizationStatus::STATUS_FOLLOWING_DR,
            "a stalled lidar worker reports degraded P1 without inventing tracking loss");
    telemetry.ObserveLocalizationStale(false);
    Require(telemetry.MakeFaultStatus(telemetry_stamp).level == lightning::msg::FaultStatus::LEVEL_NO_FAULT,
            "freshness recovery clears only the freshness degradation");
    telemetry.ObserveLocalization(loc::LocalizationStatus::FAIL, 0);
    fault = telemetry.MakeFaultStatus(telemetry_stamp);
    Require(fault.level == lightning::msg::FaultStatus::LEVEL_P0,
            "an explicit localization failure still reports P0");
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
