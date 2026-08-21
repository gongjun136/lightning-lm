//
// Created by xiang on 25-9-12.
//

#include "core/system/loc_system.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <utility>
#include <vector>

#include "core/localization/localization.h"
#include "io/yaml_io.h"
#include "wrapper/ros_utils.h"
#include "yaml-cpp/yaml.h"

namespace lightning {
namespace {

double CloudStampSec(const CloudPtr& cloud) {
    if (!cloud) return 0.0;
    return static_cast<double>(cloud->header.stamp) * 1e-9;
}

}  // namespace

LocSystem::LocSystem(LocSystem::Options options) : options_(options) {
    /// handle ctrl-c
    signal(SIGINT, lightning::debug::SigHandle);
}

LocSystem::~LocSystem() {
    Finish();
}

bool LocSystem::Init(const std::string &yaml_path, const std::string &map_path_override) {
    finished_ = false;
    lidar_input_timestamp_gate_.Reset();
    posres_timestamp_gate_.Reset();
    last_localization_stamp_ = 0.0;
    last_posres_stamp_ = 0.0;
    loc::Localization::Options opt;
    opt.online_mode_ = true;
    loc_ = std::make_shared<loc::Localization>(opt);

    YAML_IO yaml(yaml_path);

    const YAML::Node root = YAML::LoadFile(yaml_path);
    std::string map_path = map_path_override;
    if (map_path.empty()) {
        const YAML::Node configured_map = root["system"] ? root["system"]["map_path"] : YAML::Node();
        if (!configured_map || configured_map.IsNull()) {
            LOG(ERROR) << "online localization requires --map or system.map_path";
            return false;
        }
        map_path = configured_map.as<std::string>();
    }
    map_frame_ = root["output"] && root["output"]["map_frame"] ? root["output"]["map_frame"].as<std::string>() : "map";
    std::string fixed_transform_error;
    if (!sany_output::LoadFixedMapTransform(root, fixed_map_transform_, fixed_transform_error)) {
        LOG(ERROR) << fixed_transform_error;
        return false;
    }
    if (fixed_map_transform_.enabled && fixed_map_transform_.target_frame != map_frame_) {
        LOG(ERROR) << "output.fixed_map_transform.target_frame must match output.map_frame ("
                   << map_frame_ << ")";
        return false;
    }
    primary_lidar_position_in_body_ = Vec3d(2.199, 0.0, 2.740);
    if (root["output"] && root["output"]["primary_lidar_position_in_body"]) {
        const auto values = root["output"]["primary_lidar_position_in_body"].as<std::vector<double>>();
        if (values.size() != 3) {
            LOG(ERROR) << "output.primary_lidar_position_in_body must have 3 values";
            return false;
        }
        primary_lidar_position_in_body_ = Vec3d(values[0], values[1], values[2]);
    }
    const int lost_frame_threshold =
        root["relocalization"] && root["relocalization"]["lost_frame_threshold"]
            ? root["relocalization"]["lost_frame_threshold"].as<int>()
            : 5;
    if (lost_frame_threshold <= 0) {
        LOG(ERROR) << "relocalization.lost_frame_threshold must be positive";
        return false;
    }
    publication_gate_.SetLostFrameThreshold(static_cast<std::size_t>(lost_frame_threshold));
    const double max_lidar_match_age_sec =
        root["system"] && root["system"]["localization_output_max_lidar_age_sec"]
            ? root["system"]["localization_output_max_lidar_age_sec"].as<double>()
            : 0.6;
    if (max_lidar_match_age_sec <= 0.0) {
        LOG(ERROR) << "system.localization_output_max_lidar_age_sec must be positive";
        return false;
    }
    publication_gate_.SetMaxLidarMatchAge(max_lidar_match_age_sec);
    online_lidar_input_max_timestamp_lag_sec_ =
        root["system"] && root["system"]["online_lidar_input_max_timestamp_lag_sec"]
            ? root["system"]["online_lidar_input_max_timestamp_lag_sec"].as<double>()
            : 0.0;
    if (!std::isfinite(online_lidar_input_max_timestamp_lag_sec_) ||
        online_lidar_input_max_timestamp_lag_sec_ < 0.0) {
        LOG(ERROR) << "system.online_lidar_input_max_timestamp_lag_sec must be finite and non-negative";
        return false;
    }
    lidar_input_timestamp_gate_.SetMaximumLag(
        online_lidar_input_max_timestamp_lag_sec_);
    wheel_speed_observation_enabled_ =
        root["system"] && root["system"]["enable_wheel_speed_observation"]
            ? root["system"]["enable_wheel_speed_observation"].as<bool>()
            : true;
    if (const char* override_value = std::getenv("SANY_ENABLE_CAN_OBSERVATION")) {
        const std::string value(override_value);
        if (value == "0") {
            wheel_speed_observation_enabled_ = false;
        } else if (value == "1") {
            wheel_speed_observation_enabled_ = true;
        } else {
            LOG(ERROR) << "SANY_ENABLE_CAN_OBSERVATION must be 0 or 1";
            return false;
        }
    }
    wheel_speed_topic_ =
        root["system"] && root["system"]["wheel_speed_topic"]
            ? root["system"]["wheel_speed_topic"].as<std::string>()
            : "/SpeThrCAN4_topic";
    if (const char* override_topic = std::getenv("SANY_WHEEL_SPEED_TOPIC")) {
        wheel_speed_topic_ = override_topic;
    }
    wheel_speed_scale_mps_per_rpm_ =
        root["system"] && root["system"]["wheel_speed_scale_mps_per_rpm"]
            ? root["system"]["wheel_speed_scale_mps_per_rpm"].as<double>()
            : 0.00120639253574024;
    if (wheel_speed_topic_.empty() || wheel_speed_topic_.front() != '/' ||
        !std::isfinite(wheel_speed_scale_mps_per_rpm_) ||
        wheel_speed_scale_mps_per_rpm_ <= 0.0) {
        LOG(ERROR) << "invalid system wheel-speed topic or rpm conversion scale";
        return false;
    }
    telemetry_ = std::make_unique<sany_output::LocalizationTelemetryState>(
        static_cast<std::size_t>(lost_frame_threshold));

    if (!loc_->Init(yaml_path, map_path)) {
        LOG(ERROR) << "failed to initialize online localization";
        return false;
    }

    LOG(INFO) << "online mode, creating ros2 node ... ";

    /// subscribers
    node_ = std::make_shared<rclcpp::Node>("lightning_slam");

    imu_topic_ = yaml.GetValue<std::string>("common", "imu_topic");
    cloud_topic_ = yaml.GetValue<std::string>("common", "lidar_topic");
    livox_topic_ = yaml.GetValue<std::string>("common", "livox_lidar_topic");

    const auto imu_qos = rclcpp::SensorDataQoS();
    const YAML::Node lidar_qos_config = root["system"] ? root["system"]["lidar_qos"] : YAML::Node();
    const size_t lidar_qos_depth =
        lidar_qos_config && lidar_qos_config["depth"]
            ? lidar_qos_config["depth"].as<size_t>()
            : 5;
    if (lidar_qos_depth == 0) {
        LOG(ERROR) << "system.lidar_qos.depth must be positive";
        return false;
    }
    const std::string lidar_qos_reliability =
        lidar_qos_config && lidar_qos_config["reliability"]
            ? lidar_qos_config["reliability"].as<std::string>()
            : "best_effort";
    rclcpp::QoS lidar_qos{rclcpp::KeepLast(lidar_qos_depth)};
    lidar_qos.durability_volatile();
    if (lidar_qos_reliability == "reliable") {
        lidar_qos.reliable();
    } else if (lidar_qos_reliability == "best_effort") {
        lidar_qos.best_effort();
    } else {
        LOG(ERROR) << "system.lidar_qos.reliability must be reliable or best_effort";
        return false;
    }
    LOG(INFO) << "lidar input QoS: reliability=" << lidar_qos_reliability
              << ", depth=" << lidar_qos_depth;

    auto subscribe_imu = [this, &imu_qos](const std::string& topic) {
        {
            std::lock_guard<std::mutex> lock(input_stats_mutex_);
            imu_input_stats_.topic = topic;
        }
        imu_sub_ = node_->create_subscription<sensor_msgs::msg::Imu>(
            topic, imu_qos, [this](sensor_msgs::msg::Imu::SharedPtr msg) {
            IMUPtr imu = std::make_shared<IMU>();
            imu->timestamp = ToSec(msg->header.stamp);
            imu->linear_acceleration =
                Vec3d(msg->linear_acceleration.x, msg->linear_acceleration.y, msg->linear_acceleration.z);
            imu->angular_velocity = Vec3d(msg->angular_velocity.x, msg->angular_velocity.y, msg->angular_velocity.z);

            ProcessIMU(imu);
        });
    };

    if (loc_->IsMultiLidarEnabled()) {
        const auto* primary = loc_->GetMultiLidarConfig().PrimaryLidar();
        if (!primary) {
            LOG(ERROR) << "multi-lidar primary sensor is missing";
            return false;
        }
        subscribe_imu(primary->imu_topic);
        for (const auto& sensor : loc_->GetMultiLidarConfig().lidars) {
            RegisterLidarInput(sensor.id, sensor.lidar_topic);
            cloud_subs_.push_back(node_->create_subscription<sensor_msgs::msg::PointCloud2>(
                sensor.lidar_topic, lidar_qos,
                [this, id = sensor.id](sensor_msgs::msg::PointCloud2::SharedPtr cloud) {
                    Timer::Evaluate([&]() { ProcessLidar(cloud, id); }, "Proc Lidar", true);
                }));
        }
        LOG(INFO) << "online localization subscribed to " << cloud_subs_.size()
                  << " lidar topics and primary IMU " << primary->imu_topic;
    } else {
        if (imu_topic_.empty() || (cloud_topic_.empty() && livox_topic_.empty())) {
            LOG(ERROR) << "single-lidar online topics are incomplete";
            return false;
        }
        subscribe_imu(imu_topic_);
        if (!cloud_topic_.empty()) {
            RegisterLidarInput(0, cloud_topic_);
            cloud_subs_.push_back(node_->create_subscription<sensor_msgs::msg::PointCloud2>(
                cloud_topic_, lidar_qos, [this](sensor_msgs::msg::PointCloud2::SharedPtr cloud) {
                    Timer::Evaluate([&]() { ProcessLidar(cloud); }, "Proc Lidar", true);
                }));
        }
        if (!livox_topic_.empty()) {
            RegisterLidarInput(0, livox_topic_);
            livox_sub_ = node_->create_subscription<livox_ros_driver2::msg::CustomMsg>(
                livox_topic_, lidar_qos, [this](livox_ros_driver2::msg::CustomMsg::SharedPtr cloud) {
                    Timer::Evaluate([&]() { ProcessLidar(cloud); }, "Proc Lidar", true);
                });
        }
    }

    {
        std::lock_guard<std::mutex> lock(input_stats_mutex_);
        wheel_speed_input_stats_.topic = wheel_speed_topic_;
    }
    if (wheel_speed_observation_enabled_) {
        wheel_speed_sub_ = node_->create_subscription<geosun_msgs::msg::SpeThrCAN4>(
            wheel_speed_topic_, rclcpp::QoS(rclcpp::KeepLast(50)).reliable().durability_volatile(),
            [this](geosun_msgs::msg::SpeThrCAN4::SharedPtr message) {
                const double sensor_stamp = ToSec(message->header.stamp);
                ObserveWheelSpeedInput(sensor_stamp, message->x, message->y);
            });
        LOG(INFO) << "wheel-speed observation enabled: topic=" << wheel_speed_topic_
                  << ", scale=" << std::setprecision(16)
                  << wheel_speed_scale_mps_per_rpm_ << " m/s/rpm";
    } else {
        LOG(WARNING) << "wheel-speed observation disabled; using lidar/IMU-only fallback";
    }

    const auto pose_qos =
        rclcpp::QoS(rclcpp::KeepLast(1000)).reliable().durability_volatile();
    const auto cloud_qos = rclcpp::QoS(rclcpp::KeepLast(1)).reliable().durability_volatile();
    pos_res_pub_ = node_->create_publisher<geosun_msgs::msg::PosRes>("/PosRes", pose_qos);
    pose_pub_ = node_->create_publisher<geometry_msgs::msg::PoseStamped>("/slamPoseRaw_topic", pose_qos);
    inv_cloud_pub_ = node_->create_publisher<sensor_msgs::msg::PointCloud2>("/LidarDataInv", cloud_qos);
    map_cloud_pub_ = node_->create_publisher<sensor_msgs::msg::PointCloud2>("/LidarDataInL", cloud_qos);
    const auto health_qos = rclcpp::QoS(rclcpp::KeepLast(10)).reliable().durability_volatile();
    const auto path_qos = rclcpp::QoS(rclcpp::KeepLast(1)).reliable().durability_volatile();
    fault_status_pub_ =
        node_->create_publisher<lightning::msg::FaultStatus>("/localization/fault_status", health_qos);
    loc_status_pub_ =
        node_->create_publisher<lightning::msg::LocalizationStatus>("/localization/loc_status", health_qos);
    pipeline_diagnostics_pub_ = node_->create_publisher<lightning::msg::PipelineDiagnostics>(
        "/localization/pipeline_diagnostics", health_qos);
    path_pub_ = node_->create_publisher<nav_msgs::msg::Path>("/localization/path", path_qos);
    health_timer_ = node_->create_wall_timer(std::chrono::milliseconds(100),
                                             [this]() { PublishHealthStatus(); });
    path_timer_ = node_->create_wall_timer(std::chrono::seconds(2), [this]() { PublishPath(); });

    if (options_.pub_tf_) {
        tf_broadcaster_ = std::make_shared<tf2_ros::TransformBroadcaster>(node_);
        loc_->SetTFCallback([this](const geometry_msgs::msg::TransformStamped &pose) {
            if (!publication_gate_.MapOutputsEnabled(ToSec(pose.header.stamp))) return;
            tf_broadcaster_->sendTransform(
                sany_output::TransformTfForOutput(pose, fixed_map_transform_));
        });
    }
    loc_->SetLocalizationResultCallback([this](const loc::LocalizationResult& result) {
        PublishLocalizationResult(result);
    });
    loc_->SetGlobalLocalizationResultCallback([this](const loc::LocalizationResult& result) {
        CaptureGlobalLocalizationResult(result);
    });
    loc_->SetProcessedCloudCallback(
        [this](const CloudPtr& cloud, const loc::LocalizationResult& result) { PublishProcessedCloud(cloud, result); });

    LOG(INFO) << "online loc node has been created.";
    return true;
}

void LocSystem::SetInitPose(const SE3 &pose) {
    LOG(INFO) << "set init pose: " << pose.translation().transpose() << ", "
              << pose.unit_quaternion().coeffs().transpose();

    loc_->SetExternalPose(pose.unit_quaternion(), pose.translation());
    loc_started_ = true;
    if (telemetry_) telemetry_->Start();
}

void LocSystem::ProcessIMU(const IMUPtr &imu) {
    ObserveImuInput(imu ? imu->timestamp : 0.0);
    if (loc_started_) {
        loc_->ProcessIMUMsg(imu);
    }
}

void LocSystem::ProcessLidar(const sensor_msgs::msg::PointCloud2::SharedPtr &cloud) {
    const int lidar_id = loc_ && loc_->IsMultiLidarEnabled()
                             ? loc_->GetMultiLidarConfig().primary_lidar_id
                             : 0;
    if (!ObserveLidarInput(lidar_id, cloud ? ToSec(cloud->header.stamp) : 0.0)) return;
    if (loc_started_) {
        loc_->ProcessLidarMsg(cloud);
    }
}

void LocSystem::ProcessLidar(const sensor_msgs::msg::PointCloud2::SharedPtr& cloud, int lidar_id) {
    if (!ObserveLidarInput(lidar_id, cloud ? ToSec(cloud->header.stamp) : 0.0)) return;
    if (loc_started_) {
        loc_->ProcessLidarMsg(cloud, lidar_id);
    }
}

void LocSystem::ProcessLidar(const livox_ros_driver2::msg::CustomMsg::SharedPtr &cloud) {
    if (!ObserveLidarInput(0, cloud ? ToSec(cloud->header.stamp) : 0.0)) return;
    if (loc_started_) {
        loc_->ProcessLivoxLidarMsg(cloud);
    }
}

void LocSystem::Spin() {
    if (node_ != nullptr) {
        spin(node_);
    }
}

void LocSystem::Finish() {
    if (!loc_ || finished_) return;
    loc_->Finish();
    finished_ = true;
}

bool LocSystem::SaveTrajectoryTum(const std::string& path) const {
    std::lock_guard<std::mutex> lock(trajectory_mutex_);
    return WriteTrajectoryTum(path, global_localization_states_, "online global localization");
}

bool LocSystem::SaveHighFrequencyTrajectoryTum(const std::string& path) const {
    std::lock_guard<std::mutex> lock(trajectory_mutex_);
    return WriteTrajectoryTum(path, localization_states_, "online high-frequency localization");
}

bool LocSystem::WriteTrajectoryTum(const std::string& path, const std::vector<NavState>& states,
                                   const char* description) const {
    std::ofstream tum(path);
    if (!tum.is_open()) {
        LOG(ERROR) << "failed to open online localization trajectory: " << path;
        return false;
    }
    double last_timestamp = 0.0;
    int count = 0;
    for (const auto& state : states) {
        // Entries are captured only from valid LocalizationResult callbacks.
        // PGO does not currently promote result.status_ from IDLE to GOOD, so
        // NavState::pose_is_ok_ is not a valid additional filter here.  Match
        // the offline localization exporter, which accepts valid fused output.
        if (state.timestamp_ <= 0.0 || state.timestamp_ <= last_timestamp) continue;
        const auto pose = state.GetPose();
        const auto q = pose.unit_quaternion();
        const auto p = pose.translation();
        tum << std::fixed << std::setprecision(9) << state.timestamp_ << " " << std::setprecision(12) << p.x() << " "
            << p.y() << " " << p.z() << " " << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << "\n";
        last_timestamp = state.timestamp_;
        ++count;
    }
    LOG(INFO) << "wrote " << count << " " << description << " poses to " << path;
    return count > 0;
}

void LocSystem::CaptureGlobalLocalizationResult(const loc::LocalizationResult& result) {
    if (!result.valid_ || result.timestamp_ <= 0.0) return;
    const NavState state = result.ToNavState();
    std::lock_guard<std::mutex> lock(trajectory_mutex_);
    if (global_localization_states_.empty() || state.timestamp_ > global_localization_states_.back().timestamp_) {
        global_localization_states_.push_back(state);
    }
}

void LocSystem::PublishLocalizationResult(const loc::LocalizationResult& result) {
    if (!result.valid_ || result.timestamp_ <= 0.0) return;
    std::lock_guard<std::mutex> publish_lock(posres_publish_mutex_);
    const auto timestamp_decision = posres_timestamp_gate_.Observe(result.timestamp_);
    if (!timestamp_decision.accepted) {
        LOG(WARNING) << "drop non-monotonic PosRes candidate: timestamp="
                     << std::setprecision(16) << result.timestamp_
                     << ", last_accepted=" << timestamp_decision.reference_timestamp
                     << ", rollback_sec=" << timestamp_decision.lag_sec;
        return;
    }
    last_localization_stamp_ = result.timestamp_;
    if (!publication_gate_.MapOutputsEnabled(result.timestamp_)) return;
    const NavState state = result.ToNavState();
    {
        std::lock_guard<std::mutex> lock(trajectory_mutex_);
        if (localization_states_.empty() || state.timestamp_ > localization_states_.back().timestamp_) {
            localization_states_.push_back(state);
        }
    }
    if (!pos_res_pub_ && !pose_pub_) return;
    const SE3 localization_rear_axle_pose = sany_output::MakeMapRearAxlePose(
        result.pose_, loc_->GetInitialLidarRotation(), primary_lidar_position_in_body_);
    const SE3 map_rear_axle_pose =
        sany_output::TransformPoseForOutput(localization_rear_axle_pose, fixed_map_transform_);
    double vehicle_speed = result.vel_b_.x();
    {
        std::lock_guard<std::mutex> lock(input_stats_mutex_);
        constexpr double kMaxWheelSpeedTimestampDeltaSec = 0.25;
        if (wheel_speed_input_stats_.has_arrival &&
            wheel_speed_input_stats_.last_sensor_stamp > 0.0 &&
            std::abs(result.timestamp_ - wheel_speed_input_stats_.last_sensor_stamp) <=
                kMaxWheelSpeedTimestampDeltaSec) {
            vehicle_speed = last_wheel_speed_mps_;
        }
    }
    const auto position = sany_output::MakePosResMessage(
        map_rear_axle_pose, vehicle_speed, result.timestamp_, map_frame_);
    if (pos_res_pub_) {
        pos_res_pub_->publish(position);
        last_posres_stamp_ = result.timestamp_;
    }
    if (pose_pub_) {
        const auto pose = sany_output::MakePoseMessage(position);
        pose_pub_->publish(pose);
        if (telemetry_) telemetry_->ObservePose(pose);
    }
}

void LocSystem::PublishProcessedCloud(const CloudPtr& cloud, const loc::LocalizationResult& result) {
    publication_gate_.ObserveLidarMatch(result.lidar_loc_valid_, result.timestamp_);
    if (telemetry_) {
        telemetry_->ObserveLocalization(result.status_, publication_gate_.ConsecutiveLostFrames());
    }
    if (!inv_cloud_pub_ || !map_cloud_pub_ || !cloud || cloud->empty()) return;
    const double begin_time = CloudStampSec(cloud);
    const double end_time = result.timestamp_ > 0.0 ? result.timestamp_ : begin_time;
    if (begin_time <= 0.0 || end_time <= 0.0) return;
    const SO3 initial_lidar_rotation = loc_->GetInitialLidarRotation();
    inv_cloud_pub_->publish(sany_output::MakeCloudMessage(
        cloud, begin_time, end_time,
        sany_output::MakeRearAxleLidarTransform(initial_lidar_rotation, primary_lidar_position_in_body_),
        rear_axle_frame_));
    const bool publish_map_frame = map_cloud_decimator_.Tick();
    if (publish_map_frame && publication_gate_.MapOutputsEnabled(end_time)) {
        // LidarLoc registers this exact cloud in the map frame and result.pose_ is T_map_lidar.
        const SE3 output_lidar_pose =
            sany_output::TransformPoseForOutput(result.pose_, fixed_map_transform_);
        map_cloud_pub_->publish(
            sany_output::MakeCloudMessage(cloud, begin_time, end_time, output_lidar_pose, map_frame_));
    }
}

void LocSystem::PublishHealthStatus() {
    if (!node_ || !telemetry_ || !fault_status_pub_ || !loc_status_pub_ || !pipeline_diagnostics_pub_) return;
    const builtin_interfaces::msg::Time stamp = node_->now();

    lightning::msg::PipelineDiagnostics diagnostics;
    diagnostics.header.stamp = stamp;
    diagnostics.header.frame_id = map_frame_;
    const auto now = std::chrono::steady_clock::now();
    double latest_input_sensor_stamp = 0.0;
    {
        std::lock_guard<std::mutex> lock(input_stats_mutex_);
        for (const auto& [lidar_id, stats] : lidar_input_stats_) {
            diagnostics.lidar_ids.push_back(lidar_id);
            diagnostics.lidar_topics.push_back(stats.topic);
            diagnostics.lidar_message_counts.push_back(stats.message_count);
            diagnostics.lidar_stale_drop_counts.push_back(stats.stale_drop_count);
            diagnostics.lidar_last_sensor_stamps.push_back(stats.last_sensor_stamp);
            diagnostics.lidar_silence_sec.push_back(
                stats.has_arrival ? std::chrono::duration<double>(now - stats.last_arrival).count() : -1.0);
            latest_input_sensor_stamp = std::max(latest_input_sensor_stamp, stats.last_sensor_stamp);
        }
        diagnostics.imu_topic = imu_input_stats_.topic;
        diagnostics.imu_message_count = imu_input_stats_.message_count;
        diagnostics.imu_last_sensor_stamp = imu_input_stats_.last_sensor_stamp;
        diagnostics.imu_silence_sec =
            imu_input_stats_.has_arrival
                ? std::chrono::duration<double>(now - imu_input_stats_.last_arrival).count()
                : -1.0;
        latest_input_sensor_stamp =
            std::max(latest_input_sensor_stamp, imu_input_stats_.last_sensor_stamp);
        diagnostics.wheel_speed_topic = wheel_speed_input_stats_.topic;
        diagnostics.wheel_speed_message_count = wheel_speed_input_stats_.message_count;
        diagnostics.wheel_speed_last_sensor_stamp = wheel_speed_input_stats_.last_sensor_stamp;
        diagnostics.wheel_speed_silence_sec =
            wheel_speed_input_stats_.has_arrival
                ? std::chrono::duration<double>(now - wheel_speed_input_stats_.last_arrival).count()
                : -1.0;
        diagnostics.wheel_speed_mps = last_wheel_speed_mps_;
        diagnostics.motor_speed_rpm = last_motor_rpm_;
        diagnostics.motor_torque_nm = last_motor_torque_;
        constexpr double kMaxWheelSpeedTimestampDeltaSec = 0.25;
        if (imu_input_stats_.last_sensor_stamp > 0.0 &&
            wheel_speed_input_stats_.last_sensor_stamp > 0.0) {
            diagnostics.wheel_speed_imu_stamp_delta_sec =
                wheel_speed_input_stats_.last_sensor_stamp -
                imu_input_stats_.last_sensor_stamp;
            diagnostics.wheel_speed_timestamp_aligned =
                std::abs(diagnostics.wheel_speed_imu_stamp_delta_sec) <=
                kMaxWheelSpeedTimestampDeltaSec;
        }
    }
    const bool lidar_match_stale = publication_gate_.LidarMatchStale(latest_input_sensor_stamp);
    telemetry_->ObserveLocalizationStale(lidar_match_stale);
    fault_status_pub_->publish(telemetry_->MakeFaultStatus(stamp));
    loc_status_pub_->publish(telemetry_->MakeLocalizationStatus(stamp));
    const auto runtime = loc_->GetRuntimeStats();
    diagnostics.sensor_queue_pending = runtime.sensor_queue_pending;
    diagnostics.sensor_queue_dropped = runtime.sensor_queue_dropped;
    diagnostics.sensor_queue_processed = runtime.sensor_queue_processed;
    diagnostics.localization_queue_pending = runtime.localization_queue_pending;
    diagnostics.localization_queue_dropped = runtime.localization_queue_dropped;
    diagnostics.localization_queue_processed = runtime.localization_queue_processed;
    diagnostics.latest_enqueued_sensor_stamp = runtime.latest_enqueued_sensor_stamp;
    diagnostics.latest_processed_sensor_stamp = runtime.latest_processed_sensor_stamp;
    diagnostics.current_sensor_lag_sec = runtime.current_sensor_lag_sec;
    diagnostics.max_sensor_lag_sec = runtime.max_sensor_lag_sec;
    diagnostics.severe_timestamp_rollback_count = runtime.severe_timestamp_rollback_count;
    diagnostics.worst_timestamp_rollback_sec = runtime.worst_timestamp_rollback_sec;
    diagnostics.live_output_non_monotonic_drop_count =
        runtime.live_output_non_monotonic_drop_count;
    diagnostics.worst_live_output_timestamp_rollback_sec =
        runtime.worst_live_output_timestamp_rollback_sec;
    diagnostics.relocalization_attempt_count = runtime.relocalization_attempt_count;
    diagnostics.relocalization_accept_count = runtime.relocalization_accept_count;
    diagnostics.last_relocalization_candidate_found = runtime.last_relocalization_candidate_found;
    diagnostics.last_relocalization_accepted = runtime.last_relocalization_accepted;
    diagnostics.last_relocalization_candidate_id = runtime.last_relocalization_candidate_id;
    diagnostics.last_relocalization_score = runtime.last_relocalization_score;
    diagnostics.last_relocalization_search_time_ms = runtime.last_relocalization_search_time_ms;
    diagnostics.last_relocalization_reason = runtime.last_relocalization_reason;
    diagnostics.consecutive_lost_frames = publication_gate_.ConsecutiveLostFrames();
    const bool map_outputs_enabled =
        publication_gate_.MapOutputsEnabled(latest_input_sensor_stamp);
    const bool map_outputs_were_enabled =
        map_outputs_enabled_last_.exchange(map_outputs_enabled);
    if (map_outputs_enabled) {
        const bool had_enabled_output = map_outputs_ever_enabled_.exchange(true);
        if (had_enabled_output && !map_outputs_were_enabled) {
            LOG(WARNING) << "localization map outputs recovered after a fresh valid lidar match";
        }
    } else if (map_outputs_were_enabled) {
        LOG(ERROR) << "localization map outputs disabled: lidar_match_stale="
                   << lidar_match_stale << ", lidar_match_age_sec="
                   << publication_gate_.LidarMatchAgeSec(latest_input_sensor_stamp)
                   << ", consecutive_lost_frames="
                   << publication_gate_.ConsecutiveLostFrames();
    }
    diagnostics.map_outputs_enabled = map_outputs_enabled;
    diagnostics.lidar_match_stale = lidar_match_stale;
    diagnostics.lidar_match_age_sec =
        publication_gate_.LidarMatchAgeSec(latest_input_sensor_stamp);
    diagnostics.last_lidar_match_stamp = publication_gate_.LastLidarMatchStamp();
    diagnostics.last_valid_lidar_match_stamp =
        publication_gate_.LastValidLidarMatchStamp();
    diagnostics.last_localization_stamp = last_localization_stamp_.load();
    diagnostics.last_posres_stamp = last_posres_stamp_.load();
    diagnostics.posres_non_monotonic_drop_count =
        posres_timestamp_gate_.RejectedCount();
    diagnostics.worst_posres_timestamp_rollback_sec =
        posres_timestamp_gate_.WorstRollbackSec();
    pipeline_diagnostics_pub_->publish(diagnostics);
}

void LocSystem::RegisterLidarInput(int lidar_id, const std::string& topic) {
    std::lock_guard<std::mutex> lock(input_stats_mutex_);
    lidar_input_stats_[lidar_id].topic = topic;
}

bool LocSystem::ObserveLidarInput(int lidar_id, double sensor_stamp) {
    const auto timestamp_decision = lidar_input_timestamp_gate_.Observe(sensor_stamp);
    {
        std::lock_guard<std::mutex> lock(input_stats_mutex_);
        auto& stats = lidar_input_stats_[lidar_id];
        ++stats.message_count;
        stats.last_arrival = std::chrono::steady_clock::now();
        stats.has_arrival = true;
        if (timestamp_decision.accepted) {
            stats.last_sensor_stamp = std::max(stats.last_sensor_stamp, sensor_stamp);
        } else {
            ++stats.stale_drop_count;
        }
    }
    if (!timestamp_decision.accepted) {
        LOG_EVERY_N(WARNING, 10)
            << "drop stale/invalid lidar input before preprocessing: lidar_id="
            << lidar_id << ", timestamp=" << std::setprecision(16) << sensor_stamp
            << ", newest_timestamp=" << timestamp_decision.reference_timestamp
            << ", lag_sec=" << timestamp_decision.lag_sec
            << ", configured_max_lag_sec="
            << online_lidar_input_max_timestamp_lag_sec_;
    }
    return timestamp_decision.accepted;
}

void LocSystem::ObserveImuInput(double sensor_stamp) {
    std::lock_guard<std::mutex> lock(input_stats_mutex_);
    ++imu_input_stats_.message_count;
    imu_input_stats_.last_sensor_stamp = sensor_stamp;
    imu_input_stats_.last_arrival = std::chrono::steady_clock::now();
    imu_input_stats_.has_arrival = true;
}

void LocSystem::ObserveWheelSpeedInput(double sensor_stamp, double motor_rpm,
                                       double motor_torque) {
    const double speed_mps = motor_rpm * wheel_speed_scale_mps_per_rpm_;
    {
        std::lock_guard<std::mutex> lock(input_stats_mutex_);
        ++wheel_speed_input_stats_.message_count;
        wheel_speed_input_stats_.last_sensor_stamp = sensor_stamp;
        wheel_speed_input_stats_.last_arrival = std::chrono::steady_clock::now();
        wheel_speed_input_stats_.has_arrival = true;
        last_motor_rpm_ = motor_rpm;
        last_motor_torque_ = motor_torque;
        last_wheel_speed_mps_ = speed_mps;
    }
    if (loc_started_ && loc_) loc_->ProcessWheelSpeed(sensor_stamp, speed_mps);
}

void LocSystem::PublishPath() {
    if (!node_ || !telemetry_ || !path_pub_ || telemetry_->PathSize() == 0) return;
    if (!publication_gate_.MapOutputsEnabled(last_localization_stamp_.load())) return;
    const builtin_interfaces::msg::Time stamp = node_->now();
    path_pub_->publish(telemetry_->MakePath(stamp));
}

}  // namespace lightning
