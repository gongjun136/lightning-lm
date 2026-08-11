//
// Created by xiang on 25-9-12.
//

#include "core/system/loc_system.h"

#include <chrono>
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

    const auto qos = rclcpp::SensorDataQoS();
    auto subscribe_imu = [this, &qos](const std::string& topic) {
        {
            std::lock_guard<std::mutex> lock(input_stats_mutex_);
            imu_input_stats_.topic = topic;
        }
        imu_sub_ = node_->create_subscription<sensor_msgs::msg::Imu>(
            topic, qos, [this](sensor_msgs::msg::Imu::SharedPtr msg) {
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
                sensor.lidar_topic, qos,
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
                cloud_topic_, qos, [this](sensor_msgs::msg::PointCloud2::SharedPtr cloud) {
                    Timer::Evaluate([&]() { ProcessLidar(cloud); }, "Proc Lidar", true);
                }));
        }
        if (!livox_topic_.empty()) {
            RegisterLidarInput(0, livox_topic_);
            livox_sub_ = node_->create_subscription<livox_ros_driver2::msg::CustomMsg>(
                livox_topic_, qos, [this](livox_ros_driver2::msg::CustomMsg::SharedPtr cloud) {
                    Timer::Evaluate([&]() { ProcessLidar(cloud); }, "Proc Lidar", true);
                });
        }
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
        loc_->SetTFCallback(
            [this](const geometry_msgs::msg::TransformStamped &pose) { tf_broadcaster_->sendTransform(pose); });
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
    ObserveLidarInput(lidar_id, cloud ? ToSec(cloud->header.stamp) : 0.0);
    if (loc_started_) {
        loc_->ProcessLidarMsg(cloud);
    }
}

void LocSystem::ProcessLidar(const sensor_msgs::msg::PointCloud2::SharedPtr& cloud, int lidar_id) {
    ObserveLidarInput(lidar_id, cloud ? ToSec(cloud->header.stamp) : 0.0);
    if (loc_started_) {
        loc_->ProcessLidarMsg(cloud, lidar_id);
    }
}

void LocSystem::ProcessLidar(const livox_ros_driver2::msg::CustomMsg::SharedPtr &cloud) {
    ObserveLidarInput(0, cloud ? ToSec(cloud->header.stamp) : 0.0);
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
    last_localization_stamp_ = result.timestamp_;
    const NavState state = result.ToNavState();
    {
        std::lock_guard<std::mutex> lock(trajectory_mutex_);
        if (localization_states_.empty() || state.timestamp_ > localization_states_.back().timestamp_) {
            localization_states_.push_back(state);
        }
    }
    if (!publication_gate_.MapOutputsEnabled()) return;
    if (!pos_res_pub_ || !pose_pub_) return;
    const SE3 map_rear_axle_pose = sany_output::MakeMapRearAxlePose(
        result.pose_, loc_->GetInitialLidarRotation(), primary_lidar_position_in_body_);
    const auto position =
        sany_output::MakePosResMessage(map_rear_axle_pose, result.vel_b_.x(), result.timestamp_, map_frame_);
    pos_res_pub_->publish(position);
    last_posres_stamp_ = result.timestamp_;
    const auto pose = sany_output::MakePoseMessage(position);
    pose_pub_->publish(pose);
    if (telemetry_) telemetry_->ObservePose(pose);
}

void LocSystem::PublishProcessedCloud(const CloudPtr& cloud, const loc::LocalizationResult& result) {
    publication_gate_.ObserveLidarMatch(result.lidar_loc_valid_);
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
    if (publish_map_frame && publication_gate_.MapOutputsEnabled()) {
        // LidarLoc registers this exact cloud in the map frame and result.pose_ is T_map_lidar.
        map_cloud_pub_->publish(
            sany_output::MakeCloudMessage(cloud, begin_time, end_time, result.pose_, map_frame_));
    }
}

void LocSystem::PublishHealthStatus() {
    if (!node_ || !telemetry_ || !fault_status_pub_ || !loc_status_pub_ || !pipeline_diagnostics_pub_) return;
    const builtin_interfaces::msg::Time stamp = node_->now();
    fault_status_pub_->publish(telemetry_->MakeFaultStatus(stamp));
    loc_status_pub_->publish(telemetry_->MakeLocalizationStatus(stamp));

    lightning::msg::PipelineDiagnostics diagnostics;
    diagnostics.header.stamp = stamp;
    diagnostics.header.frame_id = map_frame_;
    const auto now = std::chrono::steady_clock::now();
    {
        std::lock_guard<std::mutex> lock(input_stats_mutex_);
        for (const auto& [lidar_id, stats] : lidar_input_stats_) {
            diagnostics.lidar_ids.push_back(lidar_id);
            diagnostics.lidar_topics.push_back(stats.topic);
            diagnostics.lidar_message_counts.push_back(stats.message_count);
            diagnostics.lidar_last_sensor_stamps.push_back(stats.last_sensor_stamp);
            diagnostics.lidar_silence_sec.push_back(
                stats.has_arrival ? std::chrono::duration<double>(now - stats.last_arrival).count() : -1.0);
        }
        diagnostics.imu_topic = imu_input_stats_.topic;
        diagnostics.imu_message_count = imu_input_stats_.message_count;
        diagnostics.imu_last_sensor_stamp = imu_input_stats_.last_sensor_stamp;
        diagnostics.imu_silence_sec =
            imu_input_stats_.has_arrival
                ? std::chrono::duration<double>(now - imu_input_stats_.last_arrival).count()
                : -1.0;
    }
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
    diagnostics.relocalization_attempt_count = runtime.relocalization_attempt_count;
    diagnostics.relocalization_accept_count = runtime.relocalization_accept_count;
    diagnostics.last_relocalization_candidate_found = runtime.last_relocalization_candidate_found;
    diagnostics.last_relocalization_accepted = runtime.last_relocalization_accepted;
    diagnostics.last_relocalization_candidate_id = runtime.last_relocalization_candidate_id;
    diagnostics.last_relocalization_score = runtime.last_relocalization_score;
    diagnostics.last_relocalization_search_time_ms = runtime.last_relocalization_search_time_ms;
    diagnostics.last_relocalization_reason = runtime.last_relocalization_reason;
    diagnostics.consecutive_lost_frames = publication_gate_.ConsecutiveLostFrames();
    diagnostics.map_outputs_enabled = publication_gate_.MapOutputsEnabled();
    diagnostics.last_localization_stamp = last_localization_stamp_.load();
    diagnostics.last_posres_stamp = last_posres_stamp_.load();
    pipeline_diagnostics_pub_->publish(diagnostics);
}

void LocSystem::RegisterLidarInput(int lidar_id, const std::string& topic) {
    std::lock_guard<std::mutex> lock(input_stats_mutex_);
    lidar_input_stats_[lidar_id].topic = topic;
}

void LocSystem::ObserveLidarInput(int lidar_id, double sensor_stamp) {
    std::lock_guard<std::mutex> lock(input_stats_mutex_);
    auto& stats = lidar_input_stats_[lidar_id];
    ++stats.message_count;
    stats.last_sensor_stamp = sensor_stamp;
    stats.last_arrival = std::chrono::steady_clock::now();
    stats.has_arrival = true;
}

void LocSystem::ObserveImuInput(double sensor_stamp) {
    std::lock_guard<std::mutex> lock(input_stats_mutex_);
    ++imu_input_stats_.message_count;
    imu_input_stats_.last_sensor_stamp = sensor_stamp;
    imu_input_stats_.last_arrival = std::chrono::steady_clock::now();
    imu_input_stats_.has_arrival = true;
}

void LocSystem::PublishPath() {
    if (!node_ || !telemetry_ || !path_pub_ || telemetry_->PathSize() == 0) return;
    const builtin_interfaces::msg::Time stamp = node_->now();
    path_pub_->publish(telemetry_->MakePath(stamp));
}

}  // namespace lightning
