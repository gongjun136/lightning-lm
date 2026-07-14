//
// Created by xiang on 25-9-12.
//

#include "core/system/loc_system.h"

#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include "core/localization/localization.h"
#include "io/yaml_io.h"
#include "wrapper/ros_utils.h"
#include "sensor_msgs/point_cloud2_iterator.hpp"
#include "yaml-cpp/yaml.h"

namespace lightning {
namespace {

builtin_interfaces::msg::Time ToRosStamp(double seconds) {
    const std::int64_t nanoseconds = static_cast<std::int64_t>(std::llround(seconds * 1e9));
    builtin_interfaces::msg::Time stamp;
    stamp.sec = static_cast<std::int32_t>(nanoseconds / 1000000000LL);
    stamp.nanosec = static_cast<std::uint32_t>(nanoseconds % 1000000000LL);
    return stamp;
}

std::int64_t SteadyNowNs() {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch())
        .count();
}

double CloudStampSec(const CloudPtr& cloud) {
    if (!cloud) return 0.0;
    return static_cast<double>(cloud->header.stamp) * 1e-9;
}

bool IsUsableResult(const loc::LocalizationResult& result) {
    return result.valid_ || result.lidar_loc_valid_;
}

geometry_msgs::msg::PoseStamped MakePoseMessage(const SE3& pose, double stamp, const std::string& frame_id) {
    geometry_msgs::msg::PoseStamped message;
    message.header.stamp = ToRosStamp(stamp);
    message.header.frame_id = frame_id;
    const auto quaternion = pose.unit_quaternion().normalized();
    message.pose.position.x = pose.translation().x();
    message.pose.position.y = pose.translation().y();
    message.pose.position.z = pose.translation().z();
    message.pose.orientation.x = quaternion.x();
    message.pose.orientation.y = quaternion.y();
    message.pose.orientation.z = quaternion.z();
    message.pose.orientation.w = quaternion.w();
    return message;
}

sensor_msgs::msg::PointCloud2 MakeCloudMessage(const CloudPtr& cloud, double stamp, const std::string& frame_id) {
    sensor_msgs::msg::PointCloud2 message;
    message.header.stamp = ToRosStamp(stamp);
    message.header.frame_id = frame_id;
    sensor_msgs::PointCloud2Modifier modifier(message);
    modifier.setPointCloud2Fields(6, "x", 1, sensor_msgs::msg::PointField::FLOAT32, "y", 1,
                                  sensor_msgs::msg::PointField::FLOAT32, "z", 1,
                                  sensor_msgs::msg::PointField::FLOAT32, "intensity", 1,
                                  sensor_msgs::msg::PointField::FLOAT32, "time", 1,
                                  sensor_msgs::msg::PointField::FLOAT64, "lidar_id", 1,
                                  sensor_msgs::msg::PointField::UINT8);
    const std::size_t size = cloud ? cloud->size() : 0;
    modifier.resize(size);
    sensor_msgs::PointCloud2Iterator<float> x(message, "x"), y(message, "y"), z(message, "z"),
        intensity(message, "intensity");
    sensor_msgs::PointCloud2Iterator<double> time(message, "time");
    sensor_msgs::PointCloud2Iterator<std::uint8_t> lidar_id(message, "lidar_id");
    if (cloud) {
        for (const auto& point : cloud->points) {
            *x = point.x;
            *y = point.y;
            *z = point.z;
            *intensity = point.intensity;
            *time = point.time * 1e-3;
            *lidar_id = point.lidar_id;
            ++x;
            ++y;
            ++z;
            ++intensity;
            ++time;
            ++lidar_id;
        }
    }
    return message;
}

Mat3d ReadMatrix3(const YAML::Node& node) {
    Mat3d matrix = Mat3d::Identity();
    if (!node) return matrix;
    const auto values = node.as<std::vector<double>>();
    if (values.size() != 9) {
        throw std::runtime_error("fasterlio.extrinsic_R must have 9 values");
    }
    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) {
            matrix(row, col) = values[row * 3 + col];
        }
    }
    return matrix;
}

Vec3d ReadVector3(const YAML::Node& node, const Vec3d& fallback, const char* label) {
    if (!node) return fallback;
    const auto values = node.as<std::vector<double>>();
    if (values.size() != 3) {
        throw std::runtime_error(std::string(label) + " must have 3 values");
    }
    return Vec3d(values[0], values[1], values[2]);
}

}  // namespace

LocSystem::LocSystem(LocSystem::Options options) : options_(options) {
    /// handle ctrl-c
    signal(SIGINT, lightning::debug::SigHandle);
}

LocSystem::~LocSystem() { loc_->Finish(); }

bool LocSystem::Init(const std::string &yaml_path) {
    loc::Localization::Options opt;
    opt.online_mode_ = true;
    loc_ = std::make_shared<loc::Localization>(opt);

    YAML_IO yaml(yaml_path);

    std::string map_path = yaml.GetValue<std::string>("system", "map_path");
    const YAML::Node root = YAML::LoadFile(yaml_path);
    map_frame_ = root["output"] && root["output"]["map_frame"] ? root["output"]["map_frame"].as<std::string>() : "map";
    lidar_frame_ =
        root["output"] && root["output"]["lidar_frame"] ? root["output"]["lidar_frame"].as<std::string>() : "lidar_114";
    const Vec3d primary_lidar_position =
        ReadVector3(root["output"] ? root["output"]["primary_lidar_position_in_body"] : YAML::Node(),
                    Vec3d(2.199, 0.0, 2.740), "output.primary_lidar_position_in_body");
    const Mat3d R_lidar_to_imu =
        ReadMatrix3(root["fasterlio"] ? root["fasterlio"]["extrinsic_R"] : YAML::Node());
    const Vec3d t_lidar_to_imu =
        ReadVector3(root["fasterlio"] ? root["fasterlio"]["extrinsic_T"] : YAML::Node(), Vec3d::Zero(),
                    "fasterlio.extrinsic_T");
    rear_axle_ = RearAxlePoseTransformer(R_lidar_to_imu, t_lidar_to_imu, primary_lidar_position);

    LOG(INFO) << "online mode, creating ros2 node ... ";

    /// subscribers
    node_ = std::make_shared<rclcpp::Node>("lightning_slam");

    imu_topic_ = yaml.GetValue<std::string>("common", "imu_topic");
    cloud_topic_ = yaml.GetValue<std::string>("common", "lidar_topic");
    livox_topic_ = yaml.GetValue<std::string>("common", "livox_lidar_topic");

    rclcpp::QoS qos(10);

    imu_sub_ = node_->create_subscription<sensor_msgs::msg::Imu>(
        imu_topic_, qos, [this](sensor_msgs::msg::Imu::SharedPtr msg) {
            IMUPtr imu = std::make_shared<IMU>();
            imu->timestamp = ToSec(msg->header.stamp);
            imu->linear_acceleration =
                Vec3d(msg->linear_acceleration.x, msg->linear_acceleration.y, msg->linear_acceleration.z);
            imu->angular_velocity = Vec3d(msg->angular_velocity.x, msg->angular_velocity.y, msg->angular_velocity.z);

            ProcessIMU(imu);
        });

    cloud_sub_ = node_->create_subscription<sensor_msgs::msg::PointCloud2>(
        cloud_topic_, qos, [this](sensor_msgs::msg::PointCloud2::SharedPtr cloud) {
            Timer::Evaluate([&]() { ProcessLidar(cloud); }, "Proc Lidar", true);
        });

    livox_sub_ = node_->create_subscription<livox_ros_driver2::msg::CustomMsg>(
        livox_topic_, qos, [this](livox_ros_driver2::msg::CustomMsg ::SharedPtr cloud) {
            Timer::Evaluate([&]() { ProcessLidar(cloud); }, "Proc Lidar", true);
        });

    pose_pub_ = node_->create_publisher<geometry_msgs::msg::PoseStamped>("/slamPoseRaw_topic", 10);
    cloud_pub_ = node_->create_publisher<sensor_msgs::msg::PointCloud2>("/final_points_topic", 10);
    safety_pub_ = node_->create_publisher<std_msgs::msg::Float64>("/slamSafety_topic", 10);
    state_pub_ = node_->create_publisher<std_msgs::msg::Float64>("/slamState_topic", 10);
    system_pub_ = node_->create_publisher<std_msgs::msg::Int32>("/SystemState", 10);
    state_timer_ = node_->create_wall_timer(std::chrono::milliseconds(100), [this]() { PublishStateTopics(); });

    if (options_.pub_tf_) {
        tf_broadcaster_ = std::make_shared<tf2_ros::TransformBroadcaster>(node_);
        loc_->SetTFCallback(
            [this](const geometry_msgs::msg::TransformStamped &pose) { tf_broadcaster_->sendTransform(pose); });
    }
    loc_->SetLocalizationResultCallback([this](const loc::LocalizationResult& result) {
        PublishLocalizationResult(result);
    });
    loc_->SetProcessedCloudCallback(
        [this](const CloudPtr& cloud, const loc::LocalizationResult& result) { PublishProcessedCloud(cloud, result); });

    bool ret = loc_->Init(yaml_path, map_path);
    if (ret) {
        LOG(INFO) << "online loc node has been created.";
    }

    return ret;
}

void LocSystem::SetInitPose(const SE3 &pose) {
    LOG(INFO) << "set init pose: " << pose.translation().transpose() << ", "
              << pose.unit_quaternion().coeffs().transpose();

    loc_->SetExternalPose(pose.unit_quaternion(), pose.translation());
    loc_started_ = true;
}

void LocSystem::ProcessIMU(const IMUPtr &imu) {
    if (loc_started_) {
        loc_->ProcessIMUMsg(imu);
    }
}

void LocSystem::ProcessLidar(const sensor_msgs::msg::PointCloud2::SharedPtr &cloud) {
    if (loc_started_) {
        loc_->ProcessLidarMsg(cloud);
    }
}

void LocSystem::ProcessLidar(const livox_ros_driver2::msg::CustomMsg::SharedPtr &cloud) {
    if (loc_started_) {
        loc_->ProcessLivoxLidarMsg(cloud);
    }
}

void LocSystem::Spin() {
    if (node_ != nullptr) {
        spin(node_);
    }
}

void LocSystem::PublishLocalizationResult(const loc::LocalizationResult& result) {
    if (!pose_pub_ || !result.valid_ || result.timestamp_ <= 0.0) return;
    const NavState state = result.ToNavState();
    pose_pub_->publish(MakePoseMessage(rear_axle_.Transform(state), result.timestamp_, map_frame_));
    system_initialized_.store(true);
    tracking_normal_.store(result.status_ == loc::LocalizationStatus::GOOD);
    has_output_.store(true);
    last_output_wall_ns_.store(SteadyNowNs());
}

void LocSystem::PublishProcessedCloud(const CloudPtr& cloud, const loc::LocalizationResult& result) {
    if (!cloud_pub_ || !cloud || cloud->empty() || !IsUsableResult(result)) return;
    const double stamp = result.timestamp_ > 0.0 ? result.timestamp_ : CloudStampSec(cloud);
    if (stamp <= 0.0) return;
    cloud_pub_->publish(MakeCloudMessage(cloud, stamp, lidar_frame_));
}

void LocSystem::PublishStateTopics() {
    heartbeat_ = !heartbeat_;
    std_msgs::msg::Float64 safety;
    safety.data = heartbeat_ ? 1.0 : 0.0;
    safety_pub_->publish(safety);

    bool tracking_normal = tracking_normal_.load();
    if (has_output_.load() && SteadyNowNs() - last_output_wall_ns_.load() > 500000000LL) {
        tracking_normal = false;
    }
    std_msgs::msg::Float64 state;
    state.data = tracking_normal ? 1.0 : 0.0;
    state_pub_->publish(state);

    std_msgs::msg::Int32 system_state;
    system_state.data = system_initialized_.load() ? 1 : 0;
    system_pub_->publish(system_state);
}

}  // namespace lightning
