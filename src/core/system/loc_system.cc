//
// Created by xiang on 25-9-12.
//

#include "core/system/loc_system.h"

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

bool IsUsableResult(const loc::LocalizationResult& result) {
    return result.valid_ || result.lidar_loc_valid_;
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
    rear_axle_frame_ = root["output"] && root["output"]["rear_axle_frame"]
                           ? root["output"]["rear_axle_frame"].as<std::string>()
                           : "rear_axle";
    const Vec3d primary_lidar_position =
        ReadVector3(root["output"] ? root["output"]["primary_lidar_position_in_body"] : YAML::Node(),
                    Vec3d(2.199, 0.0, 2.740), "output.primary_lidar_position_in_body");
    const Mat3d R_lidar_to_imu =
        ReadMatrix3(root["fasterlio"] ? root["fasterlio"]["extrinsic_R"] : YAML::Node());
    const Vec3d t_lidar_to_imu =
        ReadVector3(root["fasterlio"] ? root["fasterlio"]["extrinsic_T"] : YAML::Node(), Vec3d::Zero(),
                    "fasterlio.extrinsic_T");
    rear_axle_ = RearAxlePoseTransformer(R_lidar_to_imu, t_lidar_to_imu, primary_lidar_position);
    T_rear_lidar_ = SE3(SO3(), primary_lidar_position);

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

    const auto pose_qos = rclcpp::QoS(rclcpp::KeepLast(1000));
    const auto cloud_qos = rclcpp::SensorDataQoS().keep_last(1);
    pos_res_pub_ = node_->create_publisher<geosun_msgs::msg::PosRes>("/PosRes", pose_qos);
    pose_pub_ = node_->create_publisher<geometry_msgs::msg::PoseStamped>("/slamPoseRaw_topic", pose_qos);
    inv_cloud_pub_ = node_->create_publisher<sensor_msgs::msg::PointCloud2>("/LidarDataInv", cloud_qos);
    map_cloud_pub_ = node_->create_publisher<sensor_msgs::msg::PointCloud2>("/LidarDataInL", cloud_qos);

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
    if (!pos_res_pub_ || !pose_pub_ || !result.valid_ || result.timestamp_ <= 0.0) return;
    const NavState state = result.ToNavState();
    const auto position =
        sany_output::MakePosResMessage(rear_axle_.Transform(state), result.vel_b_.x(), result.timestamp_, map_frame_);
    pos_res_pub_->publish(position);
    pose_pub_->publish(sany_output::MakePoseMessage(position));
}

void LocSystem::PublishProcessedCloud(const CloudPtr& cloud, const loc::LocalizationResult& result) {
    if (!inv_cloud_pub_ || !map_cloud_pub_ || !cloud || cloud->empty() || !IsUsableResult(result)) return;
    const double begin_time = CloudStampSec(cloud);
    const double end_time = result.timestamp_ > 0.0 ? result.timestamp_ : begin_time;
    if (begin_time <= 0.0 || end_time <= 0.0) return;
    const SE3 rear_axle_pose = rear_axle_.Transform(result.ToNavState());
    inv_cloud_pub_->publish(
        sany_output::MakeCloudMessage(cloud, begin_time, end_time, T_rear_lidar_, rear_axle_frame_));
    if (map_cloud_decimator_.Tick()) {
        map_cloud_pub_->publish(sany_output::MakeCloudMessage(cloud, begin_time, end_time,
                                                              rear_axle_pose * T_rear_lidar_, map_frame_));
    }
}

}  // namespace lightning
