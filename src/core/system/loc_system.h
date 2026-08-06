//
// Created by xiang on 25-9-8.
//

#ifndef LIGHTNING_LOC_SYSTEM_H
#define LIGHTNING_LOC_SYSTEM_H

#include <tf2_ros/transform_broadcaster.h>
#include <geometry_msgs/msg/pose_stamped.hpp>
#include <geosun_msgs/msg/pos_res.hpp>
#include <nav_msgs/msg/path.hpp>
#include <rclcpp/rclcpp.hpp>
#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <memory>
#include <mutex>
#include <vector>

#include "livox_ros_driver2/msg/custom_msg.hpp"

#include "common/eigen_types.h"
#include "common/imu.h"
#include "common/keyframe.h"
#include "core/localization/localization_result.h"
#include "core/system/sany_localization_output.h"

namespace lightning {

namespace loc {
class Localization;
}

class LocSystem {
   public:
    struct Options {
        bool pub_tf_ = true;  // 是否发布tf
    };

    explicit LocSystem(Options options);
    ~LocSystem();

    /// 初始化，地图路径在yaml里配置
    bool Init(const std::string& yaml_path, const std::string& map_path_override = "");

    /// 设置初始化位姿
    void SetInitPose(const SE3& pose);

    /// 处理IMU
    void ProcessIMU(const lightning::IMUPtr& imu);

    /// 处理点云
    void ProcessLidar(const sensor_msgs::msg::PointCloud2::SharedPtr& cloud);
    void ProcessLidar(const sensor_msgs::msg::PointCloud2::SharedPtr& cloud, int lidar_id);
    void ProcessLidar(const livox_ros_driver2::msg::CustomMsg::SharedPtr& cloud);

    /// 实时模式下的spin
    void Spin();
    /// Stop workers after draining all accepted sensor and localization data.
    void Finish();
    /// PGO global corrections, matching the offline trajectory definition.
    bool SaveTrajectoryTum(const std::string& path) const;
    /// Live high-frequency extrapolated output published to ROS.
    bool SaveHighFrequencyTrajectoryTum(const std::string& path) const;

   private:
    void PublishLocalizationResult(const loc::LocalizationResult& result);
    void CaptureGlobalLocalizationResult(const loc::LocalizationResult& result);
    bool WriteTrajectoryTum(const std::string& path, const std::vector<NavState>& states,
                            const char* description) const;
    void PublishProcessedCloud(const CloudPtr& cloud, const loc::LocalizationResult& result);
    void PublishHealthStatus();
    void PublishPath();

    Options options_;

    std::shared_ptr<loc::Localization> loc_ = nullptr;  // 定位接口

    std::atomic_bool loc_started_ = false;  // 是否开启定位
    std::atomic_bool map_loaded_ = false;   // 地图是否已载入

    /// 实时模式下的ros2 node, subscribers
    rclcpp::Node::SharedPtr node_;
    std::shared_ptr<tf2_ros::TransformBroadcaster> tf_broadcaster_ = nullptr;

    std::string imu_topic_;
    std::string cloud_topic_;
    std::string livox_topic_;
    std::string map_frame_ = "map";
    std::string rear_axle_frame_ = "rear_axle";
    Vec3d primary_lidar_position_in_body_ = Vec3d::Zero();
    sany_output::LocalizationPublicationGate publication_gate_;
    sany_output::FrameDecimator map_cloud_decimator_{10};
    std::unique_ptr<sany_output::LocalizationTelemetryState> telemetry_;
    mutable std::mutex trajectory_mutex_;
    std::vector<NavState> localization_states_;
    std::vector<NavState> global_localization_states_;
    bool finished_ = false;

    rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr imu_sub_ = nullptr;
    std::vector<rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr> cloud_subs_;
    rclcpp::Subscription<livox_ros_driver2::msg::CustomMsg>::SharedPtr livox_sub_ = nullptr;

    rclcpp::Publisher<geosun_msgs::msg::PosRes>::SharedPtr pos_res_pub_ = nullptr;
    rclcpp::Publisher<geometry_msgs::msg::PoseStamped>::SharedPtr pose_pub_ = nullptr;
    rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr inv_cloud_pub_ = nullptr;
    rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr map_cloud_pub_ = nullptr;
    rclcpp::Publisher<lightning::msg::FaultStatus>::SharedPtr fault_status_pub_ = nullptr;
    rclcpp::Publisher<lightning::msg::LocalizationStatus>::SharedPtr loc_status_pub_ = nullptr;
    rclcpp::Publisher<nav_msgs::msg::Path>::SharedPtr path_pub_ = nullptr;
    rclcpp::TimerBase::SharedPtr health_timer_ = nullptr;
    rclcpp::TimerBase::SharedPtr path_timer_ = nullptr;
};

};  // namespace lightning

#endif  // LIGHTNING_LOC_SYSTEM_H
