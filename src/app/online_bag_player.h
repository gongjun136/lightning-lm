#ifndef LIGHTNING_APP_ONLINE_BAG_PLAYER_H
#define LIGHTNING_APP_ONLINE_BAG_PLAYER_H

#include <glog/logging.h>
#include <rclcpp/rclcpp.hpp>
#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <yaml-cpp/yaml.h>

#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <thread>

#include "geosun_msgs/msg/spe_thr_can4.hpp"
#include "livox_ros_driver2/msg/custom_msg.hpp"
#include "wrapper/bag_io.h"

namespace lightning {

/// A deterministic 1x bag publisher for online regression. It shares the
/// process-wide RMW participant with the system under test, while messages
/// still traverse normal ROS 2 publisher/subscription callbacks.
class OnlineBagPlayer {
   public:
    bool Init(const std::string& config_path, const std::string& bag_path, double playback_rate) {
        if (bag_path.empty()) return false;
        playback_rate_ = playback_rate;
        bag_ = std::make_unique<RosbagIO>(bag_path);
        node_ = std::make_shared<rclcpp::Node>("lightning_online_bag_player");
        const YAML::Node root = YAML::LoadFile(config_path);
        const auto imu_qos = rclcpp::SensorDataQoS().keep_last(1000);
        const YAML::Node lidar_qos_config = root["system"] ? root["system"]["lidar_qos"] : YAML::Node();
        const size_t lidar_qos_depth =
            lidar_qos_config && lidar_qos_config["depth"]
                ? lidar_qos_config["depth"].as<size_t>()
                : 1000;
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

        const std::string imu_topic = root["common"]["imu_topic"].as<std::string>();
        imu_pub_ = node_->create_publisher<sensor_msgs::msg::Imu>(imu_topic, imu_qos);
        bag_->AddRosImuHandle(imu_topic, [this](const sensor_msgs::msg::Imu::SharedPtr msg) {
            imu_pub_->publish(*msg);
            return rclcpp::ok();
        });

        const YAML::Node system = root["system"];
        bool wheel_speed_observation_enabled =
            system && system["enable_wheel_speed_observation"]
                ? system["enable_wheel_speed_observation"].as<bool>()
                : true;
        if (const char* override_value = std::getenv("SANY_ENABLE_CAN_OBSERVATION")) {
            const std::string value(override_value);
            if (value == "0") {
                wheel_speed_observation_enabled = false;
            } else if (value == "1") {
                wheel_speed_observation_enabled = true;
            } else {
                LOG(ERROR) << "SANY_ENABLE_CAN_OBSERVATION must be 0 or 1";
                return false;
            }
        }
        const std::string wheel_speed_topic =
            system && system["wheel_speed_topic"]
                ? system["wheel_speed_topic"].as<std::string>()
                : "/SpeThrCAN4_topic";
        const std::string selected_wheel_speed_topic =
            std::getenv("SANY_WHEEL_SPEED_TOPIC")
                ? std::getenv("SANY_WHEEL_SPEED_TOPIC")
                : wheel_speed_topic;
        if (wheel_speed_observation_enabled && !selected_wheel_speed_topic.empty()) {
            wheel_speed_pub_ = node_->create_publisher<geosun_msgs::msg::SpeThrCAN4>(
                selected_wheel_speed_topic, rclcpp::QoS(rclcpp::KeepLast(50)).reliable());
            bag_->AddHandle(selected_wheel_speed_topic, [this](const RosbagIO::MsgType& serialized) {
                auto msg = std::make_shared<geosun_msgs::msg::SpeThrCAN4>();
                rclcpp::SerializedMessage data(*serialized->serialized_data);
                wheel_speed_serialization_.deserialize_message(&data, msg.get());
                wheel_speed_pub_->publish(*msg);
                return rclcpp::ok();
            });
        }

        const bool multi_lidar = root["multi_lidar"] && root["multi_lidar"]["enabled"] &&
                                 root["multi_lidar"]["enabled"].as<bool>();
        if (multi_lidar) {
            const YAML::Node topics = root["multi_lidar"]["topics"];
            if (!topics || !topics.IsMap()) {
                LOG(ERROR) << "multi_lidar.topics is required for embedded online playback";
                return false;
            }
            for (const auto& entry : topics) {
                const std::string key = entry.first.as<std::string>();
                if (key.rfind("lidar_", 0) != 0) continue;
                const std::string topic = entry.second.as<std::string>();
                auto publisher = node_->create_publisher<sensor_msgs::msg::PointCloud2>(topic, lidar_qos);
                cloud_pubs_[topic] = publisher;
                bag_->AddPointCloud2Handle(topic, [publisher](const sensor_msgs::msg::PointCloud2::SharedPtr msg) {
                    publisher->publish(*msg);
                    return rclcpp::ok();
                });
            }
            if (cloud_pubs_.empty()) {
                LOG(ERROR) << "no multi-lidar PointCloud2 topics configured";
                return false;
            }
        } else {
            const std::string cloud_topic = root["common"]["lidar_topic"].as<std::string>();
            const std::string livox_topic = root["common"]["livox_lidar_topic"].as<std::string>();
            if (!cloud_topic.empty()) {
                cloud_pub_ = node_->create_publisher<sensor_msgs::msg::PointCloud2>(cloud_topic, lidar_qos);
                bag_->AddPointCloud2Handle(cloud_topic, [this](const sensor_msgs::msg::PointCloud2::SharedPtr msg) {
                    cloud_pub_->publish(*msg);
                    return rclcpp::ok();
                });
            }
            if (!livox_topic.empty()) {
                livox_pub_ = node_->create_publisher<livox_ros_driver2::msg::CustomMsg>(livox_topic, lidar_qos);
                bag_->AddLivoxCloudHandle(livox_topic,
                                          [this](const livox_ros_driver2::msg::CustomMsg::SharedPtr msg) {
                                              livox_pub_->publish(*msg);
                                              return rclcpp::ok();
                                          });
            }
            if (!cloud_pub_ && !livox_pub_) {
                LOG(ERROR) << "no single-lidar topic configured for embedded online playback";
                return false;
            }
        }
        return true;
    }

    void Run() {
        // Allow the executor to spin once before the first sensor message.
        std::this_thread::sleep_for(std::chrono::seconds(1));
        bag_->GoRealtime(playback_rate_);
    }

   private:
    double playback_rate_ = 1.0;
    std::unique_ptr<RosbagIO> bag_;
    rclcpp::Node::SharedPtr node_;
    rclcpp::Publisher<sensor_msgs::msg::Imu>::SharedPtr imu_pub_;
    rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr cloud_pub_;
    rclcpp::Publisher<livox_ros_driver2::msg::CustomMsg>::SharedPtr livox_pub_;
    rclcpp::Publisher<geosun_msgs::msg::SpeThrCAN4>::SharedPtr wheel_speed_pub_;
    rclcpp::Serialization<geosun_msgs::msg::SpeThrCAN4> wheel_speed_serialization_;
    std::map<std::string, rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr> cloud_pubs_;
};

}  // namespace lightning

#endif  // LIGHTNING_APP_ONLINE_BAG_PLAYER_H
