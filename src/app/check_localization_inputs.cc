#include <rclcpp/rclcpp.hpp>
#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>

#include <chrono>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "common/sany_wheel_speed_wire.h"

// Check deliveries, not cached ROS daemon topic names. No output is published.
int main(int argc, char** argv) {
    rclcpp::init(argc, argv);
    try {
        auto node = std::make_shared<rclcpp::Node>("localization_input_check");
        struct Input { size_t count = 0; double stamp = 0; double received = -1; };
        std::map<std::string, Input> inputs;
        std::vector<rclcpp::SubscriptionBase::SharedPtr> subscriptions;
        std::string primary_imu, wheel_topic;
        double timeout = 60;
        size_t malformed_can = 0;
        const auto start = std::chrono::steady_clock::now();
        const auto elapsed = [&]() {
            return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        };
        const auto observe = [&](const std::string& topic, double stamp) {
            auto& input = inputs.at(topic);
            if (!std::isfinite(stamp) || stamp <= input.stamp) return;
            ++input.count;
            input.stamp = stamp;
            input.received = elapsed();
        };
        const auto qos = rclcpp::QoS(1).best_effort().durability_volatile();
        for (int i = 1; i < argc; i += 2) {
            if (i + 1 == argc) throw std::runtime_error("missing input-check argument value");
            const std::string option = argv[i], topic = argv[i + 1];
            if (option == "--timeout") { timeout = std::stod(topic); continue; }
            if (topic.empty() || topic[0] != '/' || inputs.count(topic))
                throw std::runtime_error("invalid or duplicate topic: " + topic);
            inputs[topic] = {};
            if (option == "--lidar") {
                subscriptions.push_back(node->create_subscription<sensor_msgs::msg::PointCloud2>(
                    topic, qos, [&, topic](sensor_msgs::msg::PointCloud2::ConstSharedPtr msg) {
                        if (msg->width && msg->height && !msg->data.empty())
                            observe(topic, double(msg->header.stamp.sec) + msg->header.stamp.nanosec * 1e-9);
                    }));
            } else if (option == "--imu") {
                if (primary_imu.empty()) primary_imu = topic;
                subscriptions.push_back(node->create_subscription<sensor_msgs::msg::Imu>(
                    topic, qos, [&, topic](sensor_msgs::msg::Imu::ConstSharedPtr msg) {
                        observe(topic, double(msg->header.stamp.sec) + msg->header.stamp.nanosec * 1e-9);
                    }));
            } else if (option == "--wheel-speed") {
                wheel_topic = topic;
                subscriptions.push_back(node->create_generic_subscription(
                    topic, "geosun_msgs/msg/SpeThrCAN4", qos,
                    [&, topic](std::shared_ptr<rclcpp::SerializedMessage> msg) {
                        const auto& raw = msg->get_rcl_serialized_message();
                        lightning::SanyWheelSpeedWire sample;
                        if (!lightning::DecodeSanyWheelSpeed(raw.buffer, raw.buffer_length, sample)) {
                            ++malformed_can;
                            return;
                        }
                        if (inputs.at(topic).count == 0)
                            std::cout << "CAN wire layout: " << (sample.has_comm_header ? "comm_header+Header+x+y" : "Header+x+y") << std::endl;
                        observe(topic, sample.stamp);
                    }));
            } else throw std::runtime_error("unknown option: " + option);
        }
        if (inputs.empty() || primary_imu.empty() || !std::isfinite(timeout) || timeout <= 0)
            throw std::runtime_error("input check needs topics, a primary IMU and a positive timeout");
        auto report = [&]() {
            for (const auto& entry : inputs)
                std::cerr << "  " << entry.first << ": messages=" << entry.second.count
                          << ", arrival_age_sec=" << (entry.second.received < 0 ? -1 : elapsed() - entry.second.received) << '\n';
            if (!wheel_topic.empty())
                std::cerr << "  CAN-IMU stamp delta_sec=" << inputs.at(wheel_topic).stamp - inputs.at(primary_imu).stamp
                          << ", malformed_CAN=" << malformed_can << '\n';
        };
        double next_report = 5;
        while (rclcpp::ok() && elapsed() < timeout) {
            rclcpp::spin_some(node);
            bool ready = true;
            for (const auto& entry : inputs)
                ready = ready && entry.second.count >= 3 && elapsed() - entry.second.received < 1;
            if (!wheel_topic.empty())
                ready = ready && std::abs(inputs.at(wheel_topic).stamp - inputs.at(primary_imu).stamp) < 0.25;
            if (ready) {
                std::cout << "Sensor delivery check passed for " << inputs.size() << " topics." << std::endl;
                report();
                rclcpp::shutdown();
                return 0;
            }
            if (elapsed() >= next_report) { report(); next_report += 5; }
            rclcpp::sleep_for(std::chrono::milliseconds(10));
        }
        std::cerr << "ERROR: sensor delivery check timed out; check publishers, ROS_DOMAIN_ID, DDS profile and CAN message layout/timestamps.\n";
        report();
    } catch (const std::exception& error) {
        std::cerr << "ERROR: sensor delivery check: " << error.what() << '\n';
    }
    rclcpp::shutdown();
    return 1;
}
