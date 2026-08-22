#include <gflags/gflags.h>
#include <geometry_msgs/msg/pose_stamped.hpp>
#include <glog/logging.h>
#include <rclcpp/rclcpp.hpp>
#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <sensor_msgs/point_cloud2_iterator.hpp>
#include <std_msgs/msg/float64.hpp>
#include <std_msgs/msg/int32.hpp>
#include <yaml-cpp/yaml.h>

#include <cmath>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "core/lio/laser_mapping.h"
#include "core/lio/rear_axle_pose.h"
#include "wrapper/ros_utils.h"

DEFINE_string(config, "./config/default.yaml", "frontend configuration yaml");

namespace {

builtin_interfaces::msg::Time ToRosStamp(double seconds) {
    const std::int64_t nanoseconds = static_cast<std::int64_t>(std::llround(seconds * 1e9));
    builtin_interfaces::msg::Time stamp;
    stamp.sec = static_cast<std::int32_t>(nanoseconds / 1000000000LL);
    stamp.nanosec = static_cast<std::uint32_t>(nanoseconds % 1000000000LL);
    return stamp;
}

sensor_msgs::msg::PointCloud2 MakeCloudMessage(const lightning::CloudPtr& cloud, double begin_time, double end_time,
                                                const std::string& frame_id) {
    sensor_msgs::msg::PointCloud2 message;
    message.header.stamp = ToRosStamp(end_time);
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
            // Public time is seconds relative to the message's scan-end stamp: earlier points are negative.
            *time = point.time * 1e-3 - (end_time - begin_time);
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

geometry_msgs::msg::PoseStamped MakePoseMessage(const lightning::SE3& pose, double stamp,
                                                 const std::string& frame_id) {
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

class FrontendNode : public rclcpp::Node {
   public:
    explicit FrontendNode(const std::string& config_path) : Node("multi_lidar_frontend") {
        if (!lio_.Init(config_path)) throw std::runtime_error("failed to initialize LaserMapping");
        const YAML::Node root = YAML::LoadFile(config_path);
        map_frame_ = root["output"] && root["output"]["map_frame"]
                         ? root["output"]["map_frame"].as<std::string>()
                         : "map";
        lidar_frame_ = root["output"] && root["output"]["lidar_frame"]
                           ? root["output"]["lidar_frame"].as<std::string>()
                           : "lidar_114";
        lightning::Vec3d lidar_position(2.199, 0.0, 2.740);
        if (root["output"] && root["output"]["primary_lidar_position_in_body"]) {
            const auto values = root["output"]["primary_lidar_position_in_body"].as<std::vector<double>>();
            if (values.size() != 3) throw std::runtime_error("primary_lidar_position_in_body must have 3 values");
            lidar_position = lightning::Vec3d(values[0], values[1], values[2]);
        }
        rear_axle_ = lightning::RearAxlePoseTransformer(lio_.GetLidarToImuRotation(),
                                                        lio_.GetLidarToImuTranslation(), lidar_position);

        const auto output_qos =
            rclcpp::QoS(rclcpp::KeepLast(10)).best_effort().durability_volatile();
        pose_pub_ = create_publisher<geometry_msgs::msg::PoseStamped>("/slamPoseRaw_topic", output_qos);
        cloud_pub_ = create_publisher<sensor_msgs::msg::PointCloud2>("/final_points_topic", output_qos);
        safety_pub_ = create_publisher<std_msgs::msg::Float64>("/slamSafety_topic", output_qos);
        state_pub_ = create_publisher<std_msgs::msg::Float64>("/slamState_topic", output_qos);
        system_pub_ = create_publisher<std_msgs::msg::Int32>("/SystemState", output_qos);

        const auto qos = rclcpp::SensorDataQoS();
        if (lio_.IsMultiLidarEnabled()) {
            const auto* primary = lio_.GetMultiLidarConfig().PrimaryLidar();
            if (!primary) throw std::runtime_error("primary lidar configuration is missing");
            imu_sub_ = create_subscription<sensor_msgs::msg::Imu>(
                primary->imu_topic, qos, [this](sensor_msgs::msg::Imu::SharedPtr message) { OnImu(message); });
            for (const auto& sensor : lio_.GetMultiLidarConfig().lidars) {
                lidar_subs_.push_back(create_subscription<sensor_msgs::msg::PointCloud2>(
                    sensor.lidar_topic, qos,
                    [this, id = sensor.id](sensor_msgs::msg::PointCloud2::SharedPtr message) {
                        lio_.ProcessPointCloud2(message, id);
                        Drain();
                    }));
            }
        } else {
            const std::string lidar_topic = root["common"]["lidar_topic"].as<std::string>();
            const std::string imu_topic = root["common"]["imu_topic"].as<std::string>();
            imu_sub_ = create_subscription<sensor_msgs::msg::Imu>(
                imu_topic, qos, [this](sensor_msgs::msg::Imu::SharedPtr message) { OnImu(message); });
            lidar_subs_.push_back(create_subscription<sensor_msgs::msg::PointCloud2>(
                lidar_topic, qos, [this](sensor_msgs::msg::PointCloud2::SharedPtr message) {
                    lio_.ProcessPointCloud2(message);
                    Drain();
                }));
        }
        publishing_.store(true);
        output_publish_thread_ = std::thread([this]() { OutputPublishLoop(); });
        state_publish_thread_ = std::thread([this]() { StatePublishLoop(); });
    }

    ~FrontendNode() override {
        publishing_.store(false);
        {
            std::lock_guard<std::mutex> lock(state_mutex_);
            output_queue_.clear();
        }
        output_ready_.notify_all();
        if (output_publish_thread_.joinable()) output_publish_thread_.join();
        if (state_publish_thread_.joinable()) state_publish_thread_.join();
    }

   private:
    struct OutputSnapshot {
        lightning::SE3 rear_pose;
        lightning::CloudPtr cloud;
        double begin_time = 0.0;
        double end_time = 0.0;
    };

    void OnImu(const sensor_msgs::msg::Imu::SharedPtr& message) {
        auto imu = std::make_shared<lightning::IMU>();
        imu->timestamp = lightning::ToSec(message->header.stamp);
        imu->linear_acceleration = lightning::Vec3d(message->linear_acceleration.x, message->linear_acceleration.y,
                                                     message->linear_acceleration.z);
        imu->angular_velocity = lightning::Vec3d(message->angular_velocity.x, message->angular_velocity.y,
                                                  message->angular_velocity.z);
        lio_.ProcessIMU(imu);
        Drain();
    }

    void Drain() {
        while (true) {
            const auto status = lio_.RunDetailed();
            if (status == lightning::LaserMapping::RunStatus::kNoData) break;
            const bool output = status == lightning::LaserMapping::RunStatus::kOutput;
            const bool initialized = lio_.IsInitialized();
            const bool healthy = lio_.IsTrackingHealthy();
            std::lock_guard<std::mutex> lock(state_mutex_);
            system_initialized_ = system_initialized_ || initialized;
            const bool good = system_initialized_ && output && healthy;
            if (good) {
                ++consecutive_good_;
                consecutive_bad_ = 0;
                if (consecutive_good_ >= 5) tracking_normal_ = true;
            } else {
                ++consecutive_bad_;
                consecutive_good_ = 0;
                if (consecutive_bad_ >= 3) tracking_normal_ = false;
            }
            if (output) {
                const auto state = lio_.GetState();
                if (state.pose_is_ok_ && !output_publish_failed_) {
                    constexpr std::size_t kMaxOutputQueue = 50;
                    if (output_queue_.size() >= kMaxOutputQueue) {
                        output_queue_.pop_front();
                        tracking_normal_ = false;
                        LOG_EVERY_N(WARNING, 100) << "online output publisher is falling behind; dropping oldest frame";
                    }
                    output_queue_.push_back(OutputSnapshot{rear_axle_.Transform(state), lio_.GetScanUndist(),
                                                           lio_.GetLastFrameBeginTime(),
                                                           lio_.GetLastFrameEndTime()});
                    last_output_wall_time_ = std::chrono::steady_clock::now();
                    has_output_ = true;
                    output_ready_.notify_one();
                }
            }
        }
    }

    void OutputPublishLoop() {
        try {
            while (publishing_.load() && rclcpp::ok()) {
                OutputSnapshot snapshot;
                {
                    std::unique_lock<std::mutex> lock(state_mutex_);
                    output_ready_.wait(lock, [this]() {
                        return !publishing_.load() || !rclcpp::ok() || !output_queue_.empty();
                    });
                    if (!publishing_.load() || !rclcpp::ok()) {
                        output_queue_.clear();
                        break;
                    }
                    snapshot = std::move(output_queue_.front());
                    output_queue_.pop_front();
                }
                pose_pub_->publish(MakePoseMessage(snapshot.rear_pose, snapshot.end_time, map_frame_));
                cloud_pub_->publish(
                    MakeCloudMessage(snapshot.cloud, snapshot.begin_time, snapshot.end_time, lidar_frame_));
            }
        } catch (const std::exception& error) {
            if (publishing_.load() && rclcpp::ok()) {
                {
                    std::lock_guard<std::mutex> lock(state_mutex_);
                    output_publish_failed_ = true;
                    tracking_normal_ = false;
                    output_queue_.clear();
                }
                LOG(ERROR) << "output publishing thread failed: " << error.what();
            }
        }
    }

    void StatePublishLoop() {
        auto next = std::chrono::steady_clock::now();
        while (publishing_.load() && rclcpp::ok()) {
            next += std::chrono::milliseconds(100);
            bool tracking_normal = false;
            bool system_initialized = false;
            {
                std::lock_guard<std::mutex> lock(state_mutex_);
                if (has_output_ && std::chrono::steady_clock::now() - last_output_wall_time_ >
                                       std::chrono::milliseconds(300)) {
                    tracking_normal_ = false;
                }
                tracking_normal = tracking_normal_;
                if (output_publish_failed_) tracking_normal = false;
                system_initialized = system_initialized_;
                heartbeat_ = !heartbeat_;
            }
            try {
                std_msgs::msg::Float64 safety;
                safety.data = heartbeat_ ? 1.0 : 0.0;
                safety_pub_->publish(safety);
                std_msgs::msg::Int32 system_message;
                system_message.data = system_initialized ? 1 : 0;
                system_pub_->publish(system_message);
                std_msgs::msg::Float64 state_message;
                state_message.data = tracking_normal ? 1.0 : 0.0;
                state_pub_->publish(state_message);
            } catch (const std::exception& error) {
                if (publishing_.load() && rclcpp::ok()) LOG(ERROR) << "state publishing thread failed: " << error.what();
                break;
            }
            std::this_thread::sleep_until(next);
        }
    }

    lightning::LaserMapping lio_;
    lightning::RearAxlePoseTransformer rear_axle_;
    std::string map_frame_;
    std::string lidar_frame_;
    std::vector<rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr> lidar_subs_;
    rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr imu_sub_;
    rclcpp::Publisher<geometry_msgs::msg::PoseStamped>::SharedPtr pose_pub_;
    rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr cloud_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr safety_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr state_pub_;
    rclcpp::Publisher<std_msgs::msg::Int32>::SharedPtr system_pub_;
    std::mutex state_mutex_;
    std::condition_variable output_ready_;
    std::deque<OutputSnapshot> output_queue_;
    std::atomic_bool publishing_{false};
    std::thread output_publish_thread_;
    std::thread state_publish_thread_;
    std::chrono::steady_clock::time_point last_output_wall_time_ = std::chrono::steady_clock::now();
    int consecutive_good_ = 0;
    int consecutive_bad_ = 0;
    bool system_initialized_ = false;
    bool tracking_normal_ = false;
    bool output_publish_failed_ = false;
    bool heartbeat_ = false;
    bool has_output_ = false;
};

}  // namespace

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_colorlogtostderr = true;
    FLAGS_stderrthreshold = google::INFO;
    google::ParseCommandLineFlags(&argc, &argv, true);
    rclcpp::init(argc, argv);
    try {
        rclcpp::spin(std::make_shared<FrontendNode>(FLAGS_config));
    } catch (const std::exception& e) {
        LOG(ERROR) << e.what();
        rclcpp::shutdown();
        return 2;
    }
    rclcpp::shutdown();
    return 0;
}
