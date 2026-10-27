//
// Created by xiang on 23-12-14.
//

#ifndef LIGHTNING_BAG_IO_H
#define LIGHTNING_BAG_IO_H

#include <functional>
#include <map>
#include <string>

#include <nav_msgs/msg/odometry.hpp>
#include <rclcpp/serialization.hpp>
#include <rclcpp/serialized_message.hpp>
#include <rosbag2_cpp/reader.hpp>
#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/laser_scan.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>

#include "livox_ros_driver2/msg/custom_msg.hpp"

#include "common/imu.h"
#include "common/odom.h"
#include "common/point_def.h"
#include "core/lightning_math.hpp"
#include "io/dataset_type.h"
#include "wrapper/ros_utils.h"

namespace lightning {

/**
 * @brief 在调用线程中顺序读取 ROS 2 bag，并按注册话题同步派发。
 * @note 回调 bool 返回值不控制 Go/GoRealtime；异常向调用者传播。
 * 构造参数 dataset_type 为兼容接口保留，当前未用于分发或消息转换。
 */
class RosbagIO {
   public:
    explicit RosbagIO(std::string bag_file, DatasetType dataset_type = DatasetType::NCLT)
        : bag_file_(std::move(bag_file)) {
        /// handle ctrl-c
        signal(SIGINT, lightning::debug::SigHandle);
    }

    using MsgType = std::shared_ptr<rosbag2_storage::SerializedBagMessage>;
    using MessageProcessFunction = std::function<bool(const MsgType &m)>;

    /// 一些方便直接使用的topics, messages
    using Scan2DHandle = std::function<bool(sensor_msgs::msg::LaserScan::SharedPtr)>;

    using PointCloud2Handle = std::function<bool(sensor_msgs::msg::PointCloud2::SharedPtr)>;
    using LivoxCloud2Handle = std::function<bool(livox_ros_driver2::msg::CustomMsg::SharedPtr)>;
    using RosImuHandle = std::function<bool(sensor_msgs::msg::Imu::SharedPtr)>;
    using FullPointCloudHandle = std::function<bool(FullCloudPtr)>;
    using ImuHandle = std::function<bool(IMUPtr)>;
    using OdomHandle = std::function<bool(const OdomPtr &)>;

    /**
     * @brief 顺序读取消息并同步派发；回调后检查全局退出标志。
     * @par 主要流程
     * 1. <b>打开 bag</b>：按文件类型选择存储插件，以 CDR 格式读取。
     * 2. <b>按注册话题派发</b>：逐条读出消息，有匹配回调时在当前线程执行。
     * 3. <b>等待并检查退出</b>：按配置固定等待，再检查退出标志；否则继续下一条。
     * @note 阶段名称对应函数体；阶段 2–3 对每条读出消息重复，回调 bool 返回值不控制回放。
     * @param sleep_usec 每条读出消息之后的固定等待微秒数；未注册话题也等待，非正值不等待。
     */
    void Go(int sleep_usec = 0);

    /// 注册话题回调；同一话题重复注册保留首次回调，不覆盖旧值。
    // [bag-topic-registration]
    RosbagIO &AddHandle(const std::string &topic_name, MessageProcessFunction func) {
        process_func_.emplace(topic_name, func);
        return *this;
    }
    // [bag-topic-registration]

    /// point cloud 2 处理
    RosbagIO &AddPointCloud2Handle(const std::string &topic_name, PointCloud2Handle f) {
        return AddHandle(topic_name, [f, this](const MsgType &m) -> bool {
            auto msg = std::make_shared<sensor_msgs::msg::PointCloud2>();
            rclcpp::SerializedMessage data(*m->serialized_data);
            seri_cloud2_.deserialize_message(&data, msg.get());

            return f(msg);
        });
    }

    /// livox 处理
    RosbagIO &AddLivoxCloudHandle(const std::string &topic_name, LivoxCloud2Handle f) {
        return AddHandle(topic_name, [f, this](const MsgType &m) -> bool {
            auto msg = std::make_shared<livox_ros_driver2::msg::CustomMsg>();
            rclcpp::SerializedMessage data(*m->serialized_data);
            seri_livox_.deserialize_message(&data, msg.get());

            return f(msg);
        });
    }

    RosbagIO &AddRosImuHandle(const std::string &topic_name, RosImuHandle f) {
        return AddHandle(topic_name, [f, this](const MsgType &m) -> bool {
            auto msg = std::make_shared<sensor_msgs::msg::Imu>();
            rclcpp::SerializedMessage data(*m->serialized_data);
            seri_imu_.deserialize_message(&data, msg.get());
            return f(msg);
        });
    }

    RosbagIO &AddImuHandle(const std::string &topic_name, ImuHandle f) {
        return AddHandle(topic_name, [f, this](const MsgType &m) -> bool {
            auto msg = std::make_shared<sensor_msgs::msg::Imu>();
            rclcpp::SerializedMessage data(*m->serialized_data);
            seri_imu_.deserialize_message(&data, msg.get());

            IMUPtr imu = std::make_shared<IMU>();
            imu->timestamp = ToSec(msg->header.stamp);

            /// NOTE: 如果需要乘重力，请修改此处
            imu->linear_acceleration =
                Vec3d(msg->linear_acceleration.x, msg->linear_acceleration.y, msg->linear_acceleration.z);
            // constant::kGRAVITY;
            imu->angular_velocity = Vec3d(msg->angular_velocity.x, msg->angular_velocity.y, msg->angular_velocity.z);

            return f(imu);
        });
    }

    /// odom 处理
    // RosbagIO &AddOdomHandle(const std::string &topic_name, OdomHandle f) {
    //     return AddHandle(topic_name, [f, this](const MsgType &m) -> bool {
    //         auto msg = std::make_shared<nav_msgs::msg::Odometry>();
    //         rclcpp::SerializedMessage data(*m->serialized_data);
    //         seri_odom_.deserialize_message(&data, msg.get());

    //         /// nav_msg 的 odometry 转 odom
    //         return f(msg);
    //     });
    // }

    /// 清除现有的处理函数
    void CleanProcessFunc() { process_func_.clear(); }

    /**
     * @brief 按 bag 存储时间与正倍率安排 steady_clock 目标时间，保留消息 header。
     * @par 主要流程
     * 1. <b>校验倍率并打开 bag</b>：拒绝非有限或非正倍率，选择存储插件打开数据。
     * 2. <b>准备回放时钟</b>：准备首条 bag 时间和对应的单调时钟基准。
     * 3. <b>等待消息目标时刻</b>：首条消息建立基准，按存储时间差除以倍率安排等待。
     * 4. <b>派发消息并检查退出</b>：同步执行已注册回调，检查退出标志后继续读取。
     * @note 阶段名称对应函数体；阶段 3–4 逐条重复，回调耗时过长时不保证实时期限。
     */
    void GoRealtime(double playback_rate = 1.0);

   private:
    std::map<std::string, MessageProcessFunction> process_func_;  ///< 按 ROS 话题注册的消息处理回调表。

    /// 序列化
    rclcpp::Serialization<nav_msgs::msg::Odometry> seri_odom_;
    rclcpp::Serialization<sensor_msgs::msg::Imu> seri_imu_;  ///< IMU bag 消息的反序列化器。
    rclcpp::Serialization<sensor_msgs::msg::PointCloud2> seri_cloud2_;  ///< PointCloud2 bag 消息的反序列化器。
    rclcpp::Serialization<livox_ros_driver2::msg::CustomMsg> seri_livox_;  ///< Livox CustomMsg bag 消息的反序列化器。

    std::string bag_file_;  ///< 输入 rosbag 文件或目录路径。
    DatasetType dataset_type_ = DatasetType::NCLT;  ///< 兼容字段；当前回放与消息转换不读取它。
};
}  // namespace lightning

#endif  // SLAM_ROS_BAG_IO_H
