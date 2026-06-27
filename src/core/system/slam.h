//
// Created by xiang on 25-5-6.
//

#ifndef LIGHTNING_SLAM_H
#define LIGHTNING_SLAM_H

#include <rclcpp/rclcpp.hpp>
#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <string>

#include "lightning/srv/save_map.hpp"
#include "livox_ros_driver2/msg/custom_msg.hpp"

#include "common/eigen_types.h"
#include "common/imu.h"
#include "common/keyframe.h"
#include "common/nav_state.h"

namespace lightning {

class LaserMapping;  //  lio 前端
class LoopClosing;   // 回环检测

namespace ui {
class PangolinWindow;
}

namespace g2p5 {
class G2P5;
}

/**
 * @brief SLAM系统的顶层调度接口。
 *
 * SlamSystem负责串联LIO前端、回环检测、2.5D栅格建图、可视化和ROS2在线输入输出。
 * 外部通常先调用Init()完成模块初始化，再通过StartSLAM()进入建图状态；
 * 在线模式下由ROS2订阅器驱动ProcessIMU()/ProcessLidar()，离线模式下可由rosbag回放代码直接调用。
 */
class SlamSystem {
   public:
    /**
     * @brief SLAM系统运行选项。
     *
     * 这些选项会在Init()中结合YAML配置进一步更新，用于决定是否启用回环、栅格建图和可视化模块。
     */
    struct Options {
        Options() {}

        bool online_mode_ = true;  ///< 在线模式；在线模式下会创建ROS2节点、订阅器和服务。

        bool with_cc_ = true;               ///< 是否启用交叉验证相关逻辑。
        bool with_gridmap_ = true;          ///< 是否启用2.5D/2D栅格地图模块。
        bool with_loop_closing_ = true;     ///< 是否启用回环检测模块。
        bool with_visualization_ = true;    ///< 是否启用3D可视化UI。
        bool with_2dvisualization_ = true;  ///< 是否启用2D栅格地图可视化窗口。

        bool step_on_kf_ = true;  ///< 是否在关键帧处暂停，主要用于离线调试。
    };

    using SaveMapService = srv::SaveMap;

    /**
     * @brief 构造SLAM系统。
     * @param options 初始运行选项，Init()会继续从YAML配置中更新部分开关。
     */
    SlamSystem(Options options);

    /**
     * @brief 析构SLAM系统并释放已创建的模块资源。
     */
    ~SlamSystem();

    /**
     * @brief 初始化SLAM系统及其子模块。
     * @param yaml_path YAML配置文件路径。
     * @return 初始化成功返回true，否则返回false。
     */
    bool Init(const std::string& yaml_path);

    /**
     * @brief 开始建图流程。
     * @param map_name 当前地图名称，后续保存地图时会作为默认目录名或标识。
     */
    void StartSLAM(std::string map_name);

    /**
     * @brief 保存当前地图。
     * @param path 地图保存路径；为空时默认保存到./data/地图名/目录下。
     */
    void SaveMap(const std::string& path = "");

    /**
     * @brief 处理一条IMU数据。
     * @param imu 已转换为项目内部格式的IMU数据指针。
     */
    void ProcessIMU(const lightning::IMUPtr& imu);

    /**
     * @brief 处理一帧点云数据。
     *
     * 该模板接口支持标准sensor_msgs::msg::PointCloud2和Livox CustomMsg。
     * 函数内部会先交给LIO前端预处理和运行，若生成新关键帧，再分发给回环、栅格建图和UI模块。
     *
     * @tparam PointCloudMsgType ROS2点云消息类型。
     * @param cloud 点云消息智能指针。
     */
    template <typename PointCloudMsgType>
    void ProcessLidar(const std::shared_ptr<PointCloudMsgType>& cloud);

    /**
     * @brief 获取当前LIO前端状态。
     *
     * 离线评估程序用该接口导出TUM轨迹；若LIO尚未初始化，返回pose_is_ok_=false的无效状态。
     */
    NavState GetLioState() const;

    /**
     * @brief 在线模式下启动ROS2事件循环。
     */
    void Spin();

   private:
    /**
     * @brief ROS2保存地图服务的回调实现。
     * @param request 保存地图服务请求。
     * @param response 保存地图服务响应。
     */
    void SaveMap(const SaveMapService::Request::SharedPtr request, SaveMapService::Response::SharedPtr response);

    Options options_;                  ///< 系统运行选项。
    std::atomic_bool running_ = false;  ///< 系统是否处于建图运行状态。

    rclcpp::Service<SaveMapService>::SharedPtr savemap_service_ = nullptr;  ///< ROS2保存地图服务。

    std::string map_name_;  ///< 当前地图名。

    std::shared_ptr<LaserMapping> lio_ = nullptr;       ///< LIO前端。
    std::shared_ptr<LoopClosing> lc_ = nullptr;         ///< 回环检测模块。
    std::shared_ptr<ui::PangolinWindow> ui_ = nullptr;  ///< 3D可视化UI。
    std::shared_ptr<g2p5::G2P5> g2p5_ = nullptr;        ///< 2.5D/2D栅格地图模块。

    Keyframe::Ptr cur_kf_ = nullptr;  ///< 最近一次已分发给后端模块的关键帧。

    rclcpp::Node::SharedPtr node_;  ///< 在线模式下使用的ROS2节点。
    std::string imu_topic_;         ///< IMU订阅话题名。
    std::string cloud_topic_;       ///< 标准PointCloud2订阅话题名。
    std::string livox_topic_;       ///< Livox CustomMsg订阅话题名。

    rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr imu_sub_ = nullptr;  ///< IMU订阅器。
    rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr cloud_sub_ = nullptr;  ///< 标准点云订阅器。
    rclcpp::Subscription<livox_ros_driver2::msg::CustomMsg>::SharedPtr livox_sub_ = nullptr;  ///< Livox点云订阅器。
};
}  // namespace lightning

#endif  // LIGHTNING_SLAM_H
