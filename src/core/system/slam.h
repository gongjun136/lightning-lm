//
// Created by xiang on 25-5-6.
//

#ifndef LIGHTNING_SLAM_H
#define LIGHTNING_SLAM_H

#include <rclcpp/rclcpp.hpp>
#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <std_srvs/srv/trigger.hpp>
#include <string>
#include <tf2_ros/transform_broadcaster.h>
#include <vector>

#include "lightning/srv/save_map.hpp"
#include "livox_ros_driver2/msg/custom_msg.hpp"

#include "common/eigen_types.h"
#include "common/imu.h"
#include "common/keyframe.h"
#include "common/nav_state.h"
#include "core/maps/map_frame.h"
#include "core/maps/map_manifest.h"

namespace lightning {

class LaserMapping;  //  lio 前端
class LoopClosing;   // 回环检测
namespace backend {
class BackendPipeline;
}

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
     *
     * @par 主要流程
     * 1. <b>初始化 LIO 前端</b>：创建并初始化 LaserMapping，失败时返回。
     * 2. <b>读取系统和导出配置</b>：解析系统开关、地面坐标与地图身份配置。
     * 3. <b>选择并初始化后端</b>：按配置启用 BA/BTC/HBA、旧回环后端或关闭回环。
     * 4. <b>按需创建三维界面</b>：界面失败时关闭可视化并继续初始化。
     * 5. <b>配置栅格地图及回调</b>：按需创建 G2P5，连接旧回环重绘和二维显示回调。
     * 6. <b>配置在线 ROS 接口</b>：仅在线模式创建节点、后端通知、传感器订阅和保存/优化服务。
     * @note 后端初始化可能已经启动线程；此处仅注册传感器回调，ROS 事件循环在 Spin 中执行。旧回环和 G2P5 的 Init 返回值未在此检查。
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
     *
     * @par 主要流程
     * 1. <b>确定目录并等待关键帧</b>：使用默认或指定目录，等待新后端优化，拒绝无关键帧的导出。
     * 2. <b>生成地图并校验导出坐标</b>：选择 LIO 或优化位姿生成全局点云，估计地面坐标并变换点云。
     * 3. <b>替换目录并写入点云</b>：坐标校验后重建目标目录，写分块地图和 global.pcd。
     * 4. <b>导出后端资料和轨迹</b>：按需保存后端诊断、BTC 库及同一地图坐标下的轨迹。
     * 5. <b>按配置导出导航栅格</b>：优先使用导航导出器，否则可从 G2P5 写出 PGM/YAML。
     * 6. <b>写入坐标元数据和清单</b>：完成元数据与地图清单后返回 true。
     * @note 目录替换和多文件写入不是原子操作，后续失败可能留下部分产物；true 表示通过本函数现有检查，分块写入和旧栅格 imwrite 并未完整检查返回值。
     */
    bool SaveMap(const std::string& path = "");

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

    template <typename PointCloudMsgType>
    void ProcessLidar(const std::shared_ptr<PointCloudMsgType>& cloud, int lidar_id);

    /**
     * @brief 获取当前LIO前端状态。
     *
     * 离线评估程序用该接口导出TUM轨迹；若LIO尚未初始化，返回pose_is_ok_=false的无效状态。
     */
    NavState GetLioState() const;

    /**
     * @brief 导出关键帧轨迹；use_lio_pose=true导出前端LIO位姿，false导出后端优化位姿。
     */
    bool SaveKeyframeTrajectoryTum(
        const std::string& path, bool use_lio_pose,
        const map_frame::Metadata* map_metadata = nullptr) const;

    /**
     * @brief 在线模式下启动ROS2事件循环。
     */
    void Spin();

    /**
     * @brief 如果启用了3D UI，则阻塞等待用户关闭窗口。
     *
     * 离线bag回放结束后使用，避免可视化窗口在处理完成后立即析构关闭。
     */
    void WaitForUIQuit() const;

   private:
    /**
     * @brief ROS2保存地图服务的回调实现。
     * @param request 保存地图服务请求。
     * @param response 保存地图服务响应。
     */
    void SaveMap(const SaveMapService::Request::SharedPtr request, SaveMapService::Response::SharedPtr response);
    void OptimizeBackend(const std_srvs::srv::Trigger::Request::SharedPtr request,
                         std_srvs::srv::Trigger::Response::SharedPtr response);
    /**
     * @brief 排空 LIO 已同步数据并分发新关键帧。
     *
     * @par 主要流程
     * 1. <b>消费已同步帧</b>：每轮调用 RunDetailed；无数据即返回，非输出帧继续下一轮。
     * 2. <b>记录有效前端状态</b>：仅保存位姿有效且时间戳递增的状态。
     * 3. <b>筛选新关键帧</b>：跳过与上次相同或为空的关键帧。
     * 4. <b>分发关键帧</b>：按配置交给新后端或旧回环，并更新栅格地图和界面。
     * @note 消费过一帧不保证产生有效位姿或新关键帧；本函数中的阶段随每轮循环重复。
     */
    void DrainLio();
    bool SaveLioTrajectoryTum(
        const std::string& path,
        const map_frame::Metadata* map_metadata = nullptr) const;
    void PublishMapToOdom();

    Options options_;                  ///< 系统运行选项。
    std::atomic_bool running_ = false;  ///< 系统是否处于建图运行状态。

    map_frame::ExportOptions map_export_options_;  ///< 地图坐标归一化与导航栅格导出配置。
    map_manifest::ExportOptions map_manifest_options_;  ///< 地图身份与业务坐标声明。
    rclcpp::Service<SaveMapService>::SharedPtr savemap_service_ = nullptr;  ///< ROS2保存地图服务。

    std::string map_name_;  ///< 当前地图名。

    std::shared_ptr<LaserMapping> lio_ = nullptr;       ///< LIO前端。
    std::shared_ptr<LoopClosing> lc_ = nullptr;         ///< 回环检测模块。
    std::shared_ptr<backend::BackendPipeline> backend_ = nullptr;  ///< 可选的新式 BA、BTC 与 HBA 后端。
    std::shared_ptr<ui::PangolinWindow> ui_ = nullptr;  ///< 3D可视化UI。
    std::shared_ptr<g2p5::G2P5> g2p5_ = nullptr;        ///< 2.5D/2D栅格地图模块。

    Keyframe::Ptr cur_kf_ = nullptr;  ///< 最近一次已分发给后端模块的关键帧。
    std::vector<NavState> lio_states_;  ///< 在线前端产生的有效全帧状态，用于回归评估。

    rclcpp::Node::SharedPtr node_;  ///< 在线模式下使用的ROS2节点。
    std::string imu_topic_;         ///< IMU订阅话题名。
    std::string cloud_topic_;       ///< 标准PointCloud2订阅话题名。
    std::string livox_topic_;       ///< Livox CustomMsg订阅话题名。

    rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr imu_sub_ = nullptr;  ///< IMU订阅器。
    std::vector<rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr> cloud_subs_;  ///< 多雷达模式下的点云订阅器集合。
    rclcpp::Subscription<livox_ros_driver2::msg::CustomMsg>::SharedPtr livox_sub_ = nullptr;  ///< Livox点云订阅器。
    rclcpp::Service<std_srvs::srv::Trigger>::SharedPtr optimize_backend_service_ = nullptr;  ///< 触发后端全局优化的 ROS 服务。
    std::shared_ptr<tf2_ros::TransformBroadcaster> tf_broadcaster_ = nullptr;  ///< 发布建图位姿坐标关系的 TF 广播器。
};
}  // namespace lightning

#endif  // LIGHTNING_SLAM_H
