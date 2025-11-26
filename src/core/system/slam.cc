//
// Created by xiang on 25-5-6.
//

#include "core/system/slam.h"
#include "core/g2p5/g2p5.h"
#include "core/lio/laser_mapping.h"
#include "core/loop_closing/loop_closing.h"
#include "core/maps/tiled_map.h"
#include "ui/pangolin_window.h"
#include "wrapper/ros_utils.h"

#include <yaml-cpp/yaml.h>
#include <filesystem>
#include <opencv2/opencv.hpp>

namespace lightning {

SlamSystem::SlamSystem(lightning::SlamSystem::Options options) : options_(options) {
    /// handle ctrl-c
    signal(SIGINT, lightning::debug::SigHandle);
}

/**
 * @brief 初始化 SLAM 系统
 * @param yaml_path 配置文件路径
 * @return 初始化是否成功，成功返回 true，失败返回 false
 *
 * 该函数完成 SLAM 系统的完整初始化流程，包括以下主要步骤：
 * 1. **初始化 LIO 前端模块**：创建并初始化激光雷达惯性里程计模块
 * 2. **加载配置参数**：从 YAML 文件读取系统开关选项
 * 3. **初始化回环检测**：根据配置启用回环检测模块
 * 4. **初始化可视化**：根据配置启用 3D 可视化界面
 * 5. **初始化 2D 栅格地图**：根据配置启用 3D 到 2D 地图转换模块
 * 6. **创建 ROS2 节点**：如果是在线模式，创建订阅器和服务器
 *
 * 支持的配置选项包括：
 * - with_loop_closing: 是否启用回环检测
 * - with_ui: 是否启用 3D 可视化
 * - with_2dui: 是否启用 2D 可视化
 * - with_g2p5: 是否启用 3D 到 2D 地图转换
 * - step_on_kf: 是否在关键帧处暂停（调试模式）
 */
bool SlamSystem::Init(const std::string& yaml_path) {
    // 初始化 LIO（激光雷达惯性里程计）前端模块
    lio_ = std::make_shared<LaserMapping>();
    if (!lio_->Init(yaml_path)) {
        LOG(ERROR) << "failed to init lio module";
        return false;
    }

    // 从配置文件加载系统选项
    auto yaml = YAML::LoadFile(yaml_path);
    options_.with_loop_closing_ = yaml["system"]["with_loop_closing"].as<bool>();
    options_.with_visualization_ = yaml["system"]["with_ui"].as<bool>();
    options_.with_2dvisualization_ = yaml["system"]["with_2dui"].as<bool>();
    options_.with_gridmap_ = yaml["system"]["with_g2p5"].as<bool>();
    options_.step_on_kf_ = yaml["system"]["step_on_kf"].as<bool>();

    // 根据配置初始化回环检测模块
    if (options_.with_loop_closing_) {
        LOG(INFO) << "slam with loop closing";
        LoopClosing::Options options;
        options.online_mode_ = options_.online_mode_;
        lc_ = std::make_shared<LoopClosing>(options);
        lc_->Init(yaml_path);
    }

    // 根据配置初始化 3D 可视化模块
    if (options_.with_visualization_) {
        LOG(INFO) << "slam with 3D UI";
        ui_ = std::make_shared<ui::PangolinWindow>();
        ui_->Init();

        // 将 UI 界面设置到 LIO 模块
        lio_->SetUI(ui_);
    }

    // 根据配置初始化 3D 到 2D 栅格地图转换模块
    if (options_.with_gridmap_) {
        g2p5::G2P5::Options opt;
        opt.online_mode_ = options_.online_mode_;

        g2p5_ = std::make_shared<g2p5::G2P5>(opt);
        g2p5_->Init(yaml_path);

        // 如果启用回环检测，设置回环回调函数
        if (options_.with_loop_closing_) {
            /// 当发生回环时，触发一次重绘
            lc_->SetLoopClosedCB([this]() { g2p5_->RedrawGlobalMap(); });
        }

        // 如果启用 2D 可视化，设置地图更新回调函数
        if (options_.with_2dvisualization_) {
            g2p5_->SetMapUpdateCallback([this](g2p5::G2P5MapPtr map) {
                cv::Mat image = map->ToCV();
                cv::imshow("map", image);

                // 根据配置决定是否在关键帧处暂停
                if (options_.step_on_kf_) {
                    cv::waitKey(0);  // 等待按键继续（调试模式）
                } else {
                    cv::waitKey(10);  // 短暂延迟
                }
            });
        }
    }

    // 如果是在线模式，创建 ROS2 节点和相关订阅器
    if (options_.online_mode_) {
        LOG(INFO) << "online mode, creating ros2 node ... ";

        /// 创建 ROS2 节点
        node_ = std::make_shared<rclcpp::Node>("lightning_slam");

        // 从配置文件读取话题名称
        imu_topic_ = yaml["common"]["imu_topic"].as<std::string>();
        cloud_topic_ = yaml["common"]["lidar_topic"].as<std::string>();
        livox_topic_ = yaml["common"]["livox_lidar_topic"].as<std::string>();

        // 设置 QoS 策略
        rclcpp::QoS qos(10);
        // qos.best_effort();

        // 创建 IMU 数据订阅器
        imu_sub_ = node_->create_subscription<sensor_msgs::msg::Imu>(
            imu_topic_, qos, [this](sensor_msgs::msg::Imu::SharedPtr msg) {
                // 将 ROS2 IMU 消息转换为内部格式
                IMUPtr imu = std::make_shared<IMU>();
                imu->timestamp = ToSec(msg->header.stamp);
                imu->linear_acceleration =
                    Vec3d(msg->linear_acceleration.x, msg->linear_acceleration.y, msg->linear_acceleration.z);
                imu->angular_velocity =
                    Vec3d(msg->angular_velocity.x, msg->angular_velocity.y, msg->angular_velocity.z);

                // 处理 IMU 数据
                ProcessIMU(imu);
            });

        // 创建标准点云订阅器（适用于 Velodyne、Ouster 等雷达）
        cloud_sub_ = node_->create_subscription<sensor_msgs::msg::PointCloud2>(
            cloud_topic_, qos, [this](sensor_msgs::msg::PointCloud2::SharedPtr cloud) {
                // 使用计时器统计点云处理耗时
                Timer::Evaluate([&]() { ProcessLidar(cloud); }, "Proc Lidar", true);
            });

        // 创建 Livox 自定义点云消息订阅器
        livox_sub_ = node_->create_subscription<livox_ros_driver2::msg::CustomMsg>(
            livox_topic_, qos, [this](livox_ros_driver2::msg::CustomMsg ::SharedPtr cloud) {
                Timer::Evaluate([&]() { ProcessLidar(cloud); }, "Proc Lidar", true);
            });

        // 创建地图保存服务
        savemap_service_ = node_->create_service<SaveMapService>(
            "lightning/save_map", [this](const SaveMapService::Request::SharedPtr& req,
                                         SaveMapService::Response::SharedPtr res) { SaveMap(req, res); });

        LOG(INFO) << "online slam node has been created.";
    }

    return true;
}

SlamSystem::~SlamSystem() {
    if (ui_) {
        ui_->Quit();
    }
}

void SlamSystem::StartSLAM(std::string map_name) {
    map_name_ = map_name;
    running_ = true;
}

void SlamSystem::SaveMap(const SaveMapService::Request::SharedPtr request,
                         SaveMapService::Response::SharedPtr response) {
    map_name_ = request->map_id;
    std::string save_path = "./data/" + map_name_ + "/";

    SaveMap(save_path);
    response->response = 0;
}

void SlamSystem::SaveMap(const std::string& path) {
    std::string save_path = path;
    if (save_path.empty()) {
        save_path = "./data/" + map_name_ + "/";
    }

    LOG(INFO) << "slam map saving to " << save_path;

    if (!std::filesystem::exists(save_path)) {
        std::filesystem::create_directories(save_path);
    } else {
        std::filesystem::remove_all(save_path);
        std::filesystem::create_directories(save_path);
    }

    // auto global_map_no_loop = lio_->GetGlobalMap(true);
    auto global_map = lio_->GetGlobalMap(!options_.with_loop_closing_);
    // auto global_map_raw = lio_->GetGlobalMap(!options_.with_loop_closing_, false, 0.1);

    TiledMap::Options tm_options;
    tm_options.map_path_ = save_path;

    TiledMap tm(tm_options);
    SE3 start_pose = lio_->GetAllKeyframes().front()->GetOptPose();
    tm.ConvertFromFullPCD(global_map, start_pose, save_path);

    pcl::io::savePCDFileBinaryCompressed(save_path + "/global.pcd", *global_map);
    // pcl::io::savePCDFileBinaryCompressed(save_path + "/global_no_loop.pcd", *global_map_no_loop);
    // pcl::io::savePCDFileBinaryCompressed(save_path + "/global_raw.pcd", *global_map_raw);

    if (options_.with_gridmap_) {
        /// 存为ROS兼容的模式
        auto map = g2p5_->GetNewestMap()->ToROS();
        const int width = map.info.width;
        const int height = map.info.height;

        cv::Mat nav_image(height, width, CV_8UC1);
        for (int y = 0; y < height; ++y) {
            const int rowStartIndex = y * width;
            for (int x = 0; x < width; ++x) {
                const int index = rowStartIndex + x;
                int8_t data = map.data[index];
                if (data == 0) {                                   // Free
                    nav_image.at<uchar>(height - 1 - y, x) = 255;  // White
                } else if (data == 100) {                          // Occupied
                    nav_image.at<uchar>(height - 1 - y, x) = 0;    // Black
                } else {                                           // Unknown
                    nav_image.at<uchar>(height - 1 - y, x) = 128;  // Gray
                }
            }
        }

        cv::imwrite(save_path + "/map.pgm", nav_image);

        /// yaml
        std::ofstream yamlFile(save_path + "/map.yaml");
        if (!yamlFile.is_open()) {
            LOG(ERROR) << "failed to write map.yaml";
            return;  // 文件打开失败
        }

        try {
            YAML::Emitter emitter;
            emitter << YAML::BeginMap;
            emitter << YAML::Key << "image" << YAML::Value << "map.pgm";
            emitter << YAML::Key << "mode" << YAML::Value << "trinary";
            emitter << YAML::Key << "width" << YAML::Value << map.info.width;
            emitter << YAML::Key << "height" << YAML::Value << map.info.height;
            emitter << YAML::Key << "resolution" << YAML::Value << float(0.05);
            std::vector<double> orig{map.info.origin.position.x, map.info.origin.position.y, 0};
            emitter << YAML::Key << "origin" << YAML::Value << orig;
            emitter << YAML::Key << "negate" << YAML::Value << 0;
            emitter << YAML::Key << "occupied_thresh" << YAML::Value << 0.65;
            emitter << YAML::Key << "free_thresh" << YAML::Value << 0.25;

            emitter << YAML::EndMap;

            yamlFile << emitter.c_str();
            yamlFile.close();
        } catch (...) {
            yamlFile.close();
            return;
        }
    }

    LOG(INFO) << "map saved";
}

void SlamSystem::ProcessIMU(const lightning::IMUPtr& imu) {
    if (running_ == false) {
        return;
    }
    lio_->ProcessIMU(imu);
}

// 模板化的点云处理函数实现
template <typename PointCloudMsgType>
void SlamSystem::ProcessLidar(const std::shared_ptr<PointCloudMsgType>& cloud) {
    if (running_ == false) {
        return;
    }

    lio_->ProcessPointCloud2(cloud);
    lio_->Run();

    auto kf = lio_->GetKeyframe();
    if (kf != cur_kf_) {
        cur_kf_ = kf;
    } else {
        return;
    }

    if (cur_kf_ == nullptr) {
        return;
    }

    if (options_.with_loop_closing_) {
        lc_->AddKF(cur_kf_);
    }

    if (options_.with_gridmap_) {
        g2p5_->PushKeyframe(cur_kf_);
    }

    if (ui_) {
        ui_->UpdateKF(cur_kf_);
    }
}

// 显式实例化
template void SlamSystem::ProcessLidar(const sensor_msgs::msg::PointCloud2::SharedPtr& cloud);
template void SlamSystem::ProcessLidar(const livox_ros_driver2::msg::CustomMsg::SharedPtr& cloud);

void SlamSystem::Spin() {
    if (options_.online_mode_ && node_ != nullptr) {
        spin(node_);
    }
}

}  // namespace lightning
