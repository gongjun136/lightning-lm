//
// Created by xiang on 25-5-6.
//

#include "core/system/slam.h"
#include "core/backend/backend_pipeline.h"
#include "core/g2p5/g2p5.h"
#include "core/lio/laser_mapping.h"
#include "core/loop_closing/loop_closing.h"
#include "core/maps/navigation_map_export.h"
#include "core/maps/tiled_map.h"
#include "ui/pangolin_window.h"
#include "wrapper/ros_utils.h"

#include <yaml-cpp/yaml.h>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <opencv2/opencv.hpp>
#include <thread>

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
    std::string map_frame_error;
    if (!map_frame::ReadExportOptions(yaml, map_export_options_, map_frame_error)) {
        LOG(ERROR) << map_frame_error;
        return false;
    }

    // 根据配置初始化回环检测模块
    if (options_.with_loop_closing_) {
        const backend::BackendMode mode = backend::ReadBackendMode(yaml_path);
        LOG(INFO) << "slam backend mode: " << backend::BackendModeName(mode);
        if (mode == backend::BackendMode::kBaBtcHba) {
            backend_ = std::make_shared<backend::BackendPipeline>();
            if (!backend_->Init(yaml_path, options_.online_mode_)) return false;
        } else if (mode == backend::BackendMode::kLegacy) {
            LoopClosing::Options options;
            options.online_mode_ = options_.online_mode_;
            lc_ = std::make_shared<LoopClosing>(options);
            lc_->Init(yaml_path);
        } else {
            options_.with_loop_closing_ = false;
        }
    }

    // 根据配置初始化 3D 可视化模块
    if (options_.with_visualization_) {
        LOG(INFO) << "slam with 3D UI";
        ui_ = std::make_shared<ui::PangolinWindow>();
        if (ui_->Init()) {
            // 将 UI 界面设置到 LIO 模块
            lio_->SetUI(ui_);
        } else {
            LOG(ERROR) << "failed to init 3D UI, continue without Pangolin";
            ui_.reset();
            options_.with_visualization_ = false;
        }
    }

    // 根据配置初始化 3D 到 2D 栅格地图转换模块
    if (options_.with_gridmap_) {
        g2p5::G2P5::Options opt;
        opt.online_mode_ = options_.online_mode_;

        g2p5_ = std::make_shared<g2p5::G2P5>(opt);
        g2p5_->Init(yaml_path);

        // 如果启用回环检测，设置回环回调函数
        if (lc_) {
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
        if (backend_) {
            auto tf_qos = tf2_ros::DynamicBroadcasterQoS();
            tf_qos.best_effort();
            tf_broadcaster_ = std::make_shared<tf2_ros::TransformBroadcaster>(node_, tf_qos);
            backend_->SetOptimizedCallback([this]() {
                PublishMapToOdom();
                if (g2p5_) g2p5_->RedrawGlobalMap();
            });
        }

        // 从配置文件读取话题名称
        imu_topic_ = yaml["common"]["imu_topic"].as<std::string>();
        cloud_topic_ = yaml["common"]["lidar_topic"].as<std::string>();
        livox_topic_ = yaml["common"]["livox_lidar_topic"].as<std::string>();

        const auto qos = rclcpp::SensorDataQoS();
        auto subscribe_imu = [this, &qos](const std::string& topic) {
            imu_sub_ = node_->create_subscription<sensor_msgs::msg::Imu>(
                topic, qos, [this](sensor_msgs::msg::Imu::SharedPtr msg) {
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
        };

        if (lio_->IsMultiLidarEnabled()) {
            const auto* primary = lio_->GetMultiLidarConfig().PrimaryLidar();
            if (!primary) {
                LOG(ERROR) << "multi-lidar primary sensor is missing";
                return false;
            }
            subscribe_imu(primary->imu_topic);
            for (const auto& sensor : lio_->GetMultiLidarConfig().lidars) {
                cloud_subs_.push_back(node_->create_subscription<sensor_msgs::msg::PointCloud2>(
                    sensor.lidar_topic, qos,
                    [this, id = sensor.id](sensor_msgs::msg::PointCloud2::SharedPtr cloud) {
                        Timer::Evaluate([&]() { ProcessLidar(cloud, id); }, "Proc Lidar", true);
                    }));
            }
            LOG(INFO) << "online SLAM subscribed to " << cloud_subs_.size()
                      << " lidar topics and primary IMU " << primary->imu_topic;
        } else {
            if (imu_topic_.empty() || (cloud_topic_.empty() && livox_topic_.empty())) {
                LOG(ERROR) << "single-lidar online topics are incomplete";
                return false;
            }
            subscribe_imu(imu_topic_);
            if (!cloud_topic_.empty()) {
                cloud_subs_.push_back(node_->create_subscription<sensor_msgs::msg::PointCloud2>(
                    cloud_topic_, qos, [this](sensor_msgs::msg::PointCloud2::SharedPtr cloud) {
                        Timer::Evaluate([&]() { ProcessLidar(cloud); }, "Proc Lidar", true);
                    }));
            }
            if (!livox_topic_.empty()) {
                livox_sub_ = node_->create_subscription<livox_ros_driver2::msg::CustomMsg>(
                    livox_topic_, qos, [this](livox_ros_driver2::msg::CustomMsg::SharedPtr cloud) {
                        Timer::Evaluate([&]() { ProcessLidar(cloud); }, "Proc Lidar", true);
                    });
            }
        }

        // 创建地图保存服务
        savemap_service_ = node_->create_service<SaveMapService>(
            "lightning/save_map", [this](const SaveMapService::Request::SharedPtr& req,
                                         SaveMapService::Response::SharedPtr res) { SaveMap(req, res); });
        optimize_backend_service_ = node_->create_service<std_srvs::srv::Trigger>(
            "lightning/optimize_backend",
            [this](const std_srvs::srv::Trigger::Request::SharedPtr req,
                   std_srvs::srv::Trigger::Response::SharedPtr res) { OptimizeBackend(req, res); });

        LOG(INFO) << "online slam node has been created.";
    }

    return true;
}

SlamSystem::~SlamSystem() {
    if (backend_) {
        backend_->WaitUntilIdle(true);
        backend_->Shutdown();
    }
    if (ui_) {
        ui_->Quit();
    }
}

void SlamSystem::StartSLAM(std::string map_name) {
    map_name_ = map_name;
    running_ = true;
}

/**
 * @brief ROS2保存地图服务回调。
 * @param request 保存地图请求，map_id会作为当前地图名和默认目录名。
 * @param response 保存地图响应，response为0表示已完成保存流程。
 */
void SlamSystem::SaveMap(const SaveMapService::Request::SharedPtr request,
                         SaveMapService::Response::SharedPtr response) {
    /// 服务请求中的地图ID覆盖当前地图名，确保后续默认路径与请求一致。
    map_name_ = request->map_id;
    std::string save_path = "./data/" + map_name_ + "/";

    response->response = SaveMap(save_path) ? 0 : 3;
}

/**
 * @brief 将当前SLAM地图保存到指定目录。
 * @param path 地图保存目录；为空时使用./data/{map_name_}/作为默认目录。
 *
 * 保存内容包括全局点云global.pcd、分块地图数据，以及可选的ROS导航兼容栅格地图
 * map.pgm和map.yaml。若目标目录已存在，会先清空再重新创建。
 */
bool SlamSystem::SaveMap(const std::string& path) {
    std::string save_path = path;
    if (save_path.empty()) {
        save_path = "./data/" + map_name_ + "/";
    }

    LOG(INFO) << "slam map saving to " << save_path;
    if (backend_) backend_->WaitUntilIdle(true);
    const auto keyframes = lio_->GetAllKeyframes();
    if (keyframes.empty()) {
        LOG(ERROR) << "cannot save a map without keyframes";
        return false;
    }

    // auto global_map_no_loop = lio_->GetGlobalMap(true);
    /// 根据回环配置导出优化后的全局点云，关闭回环时直接使用无回环轨迹。
    const bool use_lio_pose = !options_.with_loop_closing_;
    auto global_map = lio_->GetGlobalMap(use_lio_pose, true, 0.1F, true);
    map_frame::Metadata map_metadata;
    std::string map_frame_error;
    if (!map_frame::EstimateStartGroundFrame(
            global_map, keyframes.front()->GetOptPose(), map_export_options_,
            map_metadata, map_frame_error)) {
        LOG(ERROR) << "map export aborted: " << map_frame_error;
        return false;
    }
    map_frame::TransformCloudInPlace(map_metadata, global_map);
    const map_frame::Metadata* metadata =
        map_metadata.normalized ? &map_metadata : nullptr;
    // auto global_map_raw = lio_->GetGlobalMap(!options_.with_loop_closing_, false, 0.1);

    /// 地面估计成功后再重建目录，避免失败的导出破坏已有地图。
    if (!std::filesystem::exists(save_path)) {
        std::filesystem::create_directories(save_path);
    } else {
        std::filesystem::remove_all(save_path);
        std::filesystem::create_directories(save_path);
    }

    /// 将完整点云转换为项目内部的分块地图格式，起始关键帧位姿用于建立局部地图基准。
    TiledMap::Options tm_options;
    tm_options.map_path_ = save_path;

    TiledMap tm(tm_options);
    const SE3 start_pose =
        map_frame::TransformPose(map_metadata, keyframes.front()->GetOptPose());
    if (!tm.ConvertFromFullPCD(global_map, start_pose, save_path)) {
        LOG(ERROR) << "failed to export tiled map";
        return false;
    }

    if (pcl::io::savePCDFileBinaryCompressed(
            save_path + "/global.pcd", *global_map) != 0) {
        LOG(ERROR) << "failed to save global map";
        return false;
    }
    if (backend_) {
        if (!backend_->SaveDiagnostics(
                save_path + "/backend_diagnostics", metadata) ||
            !backend_->SaveRelocalizationDatabase(
                save_path + "/btc_relocalization", metadata)) {
            LOG(ERROR) << "failed to save map-consistent backend artifacts";
            return false;
        }
    }
    if (!SaveKeyframeTrajectoryTum(
            save_path + "/trajectory_slam_keyframes_lio.tum", true, metadata) ||
        !SaveKeyframeTrajectoryTum(
            save_path + "/trajectory_slam_keyframes_opt.tum", false, metadata) ||
        !SaveLioTrajectoryTum(save_path + "/trajectory_slam.tum", metadata)) {
        LOG(ERROR) << "failed to save map-consistent trajectories";
        return false;
    }
    // pcl::io::savePCDFileBinaryCompressed(save_path + "/global_no_loop.pcd", *global_map_no_loop);
    // pcl::io::savePCDFileBinaryCompressed(save_path + "/global_raw.pcd", *global_map_raw);

    if (map_export_options_.export_pgm) {
        const SE3 T_imu_primary(
            Quatd(lio_->GetLidarToImuRotation()).normalized(),
            lio_->GetLidarToImuTranslation());
        std::vector<navigation_map::RaycastFrame> raycast_frames;
        raycast_frames.reserve(keyframes.size());
        for (const auto& keyframe : keyframes) {
            if (!keyframe || !keyframe->GetCloud()) continue;
            const SE3 pose = use_lio_pose ? keyframe->GetLIOPose()
                                         : keyframe->GetOptPose();
            raycast_frames.push_back({
                lio_->PrepareMapExportCloud(keyframe->GetCloud()),
                map_frame::TransformPose(map_metadata, pose) * T_imu_primary});
        }
        navigation_map::SensorOrigins sensor_origins{{0, Vec3d::Zero()}};
        if (lio_->IsMultiLidarEnabled()) {
            for (const auto& sensor : lio_->GetMultiLidarConfig().lidars) {
                sensor_origins[static_cast<std::uint8_t>(sensor.id)] =
                    sensor.t_lidar_to_primary;
            }
        }
        navigation_map::ExportResult pgm_result;
        if (!navigation_map::ExportPgmAndYaml(
                global_map, raycast_frames, sensor_origins, save_path,
                map_export_options_, pgm_result, map_frame_error)) {
            LOG(ERROR) << "failed to export navigation map: " << map_frame_error;
            return false;
        }
        LOG(INFO) << "exported navigation map " << pgm_result.width << "x"
                  << pgm_result.height << ", origin=[" << pgm_result.origin_x
                  << ", " << pgm_result.origin_y << "] occupied_cells="
                  << pgm_result.occupied_cells << " free_cells="
                  << pgm_result.free_cells << " unknown_cells="
                  << pgm_result.unknown_cells << " rays="
                  << pgm_result.ray_count;
    } else if (options_.with_gridmap_) {
        /// 存为ROS导航兼容的栅格地图格式。
        auto map = g2p5_->GetNewestMap()->ToROS();
        const int width = map.info.width;
        const int height = map.info.height;

        /// ROS OccupancyGrid原点在左下，PGM图像原点在左上，因此写图时需要翻转y轴。
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

        /// 写入ROS导航地图元数据，与map.pgm组成标准可加载地图。
        std::ofstream yamlFile(save_path + "/map.yaml");
        if (!yamlFile.is_open()) {
            LOG(ERROR) << "failed to write map.yaml";
            return false;  // 文件打开失败
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
            return false;
        }
    }

    if (!map_frame::SaveMetadata(save_path, map_metadata, map_frame_error)) {
        LOG(ERROR) << "map export aborted: " << map_frame_error;
        return false;
    }
    LOG(INFO) << "map saved";
    return true;
}

void SlamSystem::ProcessIMU(const lightning::IMUPtr& imu) {
    if (running_ == false) {
        return;
    }
    lio_->ProcessIMU(imu);
    DrainLio();
}

NavState SlamSystem::GetLioState() const {
    if (!lio_) {
        NavState state;
        state.pose_is_ok_ = false;
        return state;
    }
    return lio_->GetState();
}

bool SlamSystem::SaveKeyframeTrajectoryTum(
    const std::string& path, bool use_lio_pose,
    const map_frame::Metadata* map_metadata) const {
    if (!lio_) {
        LOG(ERROR) << "lio is not initialized, skip trajectory export";
        return false;
    }

    std::ofstream tum(path);
    if (!tum.is_open()) {
        LOG(ERROR) << "failed to open keyframe trajectory: " << path;
        return false;
    }

    double last_timestamp = 0.0;
    int count = 0;
    for (const auto& kf : lio_->GetAllKeyframes()) {
        if (!kf) {
            continue;
        }

        const auto state = kf->GetState();
        if (state.timestamp_ <= 0.0 || state.timestamp_ <= last_timestamp) {
            continue;
        }

        const auto source_pose =
            use_lio_pose ? kf->GetLIOPose() : kf->GetOptPose();
        const auto pose = map_metadata
                              ? map_frame::TransformPose(*map_metadata, source_pose)
                              : source_pose;
        const auto q = pose.unit_quaternion();
        const auto p = pose.translation();
        tum << std::fixed << std::setprecision(9) << state.timestamp_ << " " << std::setprecision(12) << p.x() << " "
            << p.y() << " " << p.z() << " " << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << "\n";
        last_timestamp = state.timestamp_;
        count++;
    }

    LOG(INFO) << "wrote " << count << " keyframe poses to " << path;
    return count > 0;
}

bool SlamSystem::SaveLioTrajectoryTum(
    const std::string& path,
    const map_frame::Metadata* map_metadata) const {
    std::ofstream tum(path);
    if (!tum.is_open()) {
        LOG(ERROR) << "failed to open LIO trajectory: " << path;
        return false;
    }

    double last_timestamp = 0.0;
    int count = 0;
    for (const auto& state : lio_states_) {
        if (!state.pose_is_ok_ || state.timestamp_ <= 0.0 || state.timestamp_ <= last_timestamp) continue;
        const auto source_pose = state.GetPose();
        const auto pose = map_metadata
                              ? map_frame::TransformPose(*map_metadata, source_pose)
                              : source_pose;
        const auto q = pose.unit_quaternion();
        const auto p = pose.translation();
        tum << std::fixed << std::setprecision(9) << state.timestamp_ << " " << std::setprecision(12) << p.x() << " "
            << p.y() << " " << p.z() << " " << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << "\n";
        last_timestamp = state.timestamp_;
        ++count;
    }
    LOG(INFO) << "wrote " << count << " LIO poses to " << path;
    return count > 0;
}

// 模板化的点云处理函数实现
template <typename PointCloudMsgType>
void SlamSystem::ProcessLidar(const std::shared_ptr<PointCloudMsgType>& cloud) {
    const int lidar_id = lio_->IsMultiLidarEnabled() ? lio_->GetMultiLidarConfig().primary_lidar_id : 0;
    ProcessLidar(cloud, lidar_id);
}

template <typename PointCloudMsgType>
void SlamSystem::ProcessLidar(const std::shared_ptr<PointCloudMsgType>& cloud, int lidar_id) {
    if (running_ == false) {
        return;
    }

    // 先把不同类型的原始点云统一预处理并放入LIO缓存，再触发一次前端处理。
    // Run()内部会完成时间同步、IMU去畸变、雷达观测更新，并在满足条件时创建新关键帧。
    lio_->ProcessPointCloud2(cloud, lidar_id);
    DrainLio();
}

void SlamSystem::DrainLio() {
    while (lio_) {
        const auto status = lio_->RunDetailed();
        if (status == LaserMapping::RunStatus::kNoData) return;
        if (status != LaserMapping::RunStatus::kOutput) continue;

        const auto state = lio_->GetState();
        if (state.pose_is_ok_ && (lio_states_.empty() || state.timestamp_ > lio_states_.back().timestamp_)) {
            lio_states_.push_back(state);
        }

        // 后端只关心新产生的关键帧；没有新关键帧时继续排空其余已同步帧。
        auto kf = lio_->GetKeyframe();
        if (kf == cur_kf_) continue;
        cur_kf_ = kf;
        if (cur_kf_ == nullptr) continue;

        // 新关键帧按配置分发给回环、栅格建图和UI模块。
        if (backend_) {
            backend_->AddKeyframe(cur_kf_);
        } else if (lc_) {
            lc_->AddKF(cur_kf_);
        }

        if (options_.with_gridmap_) {
            g2p5_->PushKeyframe(cur_kf_);
        }

        if (ui_) {
            ui_->UpdateKF(cur_kf_);
        }
    }
}

// 显式实例化
template void SlamSystem::ProcessLidar(const sensor_msgs::msg::PointCloud2::SharedPtr& cloud);
template void SlamSystem::ProcessLidar(const livox_ros_driver2::msg::CustomMsg::SharedPtr& cloud);
template void SlamSystem::ProcessLidar(const sensor_msgs::msg::PointCloud2::SharedPtr& cloud, int lidar_id);
template void SlamSystem::ProcessLidar(const livox_ros_driver2::msg::CustomMsg::SharedPtr& cloud, int lidar_id);

void SlamSystem::OptimizeBackend(const std_srvs::srv::Trigger::Request::SharedPtr,
                                 std_srvs::srv::Trigger::Response::SharedPtr response) {
    if (!backend_) {
        response->success = false;
        response->message = "ba_btc_hba backend is not active";
        return;
    }
    backend_->RequestGlobalOptimization("ros_service");
    response->success = true;
    response->message = "global backend optimization requested";
}

void SlamSystem::PublishMapToOdom() {
    if (!backend_ || !node_ || !tf_broadcaster_) return;
    const SE3 transform = backend_->GetMapToOdom();
    geometry_msgs::msg::TransformStamped message;
    message.header.stamp = node_->now();
    message.header.frame_id = "map";
    message.child_frame_id = "odom";
    message.transform.translation.x = transform.translation().x();
    message.transform.translation.y = transform.translation().y();
    message.transform.translation.z = transform.translation().z();
    const auto rotation = transform.unit_quaternion();
    message.transform.rotation.x = rotation.x();
    message.transform.rotation.y = rotation.y();
    message.transform.rotation.z = rotation.z();
    message.transform.rotation.w = rotation.w();
    tf_broadcaster_->sendTransform(message);
}

void SlamSystem::Spin() {
    if (options_.online_mode_ && node_ != nullptr) {
        spin(node_);
    }
}

void SlamSystem::WaitForUIQuit() const {
    if (!ui_) {
        return;
    }

    LOG(INFO) << "waiting for 3D UI window to close";
    while (!ui_->ShouldQuit()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
}

}  // namespace lightning
