//
// Created by xiang on 25-3-18.
//

/**
 * @file run_loc_online.cc
 * @brief 在线定位程序入口，选择定位模式并管理系统的启动、运行与退出。
 * @details 由运行脚本启动；LiDAR 模式使用 lightning::LocSystem，纯组合导航模式使用
 * lightning::InsLocSystem。传感器处理和业务发布由相应系统对象完成。
 * @par 阅读入口
 * 先查看本文件 main 的主要流程和同名阶段注释；ins_only 在模式分支内运行并返回，
 * 后续阶段属于 LiDAR 分支。算法细节沿下列系统函数继续阅读。
 * @see lightning::LocSystem::Init(), lightning::LocSystem::SetInitPose(),
 *      lightning::LocSystem::Spin(), lightning::LocSystem::Finish(),
 *      @ref online_localization_flow "在线定位数据流"
 */

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <yaml-cpp/yaml.h>

#include <cmath>
#include <algorithm>
#include <chrono>
#include <memory>
#include <stdexcept>
#include <thread>
#include <vector>

#include "app/online_bag_player.h"
#include "core/system/loc_system.h"
#include "core/system/ins_loc_system.h"
#include "ui/pangolin_window.h"
#include "wrapper/ros_utils.h"
#include "utils/online_logging.h"
#include "utils/thread_scheduling.h"
#include "utils/causal_trace.h"

DEFINE_string(config, "./config/default.yaml", "配置文件");
DEFINE_string(map, "", "地图目录；非空时覆盖配置中的 system.map_path");

/// 运行定位的测试
DEFINE_bool(use_config_initial_pose, true,
            "use online_localization.initial_pose, falling back to offline_localization.initial_pose");
DEFINE_string(bag, "", "optional ROS 2 bag for in-process online regression playback");
DEFINE_double(playback_rate, 1.0, "sensor-time playback rate used with --bag");
DEFINE_int32(post_wait_seconds, 15, "drain time after embedded bag playback");
DEFINE_string(output_tum, "", "optional PGO global-correction TUM trajectory");
DEFINE_string(output_high_frequency_tum, "", "optional live high-frequency TUM trajectory");
DEFINE_string(output_published_tum, "", "optional live published rear-axle TUM trajectory");

namespace {

bool ReadInitialPose(const std::string& config_path, lightning::SE3& pose) {
    const YAML::Node root = YAML::LoadFile(config_path);
    YAML::Node node;
    if (root["online_localization"] && root["online_localization"]["initial_pose"]) {
        node = root["online_localization"]["initial_pose"];
    } else if (root["offline_localization"] && root["offline_localization"]["initial_pose"]) {
        node = root["offline_localization"]["initial_pose"];
    } else {
        return false;
    }
    if (node["enabled"] && !node["enabled"].as<bool>()) return false;

    const auto translation = node["translation"].as<std::vector<double>>();
    const auto quaternion_xyzw = node["quaternion_xyzw"].as<std::vector<double>>();
    if (translation.size() != 3 || quaternion_xyzw.size() != 4) {
        throw std::runtime_error("initial pose must contain a 3D translation and XYZW quaternion");
    }
    for (const double value : translation) {
        if (!std::isfinite(value)) throw std::runtime_error("initial pose translation contains a non-finite value");
    }
    for (const double value : quaternion_xyzw) {
        if (!std::isfinite(value)) throw std::runtime_error("initial pose quaternion contains a non-finite value");
    }

    lightning::Quatd q(quaternion_xyzw[3], quaternion_xyzw[0], quaternion_xyzw[1], quaternion_xyzw[2]);
    if (q.norm() < 1e-9) throw std::runtime_error("initial pose quaternion has zero norm");
    q.normalize();
    pose = lightning::SE3(q, lightning::Vec3d(translation[0], translation[1], translation[2]));
    return true;
}

}  // namespace

/**
 * @brief 按配置启动在线定位系统，并在 ROS 关闭后完成运行收尾。
 * @details LiDAR 分支沿 LocSystem::Init()、LocSystem::SetInitPose()、LocSystem::Spin()
 * 执行；退出时调用 LocSystem::Finish() 并按 flags 保存轨迹。
 * @return 正常结束返回 0；模式、初始化、回放或输出设置失败返回非零值。
 * @see lightning::LocSystem, lightning::InsLocSystem
 *
 * @par 主要流程
 * 1. <b>准备进程与 ROS</b>：初始化日志、解析 flags、应用线程调度策略并启动 ROS。
 * 2. <b>选择定位模式</b>：ins_only 分支独立初始化、运行和退出；LiDAR 模式进入后续阶段。
 * 3. <b>初始化 LiDAR 系统</b>：创建 LocSystem，并按需开启实际发布轨迹记录。
 * 4. <b>提供初始位姿</b>：读取配置中的初始位姿，缺省使用单位位姿，再开启定位输入。
 * 5. <b>启动回放与事件循环</b>：按需启动进程内 bag 回放线程，随后进入 LocSystem::Spin。
 * 6. <b>等待线程并结束运行</b>：Spin 返回后等待回放线程、Finish 算法、按需导出轨迹并关闭 ROS。
 * @note bag 回放后的 post_wait_seconds 是固定等待时间，不能证明全部队列已排空；结尾轨迹导出的返回值未检查，返回 0 不能证明这些文件保存成功。
 */
int main(int argc, char** argv) {
    // 阶段 1：准备进程与 ROS。
    google::InitGoogleLogging(argv[0]);
    FLAGS_colorlogtostderr = true;
    FLAGS_stderrthreshold = google::INFO;

    google::ParseCommandLineFlags(&argc, &argv, true);
    using namespace lightning;

    std::string scheduling_error;
    if (!threading::InitializeOnlineThreadPolicy(&scheduling_error)) {
        LOG(ERROR) << scheduling_error;
        return 2;
    }
    OnlineLogging online_logging;
    // Start the optional trace writer while the main thread still has the
    // background mask; lazy initialization from a critical callback would
    // otherwise inherit the critical mask.
    (void)profiling::CausalTraceEnabled();
    rclcpp::init(argc, argv);

    // 阶段 2：选择定位模式。
    try {
        const auto config = YAML::LoadFile(FLAGS_config);
        const std::string mode = config["system"] && config["system"]["localization_mode"]
            ? config["system"]["localization_mode"].as<std::string>() : "lidar";
        if (mode == "ins_only") {
            if (!FLAGS_map.empty() || !FLAGS_bag.empty() || !FLAGS_output_tum.empty() || !FLAGS_output_high_frequency_tum.empty()) {
                LOG(ERROR) << "ins_only: use scripts/replay_ins_only.sh for clocked input-only bag playback; "
                           << "--output_published_tum is supported, map/PGO/LIO flags are not";
                rclcpp::shutdown();
                return 2;
            }
            InsLocSystem ins;
            if (!ins.Init(FLAGS_config, FLAGS_output_published_tum)) { rclcpp::shutdown(); return 1; }
            ins.Spin();
            rclcpp::shutdown();
            return 0;
        }
        if (mode != "lidar") throw std::invalid_argument("system.localization_mode must be lidar or ins_only");
    } catch (const std::exception& error) {
        LOG(ERROR) << error.what(); rclcpp::shutdown(); return 2;
    }

    // 阶段 3：初始化 LiDAR 系统。
    // [online-loc-create-system]
    LocSystem::Options opt;
    LocSystem loc(opt);

    if (!loc.Init(FLAGS_config, FLAGS_map)) {
        LOG(ERROR) << "failed to init loc";
        return -1;
    }
    // [online-loc-create-system]
    if (!FLAGS_output_published_tum.empty() &&
        !loc.StartPublishedTrajectoryTum(FLAGS_output_published_tum)) {
        LOG(ERROR) << "failed to start published rear-axle TUM recorder";
        return -1;
    }

    // 阶段 4：提供初始位姿。
    // [online-loc-initial-pose]
    /// 先提供初始位姿并打开输入处理，再进入 ROS 事件循环。
    SE3 initial_pose;
    try {
        if (FLAGS_use_config_initial_pose && ReadInitialPose(FLAGS_config, initial_pose)) {
            LOG(INFO) << "online localization initial pose: " << initial_pose.translation().transpose();
        } else {
            LOG(WARNING) << "online localization initial pose is not configured; using identity";
        }
    } catch (const std::exception& e) {
        LOG(ERROR) << "invalid online localization initial pose: " << e.what();
        return -1;
    }
    loc.SetInitPose(initial_pose);
    // [online-loc-initial-pose]

    // 阶段 5：启动回放与事件循环。
    std::shared_ptr<OnlineBagPlayer> bag_player;
    std::thread bag_thread;
    if (!FLAGS_bag.empty()) {
        bag_player = std::make_shared<OnlineBagPlayer>();
        if (!bag_player->Init(FLAGS_config, FLAGS_bag, FLAGS_playback_rate)) {
            LOG(ERROR) << "failed to initialize embedded online bag player";
            return -1;
        }
        bag_thread = std::thread([bag_player]() {
            bag_player->Run();
            std::this_thread::sleep_for(std::chrono::seconds(std::max(0, FLAGS_post_wait_seconds)));
            rclcpp::shutdown();
        });
    }
    // [online-loc-spin-finish]
    loc.Spin();
    // 阶段 6：等待线程并结束运行。
    if (bag_thread.joinable()) bag_thread.join();
    loc.Finish();
    // [online-loc-spin-finish]
    // 这里未检查轨迹导出的返回值，正常退出不能作为轨迹保存成功的判据。
    if (!FLAGS_output_tum.empty()) loc.SaveTrajectoryTum(FLAGS_output_tum);
    if (!FLAGS_output_high_frequency_tum.empty()) {
        loc.SaveHighFrequencyTrajectoryTum(FLAGS_output_high_frequency_tum);
    }

    rclcpp::shutdown();

    return 0;
}
