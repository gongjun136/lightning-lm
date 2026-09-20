//
// Created by xiang on 25-3-18.
//

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
#include "ui/pangolin_window.h"
#include "wrapper/ros_utils.h"

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

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_colorlogtostderr = true;
    FLAGS_stderrthreshold = google::INFO;

    google::ParseCommandLineFlags(&argc, &argv, true);
    using namespace lightning;

    rclcpp::init(argc, argv);

    LocSystem::Options opt;
    LocSystem loc(opt);

    if (!loc.Init(FLAGS_config, FLAGS_map)) {
        LOG(ERROR) << "failed to init loc";
        return -1;
    }
    if (!FLAGS_output_published_tum.empty() &&
        !loc.StartPublishedTrajectoryTum(FLAGS_output_published_tum)) {
        LOG(ERROR) << "failed to start published rear-axle TUM recorder";
        return -1;
    }

    /// 默认起点开始定位
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
    loc.Spin();
    if (bag_thread.joinable()) bag_thread.join();
    loc.Finish();
    if (!FLAGS_output_tum.empty()) loc.SaveTrajectoryTum(FLAGS_output_tum);
    if (!FLAGS_output_high_frequency_tum.empty()) {
        loc.SaveHighFrequencyTrajectoryTum(FLAGS_output_high_frequency_tum);
    }

    rclcpp::shutdown();

    return 0;
}
