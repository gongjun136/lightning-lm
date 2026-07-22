//
// Created by xiang on 25-3-18.
//

#include <gflags/gflags.h>
#include <glog/logging.h>

#include <algorithm>
#include <chrono>
#include <memory>
#include <thread>

#include "app/online_bag_player.h"
#include "core/system/slam.h"
#include "utils/timer.h"
#include "wrapper/bag_io.h"
#include "wrapper/ros_utils.h"

DEFINE_string(config, "./config/default.yaml", "配置文件");

/// 运行一个LIO前端，带可视化
DEFINE_string(bag, "", "optional ROS 2 bag for in-process online regression playback");
DEFINE_double(playback_rate, 1.0, "sensor-time playback rate used with --bag");
DEFINE_int32(post_wait_seconds, 15, "drain time after embedded bag playback");
DEFINE_string(save_map, "", "optional map directory saved after online spin exits");

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_colorlogtostderr = true;
    FLAGS_stderrthreshold = google::INFO;
    google::ParseCommandLineFlags(&argc, &argv, true);

    using namespace lightning;

    /// 需要rclcpp::init
    rclcpp::init(argc, argv);

    SlamSystem::Options options;
    options.online_mode_ = true;

    SlamSystem slam(options);
    if (!slam.Init(FLAGS_config)) {
        LOG(ERROR) << "failed to init slam";
        return -1;
    }

    slam.StartSLAM("new_map");

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
    slam.Spin();
    if (bag_thread.joinable()) bag_thread.join();

    if (!FLAGS_save_map.empty()) slam.SaveMap(FLAGS_save_map);

    Timer::PrintAll();

    rclcpp::shutdown();

    LOG(INFO) << "done";

    return 0;
}
