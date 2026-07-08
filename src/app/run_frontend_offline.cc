//
// Created by xiang on 25-3-18.
//

#include <gflags/gflags.h>
#include <glog/logging.h>

#include <chrono>
#include <fstream>
#include <iomanip>
#include <memory>
#include <thread>

#include "common/options.h"
#include "core/lio/laser_mapping.h"
#include "io/yaml_io.h"
#include "ui/pangolin_window.h"
#include "utils/timer.h"
#include "wrapper/bag_io.h"

DEFINE_string(input_bag, "", "input ROS2 bag");
DEFINE_string(config, "./config/default.yaml", "config yaml");
DEFINE_string(output_tum, "", "output front-end LIO TUM trajectory; disabled when empty");
DEFINE_bool(wait_ui, true, "wait for the 3D UI window to close after offline processing");
DEFINE_int32(max_lidar_frames, 0, "stop after processing this many lidar frames; disabled when <= 0");

namespace {
void WriteTumState(std::ofstream& tum, const lightning::NavState& state, double& last_timestamp) {
    if (!tum.is_open() || !state.pose_is_ok_ || state.timestamp_ <= 0.0 || state.timestamp_ <= last_timestamp) {
        return;
    }

    const auto q = state.rot_.unit_quaternion();
    tum << std::fixed << std::setprecision(9) << state.timestamp_ << " " << std::setprecision(12) << state.pos_.x()
        << " " << state.pos_.y() << " " << state.pos_.z() << " " << q.x() << " " << q.y() << " " << q.z() << " "
        << q.w() << "\n";
    last_timestamp = state.timestamp_;
}
}  // namespace

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_colorlogtostderr = true;
    FLAGS_stderrthreshold = google::INFO;

    google::ParseCommandLineFlags(&argc, &argv, true);
    if (FLAGS_input_bag.empty()) {
        LOG(ERROR) << "input_bag is required";
        return -1;
    }

    using namespace lightning;

    LaserMapping lio;
    if (!lio.Init(FLAGS_config)) {
        LOG(ERROR) << "failed to init lio";
        return -1;
    }

    YAML_IO yaml(FLAGS_config);
    const bool with_ui = yaml.GetValue<bool>("system", "with_ui");
    const std::string lidar_topic = yaml.GetValue<std::string>("common", "lidar_topic");
    const std::string livox_lidar_topic = yaml.GetValue<std::string>("common", "livox_lidar_topic");
    const std::string imu_topic = yaml.GetValue<std::string>("common", "imu_topic");

    std::shared_ptr<ui::PangolinWindow> ui;
    if (with_ui) {
        LOG(INFO) << "frontend with 3D UI";
        ui = std::make_shared<ui::PangolinWindow>();
        if (ui->Init()) {
            lio.SetUI(ui);
        } else {
            LOG(ERROR) << "failed to init 3D UI, continue without Pangolin";
            ui.reset();
        }
    }

    std::ofstream tum;
    double last_tum_timestamp = 0.0;
    int processed_lidar_frames = 0;
    auto finish_lidar_frame = [&processed_lidar_frames]() {
        if (FLAGS_max_lidar_frames <= 0) {
            return;
        }
        processed_lidar_frames++;
        if (processed_lidar_frames >= FLAGS_max_lidar_frames) {
            LOG(INFO) << "reached max_lidar_frames=" << FLAGS_max_lidar_frames << ", stopping offline frontend";
            lightning::debug::flg_exit = true;
        }
    };
    if (!FLAGS_output_tum.empty()) {
        tum.open(FLAGS_output_tum);
        if (!tum.is_open()) {
            LOG(ERROR) << "failed to open output_tum: " << FLAGS_output_tum;
            return -1;
        }
        LOG(INFO) << "writing front-end TUM trajectory to " << FLAGS_output_tum;
    }

    RosbagIO rosbag(FLAGS_input_bag);
    rosbag.AddImuHandle(imu_topic,
                        [&lio](IMUPtr imu) {
                            lio.ProcessIMU(imu);
                            return true;
                        });

    if (!lidar_topic.empty() && lidar_topic != livox_lidar_topic) {
        rosbag.AddPointCloud2Handle(lidar_topic,
                                    [&lio, &tum, &last_tum_timestamp, &finish_lidar_frame](
                                        sensor_msgs::msg::PointCloud2::SharedPtr cloud) {
                                        lio.ProcessPointCloud2(cloud);
                                        lio.Run();
                                        WriteTumState(tum, lio.GetState(), last_tum_timestamp);
                                        finish_lidar_frame();
                                        return true;
                                    });
    }

    if (!livox_lidar_topic.empty()) {
        rosbag.AddLivoxCloudHandle(
            livox_lidar_topic,
            [&lio, &tum, &last_tum_timestamp, &finish_lidar_frame](livox_ros_driver2::msg::CustomMsg::SharedPtr cloud) {
                lio.ProcessPointCloud2(cloud);
                lio.Run();
                WriteTumState(tum, lio.GetState(), last_tum_timestamp);
                finish_lidar_frame();
                return true;
            });
    }

    rosbag.Go();

    if (tum.is_open()) {
        tum.close();
    }
    Timer::PrintAll();

    if (ui && FLAGS_wait_ui) {
        LOG(INFO) << "waiting for 3D UI window to close";
        while (!ui->ShouldQuit()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }
    if (ui) {
        ui->Quit();
    }

    LOG(INFO) << "done";
    return 0;
}
