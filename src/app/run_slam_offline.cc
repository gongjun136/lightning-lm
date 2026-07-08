//
// Created by xiang on 25-3-18.
//

#include <gflags/gflags.h>
#include <glog/logging.h>

#include <fstream>
#include <iomanip>
#include <filesystem>

#include "core/system/slam.h"
#include "ui/pangolin_window.h"
#include "utils/timer.h"
#include "wrapper/bag_io.h"
#include "wrapper/ros_utils.h"

#include "io/yaml_io.h"

DEFINE_string(input_bag, "", "输入数据包");
DEFINE_string(config, "./config/default.yaml", "配置文件");
DEFINE_string(output_tum, "", "输出TUM轨迹文件，为空则不导出");
DEFINE_bool(wait_ui, true, "离线处理完成后是否等待3D UI窗口关闭");

namespace {
void WriteTumState(std::ofstream& tum, const lightning::NavState& state, double& last_timestamp) {
    if (!tum.is_open() || !state.pose_is_ok_ || state.timestamp_ <= 0.0 || state.timestamp_ <= last_timestamp) {
        return;
    }

    const auto q = state.rot_.unit_quaternion();
    tum << std::fixed << std::setprecision(9) << state.timestamp_ << " " << std::setprecision(12)
        << state.pos_.x() << " " << state.pos_.y() << " " << state.pos_.z() << " "
        << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << "\n";
    last_timestamp = state.timestamp_;
}
}  // namespace

/// 运行一个LIO前端，带可视化
int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_colorlogtostderr = true;
    FLAGS_stderrthreshold = google::INFO;
    std::string logfile = std::string(ROOT_DIR) + "/data/MapOffline_";
	google::SetLogFilenameExtension(".log");
	google::SetLogDestination(google::GLOG_INFO, logfile.c_str());

    google::ParseCommandLineFlags(&argc, &argv, true);
    if (FLAGS_input_bag.empty()) {
        LOG(ERROR) << "未指定输入数据";
        return -1;
    }

    using namespace lightning;

    RosbagIO rosbag(FLAGS_input_bag);

    SlamSystem::Options options;
    options.online_mode_ = false;

    SlamSystem slam(options);

    /// 实时模式好像掉帧掉的比较厉害？

    if (!slam.Init(FLAGS_config)) {
        LOG(ERROR) << "failed to init slam";
        return -1;
    }

    slam.StartSLAM("new_map");

    lightning::YAML_IO yaml(FLAGS_config);
    std::string lidar_topic = yaml.GetValue<std::string>("common", "lidar_topic");
    std::string livox_lidar_topic = yaml.GetValue<std::string>("common", "livox_lidar_topic");
    std::string imu_topic = yaml.GetValue<std::string>("common", "imu_topic");

    std::ofstream tum;
    double last_tum_timestamp = 0.0;
    if (!FLAGS_output_tum.empty()) {
        tum.open(FLAGS_output_tum);
        if (!tum.is_open()) {
            LOG(ERROR) << "failed to open output_tum: " << FLAGS_output_tum;
            return -1;
        }
        LOG(INFO) << "writing TUM trajectory to " << FLAGS_output_tum;
    }

    /// IMU 的处理
    rosbag.AddImuHandle(imu_topic,
                        [&slam](IMUPtr imu) {
                            slam.ProcessIMU(imu);
                            return true;
                        });

    /// PointCloud2 lidar 的处理
    if (!lidar_topic.empty() && lidar_topic != livox_lidar_topic) {
        rosbag.AddPointCloud2Handle(lidar_topic,
                                    [&slam, &tum, &last_tum_timestamp](sensor_msgs::msg::PointCloud2::SharedPtr msg) {
                                        slam.ProcessLidar(msg);
                                        WriteTumState(tum, slam.GetLioState(), last_tum_timestamp);
                                        return true;
                                    });
    }

    /// Livox CustomMsg lidar 的处理
    rosbag.AddLivoxCloudHandle(livox_lidar_topic,
                               [&slam, &tum, &last_tum_timestamp](livox_ros_driver2::msg::CustomMsg::SharedPtr cloud) {
                                   slam.ProcessLidar(cloud);
                                   WriteTumState(tum, slam.GetLioState(), last_tum_timestamp);
                                   return true;
                               });

    rosbag.Go();

    slam.SaveMap("");
    if (!FLAGS_output_tum.empty()) {
        const std::filesystem::path tum_path(FLAGS_output_tum);
        const auto parent = tum_path.parent_path();
        const auto stem = tum_path.stem().string();
        slam.SaveKeyframeTrajectoryTum((parent / (stem + "_keyframes_lio.tum")).string(), true);
        slam.SaveKeyframeTrajectoryTum((parent / (stem + "_keyframes_opt.tum")).string(), false);
    }
    if (tum.is_open()) {
        tum.close();
    }
    Timer::PrintAll();

    if (FLAGS_wait_ui) {
        slam.WaitForUIQuit();
    }

    LOG(INFO) << "done";

    return 0;
}
