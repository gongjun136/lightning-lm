#include <gflags/gflags.h>
#include <glog/logging.h>
#include <pcl/io/pcd_io.h>
#include <yaml-cpp/yaml.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <thread>

#include "common/options.h"
#include "core/lio/laser_mapping.h"
#include "core/lio/rear_axle_pose.h"
#include "io/yaml_io.h"
#include "ui/pangolin_window.h"
#include "utils/timer.h"
#include "wrapper/bag_io.h"

DEFINE_string(input_bag, "", "input ROS2 bag");
DEFINE_string(config, "./config/default.yaml", "config yaml");
DEFINE_string(output_tum, "", "output IMU-state TUM trajectory; disabled when empty");
DEFINE_string(output_lidar_tum, "", "output primary lidar TUM trajectory; disabled when empty");
DEFINE_string(output_rear_axle_tum, "", "output rebased rear-axle TUM trajectory; disabled when empty");
DEFINE_string(output_map, "", "output LIO map PCD; disabled when empty");
DEFINE_string(output_frame_stats_csv, "", "output per-fused-frame multi-lidar statistics; disabled when empty");
DEFINE_bool(wait_ui, true, "wait for the 3D UI window to close after offline processing");
DEFINE_int32(max_lidar_frames, 0, "stop after consuming this many fused lidar frames; disabled when <= 0");

namespace {

std::map<std::string, std::string> ReadBagTopicTypes(const std::string& bag_path) {
    std::map<std::string, std::string> types;
    const auto metadata_path = std::filesystem::path(bag_path) / "metadata.yaml";
    if (!std::filesystem::exists(metadata_path)) return types;
    try {
        const YAML::Node root = YAML::LoadFile(metadata_path.string());
        const YAML::Node topics = root["rosbag2_bagfile_information"]["topics_with_message_count"];
        if (!topics || !topics.IsSequence()) return types;
        for (const auto& item : topics) {
            const YAML::Node metadata = item["topic_metadata"];
            types[metadata["name"].as<std::string>()] = metadata["type"].as<std::string>();
        }
    } catch (const YAML::Exception& e) {
        LOG(WARNING) << "failed to parse rosbag metadata: " << e.what();
    }
    return types;
}

bool IsLivoxCustomMsg(const std::map<std::string, std::string>& types, const std::string& topic) {
    const auto it = types.find(topic);
    return it != types.end() && it->second == "livox_ros_driver2/msg/CustomMsg";
}

void OpenOutput(std::ofstream& stream, const std::string& path, const char* label) {
    if (path.empty()) return;
    stream.open(path);
    if (!stream.is_open()) throw std::runtime_error(std::string("failed to open ") + label + ": " + path);
}

void WriteTumPose(std::ofstream& stream, double timestamp, const lightning::SE3& pose, double& last_timestamp) {
    if (!stream.is_open() || timestamp <= 0.0 || timestamp <= last_timestamp) return;
    const auto q = pose.unit_quaternion();
    stream << std::fixed << std::setprecision(9) << timestamp << " " << std::setprecision(12)
           << pose.translation().x() << " " << pose.translation().y() << " " << pose.translation().z() << " "
           << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << "\n";
    last_timestamp = timestamp;
}

std::string JoinIds(const std::vector<int>& ids) {
    std::ostringstream output;
    for (std::size_t i = 0; i < ids.size(); ++i) {
        if (i != 0) output << ';';
        output << ids[i];
    }
    return output.str();
}

}  // namespace

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_colorlogtostderr = true;
    FLAGS_stderrthreshold = google::INFO;
    google::ParseCommandLineFlags(&argc, &argv, true);
    if (FLAGS_input_bag.empty()) {
        LOG(ERROR) << "input_bag is required";
        return 2;
    }

    using namespace lightning;
    LaserMapping lio;
    if (!lio.Init(FLAGS_config)) {
        LOG(ERROR) << "failed to init lio";
        return 2;
    }

    YAML_IO yaml(FLAGS_config);
    const bool with_ui = yaml.GetValue<bool>("system", "with_ui");
    const std::string lidar_topic = yaml.GetValue<std::string>("common", "lidar_topic");
    const std::string livox_lidar_topic = yaml.GetValue<std::string>("common", "livox_lidar_topic");
    std::string imu_topic = yaml.GetValue<std::string>("common", "imu_topic");
    if (lio.IsMultiLidarEnabled()) {
        const auto* primary = lio.GetMultiLidarConfig().PrimaryLidar();
        if (!primary) {
            LOG(ERROR) << "multi-lidar primary sensor is missing";
            return 2;
        }
        imu_topic = primary->imu_topic;
    }

    Vec3d primary_lidar_position_in_body(2.199, 0.0, 2.740);
    try {
        const YAML::Node root = YAML::LoadFile(FLAGS_config);
        if (root["output"] && root["output"]["primary_lidar_position_in_body"]) {
            const auto values = root["output"]["primary_lidar_position_in_body"].as<std::vector<double>>();
            if (values.size() != 3) throw std::runtime_error("primary_lidar_position_in_body must have 3 values");
            primary_lidar_position_in_body = Vec3d(values[0], values[1], values[2]);
        }
    } catch (const std::exception& e) {
        LOG(ERROR) << "invalid output pose configuration: " << e.what();
        return 2;
    }
    RearAxlePoseTransformer rear_axle(lio.GetLidarToImuRotation(), lio.GetLidarToImuTranslation(),
                                      primary_lidar_position_in_body);

    std::shared_ptr<ui::PangolinWindow> ui;
    if (with_ui) {
        ui = std::make_shared<ui::PangolinWindow>();
        if (ui->Init()) {
            lio.SetUI(ui);
        } else {
            LOG(ERROR) << "failed to init 3D UI, continue without Pangolin";
            ui.reset();
        }
    }

    std::ofstream imu_tum, lidar_tum, rear_tum, frame_stats_csv;
    try {
        OpenOutput(imu_tum, FLAGS_output_tum, "output_tum");
        OpenOutput(lidar_tum, FLAGS_output_lidar_tum, "output_lidar_tum");
        OpenOutput(rear_tum, FLAGS_output_rear_axle_tum, "output_rear_axle_tum");
        OpenOutput(frame_stats_csv, FLAGS_output_frame_stats_csv, "output_frame_stats_csv");
    } catch (const std::exception& e) {
        LOG(ERROR) << e.what();
        return 2;
    }
    double last_imu_tum = 0.0, last_lidar_tum = 0.0, last_rear_tum = 0.0;
    int consumed_frames = 0;
    int output_frames = 0;
    int complete_frames = 0;
    int partial_frames = 0;
    std::map<int, std::size_t> frames_by_lidar;
    std::map<int, std::size_t> points_by_lidar;

    if (frame_stats_csv.is_open()) {
        frame_stats_csv << "begin_time,end_time,partial,merged_points,present_lidar_ids,missing_lidar_ids";
        for (const auto& sensor : lio.GetMultiLidarConfig().lidars) {
            frame_stats_csv << ",points_lidar_" << sensor.id;
        }
        frame_stats_csv << "\n";
    }

    auto drain = [&]() {
        while (!lightning::debug::flg_exit) {
            const auto status = lio.RunDetailed();
            if (status == LaserMapping::RunStatus::kNoData) break;
            ++consumed_frames;
            if (lio.IsMultiLidarEnabled()) {
                const auto& stats = lio.GetCurrentFrameStats();
                if (stats.partial) {
                    ++partial_frames;
                } else {
                    ++complete_frames;
                }
                for (const int id : stats.present_lidar_ids) ++frames_by_lidar[id];
                for (const auto& [id, points] : stats.points_by_lidar) points_by_lidar[id] += points;
                if (frame_stats_csv.is_open()) {
                    frame_stats_csv << std::fixed << std::setprecision(9) << stats.begin_time << ',' << stats.end_time
                                    << ',' << (stats.partial ? 1 : 0) << ',' << stats.merged_points << ",\""
                                    << JoinIds(stats.present_lidar_ids) << "\",\"" << JoinIds(stats.missing_lidar_ids)
                                    << '"';
                    for (const auto& sensor : lio.GetMultiLidarConfig().lidars) {
                        const auto it = stats.points_by_lidar.find(sensor.id);
                        frame_stats_csv << ',' << (it == stats.points_by_lidar.end() ? 0 : it->second);
                    }
                    frame_stats_csv << '\n';
                }
            }
            if (status == LaserMapping::RunStatus::kOutput) {
                ++output_frames;
                const NavState state = lio.GetState();
                if (state.pose_is_ok_) {
                    WriteTumPose(imu_tum, state.timestamp_, state.GetPose(), last_imu_tum);
                    WriteTumPose(lidar_tum, state.timestamp_, rear_axle.PrimaryLidarPose(state), last_lidar_tum);
                    WriteTumPose(rear_tum, state.timestamp_, rear_axle.Transform(state), last_rear_tum);
                }
            }
            if (FLAGS_max_lidar_frames > 0 && consumed_frames >= FLAGS_max_lidar_frames) {
                lightning::debug::flg_exit = true;
            }
        }
    };

    RosbagIO rosbag(FLAGS_input_bag);
    rosbag.AddImuHandle(imu_topic, [&](IMUPtr imu) {
        lio.ProcessIMU(imu);
        drain();
        return true;
    });

    const auto topic_types = ReadBagTopicTypes(FLAGS_input_bag);
    if (lio.IsMultiLidarEnabled()) {
        for (const auto& sensor : lio.GetMultiLidarConfig().lidars) {
            if (IsLivoxCustomMsg(topic_types, sensor.lidar_topic)) {
                rosbag.AddLivoxCloudHandle(sensor.lidar_topic, [&, id = sensor.id](auto cloud) {
                    lio.ProcessPointCloud2(cloud, id);
                    drain();
                    return true;
                });
            } else {
                rosbag.AddPointCloud2Handle(sensor.lidar_topic, [&, id = sensor.id](auto cloud) {
                    lio.ProcessPointCloud2(cloud, id);
                    drain();
                    return true;
                });
            }
        }
    } else {
        if (!lidar_topic.empty() && lidar_topic != livox_lidar_topic) {
            rosbag.AddPointCloud2Handle(lidar_topic, [&](auto cloud) {
                lio.ProcessPointCloud2(cloud);
                drain();
                return true;
            });
        }
        if (!livox_lidar_topic.empty()) {
            rosbag.AddLivoxCloudHandle(livox_lidar_topic, [&](auto cloud) {
                lio.ProcessPointCloud2(cloud);
                drain();
                return true;
            });
        }
    }

    rosbag.Go();
    lio.FlushMultiLidar();
    drain();

    if (!FLAGS_output_map.empty() && !lio.GetAllKeyframes().empty()) {
        const auto map = lio.GetGlobalMap(true);
        if (pcl::io::savePCDFileBinaryCompressed(FLAGS_output_map, *map) != 0) {
            LOG(ERROR) << "failed to save output map: " << FLAGS_output_map;
            return 3;
        }
    }
    Timer::PrintAll();

    if (ui && FLAGS_wait_ui) {
        while (!ui->ShouldQuit()) std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    if (ui) ui->Quit();
    LOG(INFO) << "done, consumed fused frames=" << consumed_frames << ", output frames=" << output_frames;
    if (lio.IsMultiLidarEnabled()) {
        LOG(INFO) << "multi-lidar synchronization: complete=" << complete_frames << ", partial=" << partial_frames
                  << ", late_drops=" << lio.GetMultiLidarLateDropCount()
                  << ", duplicate_drops=" << lio.GetMultiLidarDuplicateDropCount()
                  << ", tolerance_drops=" << lio.GetMultiLidarToleranceDropCount()
                  << ", invalid_drops=" << lio.GetMultiLidarInvalidDropCount();
        LOG(INFO) << "multi-lidar lifecycle: assembled=" << lio.GetMultiLidarAssembledFrameCount()
                  << ", popped=" << consumed_frames << ", pre_imu_dropped=" << lio.GetPreImuDropCount()
                  << ", lio_with_imu=" << (consumed_frames - static_cast<int>(lio.GetPreImuDropCount()))
                  << ", lio_output=" << output_frames
                  << ", no_imu_tail_pending=" << lio.GetPendingLidarFrameCount()
                  << ", insufficient_lidar_drops=" << lio.GetMultiLidarInsufficientDropCount();
        for (const auto& sensor : lio.GetMultiLidarConfig().lidars) {
            LOG(INFO) << "lidar " << sensor.id << ": present_frames=" << frames_by_lidar[sensor.id]
                      << ", input_points=" << points_by_lidar[sensor.id];
        }
    }
    return 0;
}
