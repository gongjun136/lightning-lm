#include <gflags/gflags.h>
#include <glog/logging.h>
#include <pcl/io/pcd_io.h>
#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <thread>

#include "common/options.h"
#include "core/backend/backend_pipeline.h"
#include "core/lio/laser_mapping.h"
#include "core/loop_closing/loop_closing.h"
#include "core/maps/tiled_map.h"
#include "io/yaml_io.h"
#include "ui/pangolin_window.h"
#include "utils/timer.h"
#include "wrapper/bag_io.h"
#include "wrapper/ros_utils.h"

DEFINE_string(input_bag, "", "input ROS2 bag");
DEFINE_string(config, "./config/default.yaml", "config yaml");
DEFINE_string(output_tum, "", "output LIO TUM trajectory; disabled when empty");
DEFINE_string(output_map_dir, "./data/new_map", "output tiled map directory");
DEFINE_string(output_global_map, "", "output global map PCD; default is output_map_dir/global.pcd");
DEFINE_string(output_frame_stats_csv, "", "output per-fused-frame multi-lidar statistics; disabled when empty");
DEFINE_string(output_backend_diagnostics, "",
              "backend diagnostics directory; default is output_map_dir/backend_diagnostics");
DEFINE_bool(backend_evaluation_only, false,
            "save trajectories and backend diagnostics but skip map/relocalization exports");
DEFINE_bool(wait_ui, true, "wait for the 3D UI window to close after offline processing");
DEFINE_int32(max_lidar_frames, 0, "stop after consuming this many fused lidar frames; disabled when <= 0");
DEFINE_double(playback_rate, 0.0,
              "pace bag callbacks by sensor time at this multiple of real time; disabled when <= 0");

namespace {

class InputPacer {
   public:
    explicit InputPacer(double playback_rate) : playback_rate_(playback_rate) {}

    void Wait(double sensor_time) {
        if (playback_rate_ <= 0.0 || sensor_time <= 0.0) return;
        if (!initialized_) {
            first_sensor_time_ = sensor_time;
            first_wall_time_ = std::chrono::steady_clock::now();
            initialized_ = true;
            return;
        }
        const double elapsed_sensor_time = sensor_time - first_sensor_time_;
        if (elapsed_sensor_time <= 0.0) return;
        const auto target = first_wall_time_ + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                                                   std::chrono::duration<double>(elapsed_sensor_time / playback_rate_));
        std::this_thread::sleep_until(target);
    }

   private:
    double playback_rate_ = 0.0;
    bool initialized_ = false;
    double first_sensor_time_ = 0.0;
    std::chrono::steady_clock::time_point first_wall_time_;
};

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
    const std::filesystem::path output_path(path);
    if (!output_path.parent_path().empty()) {
        std::filesystem::create_directories(output_path.parent_path());
    }
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

bool SaveKeyframeTrajectoryTum(const std::vector<lightning::Keyframe::Ptr>& keyframes, const std::string& path,
                               bool use_lio_pose) {
    if (path.empty()) return true;
    std::ofstream tum(path);
    if (!tum.is_open()) {
        LOG(ERROR) << "failed to open keyframe trajectory: " << path;
        return false;
    }
    double last_timestamp = 0.0;
    int count = 0;
    for (const auto& kf : keyframes) {
        if (!kf) continue;
        const auto state = kf->GetState();
        if (state.timestamp_ <= 0.0 || state.timestamp_ <= last_timestamp) continue;
        const auto pose = use_lio_pose ? kf->GetLIOPose() : kf->GetOptPose();
        WriteTumPose(tum, state.timestamp_, pose, last_timestamp);
        ++count;
    }
    LOG(INFO) << "wrote " << count << " keyframe poses to " << path;
    return count > 0;
}

struct TimedPose {
    double timestamp = 0.0;
    lightning::SE3 pose;
};

bool SaveCorrectedFrameTrajectoryTum(const std::vector<TimedPose>& frames,
                                     const std::vector<lightning::Keyframe::Ptr>& keyframes,
                                     const std::string& path) {
    if (path.empty() || frames.empty() || keyframes.empty()) return false;
    std::vector<double> keyframe_times;
    std::vector<lightning::SE3> corrections;
    keyframe_times.reserve(keyframes.size());
    corrections.reserve(keyframes.size());
    for (const auto& keyframe : keyframes) {
        if (!keyframe) continue;
        keyframe_times.push_back(keyframe->GetState().timestamp_);
        corrections.push_back(keyframe->GetOptPose() * keyframe->GetLIOPose().inverse());
    }
    if (corrections.empty()) return false;

    std::ofstream output(path);
    if (!output.is_open()) return false;
    double last_timestamp = 0.0;
    std::size_t upper = 0;
    for (const auto& frame : frames) {
        while (upper < keyframe_times.size() && keyframe_times[upper] < frame.timestamp) ++upper;
        lightning::SE3 correction;
        if (upper == 0) {
            correction = corrections.front();
        } else if (upper >= keyframe_times.size()) {
            correction = corrections.back();
        } else {
            const std::size_t lower = upper - 1;
            const double duration = keyframe_times[upper] - keyframe_times[lower];
            const double ratio = duration > 1e-9
                                     ? std::clamp((frame.timestamp - keyframe_times[lower]) / duration, 0.0, 1.0)
                                     : 0.0;
            const lightning::Vec3d translation =
                (1.0 - ratio) * corrections[lower].translation() + ratio * corrections[upper].translation();
            const lightning::Quatd rotation = corrections[lower].unit_quaternion().slerp(
                ratio, corrections[upper].unit_quaternion());
            correction = lightning::SE3(rotation.normalized(), translation);
        }
        WriteTumPose(output, frame.timestamp, correction * frame.pose, last_timestamp);
    }
    LOG(INFO) << "wrote " << frames.size() << " corrected frame poses to " << path;
    return true;
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
    if (FLAGS_output_map_dir.empty()) {
        LOG(ERROR) << "output_map_dir is required";
        return 2;
    }
    if (FLAGS_playback_rate < 0.0) {
        LOG(ERROR) << "playback_rate must be non-negative";
        return 2;
    }

    using namespace lightning;

    LaserMapping lio;
    bool lio_initialized = false;
    Timer::Evaluate([&]() { lio_initialized = lio.Init(FLAGS_config); }, "Offline SLAM Initialization");
    if (!lio_initialized) {
        LOG(ERROR) << "failed to init lio";
        return 2;
    }

    YAML_IO yaml(FLAGS_config);
    const bool with_ui = yaml.GetValue<bool>("system", "with_ui");
    const bool with_backend = yaml.GetValue<bool>("system", "with_loop_closing");
    const backend::BackendMode backend_mode =
        with_backend ? backend::ReadBackendMode(FLAGS_config) : backend::BackendMode::kDisabled;
    std::shared_ptr<backend::BackendPipeline> new_backend;
    std::shared_ptr<LoopClosing> legacy_backend;
    if (backend_mode == backend::BackendMode::kBaBtcHba) {
        new_backend = std::make_shared<backend::BackendPipeline>();
        if (!new_backend->Init(FLAGS_config, false)) {
            LOG(ERROR) << "failed to initialize ba_btc_hba backend";
            return 2;
        }
    } else if (backend_mode == backend::BackendMode::kLegacy) {
        LoopClosing::Options loop_options;
        loop_options.online_mode_ = false;
        legacy_backend = std::make_shared<LoopClosing>(loop_options);
        legacy_backend->Init(FLAGS_config);
    }
    LOG(INFO) << "offline backend mode: " << backend::BackendModeName(backend_mode);
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

    std::ofstream tum, frame_stats_csv;
    try {
        OpenOutput(tum, FLAGS_output_tum, "output_tum");
        OpenOutput(frame_stats_csv, FLAGS_output_frame_stats_csv, "output_frame_stats_csv");
    } catch (const std::exception& e) {
        LOG(ERROR) << e.what();
        return 2;
    }
    double last_tum_timestamp = 0.0;

    if (frame_stats_csv.is_open()) {
        frame_stats_csv << "begin_time,end_time,partial,merged_points,present_lidar_ids,missing_lidar_ids";
        for (const auto& sensor : lio.GetMultiLidarConfig().lidars) {
            frame_stats_csv << ",points_lidar_" << sensor.id;
        }
        frame_stats_csv << "\n";
    }

    int consumed_frames = 0;
    int output_frames = 0;
    int complete_frames = 0;
    int partial_frames = 0;
    std::size_t dispatched_keyframes = 0;
    std::vector<TimedPose> frame_poses;
    std::map<int, std::size_t> frames_by_lidar;
    std::map<int, std::size_t> points_by_lidar;

    auto write_frame_stats = [&]() {
        if (!lio.IsMultiLidarEnabled()) return;
        const auto& stats = lio.GetCurrentFrameStats();
        if (stats.partial) {
            ++partial_frames;
        } else {
            ++complete_frames;
        }
        for (const int id : stats.present_lidar_ids) ++frames_by_lidar[id];
        for (const auto& [id, points] : stats.points_by_lidar) points_by_lidar[id] += points;
        if (!frame_stats_csv.is_open()) return;
        frame_stats_csv << std::fixed << std::setprecision(9) << stats.begin_time << ',' << stats.end_time << ','
                        << (stats.partial ? 1 : 0) << ',' << stats.merged_points << ",\""
                        << JoinIds(stats.present_lidar_ids) << "\",\"" << JoinIds(stats.missing_lidar_ids) << '"';
        for (const auto& sensor : lio.GetMultiLidarConfig().lidars) {
            const auto it = stats.points_by_lidar.find(sensor.id);
            frame_stats_csv << ',' << (it == stats.points_by_lidar.end() ? 0 : it->second);
        }
        frame_stats_csv << '\n';
    };

    auto drain = [&]() {
        while (!debug::flg_exit) {
            const auto status = lio.RunDetailed();
            if (status == LaserMapping::RunStatus::kNoData) break;
            ++consumed_frames;
            write_frame_stats();
            if (status == LaserMapping::RunStatus::kOutput) {
                ++output_frames;
                const NavState state = lio.GetState();
                if (state.pose_is_ok_) {
                    WriteTumPose(tum, state.timestamp_, state.GetPose(), last_tum_timestamp);
                    frame_poses.push_back({state.timestamp_, state.GetPose()});
                }
            }
            const auto current_keyframes = lio.GetAllKeyframes();
            while (dispatched_keyframes < current_keyframes.size()) {
                const auto& keyframe = current_keyframes[dispatched_keyframes++];
                if (new_backend) new_backend->AddKeyframe(keyframe);
                if (legacy_backend) legacy_backend->AddKF(keyframe);
            }
            if (FLAGS_max_lidar_frames > 0 && consumed_frames >= FLAGS_max_lidar_frames) {
                debug::flg_exit = true;
            }
        }
    };

    InputPacer input_pacer(FLAGS_playback_rate);
    RosbagIO rosbag(FLAGS_input_bag);
    rosbag.AddImuHandle(imu_topic, [&](IMUPtr imu) {
        input_pacer.Wait(imu->timestamp);
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

    Timer::Evaluate([&]() { rosbag.Go(); }, "Offline SLAM Bag Playback");
    Timer::Evaluate(
        [&]() {
            lio.FlushMultiLidar();
            drain();
        },
        "Offline SLAM Final Flush");

    if (new_backend) {
        Timer::Evaluate([&]() { new_backend->WaitUntilIdle(true); }, "Offline Backend Final Optimization");
    }

    const auto keyframes = lio.GetAllKeyframes();
    if (keyframes.empty()) {
        LOG(ERROR) << "no keyframes were generated; cannot export tiled map";
        Timer::PrintAll();
        return 4;
    }

    const std::filesystem::path map_dir(FLAGS_output_map_dir);
    const std::string global_map_path =
        FLAGS_output_global_map.empty() ? (map_dir / "global.pcd").string() : FLAGS_output_global_map;
    if (!FLAGS_backend_evaluation_only) {
        bool map_saved = false;
        Timer::Evaluate(
            [&]() {
                if (std::filesystem::exists(map_dir)) {
                    std::filesystem::remove_all(map_dir);
                }
                std::filesystem::create_directories(map_dir);
                const auto global_map = lio.GetGlobalMap(backend_mode == backend::BackendMode::kDisabled);
                TiledMap::Options tm_options;
                tm_options.map_path_ = map_dir.string();
                TiledMap tiled_map(tm_options);
                tiled_map.ConvertFromFullPCD(global_map, keyframes.front()->GetOptPose(), map_dir.string());
                map_saved = pcl::io::savePCDFileBinaryCompressed(global_map_path, *global_map) == 0;
            },
            "Offline Tiled Map Export");
        if (!map_saved) {
            LOG(ERROR) << "failed to save global map: " << global_map_path;
            Timer::PrintAll();
            return 4;
        }
    }

    if (!FLAGS_output_tum.empty()) {
        const std::filesystem::path tum_path(FLAGS_output_tum);
        const auto parent = tum_path.parent_path();
        const auto stem = tum_path.stem().string();
        SaveKeyframeTrajectoryTum(keyframes, (parent / (stem + "_keyframes_lio.tum")).string(), true);
        SaveKeyframeTrajectoryTum(keyframes, (parent / (stem + "_keyframes_opt.tum")).string(), false);
        SaveCorrectedFrameTrajectoryTum(frame_poses, keyframes, (parent / (stem + "_opt.tum")).string());
    }

    if (new_backend) {
        const std::string diagnostics = FLAGS_output_backend_diagnostics.empty()
                                            ? (map_dir / "backend_diagnostics").string()
                                            : FLAGS_output_backend_diagnostics;
        if (!new_backend->SaveDiagnostics(diagnostics)) {
            LOG(WARNING) << "failed to save backend diagnostics to " << diagnostics;
        }
        if (!FLAGS_backend_evaluation_only && new_backend->GetOptions().btc.enabled &&
            !new_backend->SaveRelocalizationDatabase((map_dir / "btc_relocalization").string())) {
            LOG(ERROR) << "failed to save required BTC relocalization database under " << map_dir;
            Timer::PrintAll();
            return 4;
        }
    }

    Timer::PrintAll();
    if (ui && FLAGS_wait_ui) {
        while (!ui->ShouldQuit()) std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    if (ui) ui->Quit();

    LOG(INFO) << "done, consumed fused frames=" << consumed_frames << ", output frames=" << output_frames
              << ", keyframes=" << keyframes.size() << ", tiled map=" << FLAGS_output_map_dir;
    if (lio.IsMultiLidarEnabled()) {
        LOG(INFO) << "multi-lidar SLAM synchronization: complete=" << complete_frames
                  << ", partial=" << partial_frames << ", late_drops=" << lio.GetMultiLidarLateDropCount()
                  << ", duplicate_drops=" << lio.GetMultiLidarDuplicateDropCount()
                  << ", tolerance_drops=" << lio.GetMultiLidarToleranceDropCount()
                  << ", invalid_drops=" << lio.GetMultiLidarInvalidDropCount();
        for (const auto& sensor : lio.GetMultiLidarConfig().lidars) {
            LOG(INFO) << "lidar " << sensor.id << ": present_frames=" << frames_by_lidar[sensor.id]
                      << ", input_points=" << points_by_lidar[sensor.id];
        }
    }
    return 0;
}
