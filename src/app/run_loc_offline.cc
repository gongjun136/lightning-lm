#include <gflags/gflags.h>
#include <geometry_msgs/msg/pose_stamped.hpp>
#include <glog/logging.h>
#include <rclcpp/rclcpp.hpp>
#include <rosbag2_cpp/writer.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <yaml-cpp/yaml.h>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <thread>
#include <vector>

#include "common/options.h"
#include "core/lio/laser_mapping.h"
#include "core/lio/multi_lidar_fusion.h"
#include "core/localization/lidar_loc/lidar_loc.h"
#include "core/maps/map_frame.h"
#include "core/localization/pose_graph/pgo.h"
#include "core/system/sany_localization_output.h"
#include "io/yaml_io.h"
#include "ui/pangolin_window.h"
#include "utils/timer.h"
#include "wrapper/bag_io.h"
#include "wrapper/ros_utils.h"

DEFINE_string(input_bag, "", "input ROS2 bag");
DEFINE_string(config, "./config/default.yaml", "config yaml");
DEFINE_string(map_path, "./data/new_map/", "tiled map directory containing index.txt and chunk PCD files");
DEFINE_string(output_tum, "", "output fused localization TUM trajectory; disabled when empty");
DEFINE_string(output_lidar_loc_tum, "", "output raw lidar-map localization TUM trajectory; disabled when empty");
DEFINE_string(output_csv, "", "output per-localization-frame CSV; disabled when empty");
DEFINE_string(output_frame_stats_csv, "", "output per-fused-frame multi-lidar statistics; disabled when empty");
DEFINE_bool(wait_ui, false, "wait for the 3D UI window to close after offline processing");
DEFINE_bool(publish_topics, true,
            "publish SANY localization outputs on /PosRes, /slamPoseRaw_topic, /LidarDataInv, and /LidarDataInL");
DEFINE_string(output_bag, "", "write the four SANY localization output topics directly to a ROS 2 bag");
DEFINE_bool(use_config_initial_pose, true,
            "use offline_localization.initial_pose from YAML; disable to require global initialization");
DEFINE_int32(max_lidar_frames, 0, "stop after consuming this many fused lidar frames; disabled when <= 0");
DEFINE_int32(force_relocalization_frame, 0,
             "offline test hook: request global relocalization immediately before this localization frame");
DEFINE_double(playback_rate, 0.0,
              "pace bag callbacks by sensor time at this multiple of real time; disabled when <= 0");

namespace {

class RosContextGuard {
   public:
    RosContextGuard(int argc, char** argv) { rclcpp::init(argc, argv); }
    ~RosContextGuard() {
        if (rclcpp::ok()) rclcpp::shutdown();
    }
};

class OfflineLocalizationPublisher {
   public:
    OfflineLocalizationPublisher(std::string map_frame, std::size_t lost_frame_threshold,
                                 const lightning::Vec3d& primary_lidar_position_in_body,
                                 double rear_axle_yaw_compensation_deg,
                                 bool publish_topics, const std::string& output_bag)
        : map_frame_(std::move(map_frame)),
          primary_lidar_position_in_body_(primary_lidar_position_in_body),
          rear_axle_yaw_compensation_deg_(rear_axle_yaw_compensation_deg),
          publication_gate_(lost_frame_threshold),
          map_cloud_decimator_(10) {
        if (publish_topics) {
            node_ = std::make_shared<rclcpp::Node>("offline_multi_lidar_localization");
            const auto pose_qos =
                rclcpp::QoS(rclcpp::KeepLast(1000)).reliable().durability_volatile();
            const auto cloud_qos =
                rclcpp::QoS(rclcpp::KeepLast(1)).reliable().durability_volatile();
            pos_res_pub_ = node_->create_publisher<geosun_msgs::msg::PosRes>("/PosRes", pose_qos);
            pose_pub_ =
                node_->create_publisher<geometry_msgs::msg::PoseStamped>("/slamPoseRaw_topic", pose_qos);
            inv_cloud_pub_ =
                node_->create_publisher<sensor_msgs::msg::PointCloud2>("/LidarDataInv", cloud_qos);
            map_cloud_pub_ =
                node_->create_publisher<sensor_msgs::msg::PointCloud2>("/LidarDataInL", cloud_qos);
        }
        if (!output_bag.empty()) {
            bag_writer_ = std::make_unique<rosbag2_cpp::Writer>();
            bag_writer_->open(output_bag);
            LOG(INFO) << "writing localization output topics directly to " << output_bag;
        }
    }

    void ObserveLidarMatch(bool valid) { publication_gate_.ObserveLidarMatch(valid); }

    void PublishPose(const lightning::loc::LocalizationResult& result,
                     const lightning::SO3& initial_lidar_rotation) {
        if (!result.valid_ || result.timestamp_ <= 0.0 ||
            !publication_gate_.MapOutputsEnabled()) {
            return;
        }
        const lightning::SO3 rear_axle_lidar_rotation =
            lightning::sany_output::ApplyRearAxleYawCompensation(
                initial_lidar_rotation, rear_axle_yaw_compensation_deg_);
        const lightning::SE3 map_rear_axle_pose = lightning::sany_output::MakeMapRearAxlePose(
            result.pose_, rear_axle_lidar_rotation, primary_lidar_position_in_body_);
        const auto position = lightning::sany_output::MakePosResMessage(
            map_rear_axle_pose, result.vel_b_.x(), result.timestamp_, map_frame_);
        if (pos_res_pub_) pos_res_pub_->publish(position);
        const auto pose = lightning::sany_output::MakePoseMessage(position);
        if (pose_pub_) pose_pub_->publish(pose);
        if (bag_writer_) {
            const rclcpp::Time stamp(position.header.stamp);
            bag_writer_->write(position, "/PosRes", stamp);
            bag_writer_->write(pose, "/slamPoseRaw_topic", stamp);
        }
    }

    void PublishCloud(const lightning::SE3& map_lidar_pose, const lightning::CloudPtr& cloud,
                      double begin_time, double end_time,
                      const lightning::SO3& initial_lidar_rotation) {
        if (!cloud || cloud->empty()) return;
        const lightning::SO3 rear_axle_lidar_rotation =
            lightning::sany_output::ApplyRearAxleYawCompensation(
                initial_lidar_rotation, rear_axle_yaw_compensation_deg_);
        const auto inv_cloud = lightning::sany_output::MakeCloudMessage(
            cloud, begin_time, end_time,
            lightning::sany_output::MakeRearAxleLidarTransform(rear_axle_lidar_rotation,
                                                               primary_lidar_position_in_body_),
            rear_axle_frame_);
        if (inv_cloud_pub_) inv_cloud_pub_->publish(inv_cloud);
        if (bag_writer_) {
            bag_writer_->write(inv_cloud, "/LidarDataInv", rclcpp::Time(inv_cloud.header.stamp));
        }
        const bool publish_map_frame = map_cloud_decimator_.Tick();
        if (publish_map_frame && publication_gate_.MapOutputsEnabled()) {
            const auto map_cloud = lightning::sany_output::MakeCloudMessage(
                cloud, begin_time, end_time, map_lidar_pose, map_frame_);
            if (map_cloud_pub_) map_cloud_pub_->publish(map_cloud);
            if (bag_writer_) {
                bag_writer_->write(map_cloud, "/LidarDataInL", rclcpp::Time(map_cloud.header.stamp));
            }
        }
    }

   private:
    rclcpp::Node::SharedPtr node_;
    rclcpp::Publisher<geosun_msgs::msg::PosRes>::SharedPtr pos_res_pub_;
    rclcpp::Publisher<geometry_msgs::msg::PoseStamped>::SharedPtr pose_pub_;
    rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr inv_cloud_pub_;
    rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr map_cloud_pub_;
    std::string map_frame_;
    std::string rear_axle_frame_ = "rear_axle";
    lightning::Vec3d primary_lidar_position_in_body_ = lightning::Vec3d::Zero();
    double rear_axle_yaw_compensation_deg_ = 0.0;
    lightning::sany_output::LocalizationPublicationGate publication_gate_;
    lightning::sany_output::FrameDecimator map_cloud_decimator_;
    std::unique_ptr<rosbag2_cpp::Writer> bag_writer_;
};

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

void WritePoseCsvFields(std::ofstream& csv, const lightning::SE3& pose) {
    const auto q = pose.unit_quaternion();
    csv << ',' << pose.translation().x() << ',' << pose.translation().y() << ',' << pose.translation().z() << ','
        << q.x() << ',' << q.y() << ',' << q.z() << ',' << q.w();
}

void WriteEmptyPoseCsvFields(std::ofstream& csv) {
    for (int i = 0; i < 7; ++i) csv << ',';
}

bool IsUsableLocResult(const lightning::loc::LocalizationResult& result) {
    return result.valid_ || result.lidar_loc_valid_;
}

bool ReadInitialPose(const std::string& config_path, const std::string& map_path,
                     lightning::SE3& pose) {
    const YAML::Node root = YAML::LoadFile(config_path);
    YAML::Node node = root["offline_localization"] ? root["offline_localization"]["initial_pose"] : YAML::Node();
    if (!node) return false;
    if (node["enabled"] && !node["enabled"].as<bool>()) return false;

    std::vector<double> translation;
    std::vector<double> quaternion_xyzw;
    if (node.IsSequence()) {
        const auto values = node.as<std::vector<double>>();
        if (values.size() != 7) {
            throw std::runtime_error("offline_localization.initial_pose sequence must be [x,y,z,qx,qy,qz,qw]");
        }
        translation = {values[0], values[1], values[2]};
        quaternion_xyzw = {values[3], values[4], values[5], values[6]};
    } else {
        if (node["translation"]) {
            translation = node["translation"].as<std::vector<double>>();
        } else if (node["position"]) {
            translation = node["position"].as<std::vector<double>>();
        } else {
            translation = {0.0, 0.0, 0.0};
        }

        if (node["quaternion_xyzw"]) {
            quaternion_xyzw = node["quaternion_xyzw"].as<std::vector<double>>();
        } else if (node["q_xyzw"]) {
            quaternion_xyzw = node["q_xyzw"].as<std::vector<double>>();
        } else if (node["rpy_deg"]) {
            const auto rpy = node["rpy_deg"].as<std::vector<double>>();
            if (rpy.size() != 3) throw std::runtime_error("offline_localization.initial_pose.rpy_deg must have 3 values");
            const double roll = rpy[0] * M_PI / 180.0;
            const double pitch = rpy[1] * M_PI / 180.0;
            const double yaw = rpy[2] * M_PI / 180.0;
            const lightning::Quatd q = Eigen::AngleAxisd(yaw, lightning::Vec3d::UnitZ()) *
                                       Eigen::AngleAxisd(pitch, lightning::Vec3d::UnitY()) *
                                       Eigen::AngleAxisd(roll, lightning::Vec3d::UnitX());
            quaternion_xyzw = {q.x(), q.y(), q.z(), q.w()};
        } else {
            quaternion_xyzw = {0.0, 0.0, 0.0, 1.0};
        }
    }

    if (translation.size() != 3 || quaternion_xyzw.size() != 4) {
        throw std::runtime_error("offline_localization.initial_pose must contain 3D translation and XYZW quaternion");
    }
    for (const double value : translation) {
        if (!std::isfinite(value)) throw std::runtime_error("initial pose translation contains non-finite value");
    }
    for (const double value : quaternion_xyzw) {
        if (!std::isfinite(value)) throw std::runtime_error("initial pose quaternion contains non-finite value");
    }

    lightning::Quatd q(quaternion_xyzw[3], quaternion_xyzw[0], quaternion_xyzw[1], quaternion_xyzw[2]);
    if (q.norm() < 1e-9) throw std::runtime_error("initial pose quaternion has zero norm");
    q.normalize();
    pose = lightning::SE3(q, lightning::Vec3d(translation[0], translation[1], translation[2]));
    const std::string frame =
        node.IsMap() && node["frame"] ? node["frame"].as<std::string>() : "map";
    if (frame == "slam") {
        lightning::map_frame::ExportOptions export_options;
        std::string map_frame_error;
        if (!lightning::map_frame::ReadExportOptions(
                root, export_options, map_frame_error)) {
            throw std::runtime_error(map_frame_error);
        }
        if (export_options.normalize_start_ground_z) {
            lightning::map_frame::Metadata metadata;
            if (!lightning::map_frame::LoadMetadata(
                    map_path, metadata, map_frame_error)) {
                throw std::runtime_error(map_frame_error);
            }
            pose = lightning::map_frame::TransformPose(metadata, pose);
        }
    } else if (frame != "map") {
        throw std::runtime_error(
            "offline_localization.initial_pose.frame must be 'map' or 'slam'");
    }
    return true;
}

void WriteLocalizationCsvHeader(std::ofstream& csv) {
    if (!csv.is_open()) return;
    csv << "frame_index,timestamp,status,valid,lidar_loc_valid,confidence,match_iterations,match_success,"
           "relocalization_attempted,relocalization_candidate_found,relocalization_accepted,"
           "relocalization_candidate_id,relocalization_score,"
           "map_consistency_evaluated,map_consistency_passed,map_consistency_points,"
           "map_inside_xy_ratio,map_inside_xyz_ratio,map_overlap_ratio,map_gravity_alignment_cos,"
           "active_map_chunks,processing_ms,loc_odom_delta,loc_odom_error_normal,smooth_flag,"
           "lo_x,lo_y,lo_z,lo_qx,lo_qy,lo_qz,lo_qw,"
           "lidar_loc_x,lidar_loc_y,lidar_loc_z,lidar_loc_qx,lidar_loc_qy,lidar_loc_qz,lidar_loc_qw,"
           "final_x,final_y,final_z,final_qx,final_qy,final_qz,final_qw,"
           "present_lidar_ids,missing_lidar_ids,merged_points\n";
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
    if (FLAGS_map_path.empty()) {
        LOG(ERROR) << "map_path is required";
        return 2;
    }
    if (FLAGS_playback_rate < 0.0) {
        LOG(ERROR) << "playback_rate must be non-negative";
        return 2;
    }
    if (!std::filesystem::exists(std::filesystem::path(FLAGS_map_path) / "index.txt")) {
        LOG(ERROR) << "tiled map index not found: " << (std::filesystem::path(FLAGS_map_path) / "index.txt");
        return 2;
    }

    using namespace lightning;

    std::unique_ptr<RosContextGuard> ros_context;
    std::unique_ptr<OfflineLocalizationPublisher> topic_publisher;
    if (FLAGS_publish_topics || !FLAGS_output_bag.empty()) {
        ros_context = std::make_unique<RosContextGuard>(argc, argv);
    }

    LaserMapping::Options lio_options;
    lio_options.is_in_slam_mode_ = false;
    LaserMapping lio(lio_options);
    bool lio_initialized = false;
    Timer::Evaluate([&]() { lio_initialized = lio.Init(FLAGS_config); }, "Offline Loc LIO Initialization");
    if (!lio_initialized) {
        LOG(ERROR) << "failed to init localization LIO frontend";
        return 2;
    }

    YAML_IO yaml(FLAGS_config);
    const YAML::Node root = YAML::LoadFile(FLAGS_config);
    const std::string map_frame =
        root["output"] && root["output"]["map_frame"] ? root["output"]["map_frame"].as<std::string>() : "map";
    Vec3d primary_lidar_position_in_body(2.199, 0.0, 2.740);
    if (root["output"] && root["output"]["primary_lidar_position_in_body"]) {
        const auto values = root["output"]["primary_lidar_position_in_body"].as<std::vector<double>>();
        if (values.size() != 3) {
            LOG(ERROR) << "output.primary_lidar_position_in_body must have 3 values";
            return 2;
        }
        primary_lidar_position_in_body = Vec3d(values[0], values[1], values[2]);
    }
    const double rear_axle_yaw_compensation_deg =
        root["output"] && root["output"]["rear_axle_yaw_compensation_deg"]
            ? root["output"]["rear_axle_yaw_compensation_deg"].as<double>()
            : 0.0;
    if (!std::isfinite(rear_axle_yaw_compensation_deg)) {
        LOG(ERROR) << "output.rear_axle_yaw_compensation_deg must be finite";
        return 2;
    }
    const int lost_frame_threshold =
        root["relocalization"] && root["relocalization"]["lost_frame_threshold"]
            ? root["relocalization"]["lost_frame_threshold"].as<int>()
            : 5;
    if (lost_frame_threshold <= 0) {
        LOG(ERROR) << "relocalization.lost_frame_threshold must be positive";
        return 2;
    }
    if (ros_context) {
        topic_publisher = std::make_unique<OfflineLocalizationPublisher>(
            map_frame, static_cast<std::size_t>(lost_frame_threshold),
            primary_lidar_position_in_body, rear_axle_yaw_compensation_deg,
            FLAGS_publish_topics, FLAGS_output_bag);
    }

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

    loc::LidarLoc::Options loc_options;
    loc_options.update_dynamic_cloud_ = yaml.GetValue<bool>("lidar_loc", "update_dynamic_cloud");
    loc_options.force_2d_ = yaml.GetValue<bool>("lidar_loc", "force_2d");
    loc_options.map_option_.enable_dynamic_polygon_ = false;
    loc_options.map_option_.map_path_ = FLAGS_map_path;
    if (!FLAGS_output_csv.empty()) {
        loc_options.relocalization_debug_dir_ =
            (std::filesystem::path(FLAGS_output_csv).parent_path() / "relocalization_debug").string();
    }

    auto lidar_loc = std::make_shared<loc::LidarLoc>(loc_options);
    std::shared_ptr<ui::PangolinWindow> ui;
    if (with_ui) {
        ui = std::make_shared<ui::PangolinWindow>();
        ui->SetCurrentScanSize(1);
        if (ui->Init()) {
            lidar_loc->SetUI(ui);
            lio.SetUI(ui);
        } else {
            LOG(ERROR) << "failed to init 3D UI, continue without Pangolin";
            ui.reset();
        }
    }
    if (!lidar_loc->Init(FLAGS_config)) {
        LOG(ERROR) << "failed to init lidar localization";
        return 2;
    }

    try {
        if (FLAGS_use_config_initial_pose) {
            SE3 initial_pose;
            if (ReadInitialPose(FLAGS_config, FLAGS_map_path, initial_pose)) {
                lidar_loc->SetInitialPose(initial_pose);
                LOG(INFO) << "offline localization initial pose: " << initial_pose.translation().transpose();
            } else {
                LOG(WARNING) << "offline localization initial pose is not configured; global initialization required";
            }
        } else {
            LOG(INFO) << "configured initial pose disabled; global BTC initialization required";
        }
    } catch (const std::exception& e) {
        LOG(ERROR) << "invalid offline localization initial pose: " << e.what();
        return 2;
    }

    loc::PGO pgo;
    pgo.SetDebug(false);
    const YAML::Node localization_pgo = YAML::LoadFile(FLAGS_config)["localization_pgo"];
    const bool enable_dr_smoothing =
        localization_pgo ? localization_pgo["enable_dr_smoothing"].as<bool>(true) : true;
    const bool enable_dr_extrapolation =
        localization_pgo ? localization_pgo["enable_dr_extrapolation"].as<bool>(true) : true;
    pgo.SetDrSmoothingEnabled(enable_dr_smoothing);
    pgo.SetDrExtrapolationEnabled(enable_dr_extrapolation);
    LOG(INFO) << "localization PGO DR smoothing: " << (enable_dr_smoothing ? "enabled" : "disabled");
    LOG(INFO) << "localization PGO DR extrapolation: "
              << (enable_dr_extrapolation ? "enabled" : "disabled");
    loc::LocalizationResult latest_final_result;
    bool latest_final_result_set = false;
    auto capture_final_result = [&](const loc::LocalizationResult& result) {
        latest_final_result = result;
        latest_final_result_set = true;
    };
    pgo.SetGlobalOutputHandleFunction(capture_final_result);
    pgo.SetHighFrequencyGlobalOutputHandleFunction([&](const loc::LocalizationResult& result) {
        capture_final_result(result);
        if (topic_publisher) topic_publisher->PublishPose(result, lio.GetInitialLidarRotation());
    });

    std::ofstream fused_tum, lidar_loc_tum, csv, frame_stats_csv;
    try {
        OpenOutput(fused_tum, FLAGS_output_tum, "output_tum");
        OpenOutput(lidar_loc_tum, FLAGS_output_lidar_loc_tum, "output_lidar_loc_tum");
        OpenOutput(csv, FLAGS_output_csv, "output_csv");
        OpenOutput(frame_stats_csv, FLAGS_output_frame_stats_csv, "output_frame_stats_csv");
    } catch (const std::exception& e) {
        LOG(ERROR) << e.what();
        return 2;
    }
    WriteLocalizationCsvHeader(csv);

    if (frame_stats_csv.is_open()) {
        frame_stats_csv << "begin_time,end_time,partial,merged_points,present_lidar_ids,missing_lidar_ids";
        for (const auto& sensor : lio.GetMultiLidarConfig().lidars) {
            frame_stats_csv << ",points_lidar_" << sensor.id;
        }
        frame_stats_csv << "\n";
    }

    double last_fused_tum = 0.0;
    double last_lidar_loc_tum = 0.0;
    int consumed_frames = 0;
    int output_frames = 0;
    int loc_frames = 0;
    int loc_valid_frames = 0;
    int complete_frames = 0;
    int partial_frames = 0;
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
                const NavState lo_state = lio.GetState();
                if (lo_state.pose_is_ok_) {
                    lidar_loc->ProcessLO(lo_state);
                    pgo.ProcessLidarOdom(lo_state);
                }

                const auto current_scan = lio.GetScanUndist();
                const auto scan = lio.GetProjCloud();
                if (FLAGS_force_relocalization_frame > 0 &&
                    loc_frames + 1 == FLAGS_force_relocalization_frame) {
                    LOG(WARNING) << "offline fault injection: forcing global relocalization before frame "
                                 << FLAGS_force_relocalization_frame;
                    lidar_loc->RequestGlobalRelocalization();
                }
                const auto start = std::chrono::steady_clock::now();
                lidar_loc->ProcessCloud(scan);
                const auto end = std::chrono::steady_clock::now();
                const double processing_ms =
                    std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start).count();

                const loc::LocalizationResult loc_result = lidar_loc->GetLocalizationResult();
                const auto match_stats = lidar_loc->GetLastMatchStats();
                if (topic_publisher) topic_publisher->ObserveLidarMatch(loc_result.lidar_loc_valid_);
                if (match_stats.relocalization_accepted) {
                    pgo.Reset();
                    latest_final_result_set = false;
                    LOG(WARNING) << "reset localization PGO after accepted BTC relocalization";
                }
                pgo.ProcessLidarLoc(loc_result);

                loc::LocalizationResult final_result = loc_result;
                if (latest_final_result_set && latest_final_result.timestamp_ >= loc_result.timestamp_ - 1e-6) {
                    final_result = latest_final_result;
                }
                if (topic_publisher) {
                    topic_publisher->PublishCloud(loc_result.pose_, current_scan ? current_scan : scan,
                                                  lio.GetLastFrameBeginTime(), lio.GetLastFrameEndTime(),
                                                  lio.GetInitialLidarRotation());
                }

                ++loc_frames;
                if (IsUsableLocResult(loc_result)) ++loc_valid_frames;
                if (IsUsableLocResult(loc_result)) {
                    WriteTumPose(lidar_loc_tum, loc_result.timestamp_, loc_result.pose_, last_lidar_loc_tum);
                }
                if (IsUsableLocResult(final_result)) {
                    WriteTumPose(fused_tum, loc_result.timestamp_, final_result.pose_, last_fused_tum);
                }

                if (IsUsableLocResult(final_result)) {
                    const NavState final_state = final_result.ToNavState();
                    if (ui) {
                        ui->UpdateNavState(final_state);
                        ui->UpdateRecentPose(final_result.pose_);
                        ui->UpdateScan(current_scan ? current_scan : scan, final_result.pose_);
                    }
                }

                if (csv.is_open()) {
                    const auto& stats = lio.GetCurrentFrameStats();
                    csv << std::fixed << std::setprecision(9) << loc_frames << ',' << loc_result.timestamp_ << ','
                        << static_cast<int>(loc_result.status_) << ',' << (IsUsableLocResult(final_result) ? 1 : 0)
                        << ',' << (loc_result.lidar_loc_valid_ ? 1 : 0) << ',' << std::setprecision(12)
                        << loc_result.confidence_ << ',' << match_stats.iterations << ','
                        << (match_stats.success ? 1 : 0) << ','
                        << (match_stats.relocalization_attempted ? 1 : 0) << ','
                        << (match_stats.relocalization_candidate_found ? 1 : 0) << ','
                        << (match_stats.relocalization_accepted ? 1 : 0) << ','
                        << match_stats.relocalization_candidate_id << ','
                        << match_stats.relocalization_score << ','
                        << (match_stats.map_consistency_evaluated ? 1 : 0) << ','
                        << (match_stats.map_consistency_passed ? 1 : 0) << ','
                        << match_stats.map_consistency_points << ','
                        << match_stats.map_inside_xy_ratio << ','
                        << match_stats.map_inside_xyz_ratio << ','
                        << match_stats.map_overlap_ratio << ','
                        << match_stats.map_gravity_alignment_cos << ','
                        << match_stats.active_map_chunks << ','
                        << processing_ms << ',' << loc_result.lidar_loc_odom_delta_ << ','
                        << (loc_result.lidar_loc_odom_error_normal_ ? 1 : 0) << ','
                        << (loc_result.lidar_loc_smooth_flag_ ? 1 : 0);
                    if (lo_state.pose_is_ok_) {
                        WritePoseCsvFields(csv, lo_state.GetPose());
                    } else {
                        WriteEmptyPoseCsvFields(csv);
                    }
                    if (IsUsableLocResult(loc_result)) {
                        WritePoseCsvFields(csv, loc_result.pose_);
                    } else {
                        WriteEmptyPoseCsvFields(csv);
                    }
                    if (IsUsableLocResult(final_result)) {
                        WritePoseCsvFields(csv, final_result.pose_);
                    } else {
                        WriteEmptyPoseCsvFields(csv);
                    }
                    csv << ",\"" << JoinIds(stats.present_lidar_ids) << "\",\"" << JoinIds(stats.missing_lidar_ids)
                        << "\"," << stats.merged_points << '\n';
                }
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
        const NavState dr_state = lio.GetIMUState();
        if (dr_state.pose_is_ok_) {
            lidar_loc->ProcessDR(dr_state);
            pgo.ProcessDR(dr_state);
        }
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

    Timer::Evaluate([&]() { rosbag.Go(); }, "Offline Loc Bag Playback");
    Timer::Evaluate(
        [&]() {
            lio.FlushMultiLidar();
            drain();
        },
        "Offline Loc Final Flush");

    Timer::PrintAll();
    lidar_loc->Finish();

    if (ui && FLAGS_wait_ui) {
        while (!ui->ShouldQuit()) std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    if (ui) ui->Quit();

    LOG(INFO) << "done, consumed fused frames=" << consumed_frames << ", lio output frames=" << output_frames
              << ", loc frames=" << loc_frames << ", valid loc frames=" << loc_valid_frames;
    if (lio.IsMultiLidarEnabled()) {
        LOG(INFO) << "multi-lidar localization synchronization: complete=" << complete_frames
                  << ", partial=" << partial_frames << ", late_drops=" << lio.GetMultiLidarLateDropCount()
                  << ", duplicate_drops=" << lio.GetMultiLidarDuplicateDropCount()
                  << ", tolerance_drops=" << lio.GetMultiLidarToleranceDropCount()
                  << ", invalid_drops=" << lio.GetMultiLidarInvalidDropCount();
        for (const auto& sensor : lio.GetMultiLidarConfig().lidars) {
            LOG(INFO) << "lidar " << sensor.id << ": present_frames=" << frames_by_lidar[sensor.id]
                      << ", input_points=" << points_by_lidar[sensor.id];
        }
    }
    return loc_frames > 0 && loc_valid_frames > 0 ? 0 : 4;
}
