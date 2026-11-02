//
// Created by xiang on 25-9-12.
//

#include "core/system/loc_system.h"
#include "common/imu_body_velocity.h"
#include "common/sany_wheel_speed_wire.h"
#include "diagnostic_build_info.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "core/localization/localization.h"
#include "core/localization/pose_covariance.h"
#include "utils/localization_dtc_reporting.h"
#include "utils/sany_radar_topology.h"
#include "utils/thread_scheduling.h"
#include "core/maps/map_manifest.h"
#include "common/debug_event.h"
#include "io/yaml_io.h"
#include "wrapper/ros_utils.h"
#include "yaml-cpp/yaml.h"

namespace lightning {
namespace {

double CloudStampSec(const CloudPtr& cloud) {
    if (!cloud) return 0.0;
    return static_cast<double>(cloud->header.stamp) * 1e-9;
}

std::int64_t CaptureStamp(const builtin_interfaces::msg::Time& stamp) {
    return std::int64_t(stamp.sec) * 1000000000LL + stamp.nanosec;
}
std::uint64_t CloudCaptureBytes(const sensor_msgs::msg::PointCloud2& cloud) {
    std::uint64_t bytes = sizeof(cloud) + cloud.data.capacity() + cloud.header.frame_id.capacity()
        + cloud.fields.capacity() * sizeof(sensor_msgs::msg::PointField);
    for (const auto& field : cloud.fields) bytes += field.name.capacity();
    return bytes;
}

functional_safety::HeartbeatConfig LoadHeartbeatConfig(rclcpp::Node& node) {
    const auto node_id = node.declare_parameter<std::int64_t>("functional_safety.node_id", 6);
    const auto topic = node.declare_parameter<std::string>(
        "functional_safety.heartbeat_topic", "/diagnostics/heartbeat/lightning_slam");
    const auto period_ms =
        node.declare_parameter<std::int64_t>("functional_safety.heartbeat_period_ms", 100);
    const auto qos_depth =
        node.declare_parameter<std::int64_t>("functional_safety.heartbeat_qos_depth", 1);
    const auto reliability = node.declare_parameter<std::string>(
        "functional_safety.heartbeat_qos_reliability", "best_effort");
    if (node_id <= 0 || node_id > std::numeric_limits<std::uint16_t>::max()) {
        throw std::invalid_argument("functional_safety.node_id must be in [1, 65535]");
    }
    if (period_ms <= 0 || qos_depth <= 0) {
        throw std::invalid_argument("functional safety heartbeat period and QoS depth must be positive");
    }
    if (reliability != "best_effort" && reliability != "reliable") {
        throw std::invalid_argument(
            "functional_safety.heartbeat_qos_reliability must be best_effort or reliable");
    }
    functional_safety::HeartbeatConfig config;
    config.node_id = static_cast<std::uint16_t>(node_id);
    config.topic = topic;
    config.period = std::chrono::milliseconds(period_ms);
    config.qos_depth = static_cast<std::size_t>(qos_depth);
    config.reliable = reliability == "reliable";
    return config;
}

}  // namespace

LocSystem::LocSystem(LocSystem::Options options) : options_(options) {
    /// handle ctrl-c
    signal(SIGINT, lightning::debug::SigHandle);
}

LocSystem::~LocSystem() {
    Finish();
}

bool LocSystem::Init(const std::string &yaml_path, const std::string &map_path_override) {
    // 阶段 1：重置状态并创建算法对象。
    finished_ = false;
    lidar_input_timestamp_gate_.Reset();
    pose_timestamp_gate_.Reset();
    speed_smoothing_shadow_.Reset();
    speed_smoothing_shadow_enabled_ = sany_output::SpeedSmoothingShadow::Enabled(
        std::getenv("LIGHTNING_LM_RUN_MODE"), std::getenv("LIGHTNING_LM_SPEED_SMOOTHING_SHADOW"),
        profiling::ReduceNonessentialOverhead());
    LOG(INFO) << "SPEED_SMOOTHING_SHADOW enabled=" << speed_smoothing_shadow_enabled_
              << " tau_ms=30 output=raw";
    last_localization_stamp_ = 0.0;
    last_published_pose_stamp_ = 0.0;
    localization_ever_good_ = false;
    // [loc-system-init-online-object]
    loc::Localization::Options opt;
    opt.online_mode_ = true;
    loc_ = std::make_shared<loc::Localization>(opt);
    // [loc-system-init-online-object]

    // 阶段 2：加载地图与输出配置。
    YAML_IO yaml(yaml_path);

    const YAML::Node root = YAML::LoadFile(yaml_path);
    if (root["system"] && root["system"]["pub_tf"]) {
        options_.pub_tf_ = root["system"]["pub_tf"].as<bool>();
    }
    std::string map_path = map_path_override;
    if (map_path.empty()) {
        const YAML::Node configured_map = root["system"] ? root["system"]["map_path"] : YAML::Node();
        if (!configured_map || configured_map.IsNull()) {
            LOG(ERROR) << "event=LOC_CONFIG_INVALID 在线定位必须设置 --map 或 system.map_path";
            return false;
        }
        map_path = configured_map.as<std::string>();
    }
    std::string fixed_transform_error;
    YAML::Node map_output;
    // 输出坐标与固定变换以地图包清单为依据，避免与当前地图使用不同的坐标约定。
    if (!map_manifest::LoadOutputConfig(map_path, map_output, fixed_transform_error)) {
        LOG(ERROR) << "event=LOC_MAP_CONFIG_REJECTED 地图输出配置被拒绝 detail=" << fixed_transform_error;
        return false;
    }
    YAML::Node output_root;
    output_root["output"] = map_output;
    map_frame_ = map_output["map_frame"].as<std::string>();
    if (!sany_output::LoadFixedMapTransform(output_root, fixed_map_transform_, fixed_transform_error)) {
        LOG(ERROR) << "event=LOC_MAP_FRAME_INVALID 地图坐标变换无效 detail=" << fixed_transform_error;
        return false;
    }
    if (fixed_map_transform_.enabled && fixed_map_transform_.target_frame != map_frame_) {
        LOG(ERROR) << "event=LOC_CONFIG_INVALID output.fixed_map_transform.target_frame 必须与 output.map_frame 一致 ("
                   << map_frame_ << ")";
        return false;
    }
    LOG(INFO) << "MAP_MANIFEST loaded=" << map_path << "/" << map_manifest::kFilename
              << " fixed_map_transform=" << fixed_map_transform_.enabled;
    primary_lidar_position_in_body_ = Vec3d(2.199, 0.0, 2.740);
    if (root["output"] && root["output"]["primary_lidar_position_in_body"]) {
        const auto values = root["output"]["primary_lidar_position_in_body"].as<std::vector<double>>();
        if (values.size() != 3) {
            LOG(ERROR) << "event=LOC_CONFIG_INVALID output.primary_lidar_position_in_body 必须包含 3 个值";
            return false;
        }
        primary_lidar_position_in_body_ = Vec3d(values[0], values[1], values[2]);
    }
    // 分别设置结果发布门控与输入时间戳门控，两者约束的对象不同。
    const int lost_frame_threshold =
        root["relocalization"] && root["relocalization"]["lost_frame_threshold"]
            ? root["relocalization"]["lost_frame_threshold"].as<int>()
            : 5;
    if (lost_frame_threshold <= 0) {
        LOG(ERROR) << "event=LOC_CONFIG_INVALID relocalization.lost_frame_threshold 必须为正数";
        return false;
    }
    publication_gate_.SetLostFrameThreshold(static_cast<std::size_t>(lost_frame_threshold));
    const double max_lidar_match_age_sec =
        root["system"] && root["system"]["localization_output_max_lidar_age_sec"]
            ? root["system"]["localization_output_max_lidar_age_sec"].as<double>()
            : 0.5;
    if (!std::isfinite(max_lidar_match_age_sec) || max_lidar_match_age_sec <= 0.0) {
        LOG(ERROR) << "event=LOC_CONFIG_INVALID system.localization_output_max_lidar_age_sec 必须为正数";
        return false;
    }
    publication_gate_.SetMaxLidarMatchAge(max_lidar_match_age_sec);
    online_lidar_input_max_timestamp_lag_sec_ =
        root["system"] && root["system"]["online_lidar_input_max_timestamp_lag_sec"]
            ? root["system"]["online_lidar_input_max_timestamp_lag_sec"].as<double>()
            : 0.0;
    if (!std::isfinite(online_lidar_input_max_timestamp_lag_sec_) ||
        online_lidar_input_max_timestamp_lag_sec_ < 0.0) {
        LOG(ERROR) << "event=LOC_CONFIG_INVALID system.online_lidar_input_max_timestamp_lag_sec 必须为非负有限数";
        return false;
    }
    lidar_input_timestamp_gate_.SetMaximumLag(
        online_lidar_input_max_timestamp_lag_sec_);
    // 先确定轮速开关、话题和换算参数；实际订阅在阶段 5 注册。
    wheel_speed_observation_enabled_ =
        root["system"] && root["system"]["enable_wheel_speed_observation"]
            ? root["system"]["enable_wheel_speed_observation"].as<bool>()
            : true;
    if (const char* override_value = std::getenv("SANY_ENABLE_CAN_OBSERVATION")) {
        const std::string value(override_value);
        if (value == "0") {
            wheel_speed_observation_enabled_ = false;
        } else if (value == "1") {
            wheel_speed_observation_enabled_ = true;
        } else {
            LOG(ERROR) << "event=LOC_CONFIG_INVALID SANY_ENABLE_CAN_OBSERVATION 必须为 0 或 1";
            return false;
        }
    }
    wheel_speed_topic_ =
        root["system"] && root["system"]["wheel_speed_topic"]
            ? root["system"]["wheel_speed_topic"].as<std::string>()
            : "/SpeThrCAN4_topic";
    if (const char* override_topic = std::getenv("SANY_WHEEL_SPEED_TOPIC")) {
        wheel_speed_topic_ = override_topic;
    }
    wheel_speed_scale_mps_per_rpm_ =
        root["system"] && root["system"]["wheel_speed_scale_mps_per_rpm"]
            ? root["system"]["wheel_speed_scale_mps_per_rpm"].as<double>()
            : 0.00120639253574024;
    if (wheel_speed_topic_.empty() || wheel_speed_topic_.front() != '/' ||
        !std::isfinite(wheel_speed_scale_mps_per_rpm_) ||
        wheel_speed_scale_mps_per_rpm_ <= 0.0) {
        LOG(ERROR) << "event=LOC_CONFIG_INVALID 轮速话题或转速换算比例无效";
        return false;
    }
    telemetry_ = std::make_unique<sany_output::LocalizationTelemetryState>(
        static_cast<std::size_t>(lost_frame_threshold));

    // 阶段 3：创建节点并初始化算法。
    // 心跳与算法后台线程在此启动；ROS 订阅和定时器回调由后续 Spin() 调度。
    LOG(INFO) << "event=LOC_NODE_CREATING 在线模式，正在创建 ROS 2 节点";
    node_ = std::make_shared<rclcpp::Node>("lightning_slam");
    try {
        heartbeat_ = std::make_unique<functional_safety::HeartbeatPublisher>(
            *node_, LoadHeartbeatConfig(*node_));
    } catch (const std::exception& error) {
        LOG(ERROR) << "event=LOC_HEARTBEAT_CONFIG_INVALID 功能安全心跳配置无效 detail=" << error.what();
        return false;
    }

    if (!loc_->Init(yaml_path, map_path)) {
        LOG(ERROR) << "event=LOC_INIT_FAILED 在线定位初始化失败";
        heartbeat_->SetState(functional_safety::NodeState::kFault);
        heartbeat_->PublishNow();
        return false;
    }

    // 阶段 4：建立诊断与 Safety 通道。
    const YAML::Node dtc = root["diagnostics"];
    try {
        const auto capture = diagnostics::IncidentRecorder::Options::FromYaml(dtc ? dtc["incident_capture"] : YAML::Node());
        if (capture.enabled) {
            incident_recorder_ = std::make_unique<diagnostics::IncidentRecorder>(capture, options_.incident_directory, yaml_path);
            using Kind = diagnostics::IncidentRecorder::Kind;
            incident_pose_ = incident_recorder_->Register<lightning::msg::VehiclePose>("/localization/pose_vel", Kind::kOutput);
            incident_pipeline_ = incident_recorder_->Register<lightning::msg::PipelineDiagnostics>("/localization/pipeline_diagnostics", Kind::kDiagnostic);
            incident_fault_ = incident_recorder_->Register<lightning::msg::FaultStatus>("/localization/fault_status", Kind::kDiagnostic);
            incident_status_ = incident_recorder_->Register<lightning::msg::LocalizationStatus>("/localization/loc_status", Kind::kDiagnostic);
            incident_heartbeat_ = incident_recorder_->Register<diagnostic_monitor_interfaces::msg::NodeHeartbeat>(
                node_->get_parameter("functional_safety.heartbeat_topic").as_string(), Kind::kDiagnostic);
        }
    } catch (const std::exception& e) {
        LOG(ERROR) << "INCIDENT_CAPTURE event=INCIDENT_DISABLED 事件采集已禁用 detail=" << e.what();
        incident_recorder_.reset();
    }
    loc::dtc::Config dtc_config;
    auto value = [&dtc](const char* key, double fallback) {
        return dtc && dtc[key] ? dtc[key].as<double>() : fallback;
    };
    dtc_config.window_sec = value("motion_window_sec", dtc_config.window_sec);
    dtc_config.max_input_gap_sec = value("max_input_gap_sec", dtc_config.max_input_gap_sec);
    dtc_config.longitudinal_m = value("longitudinal_tolerance_m", dtc_config.longitudinal_m);
    dtc_config.lateral_m = value("lateral_tolerance_m", dtc_config.lateral_m);
    dtc_config.yaw_rad = value("yaw_tolerance_deg", 2.0) * M_PI / 180.0;
    dtc_config.distance_ratio = value("distance_tolerance_ratio", dtc_config.distance_ratio);
    dtc_config.angle_ratio = value("angle_tolerance_ratio", dtc_config.angle_ratio);
    dtc_config.major_sigma_m = value("major_sigma_limit_m", dtc_config.major_sigma_m);
    dtc_config.lateral_sigma_m = value("lateral_sigma_limit_m", dtc_config.lateral_sigma_m);
    dtc_config.yaw_sigma_rad = value("yaw_sigma_limit_deg", 2.0) * M_PI / 180.0;
    dtc_config.covariance_scale = value("covariance_variance_scale", dtc_config.covariance_scale);
    dtc_config.position_sigma_floor_m = value("position_sigma_floor_m", dtc_config.position_sigma_floor_m);
    dtc_config.angle_sigma_floor_rad = value("angle_sigma_floor_deg", 0.5) * M_PI / 180.0;
    dtc_config.output_timeout_sec = value("output_timeout_sec", dtc_config.output_timeout_sec);
    dtc_config.evaluation_timeout_sec = value("evaluation_timeout_sec", dtc_config.evaluation_timeout_sec);
    dtc_config.consistency_unknown_grace_sec =
        value("consistency_unknown_grace_sec", dtc_config.consistency_unknown_grace_sec);
    dtc_config.startup_grace_sec = value("startup_grace_sec", dtc_config.startup_grace_sec);
    dtc_config.max_speed_mps = value("max_speed_mps", dtc_config.max_speed_mps);
    dtc_config.max_yaw_rate_radps = value("max_yaw_rate_radps", dtc_config.max_yaw_rate_radps);
    std::array<bool, 4> required_radar_roles;
    try {
        required_radar_roles = functional_safety::RequiredSanyRadarRoles(root);
        dtc_monitor_ = std::make_unique<loc::dtc::Monitor>(dtc_config);
        functional_safety::LidarFrameWatchdog::Config frame_config;
        frame_config.warning_sec = value("radar_frame_warning_sec", frame_config.warning_sec);
        frame_config.timeout_sec = value("radar_frame_timeout_sec", frame_config.timeout_sec);
        frame_config.startup_grace_sec = dtc_config.startup_grace_sec;
        const int recovery = dtc && dtc["radar_frame_recovery_frames"] ?
            dtc["radar_frame_recovery_frames"].as<int>() : static_cast<int>(frame_config.recovery_frames);
        if (recovery < 1) throw std::invalid_argument("radar frame recovery count must be positive");
        frame_config.recovery_frames = static_cast<unsigned>(recovery);
        lidar_frame_watchdog_ = std::make_unique<loc::RadarInputFrameMonitor>(loc_->GetMultiLidarConfig(), frame_config);
        safety_fault_reporting_enabled_ = dtc && dtc["safety_fault_reporting_enabled"]
                                              ? dtc["safety_fault_reporting_enabled"].as<bool>()
                                              : true;
        radar_time_sync_shadow_enabled_ = dtc && dtc["radar_time_sync_shadow_enabled"]
                                             ? dtc["radar_time_sync_shadow_enabled"].as<bool>()
                                             : true;
        radar_bounded_sync_warning_enabled_ = dtc && dtc["radar_bounded_sync_warning_enabled"]
                                                 ? dtc["radar_bounded_sync_warning_enabled"].as<bool>()
                                                 : false;
        radar_fault_mask_enabled_ = dtc && dtc["radar_fault_mask_enabled"]
                                        ? dtc["radar_fault_mask_enabled"].as<bool>() : false;
        const double radar_timeout = value("radar_snapshot_timeout_sec", 0.5);
        const double hms_timeout = value("radar_hms_timeout_sec", 2.5);
        if (!std::isfinite(radar_timeout) || radar_timeout <= 0 ||
            !std::isfinite(hms_timeout) || hms_timeout < radar_timeout)
            throw std::invalid_argument("invalid radar DTC source or timeout");
        std::array<uint16_t, 4> radar_ids{{7, 12, 8, 9}};
        if (dtc && dtc["radar_role_source_ids"]) {
            const auto ids = dtc["radar_role_source_ids"].as<std::vector<int>>();
            if (ids.size() != 4 || std::set<int>(ids.begin(), ids.end()).size() != 4)
                throw std::invalid_argument("four unique radar role IDs required");
            for (size_t i = 0; i < 4; ++i) {
                if (ids[i] <= 0 || ids[i] > 65535 || ids[i] == heartbeat_->NodeId())
                    throw std::invalid_argument("invalid radar role source ID");
                radar_ids[i] = ids[i];
            }
        }
        radar_dtc_state_ = std::make_unique<functional_safety::RadarDtcState>(
            radar_timeout, hms_timeout, radar_ids, dtc_config.startup_grace_sec,
            functional_safety::RadarDtcState::Clock::now(), required_radar_roles);
        if (dtc && dtc["config_version"]) dtc_config_version_ = dtc["config_version"].as<std::string>();
        if (dtc_config_version_.empty()) throw std::invalid_argument("empty DTC config version");
    } catch (const std::exception& error) {
        LOG(ERROR) << "event=LOC_DTC_CONFIG_INVALID DTC 配置无效 detail=" << error.what(); return false;
    }
    LOG(INFO) << "LOCALIZATION_DIAGNOSTIC_BUILD schema=2 contract_id=" << diagnostics::kDiagnosticContract
              << " source_revision=" << diagnostics::kDiagnosticRevision
              << " source_dirty=" << diagnostics::kDiagnosticDirty
              << " built_at=" << diagnostics::kDiagnosticBuiltAt;
    LOG(INFO) << "LOCALIZATION_DTC_CONFIG version=" << dtc_config_version_
              << " policy_version=" << loc::dtc::kSafetyPolicyVersion
              << " safety_fault_reporting_enabled=" << safety_fault_reporting_enabled_
              << " radar_fault_mask_enabled=" << radar_fault_mask_enabled_
              << " radar_bounded_sync_warning_enabled=" << radar_bounded_sync_warning_enabled_
              << " max_input_gap_sec=" << dtc_config.max_input_gap_sec
              << " consistency_unknown_grace_sec=" << dtc_config.consistency_unknown_grace_sec
              << " output_timeout_sec=" << dtc_config.output_timeout_sec
              << " evaluation_timeout_sec=" << dtc_config.evaluation_timeout_sec
              << " startup_grace_sec=" << dtc_config.startup_grace_sec
              << " radar_frame_observation_point=raw_input_grouping"
              << " radar_frame_warning_sec=" << value("radar_frame_warning_sec", 0.2)
              << " radar_frame_timeout_sec=" << value("radar_frame_timeout_sec", 0.3)
              << " radar_frame_recovery_frames=" << (dtc && dtc["radar_frame_recovery_frames"] ?
                    dtc["radar_frame_recovery_frames"].as<unsigned>() : 3)
              << " localization_output_max_lidar_age_sec=" << max_lidar_match_age_sec
              << " lost_frame_threshold=" << lost_frame_threshold
              << " degraded_continuous_output_sec=" << dtc_config.window_sec
              << " covariance_variance_scale=" << dtc_config.covariance_scale
              << " major_sigma_limit_m=" << dtc_config.major_sigma_m
              << " lateral_sigma_limit_m=" << dtc_config.lateral_sigma_m
              << " yaw_sigma_limit_deg=" << dtc_config.yaw_sigma_rad * 180.0 / M_PI;
    // 诊断配置就绪后连接 Safety 输出、必需雷达的状态输入及周期评估定时器。
    const auto dtc_qos = rclcpp::QoS(10).reliable().durability_volatile();
    const auto topic = [&dtc](const char* key, const std::string& fallback) {
        return dtc && dtc[key] ? dtc[key].as<std::string>() : fallback;
    };
    safety_pub_ = node_->create_publisher<geosun_msgs::msg::Safety>(
        topic("safety_topic", "/safety/throstle"), dtc_qos);
    if (incident_recorder_) incident_safety_ = incident_recorder_->Register<geosun_msgs::msg::Safety>(
        topic("safety_topic", "/safety/throstle"), diagnostics::IncidentRecorder::Kind::kDiagnostic);
    const std::array<std::string, 4> radar_roles{{"front", "rear", "left", "right"}};
    for (size_t i = 0; i < 4; ++i) {
        LOG(INFO) << "LOCALIZATION_RADAR_REQUIRED role=" << radar_roles[i]
                  << " required=" << required_radar_roles[i];
        if (!required_radar_roles[i]) continue;
        const auto radar_topic = "/sensors/livox/" + radar_roles[i] + "/diagnostic_state";
        const int capture_id = incident_recorder_ ? incident_recorder_->Register<livox_ros_driver2::msg::LidarDiagnosticState>(
            radar_topic, diagnostics::IncidentRecorder::Kind::kDiagnostic) : -1;
        radar_availability_subs_[i] = node_->create_subscription<livox_ros_driver2::msg::LidarDiagnosticState>(
            "/sensors/livox/" + radar_roles[i] + "/diagnostic_state", dtc_qos,
            [this, i, capture_id](const livox_ros_driver2::msg::LidarDiagnosticState::SharedPtr msg) {
                if (incident_recorder_) incident_recorder_->Capture(capture_id, msg, sizeof(*msg) + 4096);
                if (msg->role == i) radar_dtc_state_->Observe(*msg);
            });
    }
    dtc_timer_ = node_->create_wall_timer(std::chrono::milliseconds(100), [this]() { UpdateDtcHeartbeat(); });
    safety_timer_ = node_->create_wall_timer(std::chrono::milliseconds(20), [this]() {
        PublishSafety();
        const auto snapshot = EvaluateLocalizationDtc();
        bool changed = false;
        {
            std::lock_guard<std::mutex> lock(dtc_publish_mutex_);
            for (size_t i = 0; i < 2; ++i) changed |= last_dtc_conditions_[i] != snapshot.observations[i].condition;
        }
        if (changed) UpdateDtcHeartbeat();
    });

    // 阶段 5：注册传感器与轮速输入。
    imu_topic_ = yaml.GetValue<std::string>("common", "imu_topic");
    cloud_topic_ = yaml.GetValue<std::string>("common", "lidar_topic");
    livox_topic_ = yaml.GetValue<std::string>("common", "livox_lidar_topic");

    const auto imu_qos = rclcpp::SensorDataQoS();
    const YAML::Node lidar_qos_config = root["system"] ? root["system"]["lidar_qos"] : YAML::Node();
    const size_t lidar_qos_depth =
        lidar_qos_config && lidar_qos_config["depth"]
            ? lidar_qos_config["depth"].as<size_t>()
            : 5;
    if (lidar_qos_depth == 0) {
        LOG(ERROR) << "event=LOC_CONFIG_INVALID system.lidar_qos.depth 必须为正数";
        return false;
    }
    const std::string lidar_qos_reliability =
        lidar_qos_config && lidar_qos_config["reliability"]
            ? lidar_qos_config["reliability"].as<std::string>()
            : "best_effort";
    rclcpp::QoS lidar_qos{rclcpp::KeepLast(lidar_qos_depth)};
    lidar_qos.durability_volatile();
    if (lidar_qos_reliability == "reliable") {
        lidar_qos.reliable();
    } else if (lidar_qos_reliability == "best_effort") {
        lidar_qos.best_effort();
    } else {
        LOG(ERROR) << "event=LOC_CONFIG_INVALID system.lidar_qos.reliability 必须为 reliable 或 best_effort";
        return false;
    }
    LOG(INFO) << "event=LOC_LIDAR_QOS 雷达输入 QoS reliability=" << lidar_qos_reliability
              << ", depth=" << lidar_qos_depth;

    auto subscribe_imu = [this, &imu_qos](const std::string& topic) {
        const int capture_id = incident_recorder_ ? incident_recorder_->Register<sensor_msgs::msg::Imu>(
            topic, diagnostics::IncidentRecorder::Kind::kImu) : -1;
        {
            std::lock_guard<std::mutex> lock(input_stats_mutex_);
            imu_input_stats_.topic = topic;
        }
        imu_sub_ = node_->create_subscription<sensor_msgs::msg::Imu>(
            topic, imu_qos, [this, capture_id](sensor_msgs::msg::Imu::SharedPtr msg) {
            const auto arrived = std::chrono::steady_clock::now();
            if (incident_recorder_) incident_recorder_->Capture(capture_id, msg,
                sizeof(*msg) + msg->header.frame_id.capacity(), CaptureStamp(msg->header.stamp));
            IMUPtr imu = std::make_shared<IMU>();
            imu->timestamp = ToSec(msg->header.stamp);
            imu->linear_acceleration =
                Vec3d(msg->linear_acceleration.x, msg->linear_acceleration.y, msg->linear_acceleration.z);
            imu->angular_velocity = Vec3d(msg->angular_velocity.x, msg->angular_velocity.y, msg->angular_velocity.z);

            ProcessIMU(imu, arrived);
        });
    };

    // 多雷达模式只订阅主雷达 IMU，并为各点云回调保留对应雷达 ID。
    if (loc_->IsMultiLidarEnabled()) {
        const auto* primary = loc_->GetMultiLidarConfig().PrimaryLidar();
        if (!primary) {
            LOG(ERROR) << "event=LOC_PRIMARY_LIDAR_MISSING 多雷达配置缺少主雷达";
            return false;
        }
        subscribe_imu(primary->imu_topic);
        for (const auto& sensor : loc_->GetMultiLidarConfig().lidars) {
            RegisterLidarInput(sensor.id, sensor.lidar_topic);
            const int capture_id = incident_recorder_ ? incident_recorder_->Register<sensor_msgs::msg::PointCloud2>(
                sensor.lidar_topic, diagnostics::IncidentRecorder::Kind::kLidar) : -1;
            cloud_subs_.push_back(node_->create_subscription<sensor_msgs::msg::PointCloud2>(
                sensor.lidar_topic, lidar_qos,
                [this, id = sensor.id, capture_id](sensor_msgs::msg::PointCloud2::SharedPtr cloud) {
                    const auto arrived = std::chrono::steady_clock::now();
                    lidar_frame_watchdog_->ObserveInput(id, ToSec(cloud->header.stamp),
                        cloud->width > 0 && cloud->height > 0 && !cloud->data.empty());
                    if (incident_recorder_) incident_recorder_->Capture(capture_id, cloud,
                        CloudCaptureBytes(*cloud), CaptureStamp(cloud->header.stamp));
                    Timer::Evaluate([&]() { ProcessLidar(cloud, id, arrived); }, "Proc Lidar",
                                    !profiling::ReduceNonessentialOverhead());
                }));
        }
        LOG(INFO) << "event=LOC_INPUT_SUBSCRIBED 在线定位已订阅，雷达话题数量=" << cloud_subs_.size()
                  << "，主 IMU 话题=" << primary->imu_topic;
    } else {
        if (imu_topic_.empty() || (cloud_topic_.empty() && livox_topic_.empty())) {
            LOG(ERROR) << "event=LOC_INPUT_TOPICS_INCOMPLETE 单雷达在线输入话题不完整";
            return false;
        }
        subscribe_imu(imu_topic_);
        if (!cloud_topic_.empty()) {
            RegisterLidarInput(0, cloud_topic_);
            const int capture_id = incident_recorder_ ? incident_recorder_->Register<sensor_msgs::msg::PointCloud2>(
                cloud_topic_, diagnostics::IncidentRecorder::Kind::kLidar) : -1;
            cloud_subs_.push_back(node_->create_subscription<sensor_msgs::msg::PointCloud2>(
                cloud_topic_, lidar_qos, [this, capture_id](sensor_msgs::msg::PointCloud2::SharedPtr cloud) {
                    const auto arrived = std::chrono::steady_clock::now();
                    lidar_frame_watchdog_->ObserveInput(0, ToSec(cloud->header.stamp),
                        cloud->width > 0 && cloud->height > 0 && !cloud->data.empty());
                    if (incident_recorder_) incident_recorder_->Capture(capture_id, cloud,
                        CloudCaptureBytes(*cloud), CaptureStamp(cloud->header.stamp));
                    Timer::Evaluate([&]() { ProcessLidar(cloud, arrived); }, "Proc Lidar",
                                    !profiling::ReduceNonessentialOverhead());
                }));
        }
        if (!livox_topic_.empty()) {
            RegisterLidarInput(0, livox_topic_);
            const int capture_id = incident_recorder_ ? incident_recorder_->Register<livox_ros_driver2::msg::CustomMsg>(
                livox_topic_, diagnostics::IncidentRecorder::Kind::kLidar) : -1;
            livox_sub_ = node_->create_subscription<livox_ros_driver2::msg::CustomMsg>(
                livox_topic_, lidar_qos, [this, capture_id](livox_ros_driver2::msg::CustomMsg::SharedPtr cloud) {
                    const auto arrived = std::chrono::steady_clock::now();
                    lidar_frame_watchdog_->ObserveInput(0, ToSec(cloud->header.stamp),
                        cloud->point_num > 0 && !cloud->points.empty());
                    if (incident_recorder_) incident_recorder_->Capture(capture_id, cloud,
                        sizeof(*cloud) + cloud->points.capacity() * sizeof(livox_ros_driver2::msg::CustomPoint), CaptureStamp(cloud->header.stamp));
                    Timer::Evaluate([&]() { ProcessLidar(cloud, arrived); }, "Proc Lidar",
                                    !profiling::ReduceNonessentialOverhead());
                });
        }
    }

    {
        std::lock_guard<std::mutex> lock(input_stats_mutex_);
        wheel_speed_input_stats_.topic = wheel_speed_topic_;
    }
    if (wheel_speed_observation_enabled_) {
        const int capture_id = incident_recorder_ ? incident_recorder_->Register(wheel_speed_topic_,
            "geosun_msgs/msg/SpeThrCAN4", diagnostics::IncidentRecorder::Kind::kWheelSpeed,
            [](const void* value) { return std::make_shared<rclcpp::SerializedMessage>(
                *static_cast<const rclcpp::SerializedMessage*>(value)); }) : -1;
        wheel_speed_sub_ = node_->create_generic_subscription(
            wheel_speed_topic_, "geosun_msgs/msg/SpeThrCAN4",
            rclcpp::QoS(rclcpp::KeepLast(50)).reliable().durability_volatile(),
            [this, capture_id](std::shared_ptr<rclcpp::SerializedMessage> message) {
                const auto arrived = std::chrono::steady_clock::now();
                const auto& bytes = message->get_rcl_serialized_message();
                if (incident_recorder_) incident_recorder_->Capture(capture_id, message,
                    sizeof(*message) + bytes.buffer_capacity);
                SanyWheelSpeedWire sample;
                if (!DecodeSanyWheelSpeed(bytes.buffer, bytes.buffer_length, sample)) {
                    // 无法解码源时间时仍保留真实回调，避免误诊为订阅静默。
                    {
                        std::lock_guard<std::mutex> lock(input_stats_mutex_);
                        wheel_speed_input_stats_.timing.Observe(std::numeric_limits<double>::quiet_NaN(), arrived);
                    }
                    LOG_EVERY_N(ERROR, 200) << "event=LOC_CAN_LAYOUT_REJECTED 拒绝格式错误或不支持的 SpeThrCAN4 消息布局 bytes="
                                            << bytes.buffer_length;
                    return;
                }
                LOG_FIRST_N(INFO, 1) << "event=LOC_CAN_LAYOUT CAN 消息布局 layout="
                                    << (sample.has_comm_header ? "comm_header+Header+x+y" : "Header+x+y")
                                    << ", stamp=" << std::setprecision(16) << sample.stamp;
                ObserveWheelSpeedInput(sample.stamp, sample.rpm, sample.torque, arrived);
            });
        LOG(INFO) << "event=LOC_CAN_ENABLED 已启用轮速观测 topic=" << wheel_speed_topic_
                  << ", scale=" << std::setprecision(16)
                  << wheel_speed_scale_mps_per_rpm_ << " m/s/rpm";
    } else {
        LOG(WARNING) << "event=LOC_CAN_DISABLED 已禁用轮速观测，使用激光与 IMU 定位";
    }

    // 阶段 6：建立业务输出并注册算法回调。
    const auto pose_qos =
        rclcpp::QoS(rclcpp::KeepLast(1000)).best_effort().durability_volatile();
    const auto cloud_qos = rclcpp::QoS(rclcpp::KeepLast(1)).best_effort().durability_volatile();
    pose_pub_ = node_->create_publisher<geometry_msgs::msg::PoseStamped>("/slamPoseRaw_topic", pose_qos);
    vehicle_pose_pub_ = node_->create_publisher<lightning::msg::VehiclePose>(
        "/localization/pose_vel", pose_qos);
    inv_cloud_pub_ = node_->create_publisher<sensor_msgs::msg::PointCloud2>("/LidarDataInv", cloud_qos);
    map_cloud_pub_ = node_->create_publisher<sensor_msgs::msg::PointCloud2>("/LidarDataInL", cloud_qos);
    const auto health_qos = rclcpp::QoS(rclcpp::KeepLast(10)).best_effort().durability_volatile();
    const auto path_qos = rclcpp::QoS(rclcpp::KeepLast(1)).best_effort().durability_volatile();
    fault_status_pub_ =
        node_->create_publisher<lightning::msg::FaultStatus>("/localization/fault_status", health_qos);
    loc_status_pub_ =
        node_->create_publisher<lightning::msg::LocalizationStatus>("/localization/loc_status", health_qos);
    pipeline_diagnostics_pub_ = node_->create_publisher<lightning::msg::PipelineDiagnostics>(
        "/localization/pipeline_diagnostics", health_qos);
    const auto debug_qos =
        rclcpp::QoS(rclcpp::KeepLast(100)).reliable().durability_volatile();
    debug_message_pub_ = node_->create_publisher<lightning::msg::DebugMessage>(
        "/localization/debug_message", debug_qos);
    const std::weak_ptr<rclcpp::Node> weak_node = node_;
    const std::weak_ptr<rclcpp::Publisher<lightning::msg::DebugMessage>> weak_publisher =
        debug_message_pub_;
    debug_event::SetSink([weak_node, weak_publisher](const std::string& description) {
        const auto node = weak_node.lock();
        const auto publisher = weak_publisher.lock();
        if (!node || !publisher) return;
        lightning::msg::DebugMessage message;
        message.header.stamp = node->now();
        message.node_name = node->get_fully_qualified_name();
        message.message = description;
        publisher->publish(message);
    });
    path_pub_ = node_->create_publisher<nav_msgs::msg::Path>("/localization/path", path_qos);
    health_timer_ = node_->create_wall_timer(std::chrono::milliseconds(100),
                                             [this]() { PublishHealthStatus(); });
    path_timer_ = node_->create_wall_timer(std::chrono::seconds(2), [this]() { PublishPath(); });

    // 下列回调在此注册；坐标转换、发布和轨迹记录在回调被触发时执行。
    if (options_.pub_tf_) {
        auto tf_qos = tf2_ros::DynamicBroadcasterQoS();
        tf_qos.best_effort();
        tf_broadcaster_ = std::make_shared<tf2_ros::TransformBroadcaster>(node_, tf_qos);
        loc_->SetTFCallback([this](const geometry_msgs::msg::TransformStamped &pose) {
            if (!publication_gate_.MapOutputsEnabled(ToSec(pose.header.stamp))) return;
            tf_broadcaster_->sendTransform(
                sany_output::TransformTfForOutput(pose, fixed_map_transform_));
        });
    }
    // [loc-system-init-live-result]
    loc_->SetLocalizationResultCallback([this](const loc::LocalizationResult& result) {
        PublishLocalizationResult(result);
    });
    // [loc-system-init-live-result]
    loc_->SetGlobalLocalizationResultCallback([this](const loc::LocalizationResult& result) {
        CaptureGlobalLocalizationResult(result);
    });
    loc_->SetProcessedCloudCallback(
        [this](const CloudPtr& cloud, const loc::LocalizationResult& result,
               const MultiLidarFrameStats& stats, bool eligible) {
            PublishProcessedCloud(cloud, result, stats, eligible);
    });

    // 阶段 7：完成初始化并报告 Idle。
    // 初始化成功不代表定位有效；定位输入的处理开关由 SetInitPose() 打开。
    heartbeat_->SetState(functional_safety::NodeState::kIdle);
    if (incident_recorder_ && !incident_recorder_->Start()) incident_recorder_.reset();
    if (incident_recorder_) heartbeat_->SetObserver([this](const diagnostic_monitor_interfaces::msg::NodeHeartbeat& msg) {
        incident_recorder_->Capture(incident_heartbeat_,
            std::make_shared<diagnostic_monitor_interfaces::msg::NodeHeartbeat>(msg),
            sizeof(msg) + msg.active_fault.capacity() * sizeof(std::uint32_t));
    });
    LOG(INFO) << "event=LOC_NODE_READY 在线定位节点已创建";
    debug_event::Emit("在线定位节点已就绪");
    return true;
}

void LocSystem::SetInitPose(const SE3 &pose) {
    std::lock_guard<std::mutex> request_lock(initial_pose_request_mutex_);
    LOG(INFO) << "event=LOC_INITIAL_POSE 设置初始位姿 pose=" << pose.translation().transpose() << ", "
              << pose.unit_quaternion().coeffs().transpose();

    // [online-loc-start-input]
    {
        std::lock_guard<std::mutex> publish_lock(pose_publish_mutex_);
        const bool runtime_relocalization = publication_gate_.ObserveRelocalizationRequested();
        if (telemetry_) {
            telemetry_->Start();
            if (runtime_relocalization)
                telemetry_->ObserveLocalization(loc::LocalizationStatus::INITIALIZING,
                                                publication_gate_.ConsecutiveLostFrames());
        }
    }
    // Do not hold the publication mutex while draining/resetting algorithm work.
    loc_->SetExternalPose(pose.unit_quaternion(), pose.translation());
    if (lidar_frame_watchdog_) lidar_frame_watchdog_->Start();
    {
        std::lock_guard<std::mutex> publish_lock(pose_publish_mutex_);
        publication_gate_.CompleteRelocalizationRequest();
        loc_started_ = true;
    }
    // [online-loc-start-input]
}

void LocSystem::ProcessIMU(const IMUPtr &imu, diagnostics::InputTimingMonitor::Clock::time_point arrived) {
    ObserveImuInput(imu ? imu->timestamp : 0.0, arrived);
    if (imu && dtc_monitor_)
        dtc_monitor_->ObserveGyro(imu->timestamp, (loc_->GetImuToBodyRotation() * imu->angular_velocity).z());
    // [online-loc-imu-start-gate]
    if (loc_started_) {
        loc_->ProcessIMUMsg(imu);
    }
    // [online-loc-imu-start-gate]
}

void LocSystem::ProcessLidar(const sensor_msgs::msg::PointCloud2::SharedPtr &cloud,
                            diagnostics::InputTimingMonitor::Clock::time_point arrived) {
    const int lidar_id = loc_ && loc_->IsMultiLidarEnabled()
                             ? loc_->GetMultiLidarConfig().primary_lidar_id
                             : 0;
    if (!ObserveLidarInput(lidar_id, cloud ? ToSec(cloud->header.stamp) : 0.0, arrived)) return;
    if (loc_started_) {
        loc_->ProcessLidarMsg(cloud);
    }
}

void LocSystem::ProcessLidar(const sensor_msgs::msg::PointCloud2::SharedPtr& cloud, int lidar_id,
                            diagnostics::InputTimingMonitor::Clock::time_point arrived) {
    if (!ObserveLidarInput(lidar_id, cloud ? ToSec(cloud->header.stamp) : 0.0, arrived)) return;
    if (loc_started_) {
        loc_->ProcessLidarMsg(cloud, lidar_id);
    }
}

void LocSystem::ProcessLidar(const livox_ros_driver2::msg::CustomMsg::SharedPtr &cloud,
                            diagnostics::InputTimingMonitor::Clock::time_point arrived) {
    if (!ObserveLidarInput(0, cloud ? ToSec(cloud->header.stamp) : 0.0, arrived)) return;
    if (loc_started_) {
        loc_->ProcessLivoxLidarMsg(cloud);
    }
}

void LocSystem::Spin() {
    if (node_ != nullptr) {
        threading::ApplyOnlineThreadRole(threading::OnlineThreadRole::kExecutor, "loc_executor");
        spin(node_);
    }
}

void LocSystem::Finish() {
    // [loc-system-finish-resources]
    if (finished_) return;
    if (incident_recorder_) incident_recorder_->BeginShutdown();
    if (dtc_timer_) dtc_timer_->cancel();
    if (safety_timer_) safety_timer_->cancel();
    debug_event::ClearSink();
    if (loc_) loc_->Finish();
    LogInputTiming(true);
    if (heartbeat_) heartbeat_->Stop();
    ClosePublishedTrajectoryTum();
    if (incident_recorder_) incident_recorder_->Stop();
    finished_ = true;
    // [loc-system-finish-resources]
}

bool LocSystem::StartPublishedTrajectoryTum(const std::string& path) {
    if (path.empty()) {
        LOG(ERROR) << "event=LOC_TUM_PATH_INVALID 发布轨迹 TUM 路径不能为空";
        return false;
    }
    std::lock_guard<std::mutex> lock(published_tum_mutex_);
    if (published_tum_stream_.is_open()) {
        LOG(ERROR) << "event=LOC_TUM_ALREADY_OPEN 发布轨迹 TUM 记录器已经打开 path=" << published_tum_path_;
        return false;
    }

    const std::filesystem::path output_path(path);
    std::error_code directory_error;
    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path(), directory_error);
    }
    if (directory_error) {
        LOG(ERROR) << "event=LOC_TUM_DIR_FAILED 发布轨迹 TUM 目录创建失败 path=" << path
                   << ": " << directory_error.message();
        return false;
    }

    published_tum_stream_.rdbuf()->pubsetbuf(published_tum_buffer_.data(), published_tum_buffer_.size());
    published_tum_stream_.open(output_path, std::ios::out | std::ios::trunc);
    if (!published_tum_stream_) {
        LOG(ERROR) << "event=LOC_TUM_OPEN_FAILED 后轴发布轨迹 TUM 文件打开失败 path=" << path;
        return false;
    }
    published_tum_path_ = path;
    last_published_tum_timestamp_ = 0.0;
    published_tum_pose_count_ = 0;
    published_tum_write_error_count_ = 0;
    published_tum_last_flush_ = std::chrono::steady_clock::now();
    published_tum_stream_ << "# timestamp tx ty tz qx qy qz qw\n";
    published_tum_stream_.flush();
    if (!published_tum_stream_) {
        LOG(ERROR) << "event=LOC_TUM_INIT_FAILED 后轴发布轨迹 TUM 文件初始化失败 path=" << path;
        published_tum_stream_.close();
        return false;
    }
    LOG(INFO) << "event=LOC_TUM_RECORDING 正在记录后轴发布位姿 path=" << path
              << " sample_policy=all_pose_vel flush_interval_sec=1 buffer_bytes=" << published_tum_buffer_.size();
    return true;
}

void LocSystem::RecordPublishedPoseTum(const geometry_msgs::msg::PoseStamped& pose) {
    std::lock_guard<std::mutex> lock(published_tum_mutex_);
    if (!published_tum_stream_.is_open()) return;
    if (!sany_output::WriteTumPoseLine(
            published_tum_stream_, pose, last_published_tum_timestamp_, true)) {
        ++published_tum_write_error_count_;
        LOG_EVERY_N(ERROR, 100)
            << "event=LOC_TUM_WRITE_FAILED 后轴发布轨迹 TUM 位姿追加失败 stamp="
            << std::setprecision(16) << ToSec(pose.header.stamp);
        return;
    }
    // 每条已发布位姿都写入缓冲；不在高频发布线程逐帧刷新文件。
    const auto now = std::chrono::steady_clock::now();
    if (now - published_tum_last_flush_ >= std::chrono::seconds(1)) {
        published_tum_stream_.flush();
        published_tum_last_flush_ = now;
    }
    if (!published_tum_stream_) {
        ++published_tum_write_error_count_;
        LOG(ERROR) << "event=LOC_TUM_FLUSH_FAILED 后轴发布轨迹 TUM 刷新失败 path="
                   << published_tum_path_;
        return;
    }
    ++published_tum_pose_count_;
}

void LocSystem::ClosePublishedTrajectoryTum() {
    std::lock_guard<std::mutex> lock(published_tum_mutex_);
    if (!published_tum_stream_.is_open()) return;
    published_tum_stream_.flush();
    if (!published_tum_stream_) {
        ++published_tum_write_error_count_;
        LOG(ERROR) << "event=LOC_TUM_FLUSH_FAILED 后轴发布轨迹 TUM 收尾刷新失败 path=" << published_tum_path_;
    }
    published_tum_stream_.close();
    LOG(INFO) << "event=LOC_TRAJECTORY_SAVED 轨迹已保存 count=" << published_tum_pose_count_
              << " 条后轴发布位姿 path=" << published_tum_path_
              << ", write_errors=" << published_tum_write_error_count_;
}

bool LocSystem::SaveTrajectoryTum(const std::string& path) const {
    std::lock_guard<std::mutex> lock(trajectory_mutex_);
    return WriteTrajectoryTum(path, global_localization_states_, "online global localization");
}

bool LocSystem::SaveHighFrequencyTrajectoryTum(const std::string& path) const {
    std::lock_guard<std::mutex> lock(trajectory_mutex_);
    return WriteTrajectoryTum(path, localization_states_, "online high-frequency localization");
}

bool LocSystem::WriteTrajectoryTum(const std::string& path, const std::vector<NavState>& states,
                                   const char* description) const {
    std::ofstream tum(path);
    if (!tum.is_open()) {
        LOG(ERROR) << "event=LOC_TRAJECTORY_OPEN_FAILED 在线定位轨迹文件打开失败 path=" << path;
        return false;
    }
    double last_timestamp = 0.0;
    int count = 0;
    for (const auto& state : states) {
        // Entries are captured only from valid LocalizationResult callbacks.
        // PGO does not currently promote result.status_ from IDLE to GOOD, so
        // NavState::pose_is_ok_ is not a valid additional filter here.  Match
        // the offline localization exporter, which accepts valid fused output.
        if (state.timestamp_ <= 0.0 || state.timestamp_ <= last_timestamp) continue;
        const auto pose = state.GetPose();
        const auto q = pose.unit_quaternion();
        const auto p = pose.translation();
        tum << std::fixed << std::setprecision(9) << state.timestamp_ << " " << std::setprecision(12) << p.x() << " "
            << p.y() << " " << p.z() << " " << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << "\n";
        last_timestamp = state.timestamp_;
        ++count;
    }
    LOG(INFO) << "event=LOC_TRAJECTORY_SAVED 轨迹已保存 count=" << count << " " << description << " 条位姿 path=" << path;
    return count > 0;
}

void LocSystem::CaptureGlobalLocalizationResult(const loc::LocalizationResult& result) {
    if (!options_.record_global_trajectory) return;
    if (!result.valid_ || result.timestamp_ <= 0.0) return;
    const NavState state = result.ToNavState();
    std::lock_guard<std::mutex> lock(trajectory_mutex_);
    if (global_localization_states_.empty() || state.timestamp_ > global_localization_states_.back().timestamp_) {
        global_localization_states_.push_back(state);
    }
}

void LocSystem::PublishLocalizationResult(const loc::LocalizationResult& result) {
    if (incident_recorder_) incident_recorder_->PosePrepared();
    // 阶段 1：校验结果与发布条件。
    // [output-result-publication-gates]
    // 结果有效、时间递增、地图匹配可用是三个独立条件；依次通过后才构造业务消息。
    if (!result.valid_ || result.timestamp_ <= 0.0) return;
    const bool profiling_enabled = profiling::ComputeProfilingEnabled();
    profiling::Stopwatch outer_timer(profiling_enabled);
    std::unique_lock<std::mutex> publish_lock(pose_publish_mutex_);
    const auto timestamp_decision = pose_timestamp_gate_.Observe(result.timestamp_);
    if (!timestamp_decision.accepted) {
        std::ostringstream message;
        message << "event=LOC_POSE_TIMESTAMP_ROLLBACK 丢弃时间戳非递增的 pose_vel 输出 timestamp="
                << std::setprecision(16) << result.timestamp_
                << ", last_accepted=" << timestamp_decision.reference_timestamp
                << ", rollback_sec=" << timestamp_decision.lag_sec;
        LOG(WARNING) << message.str();
        debug_event::EmitThrottled("pose_vel_timestamp_rollback", message.str(),
                                   std::chrono::seconds(1));
        return;
    }
    last_localization_stamp_ = result.timestamp_;
    if (!publication_gate_.PoseOutputsEnabled(result.timestamp_)) {
        speed_smoothing_shadow_.Reset();
        return;
    }
    // [output-result-publication-gates]
    // 阶段 2：仅显式调试时缓存内部状态，并转换后轴位姿。
    if (options_.record_internal_high_frequency_trajectory) {
        const NavState state = result.ToNavState();
        std::lock_guard<std::mutex> lock(trajectory_mutex_);
        if (localization_states_.empty() || state.timestamp_ > localization_states_.back().timestamp_) {
            localization_states_.push_back(state);
        }
    }
    if (!pose_pub_ && !vehicle_pose_pub_) return;

    profiling::Stopwatch pose_transform_timer(profiling_enabled);
    // [output-rear-pose-boundary]
    const SE3 localization_rear_axle_pose = sany_output::MakeMapRearAxlePose(
        result.pose_, loc_->GetInitialLidarRotation(), primary_lidar_position_in_body_);
    const SE3 map_rear_axle_pose =
        sany_output::TransformPoseForOutput(localization_rear_axle_pose, fixed_map_transform_);
    // [output-rear-pose-boundary]
    const auto pose_transform_timing = pose_transform_timer.Stop();

    // Publish the high-frequency signed body-forward estimator velocity. CAN
    // wheel speed remains an internal filter observation, but must not replace
    // the estimator state at the downstream output boundary.
    // 阶段 3：构造业务消息。
    profiling::Stopwatch message_build_timer(profiling_enabled);
    // PGO retains IMU-frame velocity internally; rotate only at the output boundary.
    // [output-rear-speed-message]
    const double vehicle_speed =
        ImuVelocityToBody(loc_->GetImuToBodyRotation(), result.vel_b_).x() - result.RearAxleSpeedOffset();
    auto vehicle_pose = sany_output::MakeVehiclePoseMessage(
        map_rear_axle_pose, vehicle_speed, result.timestamp_, map_frame_);
    const auto pose = sany_output::MakePoseMessage(vehicle_pose);
    // [output-rear-speed-message]
    vehicle_pose.comm_header.source_id = heartbeat_->NodeId();
    vehicle_pose.comm_header.seq =
        functional_safety::HeartbeatPublisher::Advance(pose_vel_comm_seq_);
    vehicle_pose.comm_header.stamp_us =
        functional_safety::HeartbeatPublisher::UnixMicrosecondsNow();
    const auto message_build_timing = message_build_timer.Stop();

    // 阶段 4：发布位姿并更新工作状态。
    profiling::Stopwatch ros_publish_timer(profiling_enabled);
    bool position_published = false;
    if (pose_pub_) {
        pose_pub_->publish(pose);
        if (telemetry_) telemetry_->ObservePose(pose);
    }
    if (vehicle_pose_pub_) {
        // [dtc-business-work-progress]
        if (incident_recorder_) incident_recorder_->PublishBegin();
        vehicle_pose_pub_->publish(vehicle_pose);
        if (incident_recorder_) {
            incident_recorder_->PublishEnd();
            incident_recorder_->Capture(incident_pose_, std::make_shared<lightning::msg::VehiclePose>(vehicle_pose), sizeof(vehicle_pose) + 256);
        }
        position_published = true;
        last_published_pose_stamp_ = result.timestamp_;
        // DTC output freshness measures completion of the local business publish.
        if (dtc_monitor_) {
            const SE3 right = sany_output::MakeRearAxleLidarTransform(
                loc_->GetInitialLidarRotation(), primary_lidar_position_in_body_).inverse();
            const SO3 left = fixed_map_transform_.enabled ? fixed_map_transform_.T_target_localization.so3() : SO3();
            const Mat6d covariance = loc::TransformPoseCovariance(result.pose_covariance_, result.pose_, right, left);
            const SE3 raw = result.unsmoothed_pose_valid_
                ? sany_output::TransformPoseForOutput(result.unsmoothed_pose_ * right, fixed_map_transform_)
                : map_rear_axle_pose;
            dtc_monitor_->ObservePose(result.timestamp_, map_rear_axle_pose, raw, covariance,
                result.pose_covariance_valid_, result.pose_unobservable_, result.pose_epoch_);
        }
        heartbeat_->RecordWork();
        // [dtc-business-work-progress]
        const auto heartbeat_state = heartbeat_->State();
        if (heartbeat_state == functional_safety::NodeState::kNotReady ||
            heartbeat_state == functional_safety::NodeState::kIdle) {
            heartbeat_->SetState(functional_safety::NodeState::kRunning);
        }
    }
    const auto ros_publish_timing = ros_publish_timer.Stop();

    // 阶段 5：记录审计与耗时。
    profiling::Stopwatch diagnostic_io_timer(profiling_enabled);
    if (position_published && !profiling::ReduceNonessentialOverhead()) {
        if (speed_smoothing_shadow_enabled_) {
            const double shadow = speed_smoothing_shadow_.Observe(
                result.timestamp_, vehicle_speed, result.is_parking_);
            LOG(INFO) << "SPEED_SMOOTHING_SHADOW"
                      << " stamp_sec=" << vehicle_pose.header.stamp.sec
                      << " stamp_nanosec=" << vehicle_pose.header.stamp.nanosec
                      << " raw_mps=" << std::setprecision(12) << vehicle_speed
                      << " filtered_mps=" << shadow << " tau_ms=30"
                      << " parking=" << result.is_parking_
                      << " reset=" << speed_smoothing_shadow_.ResetLast()
                      << " output=raw";
        }
        // Keep this marker and key=value layout stable: it is intentionally
        // machine-readable so an exported run log can reproduce the exact
        // downstream /localization/pose_vel X/Y/speed time series without a recorded rosbag.
        LOG(INFO) << "POSE_VEL_PUBLISH topic=/localization/pose_vel"
                  << " stamp_sec=" << vehicle_pose.header.stamp.sec
                  << " stamp_nanosec=" << vehicle_pose.header.stamp.nanosec
                  << " x_m=" << std::setprecision(12) << vehicle_pose.x
                  << " y_m=" << vehicle_pose.y
                  << " z_m=" << vehicle_pose.z
                  << " speed_mps=" << vehicle_pose.speed
                  << " speed_source=estimator_body_x"
                  << " rear_axle_offset_mps=" << result.RearAxleSpeedOffset()
                  << " offset_stamp=" << result.rear_axle_speed_timestamp_;
    }
    const auto diagnostic_io_timing = diagnostic_io_timer.Stop();

    profiling::Stopwatch record_tum_timer(profiling_enabled);
    if (position_published) {
        RecordPublishedPoseTum(pose);
    }
    const auto record_tum_timing = record_tum_timer.Stop();
    const auto outer_timing = outer_timer.Stop();
    publish_lock.unlock();

    if (profiling_enabled) {
        const auto summary = output_timing_window_.Add(
            {pose_transform_timing, message_build_timing, ros_publish_timing,
             diagnostic_io_timing, record_tum_timing, outer_timing});
        if (summary) {
            LOG(INFO) << "COMPUTE_BENCH_SUMMARY module=high_frequency_output"
                      << " timestamp_s=" << std::setprecision(16) << result.timestamp_
                      << " position_published=" << (position_published ? 1 : 0)
                      << " reduce_nonessential_overhead="
                      << (profiling::ReduceNonessentialOverhead() ? 1 : 0)
                      << " " << profiling::FormatTimingWindow(*summary);
        }
    }
}

void LocSystem::PublishProcessedCloud(const CloudPtr& cloud, const loc::LocalizationResult& result,
                                      const MultiLidarFrameStats& stats, bool eligible) {
    // 阶段 1：更新匹配门控。
    {
        // Serialize invalidation with the complete group of business pose publishes.
        std::lock_guard<std::mutex> lock(pose_publish_mutex_);
        // [output-observe-lidar-match]
        const bool match_updated = publication_gate_.ObserveLidarMatch(result.lidar_loc_valid_, result.timestamp_,
            sany_output::LocalizationPublicationGate::Clock::now(), result.status_,
            result.lidar_match_stamp_, result.lidar_match_started_steady_ns_);
        // [output-observe-lidar-match]
        if (telemetry_) {
            if (match_updated)
                telemetry_->ObserveLocalization(result.status_, publication_gate_.ConsecutiveLostFrames());
            telemetry_->ObserveLocalizationStale(publication_gate_.LidarMatchStale(result.timestamp_));
        }
    }
    if (result.status_ == loc::LocalizationStatus::GOOD) {
        if (!localization_ever_good_.exchange(true)) {
            debug_event::Emit("定位结果首次可用");
        }
    }
    if (localization_ever_good_.load()) {
        debug_event::ReportState(
            "localization_following_dr",
            result.status_ == loc::LocalizationStatus::FOLLOWING_DR,
            "定位质量降级，暂停位姿输出，等待有效匹配",
            "定位已从航位推算状态恢复",
            std::chrono::milliseconds(300));
    }
    // 阶段 2：检查点云发布条件。
    debug_event::ReportState(
        "fused_cloud_suppressed", !eligible,
        "可用雷达数量不足，暂停发布融合点云",
        "融合点云发布已恢复",
        std::chrono::milliseconds(500));
    if (!eligible) {
        ++cloud_publish_suppressed_count_;
        LOG_EVERY_N(WARNING, 10)
            << "event=LOC_CLOUD_SUPPRESSED 抑制融合点云发布 present_lidars="
            << stats.present_lidar_ids.size() << ", missing_lidars="
            << stats.missing_lidar_ids.size();
        return;
    }
    if (!inv_cloud_pub_ || !map_cloud_pub_ || !cloud || cloud->empty()) return;
    const double begin_time = CloudStampSec(cloud);
    const double end_time = result.timestamp_ > 0.0 ? result.timestamp_ : begin_time;
    if (begin_time <= 0.0 || end_time <= 0.0) return;
    const bool profiling_enabled = profiling::ComputeProfilingEnabled();
    profiling::Stopwatch outer_timer(profiling_enabled);
    const SO3 initial_lidar_rotation = loc_->GetInitialLidarRotation();

    // 阶段 3：发布后轴点云。
    profiling::Stopwatch inv_message_timer(profiling_enabled);
    const auto inv_cloud_message = sany_output::MakeCloudMessage(
        cloud, begin_time, end_time,
        sany_output::MakeRearAxleLidarTransform(initial_lidar_rotation, primary_lidar_position_in_body_),
        rear_axle_frame_);
    const auto inv_message_timing = inv_message_timer.Stop();
    profiling::Stopwatch inv_publish_timer(profiling_enabled);
    inv_cloud_pub_->publish(inv_cloud_message);
    const auto inv_publish_timing = inv_publish_timer.Stop();

    // 阶段 4：按条件发布地图点云并记录耗时。
    const bool publish_map_frame = map_cloud_decimator_.Tick();
    bool map_frame_published = false;
    profiling::TimingSample map_message_timing;
    profiling::TimingSample map_publish_timing;
    if (publish_map_frame && publication_gate_.MapOutputsEnabled(end_time)) {
        // LidarLoc registers this exact cloud in the map frame and result.pose_ is T_map_lidar.
        profiling::Stopwatch map_message_timer(profiling_enabled);
        const SE3 output_lidar_pose =
            sany_output::TransformPoseForOutput(result.pose_, fixed_map_transform_);
        const auto map_cloud_message =
            sany_output::MakeCloudMessage(cloud, begin_time, end_time, output_lidar_pose, map_frame_);
        map_message_timing = map_message_timer.Stop();
        profiling::Stopwatch map_publish_timer(profiling_enabled);
        map_cloud_pub_->publish(map_cloud_message);
        map_publish_timing = map_publish_timer.Stop();
        map_frame_published = true;
    }
    const auto outer_timing = outer_timer.Stop();
    if (profiling_enabled) {
        LOG(INFO) << "COMPUTE_BENCH_FRAME module=cloud_publish"
                  << " timestamp_s=" << std::setprecision(16) << end_time
                  << " points=" << cloud->size()
                  << " present_lidars=" << stats.present_lidar_ids.size()
                  << " missing_lidars=" << stats.missing_lidar_ids.size()
                  << " map_frame_published=" << (map_frame_published ? 1 : 0)
                  << " " << profiling::FormatTimingSample("inv_message", inv_message_timing)
                  << " " << profiling::FormatTimingSample("inv_publish", inv_publish_timing)
                  << " " << profiling::FormatTimingSample("map_message", map_message_timing)
                  << " " << profiling::FormatTimingSample("map_publish", map_publish_timing)
                  << " " << profiling::FormatTimingSample("outer", outer_timing);
    }
}

void LocSystem::PublishHealthStatus() {
    LogInputTiming();
    // 阶段 1：采集输入统计。
    if (!node_ || !telemetry_ || !fault_status_pub_ || !loc_status_pub_ || !pipeline_diagnostics_pub_) return;
    std::lock_guard<std::mutex> publish_lock(pose_publish_mutex_);
    const builtin_interfaces::msg::Time stamp = node_->now();

    lightning::msg::PipelineDiagnostics diagnostics;
    diagnostics.header.stamp = stamp;
    diagnostics.header.frame_id = map_frame_;
    const auto now = std::chrono::steady_clock::now();
    double latest_input_sensor_stamp = 0.0;
    {
        std::lock_guard<std::mutex> lock(input_stats_mutex_);
        for (const auto& [lidar_id, stats] : lidar_input_stats_) {
            diagnostics.lidar_ids.push_back(lidar_id);
            diagnostics.lidar_topics.push_back(stats.topic);
            diagnostics.lidar_message_counts.push_back(stats.message_count);
            diagnostics.lidar_stale_drop_counts.push_back(stats.stale_drop_count);
            diagnostics.lidar_last_sensor_stamps.push_back(stats.last_sensor_stamp);
            diagnostics.lidar_silence_sec.push_back(
                stats.has_arrival ? std::chrono::duration<double>(now - stats.last_arrival).count() : -1.0);
            latest_input_sensor_stamp = std::max(latest_input_sensor_stamp, stats.last_sensor_stamp);
        }
        diagnostics.imu_topic = imu_input_stats_.topic;
        diagnostics.imu_message_count = imu_input_stats_.message_count;
        diagnostics.imu_last_sensor_stamp = imu_input_stats_.last_sensor_stamp;
        diagnostics.imu_silence_sec =
            imu_input_stats_.has_arrival
                ? std::chrono::duration<double>(now - imu_input_stats_.last_arrival).count()
                : -1.0;
        latest_input_sensor_stamp =
            std::max(latest_input_sensor_stamp, imu_input_stats_.last_sensor_stamp);
        diagnostics.wheel_speed_topic = wheel_speed_input_stats_.topic;
        diagnostics.wheel_speed_message_count = wheel_speed_input_stats_.message_count;
        diagnostics.wheel_speed_last_sensor_stamp = wheel_speed_input_stats_.last_sensor_stamp;
        diagnostics.wheel_speed_silence_sec =
            wheel_speed_input_stats_.has_arrival
                ? std::chrono::duration<double>(now - wheel_speed_input_stats_.last_arrival).count()
                : -1.0;
        diagnostics.wheel_speed_mps = last_wheel_speed_mps_;
        diagnostics.motor_speed_rpm = last_motor_rpm_;
        diagnostics.motor_torque_nm = last_motor_torque_;
        constexpr double kMaxWheelSpeedTimestampDeltaSec = 0.25;
        if (imu_input_stats_.last_sensor_stamp > 0.0 &&
            wheel_speed_input_stats_.last_sensor_stamp > 0.0) {
            diagnostics.wheel_speed_imu_stamp_delta_sec =
                wheel_speed_input_stats_.last_sensor_stamp -
                imu_input_stats_.last_sensor_stamp;
            diagnostics.wheel_speed_timestamp_aligned =
                std::abs(diagnostics.wheel_speed_imu_stamp_delta_sec) <=
                kMaxWheelSpeedTimestampDeltaSec;
        }
    }
    // 阶段 2：发布健康状态并更新工作状态。
    const bool lidar_match_stale = publication_gate_.LidarMatchStale(latest_input_sensor_stamp);
    telemetry_->ObserveLocalizationStale(lidar_match_stale);
    auto fault_status = telemetry_->MakeFaultStatus(stamp);
    fault_status.comm_header.source_id = heartbeat_->NodeId();
    fault_status.comm_header.seq =
        functional_safety::HeartbeatPublisher::Advance(fault_status_comm_seq_);
    fault_status.comm_header.stamp_us =
        functional_safety::HeartbeatPublisher::UnixMicrosecondsNow();
    const auto localization_status = telemetry_->MakeLocalizationStatus(stamp);
    fault_status_pub_->publish(fault_status);
    loc_status_pub_->publish(localization_status);
    if (incident_recorder_) {
        incident_recorder_->Capture(incident_fault_, std::make_shared<lightning::msg::FaultStatus>(fault_status), sizeof(fault_status) + 4096);
        incident_recorder_->Capture(incident_status_, std::make_shared<lightning::msg::LocalizationStatus>(localization_status), sizeof(localization_status) + 4096);
    }
    // Running here reflects nonzero completed work, not a fresh DTC PASS.
    // [dtc-heartbeat-state-from-work]
    if (heartbeat_->WorkSequence() == 0) {
        heartbeat_->SetState(functional_safety::NodeState::kIdle);
    } else {
        heartbeat_->SetState(functional_safety::NodeState::kRunning);
    }
    // [dtc-heartbeat-state-from-work]
    // 阶段 3：汇总队列和定位状态。
    const auto runtime = loc_->GetRuntimeStats();
    diagnostics.sensor_queue_pending = runtime.sensor_queue_pending;
    diagnostics.sensor_queue_dropped = runtime.sensor_queue_dropped;
    diagnostics.sensor_queue_processed = runtime.sensor_queue_processed;
    diagnostics.localization_queue_pending = runtime.localization_queue_pending;
    diagnostics.localization_queue_dropped = runtime.localization_queue_dropped;
    diagnostics.localization_queue_processed = runtime.localization_queue_processed;
    diagnostics.high_frequency_queue_pending = runtime.high_frequency_queue_pending;
    diagnostics.high_frequency_queue_dropped = runtime.high_frequency_queue_dropped;
    diagnostics.high_frequency_queue_processed = runtime.high_frequency_queue_processed;
    diagnostics.latest_enqueued_sensor_stamp = runtime.latest_enqueued_sensor_stamp;
    diagnostics.latest_processed_sensor_stamp = runtime.latest_processed_sensor_stamp;
    diagnostics.current_sensor_lag_sec = runtime.current_sensor_lag_sec;
    diagnostics.max_sensor_lag_sec = runtime.max_sensor_lag_sec;
    diagnostics.severe_timestamp_rollback_count = runtime.severe_timestamp_rollback_count;
    diagnostics.worst_timestamp_rollback_sec = runtime.worst_timestamp_rollback_sec;
    diagnostics.live_output_non_monotonic_drop_count =
        runtime.live_output_non_monotonic_drop_count;
    diagnostics.worst_live_output_timestamp_rollback_sec =
        runtime.worst_live_output_timestamp_rollback_sec;
    diagnostics.relocalization_attempt_count = runtime.relocalization_attempt_count;
    diagnostics.relocalization_accept_count = runtime.relocalization_accept_count;
    diagnostics.last_relocalization_candidate_found = runtime.last_relocalization_candidate_found;
    diagnostics.last_relocalization_accepted = runtime.last_relocalization_accepted;
    diagnostics.last_relocalization_candidate_id = runtime.last_relocalization_candidate_id;
    diagnostics.last_relocalization_score = runtime.last_relocalization_score;
    diagnostics.last_relocalization_search_time_ms = runtime.last_relocalization_search_time_ms;
    diagnostics.last_relocalization_reason = runtime.last_relocalization_reason;
    diagnostics.consecutive_lost_frames = publication_gate_.ConsecutiveLostFrames();
    // 阶段 4：报告地图输出门控变化。
    const bool map_outputs_enabled =
        publication_gate_.MapOutputsEnabled(latest_input_sensor_stamp);
    const bool localization_lost =
        fault_status.level == lightning::msg::FaultStatus::LEVEL_P0;
    debug_event::ReportState(
        "localization_lost", localization_lost,
        "定位已丢失，已禁用位姿和地图坐标系输出",
        "定位已从丢失状态恢复",
        std::chrono::milliseconds(0));
    debug_event::ReportState(
        "lidar_match_stale", !localization_lost && lidar_match_stale,
        "激光匹配已过期，已停止定位位姿和地图输出",
        "新鲜激光匹配已到达，定位位姿和地图输出已恢复",
        std::chrono::milliseconds(500));
    const bool map_outputs_were_enabled =
        map_outputs_enabled_last_.exchange(map_outputs_enabled);
    if (map_outputs_enabled) {
        const bool had_enabled_output = map_outputs_ever_enabled_.exchange(true);
        if (had_enabled_output && !map_outputs_were_enabled) {
            LOG(WARNING) << "event=LOC_OUTPUT_RECOVERED 新鲜有效的激光匹配已到达，恢复定位位姿和地图输出";
        }
    } else if (map_outputs_were_enabled) {
        if (lidar_match_stale) {
            LOG(WARNING) << "event=LOC_OUTPUT_STOPPED 激光匹配超时，停止定位位姿和地图输出 "
                         << "lidar_match_age_sec="
                         << publication_gate_.LidarMatchAgeSec(latest_input_sensor_stamp)
                         << ", consecutive_lost_frames="
                         << publication_gate_.ConsecutiveLostFrames();
        } else {
            LOG(ERROR) << "event=LOC_OUTPUT_DISABLED 定位结果被拒绝，禁用定位位姿和地图输出 "
                       << "consecutive_lost_frames="
                       << publication_gate_.ConsecutiveLostFrames();
        }
    }
    diagnostics.map_outputs_enabled = map_outputs_enabled;
    diagnostics.lidar_match_stale = lidar_match_stale;
    diagnostics.lidar_match_age_sec =
        publication_gate_.LidarMatchAgeSec(latest_input_sensor_stamp);
    diagnostics.last_lidar_match_stamp = publication_gate_.LastLidarMatchStamp();
    diagnostics.last_valid_lidar_match_stamp =
        publication_gate_.LastValidLidarMatchStamp();
    diagnostics.last_localization_stamp = last_localization_stamp_.load();
    diagnostics.last_posres_stamp = last_published_pose_stamp_.load();
    LOG_EVERY_N(INFO, 10) << "LOCALIZATION_HEALTH"
                        << " status=" << static_cast<int>(localization_status.status)
                        << " fault_level=" << static_cast<int>(fault_status.level)
                        << " pose_enabled=" << publication_gate_.PoseOutputsEnabled(latest_input_sensor_stamp)
                        << " lost_frames=" << diagnostics.consecutive_lost_frames
                        << " stale=" << lidar_match_stale
                        << " match_age_sec=" << diagnostics.lidar_match_age_sec
                        << " last_valid_match=" << std::setprecision(16)
                        << diagnostics.last_valid_lidar_match_stamp
                        << " last_published_pose=" << diagnostics.last_posres_stamp;
    // 阶段 5：补充降级与轮速诊断并发布。
    diagnostics.posres_non_monotonic_drop_count =
        pose_timestamp_gate_.RejectedCount();
    diagnostics.worst_posres_timestamp_rollback_sec =
        pose_timestamp_gate_.WorstRollbackSec();
    diagnostics.adaptive_lidar_load_enabled = runtime.adaptive_lidar_load_enabled;
    diagnostics.adaptive_lidar_degradation_step = runtime.adaptive_lidar_degradation_step;
    diagnostics.selected_lidar_point_stride = runtime.selected_lidar_point_stride;
    diagnostics.selected_lidar_ids = runtime.selected_lidar_ids;
    diagnostics.current_frame_lidar_ids = runtime.current_frame_lidar_ids;
    diagnostics.current_frame_missing_lidar_ids =
        runtime.current_frame_missing_lidar_ids;
    diagnostics.lidar_correction_age_sec = runtime.lidar_correction_age_sec;
    diagnostics.last_lio_processing_ms = runtime.last_lio_processing_ms;
    diagnostics.adaptive_lidar_stale_drop_count = runtime.adaptive_lidar_stale_drop_count;
    diagnostics.cloud_publish_min_lidars = runtime.cloud_publish_min_lidars;
    diagnostics.cloud_publish_eligible = runtime.cloud_publish_eligible;
    diagnostics.cloud_publish_suppressed_count =
        cloud_publish_suppressed_count_.load();
    diagnostics.wheel_speed_dr_enabled = runtime.wheel_speed_dr_enabled;
    diagnostics.wheel_speed_dr_input_count = runtime.wheel_speed_dr_stats.input_count;
    diagnostics.wheel_speed_dr_invalid_input_count =
        runtime.wheel_speed_dr_stats.invalid_input_count;
    diagnostics.wheel_speed_dr_timestamp_reject_count =
        runtime.wheel_speed_dr_stats.timestamp_reject_count;
    diagnostics.wheel_speed_dr_acceleration_reject_count =
        runtime.wheel_speed_dr_stats.acceleration_reject_count;
    diagnostics.wheel_speed_dr_lidar_filter_accepted_count =
        runtime.wheel_speed_dr_stats.lidar_filter_accepted_count;
    diagnostics.wheel_speed_dr_lidar_filter_rejected_count =
        runtime.wheel_speed_dr_stats.lidar_filter_rejected_count;
    diagnostics.wheel_speed_dr_lidar_filter_stale_count =
        runtime.wheel_speed_dr_stats.lidar_filter_stale_count;
    diagnostics.wheel_speed_dr_imu_filter_accepted_count =
        runtime.wheel_speed_dr_stats.imu_filter_accepted_count;
    diagnostics.wheel_speed_dr_imu_filter_rejected_count =
        runtime.wheel_speed_dr_stats.imu_filter_rejected_count;
    diagnostics.wheel_speed_dr_imu_filter_stale_count =
        runtime.wheel_speed_dr_stats.imu_filter_stale_count;
    diagnostics.wheel_speed_dr_last_measurement_mps =
        runtime.wheel_speed_dr_stats.last_measurement_mps;
    diagnostics.wheel_speed_dr_last_predicted_mps =
        runtime.wheel_speed_dr_stats.last_predicted_mps;
    diagnostics.wheel_speed_dr_last_innovation_mps =
        runtime.wheel_speed_dr_stats.last_innovation_mps;
    diagnostics.wheel_speed_dr_last_standard_deviation_mps =
        runtime.wheel_speed_dr_stats.last_standard_deviation_mps;
    diagnostics.wheel_speed_dr_last_nis =
        runtime.wheel_speed_dr_stats.last_normalized_innovation_squared;
    pipeline_diagnostics_pub_->publish(diagnostics);
    if (incident_recorder_) {
        lightning::diagnostics::IncidentRecorder::Pipeline p;
        p.sensor_pending = diagnostics.sensor_queue_pending;
        p.localization_pending = diagnostics.localization_queue_pending;
        p.high_frequency_pending = diagnostics.high_frequency_queue_pending;
        p.sensor_dropped = diagnostics.sensor_queue_dropped;
        p.localization_dropped = diagnostics.localization_queue_dropped;
        p.sensor_lag = diagnostics.current_sensor_lag_sec; p.match_age = diagnostics.lidar_match_age_sec;
        p.match_stale = diagnostics.lidar_match_stale; p.output_gate = diagnostics.map_outputs_enabled;
        p.lio_ms = diagnostics.last_lio_processing_ms; p.missing_lidars = diagnostics.current_frame_missing_lidar_ids.size();
        incident_recorder_->ObservePipeline(p);
        incident_recorder_->Capture(incident_pipeline_, std::make_shared<lightning::msg::PipelineDiagnostics>(diagnostics), sizeof(diagnostics) + 16384);
    }
}

void LocSystem::RegisterLidarInput(int lidar_id, const std::string& topic) {
    std::lock_guard<std::mutex> lock(input_stats_mutex_);
    lidar_input_stats_[lidar_id].topic = topic;
    // 10 Hz 点云的正常间隔约 100 ms；输入时序告警留出两帧余量。
    lidar_input_stats_[lidar_id].timing = diagnostics::InputTimingMonitor(0.2);
}

bool LocSystem::ObserveLidarInput(int lidar_id, double sensor_stamp, diagnostics::InputTimingMonitor::Clock::time_point arrived) {
    const auto timestamp_decision = lidar_input_timestamp_gate_.Observe(sensor_stamp);
    {
        std::lock_guard<std::mutex> lock(input_stats_mutex_);
        auto& stats = lidar_input_stats_[lidar_id];
        ++stats.message_count;
        stats.last_arrival = arrived;
        stats.timing.Observe(sensor_stamp, arrived);
        stats.has_arrival = true;
        if (timestamp_decision.accepted) {
            stats.last_sensor_stamp = std::max(stats.last_sensor_stamp, sensor_stamp);
        } else {
            ++stats.stale_drop_count;
        }
    }
    if (!timestamp_decision.accepted) {
        std::ostringstream message;
        message << "event=LOC_LIDAR_INPUT_DROPPED 丢弃过期或无效的雷达输入 lidar_id="
                << lidar_id << ", timestamp=" << std::setprecision(16) << sensor_stamp
                << ", newest_timestamp=" << timestamp_decision.reference_timestamp
                << ", lag_sec=" << timestamp_decision.lag_sec
                << ", configured_max_lag_sec="
                << online_lidar_input_max_timestamp_lag_sec_;
        LOG_EVERY_N(WARNING, 10) << message.str();
        debug_event::EmitThrottled(
            "stale_lidar_input_" + std::to_string(lidar_id), message.str(),
            std::chrono::seconds(1));
    }
    return timestamp_decision.accepted;
}

void LocSystem::ObserveImuInput(double sensor_stamp, diagnostics::InputTimingMonitor::Clock::time_point arrived) {
    std::lock_guard<std::mutex> lock(input_stats_mutex_);
    ++imu_input_stats_.message_count;
    imu_input_stats_.last_sensor_stamp = sensor_stamp;
    imu_input_stats_.last_arrival = arrived;
    imu_input_stats_.timing.Observe(sensor_stamp, arrived);
    imu_input_stats_.has_arrival = true;
}

void LocSystem::ObserveWheelSpeedInput(double sensor_stamp, double motor_rpm,
                                       double motor_torque, diagnostics::InputTimingMonitor::Clock::time_point arrived) {
    const double speed_mps = motor_rpm * wheel_speed_scale_mps_per_rpm_;
    if (dtc_monitor_) dtc_monitor_->ObserveSpeed(sensor_stamp, speed_mps);
    {
        std::lock_guard<std::mutex> lock(input_stats_mutex_);
        ++wheel_speed_input_stats_.message_count;
        wheel_speed_input_stats_.last_sensor_stamp = sensor_stamp;
        wheel_speed_input_stats_.last_arrival = arrived;
        wheel_speed_input_stats_.timing.Observe(sensor_stamp, arrived);
        wheel_speed_input_stats_.has_arrival = true;
        last_motor_rpm_ = motor_rpm;
        last_motor_torque_ = motor_torque;
        last_wheel_speed_mps_ = speed_mps;
    }
    if (loc_started_ && loc_) {
        loc_->ProcessWheelSpeed(sensor_stamp, speed_mps, motor_torque);
    }
}

void LocSystem::PublishPath() {
    if (!node_ || !telemetry_ || !path_pub_ || telemetry_->PathSize() == 0) return;
    if (!publication_gate_.MapOutputsEnabled(last_localization_stamp_.load())) return;
    const builtin_interfaces::msg::Time stamp = node_->now();
    path_pub_->publish(telemetry_->MakePath(stamp));
}

// [dtc-publication-gate-loss]
loc::dtc::Snapshot LocSystem::EvaluateLocalizationDtc() const {
    const auto now = loc::dtc::Clock::now();
    return dtc_monitor_->Evaluate(publication_gate_.TrackingState(last_localization_stamp_.load(), now), now);
}
// [dtc-publication-gate-loss]
void LocSystem::UpdateDtcHeartbeat() {
    if (!dtc_monitor_ || !heartbeat_) return;
    const auto snapshot = EvaluateLocalizationDtc();
    std::lock_guard<std::mutex> lock(dtc_publish_mutex_);
    functional_safety::UpdateLocalizationFaults(*heartbeat_, snapshot);
    for (size_t i = 0; i < 2; ++i) {
        const auto& check = snapshot.observations[i];
        last_dtc_conditions_[i] = check.condition;
    }
    LOG_EVERY_N(INFO, 10) << "LOCALIZATION_DTC consistency=" << int(snapshot.observations[0].condition)
        << " uncertainty=" << int(snapshot.observations[1].condition) << " severe=" << snapshot.severe
        << " error_x=" << snapshot.error_longitudinal_m << " error_y=" << snapshot.error_lateral_m
        << " error_yaw_rad=" << snapshot.error_yaw_rad << " major_sigma=" << snapshot.major_sigma_m
        << " lateral_sigma=" << snapshot.lateral_sigma_m << " yaw_sigma_rad=" << snapshot.yaw_sigma_rad
        << " covariance_valid=" << snapshot.covariance_valid
        << " pose_unobservable=" << snapshot.pose_unobservable
        << " can_stamp_lag_sec=" << snapshot.speed_stamp_lag_sec
        << " imu_stamp_lag_sec=" << snapshot.gyro_stamp_lag_sec
        << " confirmed_failure=" << snapshot.confirmed_failure
        << " unknown_tolerated=" << snapshot.consistency_unknown_tolerated
        << " optional_can_unavailable=" << snapshot.optional_can_unavailable
        << " translation_condition=" << int(snapshot.translation_condition)
        << " yaw_condition=" << int(snapshot.yaw_condition)
        << " yaw_pass_age_sec=" << snapshot.yaw_pass_age_sec
        << " consistency_pass_age_sec=" << snapshot.consistency_pass_age_sec
        << " output_age_sec=" << snapshot.output_age_sec
        << " can_arrival_age_sec=" << snapshot.speed_arrival_age_sec
        << " imu_arrival_age_sec=" << snapshot.gyro_arrival_age_sec
        << " output_stale=" << snapshot.output_stale
        << " localization_lost=" << snapshot.localization_lost;
}
void LocSystem::PublishSafety() {
    // 阶段 1：评估定位与雷达状态。
    if (!dtc_monitor_ || !safety_pub_ || !radar_dtc_state_) return;
    const auto localization = EvaluateLocalizationDtc();
    const auto radar = radar_dtc_state_->Evaluate();
    const auto frames = lidar_frame_watchdog_->Evaluate();
    geosun_msgs::msg::Safety msg;
    msg.header.stamp = node_->now();
    const int strict = loc::dtc::SafetyValue(localization, radar.all_failed, radar.fault, radar.all_pass);
    const int bounded = loc::dtc::SafetyValue(localization, radar.all_failed,
        radar.bounded_fault, radar.bounded_all_pass);
    // 阶段 2：按配置发布 Safety。
    // [dtc-safety-output-switches]
    const auto decision = loc::dtc::EvaluateFrameSafetyDecision(localization, frames, safety_fault_reporting_enabled_);
    msg.runstate = decision.published;
    if (incident_recorder_) {
        diagnostics::IncidentRecorder::Safety s;
        s.raw = decision.raw; s.selected = decision.selected; s.published = decision.published; s.strict = strict;
        s.consistency = localization.observations[0].condition; s.uncertainty = localization.observations[1].condition;
        s.covariance_valid = localization.covariance_valid; s.unobservable = localization.pose_unobservable;
        s.major_sigma = localization.major_sigma_m; s.lateral_sigma = localization.lateral_sigma_m; s.yaw_sigma = localization.yaw_sigma_rad;
        s.translation_condition = localization.translation_condition; s.yaw_condition = localization.yaw_condition;
        s.confirmed_failure = localization.confirmed_failure;
        s.unknown_tolerated = localization.consistency_unknown_tolerated;
        s.optional_can_unavailable = localization.optional_can_unavailable;
        s.output_stale = localization.output_stale; s.localization_lost = localization.localization_lost;
        s.radar_all_failed = radar.all_failed; s.non_sync_fault = radar.fault_without_time_sync;
        s.non_sync_all_pass = radar.all_pass_without_time_sync; s.bounded_policy_enabled = radar_bounded_sync_warning_enabled_;
        s.time_sync_fault_mask = radar.time_sync_fault_mask; s.bounded_sync_warning_mask = radar.bounded_sync_warning_mask;
        s.error_x = localization.error_longitudinal_m; s.error_y = localization.error_lateral_m;
        s.error_yaw_rad = localization.error_yaw_rad; s.output_age_sec = localization.output_age_sec;
        s.imu_arrival_age_sec = localization.gyro_arrival_age_sec; s.imu_stamp_lag_sec = localization.gyro_stamp_lag_sec;
        s.can_arrival_age_sec = localization.speed_arrival_age_sec; s.can_stamp_lag_sec = localization.speed_stamp_lag_sec;
        s.yaw_pass_age_sec = localization.yaw_pass_age_sec;
        s.policy_version = loc::dtc::kSafetyPolicyVersion;
        s.startup_pending = localization.startup_pending; s.startup_expired = localization.startup_expired;
        s.startup_age_sec = localization.startup_age_sec; s.valid_match_age_sec = localization.valid_match_age_sec;
        s.output_gate_open = localization.output_gate_open; s.match_stale = localization.match_stale;
        s.latest_match_valid = localization.latest_match_valid; s.output_continuous = localization.output_continuous;
        s.localization_healthy_for_degraded = loc::dtc::HealthyLocalizationForDegradedRadar(localization);
        s.partial_radar_suppressed = decision.partial_radar_suppressed;
        s.lidar_frame_started = frames.started; s.lidar_frame_seen = frames.seen;
        s.lidar_frame_warning = frames.warning; s.lidar_frame_fault = frames.fault;
        s.lidar_frame_age_sec = frames.age_sec; s.lidar_frame_recovery_count = frames.recovery_count;
        incident_recorder_->ObserveSafety(s);
        incident_recorder_->Capture(incident_safety_, std::make_shared<geosun_msgs::msg::Safety>(msg), sizeof(msg) + 256);
    }
    safety_pub_->publish(msg);
    // [dtc-safety-output-switches]
    // Always retain the pre-override decision, even with both reporting and
    // optional time-sync shadow logging disabled. Log every observed cause
    // transition and a one-second summary without high-rate file I/O.
    const std::array<int, 27> diagnostic_state{{
        decision.raw, decision.selected, decision.published, strict,
        localization.observations[0].condition, localization.observations[1].condition,
        localization.confirmed_failure, localization.output_stale, localization.localization_lost,
        localization.translation_condition, localization.yaw_condition,
        radar.all_failed, radar.fault_without_time_sync, radar.all_pass_without_time_sync,
        radar.time_sync_fault_mask, radar.bounded_sync_warning_mask,
        localization.covariance_valid, localization.pose_unobservable,
        localization.startup_pending, localization.startup_expired, localization.output_gate_open,
        localization.match_stale, localization.latest_match_valid, localization.output_continuous,
        decision.partial_radar_suppressed, frames.warning, frames.fault}};
    const bool diagnostic_changed = !safety_diagnostic_seen_ || diagnostic_state != last_safety_diagnostic_state_;
    const bool diagnostic_periodic = ++safety_diagnostic_ticks_ % 50 == 0;
    if (diagnostic_changed || diagnostic_periodic) {
        if (diagnostic_changed && frames.warning) {
            LOG(WARNING) << "event=LOC_RADAR_FRAME_SUPPLY 雷达输入组帧连续不可用或正在确认恢复"
                         << " age_sec=" << frames.age_sec << " fault=" << frames.fault
                         << " recovery_count=" << frames.recovery_count;
        }
        const double sensor_stamp = last_localization_stamp_.load();
        LOG(INFO) << "SAFETY_DIAGNOSTIC sequence=" << ++safety_diagnostic_sequence_
                  << " policy_version=" << loc::dtc::kSafetyPolicyVersion
                  << " radar_frame_observation_point=raw_input_grouping"
                  << " startup_pending=" << localization.startup_pending
                  << " startup_expired=" << localization.startup_expired
                  << " startup_age_sec=" << localization.startup_age_sec
                  << " confirmed_tracking_lost=" << localization.localization_lost
                  << " output_gate_open=" << localization.output_gate_open
                  << " latest_match_valid=" << localization.latest_match_valid
                  << " valid_match_age_sec=" << localization.valid_match_age_sec
                  << " output_continuous=" << localization.output_continuous
                  << " localization_healthy_for_degraded=" << loc::dtc::HealthyLocalizationForDegradedRadar(localization)
                  << " partial_radar_suppressed=" << decision.partial_radar_suppressed
                  << " partial_radar_suppression_reason=" << (decision.partial_radar_suppressed ? "healthy_localization" : "none")
                  << " lidar_frame_started=" << frames.started
                  << " lidar_frame_seen=" << frames.seen
                  << " lidar_frame_warning=" << frames.warning
                  << " lidar_frame_fault=" << frames.fault
                  << " lidar_frame_age_sec=" << frames.age_sec
                  << " lidar_frame_recovery_count=" << frames.recovery_count
                  << " lidar_frame_accepted=" << frames.accepted
                  << " lidar_frame_rejected=" << frames.rejected
                  << " event=" << (diagnostic_changed ? "transition" : "periodic")
                  << " reporting_enabled=" << safety_fault_reporting_enabled_
                  << " raw_runstate=" << int(decision.raw)
                  << " selected_runstate=" << int(decision.selected)
                  << " published_runstate=" << int(decision.published)
                  << " suppressed=" << (decision.raw != decision.published)
                  << " strict_runstate=" << strict << " bounded_runstate=" << bounded
                  << " bounded_policy_enabled=" << radar_bounded_sync_warning_enabled_
                  << " radar_fault_mask_enabled=" << radar_fault_mask_enabled_
                  << " radar_all_failed=" << radar.all_failed
                  << " time_sync_fault_mask=" << int(radar.time_sync_fault_mask)
                  << " bounded_sync_warning_mask=" << int(radar.bounded_sync_warning_mask)
                  << " non_sync_fault=" << radar.fault_without_time_sync
                  << " non_sync_all_pass=" << radar.all_pass_without_time_sync
                  << " localization_consistency=" << int(localization.observations[0].condition)
                  << " localization_uncertainty=" << int(localization.observations[1].condition)
                  << " translation_condition=" << int(localization.translation_condition)
                  << " yaw_condition=" << int(localization.yaw_condition)
                  << " severe=" << localization.severe
                  << " confirmed_failure=" << localization.confirmed_failure
                  << " unknown_tolerated=" << localization.consistency_unknown_tolerated
                  << " optional_can_unavailable=" << localization.optional_can_unavailable
                  << " error_x=" << localization.error_longitudinal_m
                  << " error_y=" << localization.error_lateral_m
                  << " error_yaw_rad=" << localization.error_yaw_rad
                  << " major_sigma=" << localization.major_sigma_m
                  << " lateral_sigma=" << localization.lateral_sigma_m
                  << " yaw_sigma_rad=" << localization.yaw_sigma_rad
                  << " covariance_valid=" << localization.covariance_valid
                  << " pose_unobservable=" << localization.pose_unobservable
                  << " output_stale=" << localization.output_stale
                  << " localization_lost=" << localization.localization_lost
                  << " output_age_sec=" << localization.output_age_sec
                  << " can_arrival_age_sec=" << localization.speed_arrival_age_sec
                  << " imu_arrival_age_sec=" << localization.gyro_arrival_age_sec
                  << " can_stamp_lag_sec=" << localization.speed_stamp_lag_sec
                  << " imu_stamp_lag_sec=" << localization.gyro_stamp_lag_sec
                  << " yaw_pass_age_sec=" << localization.yaw_pass_age_sec
                  << " consistency_pass_age_sec=" << localization.consistency_pass_age_sec
                  << " match_stale=" << localization.match_stale
                  << " lost_frames=" << publication_gate_.ConsecutiveLostFrames()
                  << " match_age_sec=" << publication_gate_.LidarMatchAgeSec(sensor_stamp)
                  << std::setprecision(16)
                  << " last_valid_match=" << publication_gate_.LastValidLidarMatchStamp()
                  << " last_published_pose=" << last_published_pose_stamp_.load()
                  << " steady_sec=" << std::chrono::duration<double>(
                         loc::dtc::Clock::now().time_since_epoch()).count();
    }
    last_safety_diagnostic_state_ = diagnostic_state;
    safety_diagnostic_seen_ = true;
    // 阶段 3：记录影子评估。
    if (radar_time_sync_shadow_enabled_ || radar_fault_mask_enabled_) {
        const int shadow = loc::dtc::SafetyValue(localization, radar.all_failed,
            radar.fault_without_time_sync, radar.all_pass_without_time_sync);
        const bool changed = strict != last_safety_strict_ || shadow != last_safety_shadow_ ||
                             radar.time_sync_fault_mask != last_safety_sync_mask_ ||
                             msg.runstate != last_safety_published_ ||
                             radar.bounded_sync_warning_mask != last_safety_bounded_mask_;
        // Retain a short tolerated UNKNOWN even when Safety remains unchanged.
        if (changed || localization.consistency_unknown_tolerated || ++safety_shadow_ticks_ % 50 == 0) {
            LOG(INFO) << "SAFETY_TIME_SYNC_SHADOW published_runstate=" << int(msg.runstate)
                      << " strict_runstate=" << strict << " shadow_runstate=" << shadow
                      << " bounded_runstate=" << bounded
                      << " bounded_policy_enabled=" << radar_bounded_sync_warning_enabled_
                      << " radar_fault_mask_enabled=" << radar_fault_mask_enabled_
                      << " bounded_sync_warning_mask=" << int(radar.bounded_sync_warning_mask)
                      << " time_sync_fault_mask=" << int(radar.time_sync_fault_mask)
                      << " non_sync_fault=" << radar.fault_without_time_sync
                      << " non_sync_all_pass=" << radar.all_pass_without_time_sync
                      << " localization_consistency=" << int(localization.observations[0].condition)
                      << " localization_uncertainty=" << int(localization.observations[1].condition)
                      << " severe=" << localization.severe
                      << " confirmed_failure=" << localization.confirmed_failure
                      << " unknown_tolerated=" << localization.consistency_unknown_tolerated
                      << " optional_can_unavailable=" << localization.optional_can_unavailable
                      << " translation_condition=" << int(localization.translation_condition)
                      << " yaw_condition=" << int(localization.yaw_condition)
                      << " yaw_pass_age_sec=" << localization.yaw_pass_age_sec
                      << " consistency_pass_age_sec=" << localization.consistency_pass_age_sec
                      << " output_age_sec=" << localization.output_age_sec
                      << " can_arrival_age_sec=" << localization.speed_arrival_age_sec
                      << " imu_arrival_age_sec=" << localization.gyro_arrival_age_sec
                      << " can_stamp_lag_sec=" << localization.speed_stamp_lag_sec
                      << " imu_stamp_lag_sec=" << localization.gyro_stamp_lag_sec
                      << " output_stale=" << localization.output_stale
                      << " localization_lost=" << localization.localization_lost;
        }
        last_safety_strict_ = strict;
        last_safety_shadow_ = shadow;
        last_safety_sync_mask_ = radar.time_sync_fault_mask;
        last_safety_published_ = msg.runstate;
        last_safety_bounded_mask_ = radar.bounded_sync_warning_mask;
    }
}

void LocSystem::LogInputTiming(bool force) {
    using Monitor = diagnostics::InputTimingMonitor;
    const auto now = Monitor::Clock::now();
    if (!force && now - input_timing_last_log_ < std::chrono::seconds(1)) return;
    const bool summary = force || now - input_timing_last_summary_ >= std::chrono::seconds(10);
    input_timing_last_log_ = now;
    if (summary) input_timing_last_summary_ = now;
    struct Entry { std::string stream, topic; Monitor::Record record; };
    std::vector<Entry> entries;
    {
        std::lock_guard<std::mutex> lock(input_stats_mutex_);
        auto collect = [&](const char* stream, InputTopicStats& stats) {
            if (stats.topic.empty()) return;
            stats.timing.Poll(now);
            Monitor::Record record;
            while (stats.timing.Pop(record)) entries.push_back({stream, stats.topic, record});
            if (summary) entries.push_back({stream, stats.topic, stats.timing.Summary()});
        };
        for (auto& item : lidar_input_stats_) collect("lidar", item.second);
        collect("imu", imu_input_stats_);
        if (wheel_speed_observation_enabled_) collect("can", wheel_speed_input_stats_);
    }
    const double emitted = std::chrono::duration<double>(now.time_since_epoch()).count();
    for (const auto& entry : entries) {
        const auto& r = entry.record;
        // 英文字段为稳定机器接口；正文保持中文，输出移出传感器回调和输入锁。
        LOG(INFO) << "INPUT_TIMING 输入时序观测 schema=1 stream=" << entry.stream << " topic=" << entry.topic
                  << " kind=" << Monitor::Name(r.kind) << " sample_seq=" << r.sample << " event_seq=" << r.event
                  << std::setprecision(17) << " source_stamp_sec=" << r.stamp
                  << " previous_source_stamp_sec=" << r.previous_stamp << " source_gap_sec=" << r.source_gap
                  << " arrival_gap_sec=" << r.arrival_gap << " arrival_steady_sec=" << r.arrival_steady
                  << " emitted_steady_sec=" << emitted << " threshold_sec=" << r.threshold
                  << " events_lost=" << r.lost << " max_source_gap_sec=" << r.max_source_gap
                  << " max_arrival_gap_sec=" << r.max_arrival_gap << " source_gaps=" << r.source_gaps
                  << " arrival_gaps=" << r.arrival_gaps << " duplicates=" << r.duplicates
                  << " rollbacks=" << r.rollbacks << " invalid=" << r.invalid;
    }
}

}  // namespace lightning
