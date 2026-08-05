#include <algorithm>
#include <array>
#include <execution>
#include <filesystem>
#include <limits>
#include <stdexcept>

#include <pcl/common/transforms.h>
#include <pcl/filters/passthrough.h>
#include <pcl/io/pcd_io.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/pcl_base.h>
#include <pcl/registration/ndt.h>
#include <yaml-cpp/yaml.h>

#include "pclomp/ndt_omp_impl.hpp"
#include "pclomp/voxel_grid_covariance_omp_impl.hpp"

#include "core/localization/lidar_loc/lidar_loc.h"

#include <opencv2/highgui.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

#include "glog/logging.h"
#include "core/localization/btc_relocalizer.h"
#include "core/localization/solid_relocalizer.h"
#include "core/maps/map_frame.h"
#include "io/file_io.h"
#include "io/yaml_io.h"
#include "ui/pangolin_window.h"
#include "utils/timer.h"

namespace lightning::loc {
namespace {

SE3 ReadLidarToImu(const YAML::Node& root) {
    Mat3d rotation = Mat3d::Identity();
    Vec3d translation = Vec3d::Zero();
    const YAML::Node fasterlio = root["fasterlio"];
    if (fasterlio && fasterlio["extrinsic_R"]) {
        const auto values = fasterlio["extrinsic_R"].as<std::vector<double>>();
        if (values.size() != 9) throw std::runtime_error("fasterlio.extrinsic_R must have 9 values");
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < 3; ++col) rotation(row, col) = values[row * 3 + col];
        }
    }
    if (fasterlio && fasterlio["extrinsic_T"]) {
        const auto values = fasterlio["extrinsic_T"].as<std::vector<double>>();
        if (values.size() != 3) throw std::runtime_error("fasterlio.extrinsic_T must have 3 values");
        translation = Vec3d(values[0], values[1], values[2]);
    }
    return SE3(Quatd(rotation).normalized(), translation);
}

}  // namespace

LidarLoc::LidarLoc(LidarLoc::Options options) : options_(options) {
    pcl_ndt_.reset(new NDTType());
    pcl_ndt_->setResolution(1.0);
    pcl_ndt_->setNeighborhoodSearchMethod(pclomp::DIRECT7);
    pcl_ndt_->setOulierRatio(0.45);
    pcl_ndt_->setStepSize(0.1);
    pcl_ndt_->setTransformationEpsilon(0.01);
    pcl_ndt_->setMaximumIterations(20);
    pcl_ndt_->setNumThreads(4);

    pcl_ndt_rough_.reset(new NDTType());
    pcl_ndt_rough_->setResolution(5.0);
    pcl_ndt_rough_->setNeighborhoodSearchMethod(pclomp::DIRECT7);
    pcl_ndt_rough_->setStepSize(0.1);
    pcl_ndt_rough_->setMaximumIterations(4);
    pcl_ndt_rough_->setNumThreads(4);

    pcl_icp_.reset(new ICPType());
    pcl_icp_->setMaximumIterations(4);
    pcl_icp_->setTransformationEpsilon(0.01);

    LOG(INFO) << "match name is NDT_OMP"
              << ", MaximumIterations is: " << pcl_ndt_->getMaximumIterations();
}

LidarLoc::~LidarLoc() {
    if (update_map_thread_.joinable()) {
        update_map_quit_ = true;
        update_map_thread_.join();
    }

    recover_pose_out_.close();
}

bool LidarLoc::Init(const std::string& config_path) {
    YAML_IO yaml(config_path);
    const YAML::Node root = YAML::LoadFile(config_path);
    map_frame::ExportOptions export_options;
    std::string map_frame_error;
    if (!map_frame::ReadExportOptions(root, export_options, map_frame_error)) {
        LOG(ERROR) << map_frame_error;
        return false;
    }
    if (export_options.normalize_start_ground_z) {
        map_frame::Metadata metadata;
        if (!map_frame::LoadMetadata(
                options_.map_option_.map_path_, metadata, map_frame_error)) {
            LOG(ERROR) << map_frame_error;
            return false;
        }
    }
    options_.map_option_.enable_dynamic_polygon_ = yaml.GetValue<bool>("maps", "with_dyn_area");
    options_.map_option_.max_pts_in_dyn_chunk_ = yaml.GetValue<int>("maps", "max_pts_dyn_chunk");
    options_.map_option_.load_map_size_ = yaml.GetValue<int>("maps", "load_map_size");
    options_.map_option_.unload_map_size_ = yaml.GetValue<int>("maps", "unload_map_size");

    options_.update_kf_dis_ = yaml.GetValue<double>("lidar_loc", "update_kf_dis");
    options_.update_lidar_loc_score_ = yaml.GetValue<double>("lidar_loc", "update_lidar_loc_score");
    options_.min_init_confidence_ = yaml.GetValue<float>("lidar_loc", "min_init_confidence");
    if (root["lidar_loc"] && root["lidar_loc"]["min_tracking_confidence"]) {
        options_.min_tracking_confidence_ = root["lidar_loc"]["min_tracking_confidence"].as<float>();
    } else {
        options_.min_tracking_confidence_ = options_.min_init_confidence_;
    }
    if (root["relocalization"] && root["relocalization"]["lost_frame_threshold"]) {
        options_.relocalization_lost_frame_threshold_ =
            std::max(1, root["relocalization"]["lost_frame_threshold"].as<int>());
    }
    const YAML::Node relocalization = root["relocalization"];
    if (relocalization && relocalization["enable_map_consistency"]) {
        options_.enable_relocalization_map_consistency_ =
            relocalization["enable_map_consistency"].as<bool>();
    }
    if (relocalization && relocalization["map_bounds_margin"]) {
        options_.relocalization_bounds_margin_ = relocalization["map_bounds_margin"].as<double>();
    }
    if (relocalization && relocalization["map_nearest_neighbor_distance"]) {
        options_.relocalization_nearest_neighbor_distance_ =
            relocalization["map_nearest_neighbor_distance"].as<double>();
    }
    if (relocalization && relocalization["map_min_inside_xy_ratio"]) {
        options_.relocalization_min_inside_xy_ratio_ =
            relocalization["map_min_inside_xy_ratio"].as<double>();
    }
    if (relocalization && relocalization["map_min_overlap_ratio"]) {
        options_.relocalization_min_overlap_ratio_ =
            relocalization["map_min_overlap_ratio"].as<double>();
    }
    if (relocalization && relocalization["map_min_gravity_alignment_cos"]) {
        options_.relocalization_min_gravity_alignment_cos_ =
            relocalization["map_min_gravity_alignment_cos"].as<double>();
    }
    if (relocalization && relocalization["map_consistency_max_points"]) {
        options_.relocalization_map_consistency_max_points_ =
            relocalization["map_consistency_max_points"].as<int>();
    }
    if (relocalization && relocalization["validation_workers"]) {
        options_.relocalization_validation_workers_ =
            relocalization["validation_workers"].as<int>();
    }
    if (relocalization && relocalization["confirmation_count"]) {
        options_.relocalization_confirmation_count_ =
            relocalization["confirmation_count"].as<int>();
    }
    if (relocalization && relocalization["confirmation_max_translation"]) {
        options_.relocalization_confirmation_max_translation_ =
            relocalization["confirmation_max_translation"].as<double>();
    }
    if (relocalization && relocalization["confirmation_max_rotation_deg"]) {
        options_.relocalization_confirmation_max_rotation_deg_ =
            relocalization["confirmation_max_rotation_deg"].as<double>();
    }
    if (relocalization && relocalization["confirmation_max_interval"]) {
        options_.relocalization_confirmation_max_interval_ =
            relocalization["confirmation_max_interval"].as<double>();
    }
    if (options_.relocalization_bounds_margin_ < 0.0 ||
        options_.relocalization_nearest_neighbor_distance_ <= 0.0 ||
        options_.relocalization_min_inside_xy_ratio_ < 0.0 ||
        options_.relocalization_min_inside_xy_ratio_ > 1.0 ||
        options_.relocalization_min_overlap_ratio_ < 0.0 ||
        options_.relocalization_min_overlap_ratio_ > 1.0 ||
        options_.relocalization_min_gravity_alignment_cos_ < -1.0 ||
        options_.relocalization_min_gravity_alignment_cos_ > 1.0 ||
        options_.relocalization_map_consistency_max_points_ <= 0 ||
        options_.relocalization_validation_workers_ <= 0 ||
        options_.relocalization_confirmation_count_ <= 0 ||
        options_.relocalization_confirmation_max_translation_ < 0.0 ||
        options_.relocalization_confirmation_max_rotation_deg_ < 0.0 ||
        options_.relocalization_confirmation_max_interval_ <= 0.0) {
        LOG(ERROR) << "invalid relocalization map-consistency configuration";
        return false;
    }

    // options_.filter_z_min_ = yaml.GetValue<double>("lidar_loc", "filter_z_min");
    // options_.filter_z_max_ = yaml.GetValue<double>("lidar_loc", "filter_z_max");
    // options_.filter_intensity_min_ = yaml.GetValue<double>("lidar_loc", "filter_intensity_min");
    // options_.filter_intensity_max_ = yaml.GetValue<double>("lidar_loc", "filter_intensity_max");
    options_.lidar_loc_odom_th_ = yaml.GetValue<double>("lidar_loc", "lidar_loc_odom_th");

    options_.init_with_fp_ = yaml.GetValue<bool>("lidar_loc", "init_with_fp");
    options_.enable_parking_static_ = yaml.GetValue<bool>("lidar_loc", "enable_parking_static");
    options_.enable_icp_adjust_ = yaml.GetValue<bool>("lidar_loc", "enable_icp_adjust");
    options_.with_height_ = yaml.GetValue<bool>("loop_closing", "with_height");
    options_.try_self_extrap_ = yaml.GetValue<bool>("lidar_loc", "try_self_extrap");

    lidar_loc::grid_search_angle_step = yaml.GetValue<double>("lidar_loc", "grid_search_angle_step");
    lidar_loc::grid_search_angle_range = yaml.GetValue<double>("lidar_loc", "grid_search_angle_range");

    LOG(INFO) << "min init confidence: " << options_.min_init_confidence_
              << ", min tracking confidence: " << options_.min_tracking_confidence_;

    std::string map_policy = yaml.GetValue<std::string>("maps", "dyn_cloud_policy");
    if (map_policy == "short") {
        options_.map_option_.policy_ = TiledMap::DynamicCloudPolicy::SHORT;
    } else if (map_policy == "long") {
        options_.map_option_.policy_ = TiledMap::DynamicCloudPolicy::LONG;
    } else if (map_policy == "persistent") {
        options_.map_option_.policy_ = TiledMap::DynamicCloudPolicy::PERSISTENT;
    }

    options_.map_option_.delete_when_unload_ = yaml.GetValue<bool>("maps", "delete_when_unload");
    options_.map_option_.load_dyn_cloud_ = yaml.GetValue<bool>("maps", "load_dyn_cloud");
    options_.map_option_.save_dyn_when_quit_ = yaml.GetValue<bool>("maps", "save_dyn_when_quit");
    options_.map_option_.save_dyn_when_unload_ = yaml.GetValue<bool>("maps", "save_dyn_when_unload");

    map_ = std::make_shared<TiledMap>(options_.map_option_);
    if (!map_->LoadMapIndex()) return false;
    if (options_.enable_relocalization_map_consistency_ && !BuildRelocalizationMapCache()) {
        return false;
    }

    auto fps = map_->GetAllFP();
    if (!fps.empty()) {
        map_->LoadOnPose(fps.front().pose_);
        /// 更新一次地图，保证有初始数据
        UpdateGlobalMap();
    }

    relocalization_backend_name_ = relocalization && relocalization["backend"]
                                        ? relocalization["backend"].as<std::string>()
                                        : "btc";
    if (relocalization_backend_name_ == "btc") {
        global_relocalizer_ = std::make_unique<BtcRelocalizer>();
    } else if (relocalization_backend_name_ == "solid") {
        global_relocalizer_ = std::make_unique<SolidRelocalizer>();
    } else {
        LOG(ERROR) << "unsupported global relocalization backend: "
                   << relocalization_backend_name_;
        return false;
    }
    if (!global_relocalizer_->Init(
            config_path, options_.map_option_.map_path_, ReadLidarToImu(root))) {
        if (export_options.normalize_start_ground_z) {
            LOG(ERROR) << "normalized map requires a consistent "
                       << relocalization_backend_name_ << " relocalization database";
            return false;
        }
        LOG(WARNING) << relocalization_backend_name_
                     << " relocalization is unavailable; NDT localization remains enabled";
    }

    /// load recover pose if exist
    if (PathExists(options_.recover_pose_path_)) {
        std::ifstream fin(options_.recover_pose_path_);
        double data[7] = {0, 0, 0, 0, 0, 0, 1};
        for (int i = 0; i < 7; ++i) {
            fin >> data[i];
        }

        SE3 pose(Quatd(data[6], data[3], data[4], data[5]), Vec3d(data[0], data[1], data[2]));
        FunctionalPoint fp_recover;
        fp_recover.name_ = "recover";
        fp_recover.pose_ = pose;
        map_->AddFP(fp_recover);
    }

    update_map_thread_ = std::thread([this]() { LidarLoc::UpdateMapThread(); });

    return true;
}

bool LidarLoc::ProcessCloud(CloudPtr cloud_input) {
    assert(cloud_input != nullptr);

    if (cloud_input->empty() || cloud_input->size() < 50) {
        LOG(WARNING) << "loc input is empty or invalid, sz: " << cloud_input->size();
        return false;
    }

    // CloudPtr cloud(new PointCloudType);
    // pcl::VoxelGrid<PointType> voxel;

    // float sz = 0.1;
    // voxel.setLeafSize(sz, sz, sz);
    // voxel.setInputCloud(cloud_input);
    // voxel.filter(*cloud);

    current_scan_ = cloud_input;

    Align(cloud_input);
    return true;
}

NavState LidarLoc::GetState() {
    UL lock_res(result_mutex_);
    NavState ns;
    ns.SetPose(current_abs_pose_);
    ns.timestamp_ = current_timestamp_;
    ns.confidence_ = current_score_;

    UL lock(lo_pose_mutex_);
    if (!lo_pose_queue_.empty()) {
        auto s = lo_pose_queue_.back();
        ns.SetVel(s.GetVel());
    }

    return ns;
}

LidarLoc::MatchStats LidarLoc::GetLastMatchStats() const {
    MatchStats stats = last_match_stats_;
    stats.active_map_chunks = map_ ? map_->NumActiveChunks() : 0;
    return stats;
}

bool LidarLoc::ProcessDR(const NavState& state) {
    // 未初始化成功的数据不接收
    if (!state.pose_is_ok_) {
        return false;
    }

    // DR数据check
    UL lock(dr_pose_mutex_);
    if (!dr_pose_queue_.empty()) {
        const double last_stamp = dr_pose_queue_.back().timestamp_;
        if (state.timestamp_ < last_stamp) {
            return false;
        }
    }

    dr_pose_queue_.emplace_back(state);
    while (dr_pose_queue_.size() >= 1000) {
        dr_pose_queue_.pop_front();
    }

    return true;
}

bool LidarLoc::ProcessLO(const NavState& state) {
    /// 理论上相对定位是按时间顺序到达的
    UL lock(lo_pose_mutex_);
    if (!lo_pose_queue_.empty()) {
        const double last_stamp = lo_pose_queue_.back().timestamp_;
        if (state.timestamp_ < last_stamp) {
            LOG(WARNING) << "当前相对定位的结果的时间戳应当比上一个时间戳数值大，实际相减得"
                         << state.timestamp_ - last_stamp;
            return false;
        }
    }

    lo_pose_queue_.emplace_back(state);

    while (lo_pose_queue_.size() >= 50) {
        lo_pose_queue_.pop_front();
    }

    if (state.lidar_odom_reliable_ == false) {
        lo_reliable_ = false;
        lo_reliable_cnt_ = 10;
    } else {
        if (state.lidar_odom_reliable_ && lo_reliable_cnt_ > 0) {
            lo_reliable_cnt_--;
        }

        if (lo_reliable_cnt_ == 0) {
            lo_reliable_ = true;
        }
    }

    return true;
}

bool LidarLoc::YawSearch(SE3& pose, double& confidence, CloudPtr input, CloudPtr output) {
    SE3 init_pose = pose;
    auto RPYXYZ = math::SE3ToRollPitchYaw(init_pose);
    double init_yaw = RPYXYZ.yaw;

    confidence = 0;
    bool yaw_search_success = false;

    int step = lidar_loc::grid_search_angle_step;
    double radius = lidar_loc::grid_search_angle_range * constant::kDEG2RAD;
    double angle_search_step = 2 * radius / step;

    std::vector<double> searched_yaw;
    std::vector<double> scores(step);
    std::vector<int> index;
    std::vector<SE3> pose_opti(step);

    for (int i = 0; i < step; ++i) {
        double search_yaw = init_yaw + i * angle_search_step - radius;
        searched_yaw.emplace_back(search_yaw);
        index.emplace_back(i);
    }

    LOG(INFO) << "init yaw: " << init_yaw << ", p: " << RPYXYZ.pitch << ", ro: " << RPYXYZ.roll << ", search from "
              << searched_yaw.front() << " to " << searched_yaw.back();

    /// 粗分辨率
    std::for_each(index.begin(), index.end(), [&](int i) {
        double fitness_score = 0;
        RPYXYZ.yaw = searched_yaw[i];
        SE3 pose_esti = math::XYZRPYToSE3(RPYXYZ);

        Localize(pose_esti, fitness_score, input, output, true);

        scores[i] = fitness_score;
        pose_opti[i] = pose_esti;
    });

    // find best match
    auto best_score_idx = std::max_element(scores.begin(), scores.end()) - scores.begin();
    confidence = scores.at(best_score_idx);
    pose = pose_opti.at(best_score_idx);

    /// 高分辨率
    if (confidence > options_.min_init_confidence_) {
        Localize(pose, confidence, input, output, false);
    }

    if (confidence > options_.min_init_confidence_) {
        LOG(INFO) << "init success, score: " << confidence << ", th=" << options_.min_init_confidence_;
        Eigen::Vector3d suc_translation = pose.translation();
        Eigen::Matrix3d suc_rotation_matrix = pose.rotationMatrix();
        double suc_x = suc_translation.x();
        double suc_y = suc_translation.y();
        double suc_yaw = atan2(suc_rotation_matrix(1, 0), suc_rotation_matrix(0, 0));
        LOG(INFO) << "localization init success, pose: " << suc_x << ", " << suc_y << ", " << suc_yaw
                  << ", conf: " << confidence;
        yaw_search_success = true;
    }

    return yaw_search_success;
}

bool LidarLoc::InitWithFP(CloudPtr input, const SE3& fp_pose) {
    assert(input != nullptr && !input->empty());

    // 使用功能点的位置进行定位初始化
    double fitness_score;
    SE3 pose_esti = fp_pose;
    CloudPtr output_cloud(new PointCloudType);
    // loc_inited_ = YawSearch(pose_esti, fitness_score, input, output_cloud);
    loc_inited_ = Localize(pose_esti, fitness_score, input, output_cloud);
    if (loc_inited_ && !ValidateRelocalizationMapConsistency(input, pose_esti)) {
        loc_inited_ = false;
        last_match_stats_.success = false;
    }

    if (loc_inited_) {
        current_timestamp_ = math::ToSec(input->header.stamp);
        localization_result_.confidence_ = fitness_score;
        current_abs_pose_ = pose_esti;
        localization_result_.pose_ = pose_esti;
        localization_result_.timestamp_ = current_timestamp_;
        localization_result_.lidar_loc_valid_ = true;
        localization_result_.status_ = LocalizationStatus::GOOD;

        last_abs_pose_set_ = true;
        last_abs_pose_ = pose_esti;

        current_score_ = fitness_score;
        LOG(INFO) << "fitness_score is: " << fitness_score << ", global_pose is: " << fp_pose.translation().transpose();
        LOG(INFO) << " [Loc init pose]: " << last_abs_pose_.translation().transpose();
        map_height_ = fp_pose.translation()[2];

        if (current_lo_pose_set_) {
            // 设置上一次的相对定位结果
            last_lo_pose_ = current_lo_pose_;
            last_lo_pose_set_ = true;

            last_dr_pose_ = current_dr_pose_;
            last_dr_pose_set_ = true;
        }

        //  定位成功，则清空失败记录
        fp_init_fail_pose_vec_.clear();
    } else {
        // 添加失败历史记录
        LOG(INFO) << "init failed, score: " << fitness_score;
        fp_init_fail_pose_vec_.emplace_back(fp_pose);
        fp_last_tried_time_ = 1e-6 * static_cast<double>(input->header.stamp);
    }
    return loc_inited_;
}

bool LidarLoc::BuildRelocalizationMapCache() {
    const auto begin = std::chrono::steady_clock::now();
    const std::filesystem::path global_map_path =
        std::filesystem::path(options_.map_option_.map_path_) / "global.pcd";
    CloudPtr static_map(new PointCloudType);
    if (pcl::io::loadPCDFile(global_map_path.string(), *static_map) != 0 ||
        static_map->empty()) {
        LOG(ERROR) << "failed to load relocalization validation map: " << global_map_path;
        return false;
    }

    Vec3d minimum = Vec3d::Constant(std::numeric_limits<double>::infinity());
    Vec3d maximum = Vec3d::Constant(-std::numeric_limits<double>::infinity());
    for (const auto& point : static_map->points) {
        if (!std::isfinite(point.x) || !std::isfinite(point.y) || !std::isfinite(point.z)) continue;
        const Vec3d position(point.x, point.y, point.z);
        minimum = minimum.cwiseMin(position);
        maximum = maximum.cwiseMax(position);
    }
    if (!minimum.allFinite() || !maximum.allFinite()) {
        LOG(ERROR) << "relocalization validation map has no finite points: " << global_map_path;
        return false;
    }

    std::vector<std::unique_ptr<pcl::KdTreeFLANN<PointType>>> kdtrees;
    kdtrees.reserve(static_cast<std::size_t>(options_.relocalization_validation_workers_));
    for (int worker = 0; worker < options_.relocalization_validation_workers_; ++worker) {
        auto kdtree = std::make_unique<pcl::KdTreeFLANN<PointType>>();
        kdtree->setInputCloud(static_map);
        kdtrees.push_back(std::move(kdtree));
    }

    relocalization_static_map_ = std::move(static_map);
    relocalization_map_min_ = minimum;
    relocalization_map_max_ = maximum;
    relocalization_kdtrees_ = std::move(kdtrees);
    const double elapsed_ms = std::chrono::duration<double, std::milli>(
                                  std::chrono::steady_clock::now() - begin)
                                  .count();
    LOG(INFO) << "cached relocalization validation map: points="
              << relocalization_static_map_->size()
              << ", workers=" << relocalization_kdtrees_.size()
              << ", build_ms=" << elapsed_ms
              << ", path=" << global_map_path;
    return true;
}

bool LidarLoc::TryGlobalRelocalization(const CloudPtr& input) {
    if (!global_relocalizer_ || !global_relocalizer_->IsReady() ||
        !current_lo_pose_set_) return false;

    const bool reuse_pending = pending_relocalization_.valid &&
        current_timestamp_ - pending_relocalization_.timestamp <=
            options_.relocalization_confirmation_max_interval_;
    std::optional<RelocalizationResult> result;
    if (reuse_pending) {
        RelocalizationResult pending_result;
        pending_result.attempted = true;
        pending_result.candidate_found = true;
        pending_result.accepted = true;
        pending_result.timestamp = current_timestamp_;
        pending_result.candidate_id = pending_relocalization_.candidate_id;
        pending_result.score = pending_relocalization_.score;
        pending_result.reason = "pending_confirmation";
        RelocalizationCandidate candidate;
        candidate.candidate_id = pending_relocalization_.candidate_id;
        candidate.score = pending_relocalization_.score;
        candidate.query_submap_size =
            pending_relocalization_.query_submap_size;
        candidate.T_world_imu =
            pending_relocalization_.T_map_odom * current_lo_pose_;
        pending_result.T_world_imu = candidate.T_world_imu;
        pending_result.candidates.push_back(std::move(candidate));
        result = std::move(pending_result);
    } else {
        result = global_relocalizer_->AddFrame(
            input, current_lo_pose_, current_timestamp_);
    }
    if (!result) return false;

    MatchStats summary;
    summary.relocalization_attempted = result->attempted;
    summary.relocalization_candidate_found = result->candidate_found;
    summary.relocalization_candidate_id = result->candidate_id;
    summary.relocalization_score = result->score;
    summary.relocalization_candidate_count = static_cast<int>(result->candidates.size());
    summary.relocalization_search_time_ms = result->search_time_ms;
    summary.relocalization_reason = result->reason;
    if (!result->accepted || result->candidates.empty()) {
        last_match_stats_ = summary;
        LOG(WARNING) << "GLOBAL_RELOCALIZATION[" << relocalization_backend_name_
                     << "] rejected: reason=" << result->reason
                     << ", candidate=" << result->candidate_id << ", score=" << result->score
                     << ", points=" << result->point_count
                     << ", descriptors=" << result->descriptor_count
                     << ", search_ms=" << result->search_time_ms;
        return false;
    }

    std::vector<MapConsistencyResult> prechecks(result->candidates.size());
    if (!options_.enable_relocalization_map_consistency_) {
        for (auto& precheck : prechecks) precheck.passed = true;
    } else {
        const std::size_t worker_count = std::min(
            result->candidates.size(), relocalization_kdtrees_.size());
        std::vector<std::future<void>> workers;
        workers.reserve(worker_count);
        for (std::size_t worker = 0; worker < worker_count; ++worker) {
            workers.emplace_back(std::async(std::launch::async, [&, worker]() {
                for (std::size_t index = worker; index < result->candidates.size();
                     index += worker_count) {
                    prechecks[index] = EvaluateRelocalizationMapConsistency(
                        input, result->candidates[index].T_world_imu, worker);
                }
            }));
        }
        for (auto& worker : workers) worker.get();
    }
    summary.relocalization_candidates_prechecked = static_cast<int>(prechecks.size());

    for (std::size_t index = 0; index < result->candidates.size(); ++index) {
        const auto& candidate = result->candidates[index];
        if (!prechecks[index].passed) {
            LOG(WARNING) << "GLOBAL_RELOCALIZATION[" << relocalization_backend_name_
                         << "] candidate rejected by map precheck: candidate="
                         << candidate.candidate_id << ", query_frames="
                         << candidate.query_submap_size << ", retrieval score=" << candidate.score
                         << ", pose=" << candidate.T_world_imu.translation().transpose()
                         << ", inside_xy=" << prechecks[index].inside_xy_ratio
                         << ", overlap=" << prechecks[index].overlap_ratio
                         << ", gravity_alignment_cos="
                         << prechecks[index].gravity_alignment_cos;
            continue;
        }

        map_->LoadOnPose(candidate.T_world_imu);
        UpdateGlobalMap();
        SE3 refined_pose = candidate.T_world_imu;
        double ndt_confidence = 0.0;
        CloudPtr output(new PointCloudType);
        const bool ndt_accepted = Localize(
            refined_pose, ndt_confidence, input, output);
        if (!ndt_accepted ||
            !ValidateRelocalizationMapConsistency(input, refined_pose)) {
            LOG(WARNING) << "GLOBAL_RELOCALIZATION[" << relocalization_backend_name_
                         << "] candidate rejected by NDT/map: candidate="
                         << candidate.candidate_id << ", retrieval score=" << candidate.score
                         << ", NDT confidence=" << ndt_confidence
                         << ", refined_pose=" << refined_pose.translation().transpose()
                         << ", overlap=" << last_match_stats_.map_overlap_ratio;
            continue;
        }

        const MatchStats candidate_stats = last_match_stats_;
        summary.confidence = candidate_stats.confidence;
        summary.iterations = candidate_stats.iterations;
        summary.active_map_chunks = candidate_stats.active_map_chunks;
        summary.map_consistency_evaluated = candidate_stats.map_consistency_evaluated;
        summary.map_consistency_passed = candidate_stats.map_consistency_passed;
        summary.map_consistency_points = candidate_stats.map_consistency_points;
        summary.map_inside_xy_ratio = candidate_stats.map_inside_xy_ratio;
        summary.map_inside_xyz_ratio = candidate_stats.map_inside_xyz_ratio;
        summary.map_overlap_ratio = candidate_stats.map_overlap_ratio;
        summary.map_gravity_alignment_cos = candidate_stats.map_gravity_alignment_cos;
        summary.relocalization_candidate_found = true;
        summary.relocalization_candidate_id = candidate.candidate_id;
        summary.relocalization_score = candidate.score;
        summary.relocalization_query_submap_size = candidate.query_submap_size;

        const SE3 T_map_odom = refined_pose * current_lo_pose_.inverse();
        const bool close_in_time = pending_relocalization_.valid &&
            current_timestamp_ - pending_relocalization_.timestamp <=
                options_.relocalization_confirmation_max_interval_;
        const SE3 delta = pending_relocalization_.valid
                              ? pending_relocalization_.T_map_odom.inverse() * T_map_odom
                              : SE3();
        const double rotation_delta_deg =
            delta.so3().log().norm() * 180.0 / M_PI;
        const bool consistent = close_in_time &&
            delta.translation().norm() <=
                options_.relocalization_confirmation_max_translation_ &&
            rotation_delta_deg <=
                options_.relocalization_confirmation_max_rotation_deg_;
        const int confirmation_count = consistent
                                           ? pending_relocalization_.confirmation_count + 1
                                           : 1;
        pending_relocalization_.valid = true;
        pending_relocalization_.T_map_odom = T_map_odom;
        pending_relocalization_.pose = refined_pose;
        pending_relocalization_.candidate_id = candidate.candidate_id;
        pending_relocalization_.query_submap_size = candidate.query_submap_size;
        pending_relocalization_.score = candidate.score;
        pending_relocalization_.ndt_confidence = ndt_confidence;
        pending_relocalization_.timestamp = current_timestamp_;
        pending_relocalization_.confirmation_count = confirmation_count;
        summary.relocalization_confirmation_count = confirmation_count;

        if (confirmation_count < options_.relocalization_confirmation_count_) {
            summary.relocalization_reason = "awaiting_consistent_window";
            last_match_stats_ = summary;
            LOG(INFO) << "GLOBAL_RELOCALIZATION[" << relocalization_backend_name_
                      << "] awaiting confirmation: candidate="
                      << candidate.candidate_id << ", count=" << confirmation_count
                      << "/" << options_.relocalization_confirmation_count_
                      << ", retrieval score=" << candidate.score
                      << ", NDT confidence=" << ndt_confidence
                      << ", overlap=" << summary.map_overlap_ratio
                      << ", pose=" << refined_pose.translation().transpose();
            return false;
        }

        loc_inited_ = true;
        current_abs_pose_ = refined_pose;
        last_abs_pose_ = refined_pose;
        last_abs_pose_set_ = true;
        current_score_ = ndt_confidence;
        map_height_ = refined_pose.translation().z();
        localization_result_.timestamp_ = current_timestamp_;
        localization_result_.confidence_ = ndt_confidence;
        localization_result_.pose_ = refined_pose;
        localization_result_.lidar_loc_valid_ = true;
        localization_result_.status_ = LocalizationStatus::GOOD;
        if (current_lo_pose_set_) {
            last_lo_pose_ = current_lo_pose_;
            last_lo_pose_set_ = true;
        }
        if (current_dr_pose_set_) {
            last_dr_pose_ = current_dr_pose_;
            last_dr_pose_set_ = true;
        }
        fp_init_fail_pose_vec_.clear();
        fp_last_tried_time_ = 0.0;
        match_fail_count_ = 0;
        initial_pose_set_ = false;
        summary.success = true;
        summary.relocalization_accepted = true;
        summary.relocalization_reason = "accepted";
        last_match_stats_ = summary;
        const int accepted_confirmation_count = confirmation_count;
        pending_relocalization_ = PendingRelocalization{};
        global_relocalizer_->ResetQuery();
        LOG(INFO) << "GLOBAL_RELOCALIZATION[" << relocalization_backend_name_
                  << "] accepted: candidate=" << candidate.candidate_id
                  << ", query_frames=" << candidate.query_submap_size
                  << ", confirmations=" << accepted_confirmation_count
                  << ", retrieval score=" << candidate.score
                  << ", NDT confidence=" << ndt_confidence
                  << ", overlap=" << summary.map_overlap_ratio
                  << ", pose=" << current_abs_pose_.translation().transpose();
        return true;
    }

    summary.relocalization_reason = "all_candidates_rejected";
    last_match_stats_ = summary;
    if (reuse_pending) pending_relocalization_ = PendingRelocalization{};
    LOG(WARNING) << "GLOBAL_RELOCALIZATION[" << relocalization_backend_name_
                 << "] all Top-K candidates rejected: count="
                 << result->candidates.size() << ", search_ms=" << result->search_time_ms;
    return false;
}

LidarLoc::MapConsistencyResult LidarLoc::EvaluateRelocalizationMapConsistency(
    const CloudPtr& input, const SE3& pose, std::size_t worker_index) const {
    MapConsistencyResult result;
    result.evaluated = true;
    if (!input || input->empty() || !current_lo_pose_set_ ||
        !relocalization_static_map_ || relocalization_static_map_->empty() ||
        worker_index >= relocalization_kdtrees_.size() ||
        !relocalization_kdtrees_[worker_index]) {
        return result;
    }

    const std::size_t maximum =
        static_cast<std::size_t>(options_.relocalization_map_consistency_max_points_);
    const std::size_t stride = std::max<std::size_t>(1, (input->size() + maximum - 1) / maximum);
    std::size_t finite_index = 0;
    std::size_t evaluated_points = 0;
    std::size_t inside_xy = 0;
    std::size_t inside_xyz = 0;
    std::size_t overlap = 0;
    const double margin = options_.relocalization_bounds_margin_;
    const double maximum_distance_sq = options_.relocalization_nearest_neighbor_distance_ *
                                       options_.relocalization_nearest_neighbor_distance_;
    std::vector<int> nearest_index(1);
    std::vector<float> nearest_distance_sq(1);
    for (const auto& source : input->points) {
        if (!std::isfinite(source.x) || !std::isfinite(source.y) || !std::isfinite(source.z)) continue;
        if (finite_index++ % stride != 0) continue;

        const Vec3d position = pose * Vec3d(source.x, source.y, source.z);
        if (!position.allFinite()) continue;
        ++evaluated_points;
        PointType transformed = source;
        transformed.x = static_cast<float>(position.x());
        transformed.y = static_cast<float>(position.y());
        transformed.z = static_cast<float>(position.z());

        const bool xy_ok = position.x() >= relocalization_map_min_.x() - margin &&
                           position.x() <= relocalization_map_max_.x() + margin &&
                           position.y() >= relocalization_map_min_.y() - margin &&
                           position.y() <= relocalization_map_max_.y() + margin;
        const bool xyz_ok = xy_ok &&
                            position.z() >= relocalization_map_min_.z() - margin &&
                            position.z() <= relocalization_map_max_.z() + margin;
        if (xy_ok) ++inside_xy;
        if (xyz_ok) ++inside_xyz;
        if (relocalization_kdtrees_[worker_index]->nearestKSearch(
                transformed, 1, nearest_index, nearest_distance_sq) > 0 &&
            nearest_distance_sq[0] <= maximum_distance_sq) {
            ++overlap;
        }
    }

    result.points = evaluated_points;
    if (evaluated_points == 0) return result;
    const double denominator = static_cast<double>(evaluated_points);
    result.inside_xy_ratio = static_cast<double>(inside_xy) / denominator;
    result.inside_xyz_ratio = static_cast<double>(inside_xyz) / denominator;
    result.overlap_ratio = static_cast<double>(overlap) / denominator;
    const SE3 T_map_odom = pose * current_lo_pose_.inverse();
    result.gravity_alignment_cos = T_map_odom.rotationMatrix()(2, 2);
    result.passed =
        result.inside_xy_ratio >= options_.relocalization_min_inside_xy_ratio_ &&
        result.overlap_ratio >= options_.relocalization_min_overlap_ratio_ &&
        result.gravity_alignment_cos >= options_.relocalization_min_gravity_alignment_cos_;
    return result;
}

void LidarLoc::ApplyMapConsistencyResult(const MapConsistencyResult& result) {
    last_match_stats_.map_consistency_evaluated = result.evaluated;
    last_match_stats_.map_consistency_passed = result.passed;
    last_match_stats_.map_consistency_points = result.points;
    last_match_stats_.map_inside_xy_ratio = result.inside_xy_ratio;
    last_match_stats_.map_inside_xyz_ratio = result.inside_xyz_ratio;
    last_match_stats_.map_overlap_ratio = result.overlap_ratio;
    last_match_stats_.map_gravity_alignment_cos = result.gravity_alignment_cos;
}

bool LidarLoc::ValidateRelocalizationMapConsistency(const CloudPtr& input, const SE3& pose) {
    if (!options_.enable_relocalization_map_consistency_) return true;
    const MapConsistencyResult result =
        EvaluateRelocalizationMapConsistency(input, pose, 0);
    ApplyMapConsistencyResult(result);

    if (!options_.relocalization_debug_dir_.empty() && input && !input->empty()) {
        CloudPtr scan_world(new PointCloudType);
        const std::size_t maximum =
            static_cast<std::size_t>(options_.relocalization_map_consistency_max_points_);
        const std::size_t stride =
            std::max<std::size_t>(1, (input->size() + maximum - 1) / maximum);
        std::size_t finite_index = 0;
        scan_world->reserve(std::min(input->size(), maximum));
        for (const auto& source : input->points) {
            if (!std::isfinite(source.x) || !std::isfinite(source.y) ||
                !std::isfinite(source.z)) {
                continue;
            }
            if (finite_index++ % stride != 0) continue;
            const Vec3d position = pose * Vec3d(source.x, source.y, source.z);
            if (!position.allFinite()) continue;
            PointType transformed = source;
            transformed.x = static_cast<float>(position.x());
            transformed.y = static_cast<float>(position.y());
            transformed.z = static_cast<float>(position.z());
            scan_world->push_back(transformed);
        }
        SaveRelocalizationBirdseye(
            relocalization_static_map_, scan_world,
            relocalization_map_min_, relocalization_map_max_, last_match_stats_);
    }

    LOG(INFO) << "MAP_CONSISTENCY points=" << result.points
              << ", global_map_points="
              << (relocalization_static_map_ ? relocalization_static_map_->size() : 0)
              << ", inside_xy=" << result.inside_xy_ratio
              << ", inside_xyz=" << result.inside_xyz_ratio
              << ", overlap=" << result.overlap_ratio
              << ", gravity_alignment_cos=" << result.gravity_alignment_cos
              << ", map_z=[" << relocalization_map_min_.z() << ", "
              << relocalization_map_max_.z() << "]"
              << ", passed=" << result.passed;
    return result.passed;
}

void LidarLoc::SaveRelocalizationBirdseye(const CloudPtr& static_map, const CloudPtr& scan_world,
                                          const Vec3d& map_min, const Vec3d& map_max,
                                          const MatchStats& stats) {
    if (options_.relocalization_debug_dir_.empty() || !static_map || !scan_world) return;

    try {
        const std::filesystem::path directory(options_.relocalization_debug_dir_);
        std::filesystem::create_directories(directory);
        const int index = ++relocalization_debug_index_;
        std::ostringstream stem;
        stem << "map_consistency_" << std::setw(4) << std::setfill('0') << index;
        pcl::io::savePCDFileBinaryCompressed((directory / (stem.str() + "_scan_world.pcd")).string(),
                                             *scan_world);

        Vec3d view_min = map_min;
        Vec3d view_max = map_max;
        for (const auto& point : scan_world->points) {
            if (!std::isfinite(point.x) || !std::isfinite(point.y)) continue;
            view_min.x() = std::min(view_min.x(), static_cast<double>(point.x));
            view_min.y() = std::min(view_min.y(), static_cast<double>(point.y));
            view_max.x() = std::max(view_max.x(), static_cast<double>(point.x));
            view_max.y() = std::max(view_max.y(), static_cast<double>(point.y));
        }
        view_min.x() -= 2.0;
        view_min.y() -= 2.0;
        view_max.x() += 2.0;
        view_max.y() += 2.0;
        constexpr int image_size = 1600;
        constexpr int padding = 40;
        cv::Mat image(image_size, image_size, CV_8UC3, cv::Scalar(20, 20, 20));
        const double range_x = std::max(1e-6, view_max.x() - view_min.x());
        const double range_y = std::max(1e-6, view_max.y() - view_min.y());
        const double scale = std::min((image_size - 2.0 * padding) / range_x,
                                      (image_size - 2.0 * padding) / range_y);
        auto pixel = [&](double x, double y) {
            return cv::Point(static_cast<int>(padding + (x - view_min.x()) * scale),
                             static_cast<int>(image_size - padding - (y - view_min.y()) * scale));
        };

        for (const auto& point : static_map->points) {
            if (!std::isfinite(point.x) || !std::isfinite(point.y)) continue;
            const cv::Point p = pixel(point.x, point.y);
            if (p.x >= 0 && p.x < image.cols && p.y >= 0 && p.y < image.rows) {
                image.at<cv::Vec3b>(p) = cv::Vec3b(95, 95, 95);
            }
        }
        cv::rectangle(image, pixel(map_min.x(), map_max.y()), pixel(map_max.x(), map_min.y()),
                      cv::Scalar(0, 210, 255), 2);
        const std::array<cv::Scalar, 4> lidar_colors = {
            cv::Scalar(255, 120, 40), cv::Scalar(80, 220, 80),
            cv::Scalar(40, 200, 255), cv::Scalar(220, 80, 220)};
        const double margin = options_.relocalization_bounds_margin_;
        for (const auto& point : scan_world->points) {
            if (!std::isfinite(point.x) || !std::isfinite(point.y)) continue;
            const bool inside = point.x >= map_min.x() - margin && point.x <= map_max.x() + margin &&
                                point.y >= map_min.y() - margin && point.y <= map_max.y() + margin;
            const cv::Scalar color = inside ? lidar_colors[static_cast<std::size_t>(point.lidar_id) % 4]
                                            : cv::Scalar(30, 30, 255);
            cv::circle(image, pixel(point.x, point.y), 1, color, -1, cv::LINE_AA);
        }
        std::ostringstream label;
        label << std::fixed << std::setprecision(3)
              << "inside_xy=" << stats.map_inside_xy_ratio
              << " inside_xyz=" << stats.map_inside_xyz_ratio
              << " overlap=" << stats.map_overlap_ratio
              << " gravity=" << stats.map_gravity_alignment_cos
              << " pass=" << (stats.map_consistency_passed ? "yes" : "no");
        cv::putText(image, label.str(), cv::Point(45, 32), cv::FONT_HERSHEY_SIMPLEX, 0.75,
                    stats.map_consistency_passed ? cv::Scalar(80, 230, 80) : cv::Scalar(50, 80, 255), 2,
                    cv::LINE_AA);
        cv::imwrite((directory / (stem.str() + "_birdseye.png")).string(), image);
    } catch (const std::exception& error) {
        LOG(WARNING) << "failed to save relocalization birdseye debug output: " << error.what();
    }
}

void LidarLoc::ResetLastPose(const SE3& last_pose) {
    last_abs_pose_ = last_pose;

    // TODO：清空动态图层

    return;
}

bool LidarLoc::TryOtherSolution(CloudPtr input, SE3& pose) {
    double fitness_score;
    SE3 pose_esti = pose;
    CloudPtr output_cloud(new PointCloudType);

    bool loc_success = Localize(pose_esti, fitness_score, input, output_cloud);

    if (loc_success) {
        // 激光重置逻辑
        float score_th = std::min(1.5 * current_score_, current_score_ + 0.3);
        if (fitness_score > score_th && fitness_score > 1.0) {
            // 显著好于现在的估计
            LOG(WARNING) << "rtk solution is significantly better: " << fitness_score << " " << current_score_;
            pose = pose_esti;
            localization_result_.lidar_loc_smooth_flag_ = false;
            return true;
        } else {
            LOG(INFO) << "not using rtk solution: " << fitness_score << " " << current_score_;
            return false;
        }
    }
    return false;
}

bool LidarLoc::UpdateGlobalMap() {
    NDTType::Ptr ndt(new NDTType());
    ndt->setResolution(1.0);
    ndt->setNeighborhoodSearchMethod(pclomp::DIRECT7);
    ndt->setStepSize(0.1);
    ndt->setMaximumIterations(4);
    ndt->setNumThreads(4);

    map_->SetNewTargetForNDT(ndt);
    ndt->initCompute();

    UL lock(match_mutex_);
    pcl_ndt_ = ndt;

    if (!loc_inited_) {
        NDTType::Ptr ndt_rough(new NDTType());
        ndt_rough->setResolution(5.0);
        ndt_rough->setNeighborhoodSearchMethod(pclomp::DIRECT7);
        ndt_rough->setStepSize(0.1);
        ndt_rough->setMaximumIterations(4);
        ndt_rough->setNumThreads(4);

        map_->SetNewTargetForNDT(ndt_rough);
        // ndt_rough->initCompute();

        pcl_ndt_rough_ = ndt_rough;
    }

    if (options_.enable_icp_adjust_) {
        ICPType::Ptr icp(new ICPType());
        CloudPtr map_cloud(new PointCloudType);
        pcl::VoxelGrid<PointType> voxel;
        auto sz = 0.5;
        voxel.setLeafSize(sz, sz, sz);
        voxel.setInputCloud(map_->GetAllMap());
        voxel.filter(*map_cloud);
        icp->setInputTarget(map_cloud);
        icp->setMaximumIterations(4);
        icp->setTransformationEpsilon(0.01);
        pcl_icp_ = icp;
    }

    return true;
}

void LidarLoc::UpdateMapThread() {
    LOG(INFO) << "UpdateMapThread thread is running";
    while (!update_map_quit_) {
        if (map_->MapUpdated() || map_->DynamicMapUpdated()) {
            UpdateGlobalMap();

            if (ui_) {
                ui_->UpdatePointCloudGlobal(map_->GetStaticCloud());
                ui_->UpdatePointCloudDynamic(map_->GetDynamicCloud());
            }

            map_->CleanMapUpdate();
        }
        usleep(10000);
    }
}

void LidarLoc::SetInitialPose(SE3 init_pose) {
    UL lock(initial_pose_mutex_);
    loc_inited_ = false;
    // map_->ClearMap();

    initial_pose_set_ = true;
    initial_pose_ = init_pose;
    LOG(INFO) << "Set initial pose is: " << initial_pose_.translation().transpose();
}

void LidarLoc::RequestGlobalRelocalization() {
    UL lock(initial_pose_mutex_);
    loc_inited_ = false;
    initial_pose_set_ = false;
    match_fail_count_ = 0;
    last_match_stats_ = MatchStats{};
    pending_relocalization_ = PendingRelocalization{};
    if (global_relocalizer_) global_relocalizer_->ResetQuery();
    LOG(WARNING) << "GLOBAL_RELOCALIZATION[" << relocalization_backend_name_
                 << "] requested";
}

void LidarLoc::Align(const CloudPtr& input) {
    // 输入必须非空
    assert(input != nullptr);

    // 点云去畸变定到了结束时间，所以该点云的定位也是到结束时间的
    double current_time = math::ToSec(input->header.stamp) + lo::lidar_time_interval;
    current_timestamp_ = current_time;

    LOG(INFO) << "current time: " << std::fixed << std::setprecision(12) << current_timestamp_;

    /// 设置当前帧对应的rel_pose
    if (!AssignLOPose(current_time)) {
        LOG(WARNING) << "assign LO pose failed";
    }

    if (!AssignDRPose(current_time)) {
        LOG(WARNING) << "assign DR pose failed";
    }

    /// 1. 车辆静止处理
    if (parking_ && loc_inited_) {
        LOG(INFO) << "车辆静止，不做匹配";

        UpdateState(input);
        current_abs_pose_ = last_abs_pose_;
        lidar_loc_pose_queue_.emplace_back(current_time, current_abs_pose_);

        UL lock(result_mutex_);
        localization_result_.timestamp_ = current_time;
        localization_result_.pose_ = current_abs_pose_;

        if (options_.enable_parking_static_) {
            localization_result_.is_parking_ = true;
            localization_result_.valid_ = true;
        }
        return;
    }

    /// 2. 初始化处理
    if (!loc_inited_) {
        UL lock_init(initial_pose_mutex_);
        LOG(INFO) << "initing lidarloc";
        SetInitRltState();
        last_match_stats_ = MatchStats{};

        if (initial_pose_set_) {
            /// 尝试在给定点初始化
            if (InitWithFP(input, initial_pose_)) {
                LOG(INFO) << "init with external pose: " << initial_pose_.translation().transpose();
                initial_pose_set_ = false;
                return;
            }
        }

        if (options_.init_with_fp_) {
            /// 从功能点初始化
            /// 如果之前尝试过，那么需要间隔一段时间再进行搜索
            if (!fp_init_fail_pose_vec_.empty() && current_dr_pose_set_) {
                SE3 last_tried_pose = fp_init_fail_pose_vec_.back();
                bool should_try =
                    (current_time - fp_last_tried_time_) > 2.0 ||
                    (current_dr_pose_.translation() - last_tried_pose.translation()).norm() > 0.3 ||
                    (current_dr_pose_.so3().inverse() * last_tried_pose.so3()).log().norm() > 10 * M_PI / 180.0;
                if (!should_try) {
                    LOG(INFO) << "skip trying init, please move to another place.";
                    return;
                }
            } else {
                LOG(INFO) << "fp tried pose: " << fp_init_fail_pose_vec_.size()
                          << ", dr pose set: " << current_dr_pose_set_;
            }

            auto all_fps = map_->GetAllFP();
            bool fp_init_success = false;
            for (const auto& fp : all_fps) {
                map_->LoadOnPose(fp.pose_);
                if (InitWithFP(input, fp.pose_)) {
                    LOG(INFO) << "init with fp: " << fp.name_;
                    fp_init_success = true;
                    break;
                }
            }

            if (!fp_init_success) {
                LOG(INFO) << "FP init failed.";
                if (current_dr_pose_set_) {
                    LOG(INFO) << "record fp failed time: " << std::setprecision(12) << current_time
                              << ", pose: " << current_dr_pose_.translation().transpose();
                    fp_last_tried_time_ = current_time;
                    fp_init_fail_pose_vec_.emplace_back(current_dr_pose_);
                }
            } else {
                fp_last_tried_time_ = 0;
                fp_init_fail_pose_vec_.clear();
            }
        }

        if (TryGlobalRelocalization(input)) return;

        /// 初始化未成功时，不往下走流程
        return;
    }

    /// 4. 设置当前帧对应的 pose guess
    /// NOTE: LO设置预测的位置和LidarLoc自身递推设置预测的方法并不完全一致，自身外推容易受噪声影响

    SE3 guess_from_lo = last_abs_pose_;
    if (last_lo_pose_set_ && current_lo_pose_set_) {
        // 如果有里程计，则用两个时刻的相对定位来递推，估计一个当前pose的初值
        const SE3 delta = last_lo_pose_.inverse() * current_lo_pose_;
        guess_from_lo = last_abs_pose_ * delta;

        LOG(INFO) << "current lo pose: " << current_lo_pose_.translation().transpose();
        LOG(INFO) << "last lo pose: " << last_lo_pose_.translation().transpose();
        LOG(INFO) << "lo motion: " << delta.translation().transpose();
        LOG(INFO) << "last abs pose: " << last_abs_pose_.translation().transpose();
        // guess_from_lo.translation()[2] = 0;
        LOG(INFO) << "loc using lo guess: " << guess_from_lo.translation().transpose();
    }

    SE3 guess_from_self = guess_from_lo;
    if (lidar_loc_pose_queue_.size() >= 2) {
        SE3 pred;
        TimedPose match;
        if (math::PoseInterp<TimedPose>(
                current_time, lidar_loc_pose_queue_, [](const TimedPose& p) { return p.timestamp_; },
                [](const TimedPose& p) { return p.pose_; }, pred, match, 2.0)) {
            guess_from_self = pred;
        }
    }

    // SE3 guess_from_dr = guess_from_lo;
    // if (last_dr_pose_set_ && current_dr_pose_set_) {
    //     const SE3 delta = last_dr_pose_.inverse() * current_dr_pose_;
    //     guess_from_dr = last_abs_pose_ * delta;
    //     // guess_from_dr.translation()[2] = 0;
    // }

    // bool try_dr = false;
    // if (((guess_from_dr.translation() - guess_from_lo.translation()).norm() >= try_other_guess_trans_th_ ||
    //      (guess_from_dr.so3().inverse() * guess_from_lo.so3()).log().norm() >= try_other_guess_rot_th_)) {
    //     LOG(INFO) << "trying dr pose: " << guess_from_dr.translation().transpose() << ", "
    //               << (guess_from_dr.so3().inverse() * guess_from_lo.so3()).log().norm()
    //               << ", vel_norm: " << current_vel_b_.norm();
    //     try_dr = true;
    // }

    bool try_self = false;
    // if (options_.try_self_extrap_) {
    //     if (((guess_from_self.translation() - guess_from_lo.translation()).norm() >= try_other_guess_trans_th_ ||
    //          (guess_from_self.so3().inverse() * guess_from_lo.so3()).log().norm() >= try_other_guess_rot_th_) &&
    //         ((guess_from_dr.translation() - guess_from_self.translation()).norm() >= try_other_guess_trans_th_ ||
    //          (guess_from_dr.so3().inverse() * guess_from_self.so3()).log().norm() >= try_other_guess_rot_th_)) {
    //         LOG(INFO) << "trying self extrap pose: " << guess_from_self.translation().transpose() << ", "
    //                   << (guess_from_self.so3().inverse() * guess_from_lo.so3()).log().norm();
    //         try_self = true;
    //     }
    // }

    /// 5. 载入地图, 与地图匹配定位
    /// 尝试各种初始估计
    CloudPtr output_cloud(new PointCloudType);
    double fitness_score = 0;
    SE3 current_pose_esti = guess_from_lo;
    bool loc_success_lo, loc_success_self, loc_success_dr;
    loc_success_lo = loc_success_self = loc_success_dr = false;
    bool loc_success = false;

    /// 注意load on pose存在滞后，优先load on DR
    map_->LoadOnPose(guess_from_lo);

    loc_success_lo = Localize(current_pose_esti, fitness_score, input, output_cloud);  // LO 那个肯定会算
    double score_lo = fitness_score;

    SE3 res_of_lo = current_pose_esti;
    SE3 res_of_dr = current_pose_esti;
    SE3 res_of_self = current_pose_esti;
    double score_dr = 0;

    // 先尝试外部预测，最后用自身
    // if (try_dr) {
    //     /// 尝试DR外推的pose
    //     res_of_dr = guess_from_dr;
    //     loc_success_dr = Localize(res_of_dr, score_dr, input, output_cloud);
    //     if (score_dr > (fitness_score - 0.1)) {
    //         current_pose_esti = res_of_dr;
    //         fitness_score = score_dr;
    //         LOG(INFO) << "take dr guess: " << current_pose_esti.translation().transpose()
    //                   << " , confidence: " << score_dr << ", v_norm: " << current_vel_b_.norm();
    //     }
    // }

    // 用纯激光定位有点太抖了，加一些权重
    Vec6d delta = (guess_from_lo.inverse() * current_pose_esti).log();
    SE3 esti_balanced = guess_from_lo * SE3::exp(delta * 0.1);
    current_pose_esti = esti_balanced;

    // double score_self = 0;
    // if (try_self) {
    //     /// 尝试自身外推的pose
    //     LOG(INFO) << "localize with extrap";

    //     res_of_self = guess_from_self;
    //     loc_success_self = Localize(res_of_self, score_self, input, output_cloud);

    //     // 避免分值接近但长时间采信自身预测，此处更相信外部预测源
    //     if (score_self > (fitness_score + 0.1)) {
    //         current_pose_esti = res_of_self;
    //         fitness_score = score_self;
    //         LOG(INFO) << "take self guess: " << current_pose_esti.translation().transpose()
    //                   << " , confidence: " << score_self;
    //     }
    // }

    /// NOTE 如果LO, DR出发点和收敛点不同，但分值相近，说明场景可能处在退化状态，此时使用DR预测的Pose
    // if (try_dr && (res_of_lo.translation() - res_of_dr.translation()).head<2>().norm() > 0.2 &&
    //     fabs(score_lo - score_dr) < 0.2 && score_lo < 1.2) {
    //     LOG(WARNING) << "判定激光定位进入退化状态，现在会使用DR递推pose而不是激光定位位置";
    //     current_pose_esti = guess_from_dr;
    // }

    if (options_.force_2d_) {
        PoseRPYD RPYXYZ = math::SE3ToRollPitchYaw(current_pose_esti);
        RPYXYZ.roll = 0;
        RPYXYZ.pitch = 0;
        RPYXYZ.z = 0;
        current_pose_esti = math::XYZRPYToSE3(RPYXYZ);
    }

    // if (options_.with_height_) {
    //     current_pose_esti.translation()[2] = map_height_;
    //     LOG(INFO) << "adjust current pose to : " << current_pose_esti.translation().transpose();
    // }

    current_score_ = fitness_score;
    double delta_rel_abs_pose = 0;
    bool lidar_loc_odom_valid = true;

    if (loc_success_lo || loc_success_self || loc_success_dr) {
        loc_success = true;
    } else {
        LOG(INFO) << "loc success is false.";
    }

    if (loc_success) {
        lidar_loc_odom_valid = CheckLidarOdomValid(current_pose_esti, delta_rel_abs_pose);
        match_fail_count_ = 0;
        last_timestamp_ = current_timestamp_;  // 成功时，更新上一时刻激光定位时间
    } else {
        current_score_ = fitness_score;
        LOG(WARNING) << "localization failed! score: " << current_score_;
        ++match_fail_count_;
        // Do not propagate a rejected NDT transform. Lidar odometry remains
        // the short-term motion source while global relocalization starts.
        current_pose_esti = guess_from_lo;
        if (global_relocalizer_ && global_relocalizer_->IsReady() &&
            match_fail_count_ >= options_.relocalization_lost_frame_threshold_) {
            loc_inited_ = false;
            global_relocalizer_->ResetQuery();
            pending_relocalization_ = PendingRelocalization{};
            LOG(WARNING) << "GLOBAL_RELOCALIZATION[" << relocalization_backend_name_
                         << "] tracking lost after " << match_fail_count_
                         << " consecutive rejected NDT matches";
        }
    }

    current_abs_pose_ = current_pose_esti;

    /// 确定激光定位是否满足平滑性要求
    Vec3d dpred = current_abs_pose_.translation() - guess_from_self.translation();
    if (fabs(dpred[0]) < 0.5 && fabs(dpred[1]) < 0.5 &&
        (current_abs_pose_.so3().inverse() * guess_from_self.so3()).log().norm() < 2.0 * M_PI / 180.0) {
        localization_result_.lidar_loc_smooth_flag_ = true;
    } else {
        localization_result_.lidar_loc_smooth_flag_ = false;
    }

    localization_result_.lidar_loc_odom_reliable_ = lo_reliable_;
    localization_result_.is_parking_ = false;

    // if (ui_) {
    //     ui_->UpdatePredictPose(guess_from_lo);
    // }

    /// 7. 输出结果
    {
        UL lock(result_mutex_);
        localization_result_.timestamp_ = current_timestamp_;
        localization_result_.confidence_ = fitness_score;
        if (!loc_inited_) {
            localization_result_.lidar_loc_valid_ = false;
            localization_result_.status_ = LocalizationStatus::INITIALIZING;
        } else if (loc_success && lidar_loc_odom_valid) {
            localization_result_.lidar_loc_valid_ = true;
            localization_result_.status_ = LocalizationStatus::GOOD;
        } else {
            localization_result_.lidar_loc_valid_ = false;
            localization_result_.status_ = LocalizationStatus::FOLLOWING_DR;
        }

        localization_result_.lidar_loc_odom_delta_ = delta_rel_abs_pose;
        localization_result_.lidar_loc_odom_error_normal_ = lidar_loc_odom_valid;
        localization_result_.pose_ = current_pose_esti;
    }

    UpdateState(input);

    /// 8. 更新动态图层
    /// 条件：1. 定位成功 2. 与上次更新间隔一定距离 3. RTK与激光定位横纵向误差都小于0.3，或者匹配分值大于1.0
    bool score_cond = current_score_ > options_.update_lidar_loc_score_;

    if (options_.update_dynamic_cloud_ && loc_success &&
        (((current_pose_esti.translation() - last_dyn_upd_pose_.pose_.translation()).norm() >
          options_.update_kf_dis_) ||
         fabs(current_time - last_dyn_upd_pose_.timestamp_) > options_.update_kf_time_)) {
        if (score_cond /*  || update_cache_dis_ < options_.max_update_cache_dis_ */) {
            // LOG(INFO) << "passing through z filter, input:" << input->size();
            pcl::PassThrough<PointType> pass;
            pass.setInputCloud(input);

            pass.setFilterFieldName("z");
            pass.setFilterLimits(0.5, options_.filter_z_max_);

            CloudPtr input_z_filter(new PointCloudType());
            pass.filter(*input_z_filter);

            if (!input_z_filter->empty()) {
                CloudPtr cloud_t(new PointCloudType());
                pcl::transformPointCloud(*input_z_filter, *cloud_t, current_pose_esti.matrix());

                // 以现在的scan来更新地图
                map_->UpdateDynamicCloud(cloud_t, true);

                last_dyn_upd_pose_.timestamp_ = current_time;
                last_dyn_upd_pose_.pose_ = current_pose_esti;
            }
        }
    }

    if (lidar_loc_pose_queue_.empty()) {
        lidar_loc_pose_queue_.emplace_back(current_time, current_abs_pose_);
    } else if (current_time > lidar_loc_pose_queue_.back().timestamp_) {
        lidar_loc_pose_queue_.emplace_back(current_time, current_abs_pose_);
    }

    while (lidar_loc_pose_queue_.size() > 1000) {
        lidar_loc_pose_queue_.pop_front();
    }

    ave_scores_.emplace_back(current_score_);
    while (ave_scores_.size() > 20) {
        ave_scores_.pop_front();
    }

    /// 9. save for recover pose
    recover_pose_out_.open(options_.recover_pose_path_);
    if (recover_pose_out_) {
        Vec3d t = current_pose_esti.translation();
        Quatd q = current_pose_esti.unit_quaternion();
        recover_pose_out_ << t[0] << " " << t[1] << " " << t[2] << " " << q.x() << " " << q.y() << " " << q.z() << " "
                          << q.w();
        recover_pose_out_.close();
    }
}

bool LidarLoc::CheckLidarOdomValid(const SE3& current_pose_esti, double& delta_posi) {
    delta_posi = ((last_lo_pose_.inverse() * current_lo_pose_).translation() -
                  (last_abs_pose_.inverse() * current_pose_esti).translation())
                     .head(2)
                     .norm();

    bool valid = true;

    if (delta_posi > options_.lidar_loc_odom_th_ && last_lo_pose_set_) {
        LOG(INFO) << "delta_rel_abs_pose is: " << delta_posi;
        LOG(INFO) << "LO相对pose: " << last_lo_pose_.translation().transpose() << " "
                  << current_lo_pose_.translation().transpose() << " "
                  << (last_lo_pose_.inverse() * current_lo_pose_).translation().transpose();
        LOG(INFO) << "Lidar Loc计算相对pose: " << last_abs_pose_.translation().transpose() << " "
                  << current_pose_esti.translation().transpose() << " "
                  << (last_abs_pose_.inverse() * current_pose_esti).translation().transpose();
        lo_reliable_ = false;
        lo_reliable_cnt_ = 10;
        valid = false;
    }

    last_abs_pose_ = current_pose_esti;
    last_lo_pose_ = current_lo_pose_;
    last_lo_pose_set_ = true;
    last_dr_pose_ = current_dr_pose_;
    last_dr_pose_set_ = true;

    return valid;
}

bool LidarLoc::Localize(SE3& pose, double& confidence, CloudPtr input, CloudPtr output, bool use_rough_res) {
    Eigen::Matrix4f trans;
    bool loc_success = false;
    Eigen::Matrix4f guess_pose = pose.matrix().cast<float>();
    last_match_stats_ = MatchStats{};
    last_match_stats_.active_map_chunks = map_ ? map_->NumActiveChunks() : 0;

    LOG(INFO) << "loc from: " << pose.translation().transpose();

    if (pcl_ndt_->getInputTarget() == nullptr) {
        LOG(INFO) << "lidar loc target is null, skip";
        return false;
    }

    UL lock(match_mutex_);
    NDTType::Ptr ndt = nullptr;

    if (use_rough_res) {
        ndt = pcl_ndt_rough_;
    } else {
        ndt = pcl_ndt_;
    }

    ndt->setInputSource(input);
    ndt->align(*output, guess_pose);
    trans = ndt->getFinalTransformation();
    confidence = ndt->getTransformationProbability();
    last_match_stats_.confidence = confidence;
    last_match_stats_.iterations = ndt->getFinalNumIteration();

    const double confidence_threshold =
        loc_inited_ ? options_.min_tracking_confidence_ : options_.min_init_confidence_;
    loc_success = ndt->hasConverged() && std::isfinite(confidence) && trans.allFinite() &&
                  confidence >= confidence_threshold;

    if (options_.enable_icp_adjust_ && loc_inited_ && loc_success) {
        Eigen::Matrix4f adjust_trans;
        CloudPtr input_voxel(new PointCloudType);
        pcl::VoxelGrid<PointType> voxel_icp;

        double ls = 0.2;
        voxel_icp.setLeafSize(ls, ls, ls);
        voxel_icp.setInputCloud(input);
        voxel_icp.filter(*input_voxel);
        pcl_icp_->setInputSource(input_voxel);
        Timer::Evaluate([&]() { pcl_icp_->align(*output, trans); }, "pcl_icp adjust", true);
        adjust_trans = pcl_icp_->getFinalTransformation();

        Eigen::Matrix3f rotation_diff = trans.block<3, 3>(0, 0).transpose() * adjust_trans.block<3, 3>(0, 0);
        Eigen::AngleAxisf angle_axis(rotation_diff);
        float a = angle_axis.angle();
        float d = (trans.block<3, 1>(0, 3) - adjust_trans.block<3, 1>(0, 3)).norm();
        LOG(INFO) << "icp adjust d: " << d << ", a: " << a;

        if (pcl_icp_->hasConverged() && std::fabs(d) <= 0.05 && std::fabs(a) <= 0.05) {
            LOG(INFO) << "icp ajust trans set success";
            trans = adjust_trans;
        }
    }

    Vec3d t_3d = pose.translation();
    if (trans.allFinite()) {
        Eigen::Matrix3d rot = trans.block<3, 3>(0, 0).cast<double>();
        Quatd q_3d = Quatd(rot);
        t_3d = trans.block<3, 1>(0, 3).cast<double>();
        if (q_3d.norm() > 1e-9 && q_3d.coeffs().allFinite()) {
            q_3d.normalize();
            pose = SE3(q_3d, t_3d);
        } else {
            loc_success = false;
        }
    }

    LOG(INFO) << "confidence: " << confidence << ", t: " << t_3d.transpose() << ", succ: " << loc_success;
    last_match_stats_.success = loc_success;
    last_match_stats_.active_map_chunks = map_ ? map_->NumActiveChunks() : 0;

    return loc_success;
}

bool LidarLoc::CheckStatic(double timestamp) {
    if (parking_) {
        // if (current_vel_b_.norm() < common::options::lo::parking_speed) {
        // LOG(INFO) << "car is in static mode";
        static_count_++;
        if (static_count_ >= lo::parking_count) {
            static_count_ = 0;
            return false;
        }
        return true;
    } else {
        static_count_ = 0;
        return false;
    }
}

void LidarLoc::UpdateState(const CloudPtr& input) { last_timestamp_ = current_timestamp_; }

void LidarLoc::SetInitRltState() {
    UL lock(result_mutex_);
    localization_result_.confidence_ = 0.0;
    localization_result_.timestamp_ = current_timestamp_;
    localization_result_.lidar_loc_valid_ = false;
    localization_result_.status_ = LocalizationStatus::INITIALIZING;
}

bool LidarLoc::LocInited() {
    // UL lock( data_mutex_);
    return loc_inited_;
}

bool LidarLoc::AssignLOPose(double timestamp) {
    UL lock(lo_pose_mutex_);
    SE3 interp_pose;
    NavState best_match;
    // 无法拿到最新的，插值结果都是外推出来的，和真实有时偏差会很大（╯﹏╰）
    // if (!lo_pose_queue_.empty()) {
    //     LOG(INFO) << "lo interp: " << timestamp << " " << lo_pose_queue_.back().timestamp_ << " "
    //               << lo_pose_queue_.back().pose_.translation().transpose();
    // }

    bool pose_interp_success = math::PoseInterp<NavState>(
        timestamp, lo_pose_queue_, [](const NavState& dr) { return dr.timestamp_; },
        [](const NavState& dr) { return dr.GetPose(); }, interp_pose, best_match, 5.0);

    if (pose_interp_success) {
        current_lo_pose_ = interp_pose;
        current_lo_pose_set_ = true;

        current_vel_b_ = best_match.GetRot().inverse() * best_match.GetVel();
        current_vel_ = best_match.GetVel();

        // if (options_.with_height_) {
        //     current_lo_pose_.translation()[2] = map_height_;
        // }

        return true;
    } else {
        current_lo_pose_set_ = false;
        return false;
    }
}

bool LidarLoc::AssignDRPose(double timestamp) {
    UL lock(dr_pose_mutex_);
    SE3 interp_pose;
    NavState best_match;
    bool pose_interp_success = math::PoseInterp<NavState>(
        timestamp, dr_pose_queue_, [](const NavState& dr) { return dr.timestamp_; },
        [](const NavState& dr) { return dr.GetPose(); }, interp_pose, best_match, 5.0);

    if (pose_interp_success) {
        parking_ = best_match.is_parking_;
        current_dr_pose_ = interp_pose;
        current_dr_pose_set_ = true;

        // if (options_.with_height_) {
        //     current_dr_pose_.translation()[2] = map_height_;
        // }

        return true;
    } else {
        parking_ = false;
        current_dr_pose_set_ = false;
        return false;
    }
}

void LidarLoc::Finish() {
    if (map_) {
        LOG(INFO) << "saving maps";
        update_map_quit_ = true;
        update_map_thread_.join();

        /// 永久保存时，再存储地图
        if (options_.map_option_.policy_ == TiledMap::DynamicCloudPolicy::PERSISTENT &&
            options_.map_option_.save_dyn_when_quit_ && !has_set_pose_) {
            map_->SaveToBin(true);
            LOG(INFO) << "dynamic maps saved";
        }
    }
}

}  // namespace lightning::loc
