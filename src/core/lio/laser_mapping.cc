#include <pcl/common/transforms.h>
#include <yaml-cpp/yaml.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>

#include "common/options.h"
#include "core/lightning_math.hpp"
#include "laser_mapping.h"

#include <opencv2/core/mat.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include "ui/pangolin_window.h"
#include "utils/compute_profiling.h"
#include "wrapper/ros_utils.h"

namespace lightning {

bool LaserMapping::Init(const std::string &config_yaml) {
    LOG(INFO) << "init laser mapping from " << config_yaml;
    if (!LoadParamsFromYAML(config_yaml)) {
        return false;
    }

    // localmap init (after LoadParams)
    ivox_ = std::make_shared<IVoxType>(ivox_options_);

    // esekf init
    ESKF::Options eskf_options;
    eskf_options.max_iterations_ = fasterlio::NUM_MAX_ITERATIONS;
    eskf_options.epsi_ = 1e-3 * Eigen::Matrix<double, ESKF::state_dim_, 1>::Ones();
    // 使用Lambda表达式将LaserMapping::ObsModel方法绑定到lidar_obs_func_
    eskf_options.lidar_obs_func_ = [this](NavState &s, ESKF::CustomObservationModel &obs) { ObsModel(s, obs); };
    eskf_options.propagate_velocity_ = propagate_velocity_;
    eskf_options.lidar_update_pose_only_ = lidar_update_pose_only_;
    eskf_options.lidar_update_inertial_states_ = lidar_update_inertial_states_;
    eskf_options.max_update_velocity_step_ = max_update_velocity_step_;
    eskf_options.max_update_gyro_bias_step_ = max_update_gyro_bias_step_;
    eskf_options.max_update_acc_bias_step_ = max_update_acc_bias_step_;
    eskf_options.max_update_gravity_step_ = max_update_gravity_step_;
    eskf_options.use_aa_ = use_aa_;
    kf_.Init(eskf_options);
    velocity_propagation_active_ = propagate_velocity_;
    p_imu_->SetPostPredictCallback([this](ESKF& filter, double timestamp) {
        ApplyWheelSpeedObservation(filter, timestamp,
                                   last_lidar_filter_wheel_timestamp_, false);
    });

    LOG(INFO) << "ESKF lidar velocity gate=" << max_update_velocity_step_
              << " m/s, adaptive velocity propagation=" << adaptive_velocity_propagation_
              << ", point covariance model=" << (point_noise_enabled_ ? "enabled" : "disabled");
    LOG(INFO) << "wheel-speed DR soft observation="
              << (wheel_speed_dr_config_.enabled ? "enabled" : "disabled")
              << ", base_std=" << wheel_speed_dr_config_.base_std_mps
              << " m/s, max_age=" << wheel_speed_dr_config_.max_age_sec << " sec";

    return true;
}

bool LaserMapping::LoadParamsFromYAML(const std::string &yaml_file) {
    // get params from yaml
    int lidar_type, ivox_nearby_type;
    double gyr_cov, acc_cov, b_gyr_cov, b_acc_cov;
    double filter_size_scan;

    auto yaml = YAML::LoadFile(yaml_file);
    try {
        fasterlio::NUM_MAX_ITERATIONS = yaml["fasterlio"]["max_iteration"].as<int>();
        fasterlio::ESTI_PLANE_THRESHOLD = yaml["fasterlio"]["esti_plane_threshold"].as<float>();

        filter_size_scan = yaml["fasterlio"]["filter_size_scan"].as<float>();
        filter_size_map_min_ = yaml["fasterlio"]["filter_size_map"].as<float>();
        keep_first_imu_estimation_ = yaml["fasterlio"]["keep_first_imu_estimation"].as<bool>();
        gyr_cov = yaml["fasterlio"]["gyr_cov"].as<float>();
        acc_cov = yaml["fasterlio"]["acc_cov"].as<float>();
        b_gyr_cov = yaml["fasterlio"]["b_gyr_cov"].as<float>();
        b_acc_cov = yaml["fasterlio"]["b_acc_cov"].as<float>();
        preprocess_->Blind() = yaml["fasterlio"]["blind"].as<double>();
        preprocess_->TimeScale() = yaml["fasterlio"]["time_scale"].as<double>();
        if (yaml["fasterlio"]["livox_point_time_scale"]) {
            preprocess_->LivoxPointTimeScale() =
                yaml["fasterlio"]["livox_point_time_scale"].as<double>();
        }
        if (!std::isfinite(preprocess_->LivoxPointTimeScale()) ||
            preprocess_->LivoxPointTimeScale() <= 0.0) {
            LOG(ERROR) << "invalid livox_point_time_scale";
            return false;
        }
        lidar_type = yaml["fasterlio"]["lidar_type"].as<int>();
        preprocess_->NumScans() = yaml["fasterlio"]["scan_line"].as<int>();
        preprocess_->PointFilterNum() = yaml["fasterlio"]["point_filter_num"].as<int>();

        extrinT_ = yaml["fasterlio"]["extrinsic_T"].as<std::vector<double>>();
        extrinR_ = yaml["fasterlio"]["extrinsic_R"].as<std::vector<double>>();

        ivox_options_.resolution_ = yaml["fasterlio"]["ivox_grid_resolution"].as<float>();
        ivox_nearby_type = yaml["fasterlio"]["ivox_nearby_type"].as<int>();
        use_aa_ = yaml["fasterlio"]["use_aa"].as<bool>();
        if (yaml["fasterlio"]["propagate_velocity"]) {
            propagate_velocity_ = yaml["fasterlio"]["propagate_velocity"].as<bool>();
        }
        if (yaml["fasterlio"]["lidar_update_pose_only"]) {
            lidar_update_pose_only_ = yaml["fasterlio"]["lidar_update_pose_only"].as<bool>();
        }
        if (yaml["fasterlio"]["lidar_update_inertial_states"]) {
            lidar_update_inertial_states_ = yaml["fasterlio"]["lidar_update_inertial_states"].as<bool>();
        }
        if (yaml["fasterlio"]["max_update_velocity_step"]) {
            max_update_velocity_step_ = yaml["fasterlio"]["max_update_velocity_step"].as<double>();
        }
        if (yaml["fasterlio"]["max_update_gyro_bias_step"]) {
            max_update_gyro_bias_step_ = yaml["fasterlio"]["max_update_gyro_bias_step"].as<double>();
        }
        if (yaml["fasterlio"]["max_update_acc_bias_step"]) {
            max_update_acc_bias_step_ = yaml["fasterlio"]["max_update_acc_bias_step"].as<double>();
        }
        if (yaml["fasterlio"]["max_update_gravity_step"]) {
            max_update_gravity_step_ = yaml["fasterlio"]["max_update_gravity_step"].as<double>();
        }
        if (yaml["fasterlio"]["adaptive_velocity_propagation"]) {
            adaptive_velocity_propagation_ =
                yaml["fasterlio"]["adaptive_velocity_propagation"].as<bool>();
        }
        if (yaml["fasterlio"]["velocity_innovation_ema_alpha"]) {
            velocity_innovation_ema_alpha_ =
                yaml["fasterlio"]["velocity_innovation_ema_alpha"].as<double>();
        }
        if (yaml["fasterlio"]["velocity_innovation_enable_threshold"]) {
            velocity_innovation_enable_threshold_ =
                yaml["fasterlio"]["velocity_innovation_enable_threshold"].as<double>();
        }
        if (yaml["fasterlio"]["velocity_innovation_disable_threshold"]) {
            velocity_innovation_disable_threshold_ =
                yaml["fasterlio"]["velocity_innovation_disable_threshold"].as<double>();
        }
        if (yaml["fasterlio"]["velocity_propagation_max_active_updates"]) {
            velocity_propagation_max_active_updates_ =
                yaml["fasterlio"]["velocity_propagation_max_active_updates"].as<int>();
        }
        if (yaml["fasterlio"]["velocity_propagation_cooldown_updates"]) {
            velocity_propagation_cooldown_updates_ =
                yaml["fasterlio"]["velocity_propagation_cooldown_updates"].as<int>();
        }
        if (velocity_innovation_ema_alpha_ <= 0.0 || velocity_innovation_ema_alpha_ > 1.0 ||
            velocity_innovation_disable_threshold_ < 0.0 ||
            velocity_innovation_enable_threshold_ <= velocity_innovation_disable_threshold_ ||
            velocity_propagation_max_active_updates_ < 0 ||
            velocity_propagation_cooldown_updates_ < 0) {
            LOG(ERROR) << "invalid adaptive velocity propagation thresholds";
            return false;
        }
        skip_lidar_num_ = yaml["fasterlio"]["skip_lidar_num"].as<int>();
        enable_skip_lidar_ = skip_lidar_num_ > 0;

        float height_max = yaml["roi"]["height_max"].as<float>();
        float height_min = yaml["roi"]["height_min"].as<float>();

        preprocess_->SetHeightROI(height_max, height_min);

        options_.kf_dis_th_ = yaml["fasterlio"]["kf_dis_th"].as<double>();
        options_.kf_angle_th_ = yaml["fasterlio"]["kf_angle_th"].as<double>() * M_PI / 180.0;
        options_.enable_icp_part_ = yaml["fasterlio"]["enable_icp_part"].as<bool>();
        options_.min_pts = yaml["fasterlio"]["min_pts"].as<int>();
        options_.plane_icp_weight_ = yaml["fasterlio"]["plane_icp_weight"].as<float>();

        bool use_imu_filter = yaml["fasterlio"]["imu_filter"].as<bool>();
        p_imu_->SetUseIMUFilter(use_imu_filter);
        options_.proj_kfs_ = yaml["fasterlio"]["proj_kfs"].as<bool>();

        if (yaml["imu_initialization"]) {
            const YAML::Node init = yaml["imu_initialization"];
            ImuProcess::InitializationOptions init_options;
            if (init["min_duration"]) init_options.min_duration = init["min_duration"].as<double>();
            if (init["min_samples"]) init_options.min_samples = init["min_samples"].as<int>();
            if (init["max_mean_gyro_norm"]) {
                init_options.max_mean_gyro_norm = init["max_mean_gyro_norm"].as<double>();
            }
            if (init["max_gyro_std"]) init_options.max_gyro_std = init["max_gyro_std"].as<double>();
            if (init["max_acc_std"]) init_options.max_acc_std = init["max_acc_std"].as<double>();
            if (init["min_mean_acc_norm"]) {
                init_options.min_mean_acc_norm = init["min_mean_acc_norm"].as<double>();
            }
            if (init["max_mean_acc_norm"]) {
                init_options.max_mean_acc_norm = init["max_mean_acc_norm"].as<double>();
            }
            if (init["initial_yaw_deg"]) {
                init_options.initial_yaw_deg = init["initial_yaw_deg"].as<double>();
            }
            if (init_options.min_duration <= 0.0 || init_options.min_samples < 2 ||
                init_options.min_mean_acc_norm < 0.0 ||
                init_options.max_mean_acc_norm < init_options.min_mean_acc_norm ||
                !std::isfinite(init_options.initial_yaw_deg)) {
                LOG(ERROR) << "invalid imu_initialization configuration";
                return false;
            }
            p_imu_->SetInitializationOptions(init_options);
        }

        const YAML::Node system = yaml["system"];
        if (system) {
            if (system["enable_wheel_speed_dr_observation"]) {
                wheel_speed_dr_config_.enabled =
                    system["enable_wheel_speed_dr_observation"].as<bool>();
            }
            if (system["wheel_speed_dr_base_std_mps"])
                wheel_speed_dr_config_.base_std_mps = system["wheel_speed_dr_base_std_mps"].as<double>();
            if (system["wheel_speed_dr_stationary_std_mps"])
                wheel_speed_dr_config_.stationary_std_mps = system["wheel_speed_dr_stationary_std_mps"].as<double>();
            if (system["wheel_speed_dr_stationary_threshold_mps"])
                wheel_speed_dr_config_.stationary_speed_threshold_mps = system["wheel_speed_dr_stationary_threshold_mps"].as<double>();
            if (system["wheel_speed_dr_max_age_sec"])
                wheel_speed_dr_config_.max_age_sec = system["wheel_speed_dr_max_age_sec"].as<double>();
            if (system["wheel_speed_dr_future_tolerance_sec"])
                wheel_speed_dr_config_.future_tolerance_sec = system["wheel_speed_dr_future_tolerance_sec"].as<double>();
            if (system["wheel_speed_dr_max_acceleration_mps2"])
                wheel_speed_dr_config_.max_acceleration_mps2 = system["wheel_speed_dr_max_acceleration_mps2"].as<double>();
            if (system["wheel_speed_dr_max_abs_innovation_mps"])
                wheel_speed_dr_config_.max_abs_innovation_mps = system["wheel_speed_dr_max_abs_innovation_mps"].as<double>();
            if (system["wheel_speed_dr_nis_gate"])
                wheel_speed_dr_config_.normalized_innovation_squared_gate = system["wheel_speed_dr_nis_gate"].as<double>();
            if (system["wheel_speed_dr_torque_reference_nm"])
                wheel_speed_dr_config_.torque_reference_nm = system["wheel_speed_dr_torque_reference_nm"].as<double>();
            if (system["wheel_speed_dr_torque_std_scale"])
                wheel_speed_dr_config_.torque_std_scale = system["wheel_speed_dr_torque_std_scale"].as<double>();
            if (system["wheel_speed_dr_max_velocity_step_mps"])
                wheel_speed_dr_config_.max_velocity_step_mps = system["wheel_speed_dr_max_velocity_step_mps"].as<double>();
        }
        if (wheel_speed_dr_config_.base_std_mps <= 0.0 ||
            wheel_speed_dr_config_.stationary_std_mps <= 0.0 ||
            wheel_speed_dr_config_.stationary_speed_threshold_mps < 0.0 ||
            wheel_speed_dr_config_.max_age_sec <= 0.0 ||
            wheel_speed_dr_config_.future_tolerance_sec < 0.0 ||
            wheel_speed_dr_config_.max_acceleration_mps2 <= 0.0 ||
            wheel_speed_dr_config_.max_abs_innovation_mps <= 0.0 ||
            wheel_speed_dr_config_.normalized_innovation_squared_gate <= 0.0 ||
            wheel_speed_dr_config_.torque_reference_nm <= 0.0 ||
            wheel_speed_dr_config_.torque_std_scale < 0.0 ||
            wheel_speed_dr_config_.max_velocity_step_mps <= 0.0) {
            LOG(ERROR) << "invalid wheel-speed DR observation configuration";
            return false;
        }

    } catch (...) {
        LOG(ERROR) << "bad conversion";
        return false;
    }

    LOG(INFO) << "lidar_type " << lidar_type;
    if (lidar_type == 1) {
        preprocess_->SetLidarType(LidarType::AVIA);
        LOG(INFO) << "Using AVIA Lidar, point timestamp scale "
                  << preprocess_->LivoxPointTimeScale();
    } else if (lidar_type == 2) {
        preprocess_->SetLidarType(LidarType::VELO32);
        LOG(INFO) << "Using Velodyne 32 Lidar";
    } else if (lidar_type == 3) {
        preprocess_->SetLidarType(LidarType::OUST64);
        LOG(INFO) << "Using OUST 64 Lidar";
    } else if (lidar_type == 4) {
        preprocess_->SetLidarType(LidarType::ROBOSENSE);
        LOG(INFO) << "Using RoboSense Lidar";
    } else {
        LOG(WARNING) << "unknown lidar_type";
        return false;
    }

    if (ivox_nearby_type == 0) {
        ivox_options_.nearby_type_ = IVoxType::NearbyType::CENTER;
    } else if (ivox_nearby_type == 6) {
        ivox_options_.nearby_type_ = IVoxType::NearbyType::NEARBY6;
    } else if (ivox_nearby_type == 18) {
        ivox_options_.nearby_type_ = IVoxType::NearbyType::NEARBY18;
    } else if (ivox_nearby_type == 26) {
        ivox_options_.nearby_type_ = IVoxType::NearbyType::NEARBY26;
    } else {
        LOG(WARNING) << "unknown ivox_nearby_type, use NEARBY18";
        ivox_options_.nearby_type_ = IVoxType::NearbyType::NEARBY18;
    }

    filter_size_scan_ = filter_size_scan;
    voxel_scan_.setLeafSize(filter_size_scan_, filter_size_scan_, filter_size_scan_);

    offset_t_lidar_fixed_ = math::VecFromArray<double>(extrinT_);
    offset_R_lidar_fixed_ = math::MatFromArray<double>(extrinR_);

    p_imu_->SetExtrinsic(offset_t_lidar_fixed_, offset_R_lidar_fixed_);
    p_imu_->SetGyrCov(Vec3d(gyr_cov, gyr_cov, gyr_cov));
    p_imu_->SetAccCov(Vec3d(acc_cov, acc_cov, acc_cov));
    p_imu_->SetGyrBiasCov(Vec3d(b_gyr_cov, b_gyr_cov, b_gyr_cov));
    p_imu_->SetAccBiasCov(Vec3d(b_acc_cov, b_acc_cov, b_acc_cov));

    std::string multi_lidar_error;
    if (!LoadMultiLidarConfig(yaml, multi_lidar_config_, &multi_lidar_error)) {
        LOG(ERROR) << "invalid multi-lidar configuration: " << multi_lidar_error;
        return false;
    }
    if (multi_lidar_config_.enabled) {
        multi_lidar_assembler_.Reset(multi_lidar_config_);
        LOG(INFO) << "multi-lidar frontend enabled with " << multi_lidar_config_.lidars.size()
                  << " sensors, primary id=" << multi_lidar_config_.primary_lidar_id
                  << ", reorder_window=" << multi_lidar_config_.reorder_window;
    }
    adaptive_lidar_load_controller_.Reset(multi_lidar_config_);
    if (multi_lidar_config_.adaptive_load.enabled && enable_skip_lidar_) {
        LOG(WARNING) << "ignore fasterlio.skip_lidar_num=" << skip_lidar_num_
                     << " because adaptive multi-lidar load control is enabled; fixed skipping "
                        "would lengthen dead reckoning before exhausting point/lidar degradation";
        enable_skip_lidar_ = false;
        skip_lidar_cnt_ = 0;
    }

    std::string self_filter_error;
    if (!LoadSelfPointFilterConfig(yaml, self_point_filter_config_, &self_filter_error)) {
        LOG(ERROR) << "invalid self-point filter configuration: " << self_filter_error;
        return false;
    }
    if (self_point_filter_config_.enabled) {
        LOG(INFO) << "body-aligned self-point box enabled, min="
                  << self_point_filter_config_.min_body.transpose() << ", max="
                  << self_point_filter_config_.max_body.transpose();
    }

    const YAML::Node map_export = yaml["map_export"];
    std::string map_export_filter_error;
    if (map_export && !LoadSelfPointFilterConfig(
                          map_export, map_export_self_point_filter_config_,
                          &map_export_filter_error)) {
        LOG(ERROR) << "invalid map-export self-point filter configuration: "
                   << map_export_filter_error;
        return false;
    }
    if (map_export_self_point_filter_config_.enabled) {
        LOG(INFO) << "map-export self-point box enabled, min="
                  << map_export_self_point_filter_config_.min_body.transpose()
                  << ", max="
                  << map_export_self_point_filter_config_.max_body.transpose();
    }

    const YAML::Node noise = yaml["lidar_noise_model"];
    if (noise) {
        point_noise_enabled_ = noise["enabled"] ? noise["enabled"].as<bool>() : false;
        if (noise["range_sigma"]) range_noise_sigma_ = noise["range_sigma"].as<double>();
        if (noise["angular_sigma_deg"]) {
            angular_noise_sigma_rad_ = noise["angular_sigma_deg"].as<double>() * M_PI / 180.0;
        }
        if (noise["reference_sigma"]) noise_reference_sigma_ = noise["reference_sigma"].as<double>();
        if (noise["min_information_scale"]) {
            min_information_scale_ = noise["min_information_scale"].as<double>();
        }
        if (noise["max_information_scale"]) {
            max_information_scale_ = noise["max_information_scale"].as<double>();
        }
        if (point_noise_enabled_ &&
            (range_noise_sigma_ <= 0.0 || angular_noise_sigma_rad_ <= 0.0 || noise_reference_sigma_ <= 0.0 ||
             min_information_scale_ <= 0.0 || max_information_scale_ < min_information_scale_)) {
            LOG(ERROR) << "invalid lidar_noise_model parameters";
            return false;
        }
    }
    return true;
}

LaserMapping::LaserMapping(Options options) : options_(options) {
    preprocess_.reset(new PointCloudPreprocess());
    p_imu_.reset(new ImuProcess());
}

void LaserMapping::ProcessWheelSpeed(double timestamp, double longitudinal_speed_mps,
                                     double motor_torque_nm) {
    std::lock_guard<std::mutex> lock(wheel_speed_mutex_);
    ++wheel_speed_dr_stats_.input_count;
    wheel_speed_dr_stats_.last_measurement_mps = longitudinal_speed_mps;
    if (!wheel_speed_dr_config_.enabled) return;
    if (!std::isfinite(timestamp) || timestamp <= 0.0 ||
        !std::isfinite(longitudinal_speed_mps) || !std::isfinite(motor_torque_nm)) {
        ++wheel_speed_dr_stats_.invalid_input_count;
        return;
    }
    if (timestamp <= last_raw_wheel_speed_timestamp_) {
        ++wheel_speed_dr_stats_.timestamp_reject_count;
        return;
    }

    bool acceleration_valid = true;
    if (last_raw_wheel_speed_timestamp_ > std::numeric_limits<double>::lowest()) {
        const double dt = timestamp - last_raw_wheel_speed_timestamp_;
        if (dt > 1e-3) {
            const double acceleration =
                std::abs(longitudinal_speed_mps - last_raw_wheel_speed_mps_) / dt;
            acceleration_valid =
                acceleration <= wheel_speed_dr_config_.max_acceleration_mps2;
        }
    }
    last_raw_wheel_speed_timestamp_ = timestamp;
    last_raw_wheel_speed_mps_ = longitudinal_speed_mps;
    if (!acceleration_valid) {
        ++wheel_speed_dr_stats_.acceleration_reject_count;
        return;
    }

    wheel_speed_buffer_.push_back({timestamp, longitudinal_speed_mps, motor_torque_nm});
    constexpr double kBufferRetentionSec = 2.0;
    while (!wheel_speed_buffer_.empty() &&
           timestamp - wheel_speed_buffer_.front().timestamp > kBufferRetentionSec) {
        wheel_speed_buffer_.pop_front();
    }
}

WheelSpeedDrStats LaserMapping::GetWheelSpeedDrStats() const {
    std::lock_guard<std::mutex> lock(wheel_speed_mutex_);
    return wheel_speed_dr_stats_;
}

bool LaserMapping::ApplyWheelSpeedObservation(
    ESKF& filter, double state_timestamp,
    double& last_applied_observation_timestamp, bool high_frequency_filter) {
    WheelSpeedSample observation;
    WheelSpeedDrConfig config;
    bool found = false;
    {
        std::lock_guard<std::mutex> lock(wheel_speed_mutex_);
        if (!wheel_speed_dr_config_.enabled || !std::isfinite(state_timestamp)) return false;
        config = wheel_speed_dr_config_;
        for (const auto& sample : wheel_speed_buffer_) {
            if (sample.timestamp <= last_applied_observation_timestamp + 1e-9) continue;
            if (sample.timestamp > state_timestamp + config.future_tolerance_sec) break;
            observation = sample;
            found = true;
        }
        if (!found) return false;
        last_applied_observation_timestamp = observation.timestamp;
        if (state_timestamp - observation.timestamp > config.max_age_sec) {
            if (high_frequency_filter) {
                ++wheel_speed_dr_stats_.imu_filter_stale_count;
            } else {
                ++wheel_speed_dr_stats_.lidar_filter_stale_count;
            }
            return false;
        }
    }

    const double base_std =
        std::abs(observation.speed_mps) <= config.stationary_speed_threshold_mps
            ? config.stationary_std_mps
            : config.base_std_mps;
    const double torque_ratio =
        std::min(1.0, std::abs(observation.torque_nm) / config.torque_reference_nm);
    const double standard_deviation =
        base_std * (1.0 + config.torque_std_scale * torque_ratio);
    const auto result = filter.UpdateBodyForwardSpeed(
        observation.speed_mps, standard_deviation * standard_deviation,
        config.max_abs_innovation_mps,
        config.normalized_innovation_squared_gate,
        config.max_velocity_step_mps);

    {
        std::lock_guard<std::mutex> lock(wheel_speed_mutex_);
        wheel_speed_dr_stats_.last_measurement_mps = observation.speed_mps;
        wheel_speed_dr_stats_.last_predicted_mps = result.predicted_speed_mps;
        wheel_speed_dr_stats_.last_innovation_mps = result.innovation_mps;
        wheel_speed_dr_stats_.last_standard_deviation_mps = standard_deviation;
        wheel_speed_dr_stats_.last_normalized_innovation_squared =
            result.normalized_innovation_squared;
        if (high_frequency_filter) {
            if (result.accepted) ++wheel_speed_dr_stats_.imu_filter_accepted_count;
            else ++wheel_speed_dr_stats_.imu_filter_rejected_count;
        } else {
            if (result.accepted) ++wheel_speed_dr_stats_.lidar_filter_accepted_count;
            else ++wheel_speed_dr_stats_.lidar_filter_rejected_count;
        }
    }
    return result.accepted;
}

void LaserMapping::ResetWheelSpeedIntegrationBridge(double timestamp) {
    std::lock_guard<std::mutex> lock(wheel_speed_mutex_);
    last_lidar_filter_wheel_timestamp_ = timestamp;
    last_imu_filter_wheel_timestamp_ = timestamp;
}

void LaserMapping::ProcessIMU(const lightning::IMUPtr &imu) {
    publish_count_++;

    double timestamp = imu->timestamp;

    UL lock(mtx_buffer_);
    if (timestamp < last_timestamp_imu_) {
        LOG(WARNING) << "imu loop back, clear buffer";
        imu_buffer_.clear();
        std::lock_guard<std::mutex> wheel_lock(wheel_speed_mutex_);
        wheel_speed_buffer_.clear();
        last_raw_wheel_speed_timestamp_ = std::numeric_limits<double>::lowest();
        last_lidar_filter_wheel_timestamp_ = std::numeric_limits<double>::lowest();
        last_imu_filter_wheel_timestamp_ = std::numeric_limits<double>::lowest();
    }

    if (p_imu_->IsIMUInited()) {
        /// 更新最新imu状态
        const Vec3d acc = p_imu_->ScaleAccelerationForPrediction(imu->linear_acceleration);
        kf_imu_.Predict(timestamp - last_timestamp_imu_, p_imu_->Q_, imu->angular_velocity, acc);
        ApplyWheelSpeedObservation(kf_imu_, timestamp,
                                   last_imu_filter_wheel_timestamp_, true);

        // LOG(INFO) << "newest wrt lidar: " << timestamp - kf_.GetX().timestamp_;

        /// 更新ui
        if (ui_) {
            ui_->UpdateNavState(kf_imu_.GetX());
        }
    }

    last_timestamp_imu_ = timestamp;

    imu_buffer_.emplace_back(imu);
}

/**
 * @brief 处理一帧已经缓存并完成时间同步的Lidar数据。
 *
 * Run()是LIO前端的主循环单步：
 * 1. 从缓存中同步出一帧Lidar和覆盖扫描周期的IMU；
 * 2. 用IMU预测ESKF并对点云去畸变；
 * 3. 首帧直接建初始局部地图，后续帧进入降采样和地图匹配；
 * 4. 通过Lidar观测更新ESKF；
 * 5. 根据运动幅度创建关键帧、维护局部地图和高频IMU显示状态。
 *
 * @return true表示本帧完成了有效前端处理，false表示同步失败、初始化等待、跳帧或点数不足。
 */
bool LaserMapping::Run() { return RunDetailed() == RunStatus::kOutput; }

LaserMapping::RunStatus LaserMapping::RunDetailed() {
    const bool profiling_enabled = profiling::ComputeProfilingEnabled();
    profiling::Stopwatch outer_profile_timer(profiling_enabled);
    profiling::Stopwatch sync_profile_timer(profiling_enabled);
    // SyncPackages()只在IMU已经覆盖当前Lidar扫描结束时间时才会成功。
    // 因此这里处理的可能不是最新进入缓存的点云，而是第一帧已经等到足够IMU的点云。
    if (!SyncPackages()) {
        return RunStatus::kNoData;
    }
    const profiling::TimingSample sync_timing = sync_profile_timer.Stop();

    if (measures_.imu_.empty()) {
        ++pre_imu_drop_count_;
        last_tracking_healthy_ = false;
        return RunStatus::kConsumed;
    }

    const double lidar_gap =
        last_lidar_time_ > 0.0 ? measures_.lidar_begin_time_ - last_lidar_time_ : 0.0;
    if (lidar_gap > 0.5) {
        LOG(ERROR) << "检测到雷达断流，时长：" << lidar_gap
                   << "; skip discontinuous frame and reset IMU integration bridge";
        auto safe_state = kf_.GetX();
        safe_state.timestamp_ = measures_.lidar_end_time_;
        kf_.ChangeX(safe_state);
        kf_imu_ = kf_;
        ResetWheelSpeedIntegrationBridge(measures_.lidar_end_time_);
        p_imu_->ResetIntegrationBridge(measures_.imu_.back(), measures_.lidar_end_time_, safe_state);
        last_lidar_time_ = measures_.lidar_begin_time_;
        last_tracking_healthy_ = false;
        return RunStatus::kConsumed;
    }
    last_lidar_time_ = measures_.lidar_begin_time_;

    using BenchClock = std::chrono::steady_clock;
    const auto elapsed_ms = [](const BenchClock::time_point &start) {
        return std::chrono::duration<double, std::milli>(BenchClock::now() - start).count();
    };
    const std::size_t input_points = measures_.scan_ ? measures_.scan_->size() : 0;
    std::size_t self_filter_removed = 0;
    double imu_undistort_ms = 0.0;
    double downsample_ms = 0.0;
    double match_setup_ms = 0.0;
    double scan_match_ms = 0.0;
    double map_update_ms = 0.0;
    profiling::TimingSample allocation_timing;
    profiling::TimingSample imu_timing;
    profiling::TimingSample self_filter_timing;
    profiling::TimingSample selection_timing;
    profiling::TimingSample downsample_timing;
    profiling::TimingSample match_setup_timing;
    profiling::TimingSample scan_match_timing;
    profiling::TimingSample map_update_timing;
    profiling::TimingSample adaptive_timing;
    profiling::TimingSample post_imu_timing;
    const auto emit_pipeline_benchmark = [&](const char* phase, std::size_t output_points) {
        if (!profiling_enabled) return;
        const profiling::TimingSample outer_timing = outer_profile_timer.Stop();
        LOG(INFO) << std::fixed << std::setprecision(6)
                  << "COMPUTE_BENCH_FRAME module=lio_pipeline phase=" << phase
                  << " timestamp_s=" << measures_.lidar_end_time_
                  << " input_points=" << input_points
                  << " output_points=" << output_points
                  << " imu_samples=" << measures_.imu_.size()
                  << " buffered_imu_samples=" << imu_buffer_.size()
                  << " self_filter_removed=" << self_filter_removed
                  << ' ' << profiling::FormatTimingSample("sync", sync_timing)
                  << ' ' << profiling::FormatTimingSample("allocation", allocation_timing)
                  << ' ' << profiling::FormatTimingSample("imu_undistort", imu_timing)
                  << ' ' << profiling::FormatTimingSample("self_filter", self_filter_timing)
                  << ' ' << profiling::FormatTimingSample("lidar_selection", selection_timing)
                  << ' ' << profiling::FormatTimingSample("downsample", downsample_timing)
                  << ' ' << profiling::FormatTimingSample("match_setup", match_setup_timing)
                  << ' ' << profiling::FormatTimingSample("scan_match", scan_match_timing)
                  << ' ' << profiling::FormatTimingSample("map_update", map_update_timing)
                  << ' ' << profiling::FormatTimingSample("adaptive_control", adaptive_timing)
                  << ' ' << profiling::FormatTimingSample("post_imu_catchup", post_imu_timing)
                  << ' ' << profiling::FormatTimingSample("outer", outer_timing);
    };
    const auto emit_benchmark = [&](const char *phase, std::size_t output_points) {
        const double core_update_ms =
            imu_undistort_ms + downsample_ms + match_setup_ms + scan_match_ms + map_update_ms;
        if (profiling_enabled) {
            LOG(INFO) << std::fixed << std::setprecision(6)
                      << "LIO_BENCH_FRAME method=lightning_lm phase=" << phase
                      << " timestamp_s=" << measures_.lidar_end_time_
                      << " preprocess_ms=" << current_preprocess_ms_
                      << " imu_undistort_ms=" << imu_undistort_ms
                      << " downsample_ms=" << downsample_ms
                      << " match_setup_ms=" << match_setup_ms
                      << " scan_match_ms=" << scan_match_ms
                      << " map_update_ms=" << map_update_ms
                      << " core_update_ms=" << core_update_ms
                      << " total_ms=" << (current_preprocess_ms_ + core_update_ms)
                      << " input_points=" << input_points
                      << " self_filter_removed=" << self_filter_removed
                      << " output_points=" << output_points
                      << " adaptive_step=" << adaptive_lidar_load_controller_.DegradationStep()
                      << " selected_lidars=" << current_lidar_selection_.lidar_ids.size()
                      << " point_stride=" << current_lidar_selection_.point_stride
                      << " lidar_latency_ms=" << last_lidar_latency_sec_ * 1e3;
        }
        last_frame_processing_ms_ = current_preprocess_ms_ + core_update_ms;
    };

    // IMU处理包含两种情况：
    // - 初始化未完成：继续累计IMU均值/方差，并直接返回空点云；
    // - 初始化完成：预测kf_到当前扫描结束时刻，并把点云补偿到扫描结束时刻。
    // Keyframes keep the previous cloud pointer, so allocate a new output instead of clearing it in place.
    profiling::Stopwatch allocation_profile_timer(profiling_enabled);
    scan_undistort_full_.reset(new PointCloudType());
    allocation_timing = allocation_profile_timer.Stop();
    profiling::Stopwatch imu_profile_timer(profiling_enabled);
    const auto imu_start = BenchClock::now();
    p_imu_->Process(measures_, kf_, scan_undistort_full_);
    imu_undistort_ms = elapsed_ms(imu_start);
    imu_timing = imu_profile_timer.Stop();

    profiling::Stopwatch self_filter_profile_timer(profiling_enabled);
    if (scan_undistort_full_ && !scan_undistort_full_->empty() && self_point_filter_config_.enabled) {
        self_filter_removed = FilterSelfPoints(
            *scan_undistort_full_, self_point_filter_config_, GetInitialLidarRotation().matrix());
    }
    self_filter_timing = self_filter_profile_timer.Stop();

    if (!scan_undistort_full_ || scan_undistort_full_->empty()) {
        LOG(WARNING) << "No point, skip this scan!";
        last_tracking_healthy_ = false;
        emit_pipeline_benchmark("empty_scan", 0);
        return RunStatus::kConsumed;
    }
    last_lidar_latency_sec_ = latest_input_sensor_timestamp_ > 0.0
                                  ? std::max(0.0, latest_input_sensor_timestamp_ - measures_.lidar_end_time_)
                                  : 0.0;
    if (adaptive_lidar_load_controller_.IsHardStale(last_lidar_latency_sec_)) {
        ++adaptive_stale_drop_count_;
        last_frame_processing_ms_ = current_preprocess_ms_ + imu_undistort_ms;
        last_tracking_healthy_ = false;
        LOG(WARNING) << "drop stale lidar correction before matching: age="
                     << last_lidar_latency_sec_ << " sec, hard_deadline="
                     << multi_lidar_config_.adaptive_load.hard_latency_sec;
        emit_pipeline_benchmark("hard_stale", 0);
        return RunStatus::kConsumed;
    }
    profiling::Stopwatch selection_profile_timer(profiling_enabled);
    if (multi_lidar_config_.enabled) {
        current_lidar_selection_ = adaptive_lidar_load_controller_.Select(current_lidar_stats_);
        scan_undistort_ = SelectLidarPoints(scan_undistort_full_, current_lidar_selection_);
        if (!scan_undistort_ || scan_undistort_->empty()) {
            last_tracking_healthy_ = false;
            LOG_EVERY_N(WARNING, 20)
                << "skip lidar frame: available sources do not satisfy current localization minimum";
            selection_timing = selection_profile_timer.Stop();
            emit_pipeline_benchmark("empty_selection", 0);
            return RunStatus::kConsumed;
        }
    } else {
        scan_undistort_ = scan_undistort_full_;
        current_lidar_selection_.lidar_ids = {0};
        current_lidar_selection_.point_stride = 1;
        current_lidar_selection_.degradation_step = 0;
    }
    selection_timing = selection_profile_timer.Stop();

    // 第一帧没有可匹配的局部地图，因此不做ESKF观测更新，直接把去畸变点云转到世界系作为初始地图。
    if (flg_first_scan_) {
        profiling::Stopwatch initial_map_profile_timer(profiling_enabled);
        const auto initial_map_start = BenchClock::now();
        LOG(INFO) << "first scan pts: " << scan_undistort_->size();

        state_point_ = kf_.GetX();
        scan_down_world_->resize(scan_undistort_->size());
        for (int i = 0; i < scan_undistort_->size(); i++) {
            PointBodyToWorld(scan_undistort_->points[i], scan_down_world_->points[i]);
        }
        ivox_->AddPoints(scan_down_world_->points);

        // 记录第一帧时间，后续用INIT_TIME判断地图是否足够稳定，可以开始正常观测更新/地图筛选。
        first_lidar_time_ = measures_.lidar_end_time_;
        state_point_.timestamp_ = lidar_end_time_;
        flg_first_scan_ = false;
        last_tracking_healthy_ = true;
        map_update_ms = elapsed_ms(initial_map_start);
        map_update_timing = initial_map_profile_timer.Stop();
        emit_benchmark("initialization", scan_undistort_->size());
        emit_pipeline_benchmark("initialization", scan_undistort_->size());
        return RunStatus::kOutput;
    }

    // 可选跳帧：仍然完成了IMU预测和去畸变，但直接返回，跳过后续降采样、地图匹配、ESKF观测更新、关键帧判断和地图更新。
    // 这样UI可以继续显示预测位姿，同时降低前端计算负载。
    if (enable_skip_lidar_) {
        skip_lidar_cnt_++;
        skip_lidar_cnt_ = skip_lidar_cnt_ % skip_lidar_num_;

        if (skip_lidar_cnt_ != 0) {
            /// 更新UI中的内容
            if (ui_) {
                ui_->UpdateNavState(kf_.GetX());
                ui_->UpdateScan(scan_undistort_, kf_.GetX().GetPose());
            }

            emit_pipeline_benchmark("fixed_skip", scan_undistort_->size());
            return RunStatus::kConsumed;
        }
    }

    if (!profiling::ReduceNonessentialOverhead()) {
        LOG(INFO) << "=============================";
        LOG(INFO) << "LIO get cloud at beg: " << std::setprecision(14)
                  << measures_.lidar_begin_time_ << ", end: " << measures_.lidar_end_time_;
    }

    // 初始若干秒内地图还很稀疏，ObsModel和MapIncremental会根据这个标志放宽部分逻辑。
    flg_EKF_inited_ = (measures_.lidar_begin_time_ - first_lidar_time_) >= fasterlio::INIT_TIME;

    // 对当前去畸变点云降采样，后续匹配和建图都使用scan_down_lidar_，避免逐点处理原始大点云。
    profiling::Stopwatch downsample_profile_timer(profiling_enabled);
    const auto downsample_start = BenchClock::now();
    if (multi_lidar_config_.enabled) {
        scan_down_lidar_ = DownsamplePreservingSource(scan_undistort_, filter_size_scan_);
    } else {
        voxel_scan_.setInputCloud(scan_undistort_);
        voxel_scan_.filter(*scan_down_lidar_);
    }

    // if (options_.proj_kfs_) {
    //     ProjectKFs();
    // }

    int cur_pts = scan_down_lidar_->size();

    // 如果yaml配置的体素分辨率导致点数过少，则临时用0.1m重新降采样。
    // 这是一个保护分支：保证观测模型至少有足够点参与匹配。
    if (cur_pts < (scan_undistort_->size() * 0.1) || cur_pts < options_.min_pts) {
        /// 降采样太狠了,有效点数不够，用0.1分辨率代替
        // LOG(INFO) << "too few points, using 0.1 resol";
        if (multi_lidar_config_.enabled) {
            scan_down_lidar_ = DownsamplePreservingSource(scan_undistort_, 0.1);
        } else {
            auto v = voxel_scan_;
            v.setLeafSize(0.1, 0.1, 0.1);
            v.setInputCloud(scan_undistort_);
            v.filter(*scan_down_lidar_);
        }

        // LOG(INFO) << "Now pts: " << scan_down_lidar_->size() << ", before: " << cur_pts;
        cur_pts = scan_down_lidar_->size();
    }
    downsample_ms = elapsed_ms(downsample_start);
    downsample_timing = downsample_profile_timer.Stop();

    // 极端情况下仍然点数不足，继续匹配会让最近邻和平面拟合没有意义，直接跳过。
    if (cur_pts < 5) {
        LOG(WARNING) << "Too few points, skip this scan!" << scan_undistort_->size() << ", "
                     << scan_down_lidar_->size();
        last_tracking_healthy_ = false;
        emit_pipeline_benchmark("too_few_downsampled", scan_down_lidar_->size());
        return RunStatus::kConsumed;
    }

    profiling::Stopwatch match_setup_profile_timer(profiling_enabled);
    const auto match_setup_start = BenchClock::now();
    scan_down_world_->resize(cur_pts);
    nearest_points_.resize(cur_pts);

    // 按当前帧点数预分配观测模型缓存，ObsModel()中会复用这些数组保存最近邻、残差和平面。
    residuals_.resize(cur_pts, 0);
    point_selected_surf_.resize(cur_pts, 1);
    point_selected_icp_.resize(cur_pts, 1);
    plane_coef_.resize(cur_pts, Vec4f::Zero());

    // 保存预测状态，后面用来统计Lidar观测更新带来的位姿修正量。
    auto pred_state = kf_.GetX();
    match_setup_ms = elapsed_ms(match_setup_start);
    match_setup_timing = match_setup_profile_timer.Stop();
    // pred_state.pos_ = state_point_.pos_;  // 假定位置不动行不行,防止速度漂移
    // kf_.ChangeX(pred_state);

    // Lidar观测更新：ESKF内部会多次调用ObsModel()，构造点面/点点残差的HTH和HTr。
    profiling::Stopwatch scan_match_profile_timer(profiling_enabled);
    const auto scan_match_start = BenchClock::now();
    kf_.Update(ESKF::ObsType::LIDAR, 1.0);

    // 更新当前Lidar帧结束时刻的前端状态，供建图、关键帧和外部查询使用。
    state_point_ = kf_.GetX();
    state_point_.timestamp_ = measures_.lidar_end_time_;
    last_tracking_healthy_ = kf_.LastUpdateAccepted() && effect_feat_surf_ >= 20 && current_max_imu_gap_ <= 0.25 &&
                             state_point_.pos_.allFinite() && state_point_.rot_.matrix().allFinite();

    // 统计本次观测更新相对IMU预测的修正量，主要用于日志和异常诊断。
    const double delta_translation = (pred_state.pos_ - state_point_.pos_).norm();
    const double delta_rotation_deg = (pred_state.rot_.inverse() * state_point_.rot_).log().norm() * 180.0 / M_PI;
    const double delta_velocity = (pred_state.vel_ - state_point_.vel_).norm();

    const double current_speed = state_point_.vel_.norm();

    if (adaptive_velocity_propagation_ && kf_.LastUpdateAccepted()) {
        if (!velocity_innovation_initialized_) {
            velocity_innovation_ema_ = delta_translation;
            velocity_innovation_initialized_ = true;
        } else {
            velocity_innovation_ema_ =
                (1.0 - velocity_innovation_ema_alpha_) * velocity_innovation_ema_ +
                velocity_innovation_ema_alpha_ * delta_translation;
        }
        if (velocity_propagation_cooldown_remaining_ > 0) {
            --velocity_propagation_cooldown_remaining_;
        }

        bool requested = velocity_propagation_active_;
        const char *transition_reason = "innovation threshold";
        if (requested) {
            ++velocity_propagation_active_updates_;
            if (velocity_innovation_ema_ < velocity_innovation_disable_threshold_) {
                requested = false;
            } else if (velocity_propagation_max_active_updates_ > 0 &&
                       velocity_propagation_active_updates_ >=
                           velocity_propagation_max_active_updates_) {
                requested = false;
                transition_reason = "active-update safety limit";
                velocity_propagation_cooldown_remaining_ =
                    velocity_propagation_cooldown_updates_;
                velocity_propagation_safety_lockout_ = true;
            }
        } else {
            velocity_propagation_active_updates_ = 0;
            if (velocity_propagation_safety_lockout_ &&
                velocity_innovation_ema_ < velocity_innovation_disable_threshold_) {
                velocity_propagation_safety_lockout_ = false;
            }
            if (!velocity_propagation_safety_lockout_ &&
                velocity_propagation_cooldown_remaining_ == 0 &&
                velocity_innovation_ema_ > velocity_innovation_enable_threshold_) {
                requested = true;
            }
        }
        if (requested != velocity_propagation_active_) {
            velocity_propagation_active_ = requested;
            kf_.SetPropagateVelocity(requested);
            LOG(WARNING) << "Adaptive velocity propagation " << (requested ? "enabled" : "disabled")
                         << ", reason: " << transition_reason
                         << ", lidar innovation EMA: " << velocity_innovation_ema_;
        }
    }

    scan_match_ms = elapsed_ms(scan_match_start);
    scan_match_timing = scan_match_profile_timer.Stop();

    if (!profiling::ReduceNonessentialOverhead()) {
        LOG(INFO) << "[ mapping ]: In num: " << scan_undistort_->points.size()
                  << " down " << cur_pts << " Map grid num: " << ivox_->NumValidGrids()
                  << " effect num : " << effect_feat_surf_ << ", " << effect_feat_icp_;
        LOG(INFO) << "delta trans: " << (pred_state.pos_ - state_point_.pos_).transpose()
                  << ", ang: " << delta_rotation_deg;
    }
    // LOG(INFO) << "P diag: " << kf_.GetP().diagonal().transpose();

    // Vec3d v_from_last = (state_point_.pos_ - last_state.pos_) / (state_point_.timestamp_ - last_state.timestamp_);
    // LOG(INFO) << "v from last: " << v_from_last.transpose();

    // if (delta_velocity > 1.0 || current_speed > 4.0) {
    //     LOG(ERROR) << "detected very large vel change, last: " << last_state.vel_.transpose()
    //                << ", pred: " << pred_state.vel_.transpose() << ", cur:" << state_point_.vel_.transpose();
    //     LOG(ERROR) << "please check";
    // }

    /// keyframes - 智能关键帧创建决策
    // 只有创建关键帧时才会调用MakeKF()，而MakeKF()内部会把当前帧点云增量加入IVox地图。
    // 因此关键帧阈值也间接控制了局部地图更新频率。
    profiling::Stopwatch map_update_profile_timer(profiling_enabled);
    const auto map_update_start = BenchClock::now();
    if (last_kf_ == nullptr) {
        MakeKF();  // 第一个关键帧：直接创建
    } else {
        SE3 last_pose = last_kf_->GetLIOPose();  // 上一关键帧的LIO位姿
        SE3 cur_pose = state_point_.GetPose();   // 当前帧的LIO位姿

        // 条件1：空间变化足够大（运动显著）
        if ((last_pose.translation() - cur_pose.translation()).norm() > options_.kf_dis_th_ ||    // 平移距离超过阈值
            (last_pose.so3().inverse() * cur_pose.so3()).log().norm() > options_.kf_angle_th_) {  // 旋转角度超过阈值
            MakeKF();
        } else if ((last_pose.so3().inverse() * cur_pose.so3()).log().norm() > 1.0 * M_PI / 180.0) {
            // MapIncremental();
        }
        // 条件2：时间间隔过长（定位模式下的时间保险）
        else if (!options_.is_in_slam_mode_ && (state_point_.timestamp_ - last_kf_->GetState().timestamp_) > 2.0) {
            MakeKF();  // 非SLAM模式下，超过2秒强制创建关键帧，防止长时间无关键帧
        }
    }
    map_update_ms = elapsed_ms(map_update_start);
    map_update_timing = map_update_profile_timer.Stop();
    emit_benchmark("tracking", scan_down_lidar_->size());
    profiling::Stopwatch adaptive_profile_timer(profiling_enabled);
    const int old_adaptive_step = adaptive_lidar_load_controller_.DegradationStep();
    adaptive_lidar_load_controller_.Observe(
        last_frame_processing_ms_ * 1e-3, last_lidar_latency_sec_, last_tracking_healthy_);
    const int new_adaptive_step = adaptive_lidar_load_controller_.DegradationStep();
    if (new_adaptive_step != old_adaptive_step) {
        LOG(WARNING) << "adaptive lidar load step " << old_adaptive_step << " -> "
                     << new_adaptive_step << ", processing_ms=" << last_frame_processing_ms_
                     << ", latency_ms=" << last_lidar_latency_sec_ * 1e3
                     << ", tracking_healthy=" << last_tracking_healthy_;
    }
    adaptive_timing = adaptive_profile_timer.Stop();

    // 维护一份“最新IMU时刻”的ESKF状态给UI显示。
    // kf_只到当前Lidar结束时刻；imu_buffer_中可能还有更晚的IMU，所以从kf_继续预测到最新IMU。
    profiling::Stopwatch post_imu_profile_timer(profiling_enabled);
    kf_imu_ = kf_;
    {
        std::lock_guard<std::mutex> lock(wheel_speed_mutex_);
        last_imu_filter_wheel_timestamp_ = last_lidar_filter_wheel_timestamp_;
    }
    if (!measures_.imu_.empty()) {
        double t = measures_.imu_.back()->timestamp;
        for (auto &imu : imu_buffer_) {
            double dt = imu->timestamp - t;
            // 这里做高频显示预测，不参与Lidar帧的去畸变输出。
            const Vec3d acc = p_imu_->ScaleAccelerationForPrediction(imu->linear_acceleration);
            kf_imu_.Predict(dt, p_imu_->Q_, imu->angular_velocity, acc);
            ApplyWheelSpeedObservation(kf_imu_, imu->timestamp,
                                       last_imu_filter_wheel_timestamp_, true);
            t = imu->timestamp;
        }
    }
    post_imu_timing = post_imu_profile_timer.Stop();

    if (ui_) {
        // 显示当前帧降采样点云和Lidar更新后的位姿。
        ui_->UpdateScan(scan_down_lidar_, state_point_.GetPose());
    }

    if (!profiling::ReduceNonessentialOverhead()) {
        LOG(INFO) << "LIO state: " << state_point_.pos_.transpose() << ", yaw "
                  << state_point_.rot_.angleZ<double>() * 180 / M_PI
                  << ", vel: " << state_point_.vel_.transpose()
                  << ", bg: " << state_point_.bg_.transpose()
                  << ", ba: " << state_point_.ba_.transpose()
                  << ", grav: " << state_point_.grav_.transpose()
                  << ", grav norm: " << state_point_.grav_.norm();
    }

    emit_pipeline_benchmark("tracking", scan_down_lidar_->size());

    return RunStatus::kOutput;
}

void LaserMapping::ProjectKFs(CloudPtr cloud, int size_limit) {
    auto state = kf_.GetX();
    SE3 pose_cur(state.rot_, state.pos_);
    pose_cur = pose_cur.inverse();

    for (auto kf : proj_kfs_) {
        // LOG(INFO) << "projecting kf: " << kf->GetID();
        // if (last_kf_) {
        // auto kf = last_kf_;
        SE3 pose = pose_cur * kf->GetLIOPose();

        int cnt = 0;
        for (auto &pt : kf->GetCloud()->points) {
            Vec3d p = pose * ToVec3d(pt);
            PointType pcl_pt;

            pcl_pt.x = p.x();
            pcl_pt.y = p.y();
            pcl_pt.z = p.z();
            pcl_pt.intensity = pt.intensity;

            cloud->push_back(pcl_pt);
            cnt++;

            if (cnt > size_limit) {
                break;
            }
        }
        // }
    }
}

/**
 * @brief 创建当前Lidar帧对应的关键帧，并维护关键帧列表与局部地图。
 * @details 关键帧会继承当前ESKF状态；优化位姿优先沿用上一关键帧的优化位姿递推。
 *          SLAM模式下额外保存到全局关键帧序列，并将当前帧增量加入IVox局部地图。
 */
void LaserMapping::MakeKF() {
    // 创建关键帧对象，包含ID、点云和当前状态
    Keyframe::Ptr kf = std::make_shared<Keyframe>(kf_id_++, scan_undistort_, state_point_);

    if (last_kf_) {
        /// opt pose 用之前的递推
        SE3 delta = last_kf_->GetLIOPose().inverse() * kf->GetLIOPose();
        // 基于上一关键帧的优化位姿递推当前关键帧的优化位姿
        // TODO. 这里为什么不用LIO前端里程计 kf->SetOptPose(kf->GetLIOPose())
        kf->SetOptPose(last_kf_->GetOptPose() * delta);
    } else {
        // 第一个关键帧：优化位姿直接使用LIO位姿
        kf->SetOptPose(kf->GetLIOPose());
    }

    // 设置关键帧的ESKF状态信息
    kf->SetState(state_point_);

    // 记录关键帧创建信息
    if (!profiling::ReduceNonessentialOverhead()) {
        LOG(INFO) << "LIO: create kf " << kf->GetID()
                  << ", state: " << state_point_.pos_.transpose()
                  << ", kf opt pose: " << kf->GetOptPose().translation().transpose()
                  << ", lio pose: " << kf->GetLIOPose().translation().transpose()
                  << ", time: " << std::setprecision(14) << state_point_.timestamp_;
    }

    // 只在SLAM模式下保存关键帧到列表
    if (options_.is_in_slam_mode_) {
        all_keyframes_.emplace_back(kf);
    }

    // 更新最新关键帧指针
    last_kf_ = kf;

    // 有keyframes时更新local map
    Timer::Evaluate([&, this]() { MapIncremental(); }, "    Incremental Mapping");

    /// 更新project kfs
    if (proj_kfs_.size() >= options_.max_proj_kfs_) {
        auto last = proj_kfs_.back();

        SE3 delta = last->GetLIOPose().inverse() * kf->GetLIOPose();

        if (delta.translation().norm() < 3 || delta.so3().log().norm() < 20 / 180 * M_PI) {
            // proj_kfs_.pop_back();
        } else {
            proj_kfs_.pop_front();
            proj_kfs_.emplace_back(kf);
        }
    } else {
        proj_kfs_.emplace_back(kf);
    }

    // for (auto &kf : proj_kfs_) {
    //     LOG(INFO) << "proj kf: " << kf->GetID();
    // }
}

bool LaserMapping::EnqueueCloud(double timestamp, CloudPtr cloud, const MultiLidarFrameStats *stats,
                                double preprocess_ms) {
    if (!cloud || cloud->empty()) {
        return false;
    }
    if (timestamp < last_timestamp_lidar_) {
        LOG(ERROR) << "fused lidar timestamp loop back, drop frame: " << std::setprecision(14) << timestamp
                   << " < " << last_timestamp_lidar_;
        return false;
    }
    lidar_buffer_.push_back(std::move(cloud));
    time_buffer_.push_back(timestamp);
    preprocess_time_buffer_ms_.push_back(preprocess_ms);
    if (stats) {
        lidar_stats_buffer_.push_back(*stats);
    } else {
        MultiLidarFrameStats single;
        single.begin_time = timestamp;
        single.present_lidar_ids = {multi_lidar_config_.primary_lidar_id};
        single.points_by_lidar[multi_lidar_config_.primary_lidar_id] = lidar_buffer_.back()->size();
        single.merged_points = lidar_buffer_.back()->size();
        lidar_stats_buffer_.push_back(std::move(single));
    }
    last_timestamp_lidar_ = timestamp;
    return true;
}

bool LaserMapping::DrainAssembledFrames() {
    bool enqueued = false;
    FusedLidarFrame frame;
    while (multi_lidar_assembler_.PopReady(frame)) {
        enqueued = EnqueueCloud(frame.stats.begin_time, std::move(frame.cloud), &frame.stats) || enqueued;
    }
    return enqueued;
}

bool LaserMapping::ProcessPointCloud2(const sensor_msgs::msg::PointCloud2::SharedPtr &msg) {
    const int id = multi_lidar_config_.enabled ? multi_lidar_config_.primary_lidar_id : 0;
    return ProcessPointCloud2(msg, id);
}

bool LaserMapping::ProcessPointCloud2(const sensor_msgs::msg::PointCloud2::SharedPtr &msg, int lidar_id) {
    UL lock(mtx_buffer_);
    bool accepted = false;
    double preprocess_ms = 0.0;
    Timer::Evaluate(
        [&, this]() {
            ++scan_count_;
            const double timestamp = ToSec(msg->header.stamp);
            CloudPtr cloud(new PointCloudType());
            const auto preprocess_start = std::chrono::steady_clock::now();
            preprocess_->Process(msg, cloud);
            preprocess_ms = std::chrono::duration<double, std::milli>(
                                std::chrono::steady_clock::now() - preprocess_start)
                                .count();
            // PointCloudPreprocess stores per-point relative time, but not every
            // handler preserves the ROS header in the internal PCL cloud.  The
            // localization pipeline reads this header as an absolute nanosecond
            // timestamp, so propagate it explicitly at the input boundary.
            cloud->header.stamp = static_cast<std::uint64_t>(std::llround(timestamp * 1e9));
            if (multi_lidar_config_.enabled) {
                accepted = multi_lidar_assembler_.AddCloud(lidar_id, timestamp, cloud);
                DrainAssembledFrames();
            } else {
                accepted = EnqueueCloud(timestamp, cloud, nullptr, preprocess_ms);
            }
        },
        "Preprocess (Standard)");
    return accepted;
}

bool LaserMapping::ProcessPointCloud2(const livox_ros_driver2::msg::CustomMsg::SharedPtr &msg) {
    const int id = multi_lidar_config_.enabled ? multi_lidar_config_.primary_lidar_id : 0;
    return ProcessPointCloud2(msg, id);
}

bool LaserMapping::ProcessPointCloud2(const livox_ros_driver2::msg::CustomMsg::SharedPtr &msg, int lidar_id) {
    UL lock(mtx_buffer_);
    bool accepted = false;
    double preprocess_ms = 0.0;
    Timer::Evaluate(
        [&, this]() {
            ++scan_count_;
            const double timestamp = ToSec(msg->header.stamp);
            CloudPtr cloud(new PointCloudType());
            const auto preprocess_start = std::chrono::steady_clock::now();
            preprocess_->Process(msg, cloud);
            preprocess_ms = std::chrono::duration<double, std::milli>(
                                std::chrono::steady_clock::now() - preprocess_start)
                                .count();
            cloud->header.stamp = static_cast<std::uint64_t>(std::llround(timestamp * 1e9));
            if (multi_lidar_config_.enabled) {
                accepted = multi_lidar_assembler_.AddCloud(lidar_id, timestamp, cloud);
                DrainAssembledFrames();
            } else {
                accepted = EnqueueCloud(timestamp, cloud, nullptr, preprocess_ms);
            }
        },
        "Preprocess (Livox)");
    return accepted;
}

bool LaserMapping::ProcessPointCloud2(CloudPtr cloud) {
    const int id = multi_lidar_config_.enabled ? multi_lidar_config_.primary_lidar_id : 0;
    return ProcessPointCloud2(std::move(cloud), id);
}

bool LaserMapping::ProcessPointCloud2(CloudPtr cloud, int lidar_id) {
    UL lock(mtx_buffer_);
    if (!cloud) return false;
    const double timestamp = math::ToSec(cloud->header.stamp);
    if (multi_lidar_config_.enabled) {
        const bool accepted = multi_lidar_assembler_.AddCloud(lidar_id, timestamp, cloud);
        DrainAssembledFrames();
        return accepted;
    }
    return EnqueueCloud(timestamp, std::move(cloud));
}

void LaserMapping::FlushMultiLidar() {
    UL lock(mtx_buffer_);
    if (!multi_lidar_config_.enabled) return;
    multi_lidar_assembler_.Flush();
    DrainAssembledFrames();
}

bool LaserMapping::SyncPackages() {
    if (lidar_buffer_.empty() || imu_buffer_.empty()) {
        return false;
    }

    /*** push a lidar scan ***/
    if (!lidar_pushed_) {
        measures_.scan_ = lidar_buffer_.front();
        measures_.lidar_begin_time_ = time_buffer_.front();
        current_preprocess_ms_ =
            preprocess_time_buffer_ms_.empty() ? 0.0 : preprocess_time_buffer_ms_.front();
        current_lidar_stats_ = lidar_stats_buffer_.empty() ? MultiLidarFrameStats() : lidar_stats_buffer_.front();

        if (measures_.scan_->points.size() <= 1) {
            LOG(WARNING) << "Too few input point cloud!";
            lidar_end_time_ = measures_.lidar_begin_time_ + lidar_mean_scantime_;
        } else if (measures_.scan_->points.back().time / double(1000) < 0.5 * lidar_mean_scantime_) {
            lidar_end_time_ = measures_.lidar_begin_time_ + lidar_mean_scantime_;
        } else {
            scan_num_++;
            lidar_end_time_ = measures_.lidar_begin_time_ + measures_.scan_->points.back().time / double(1000);

            lidar_mean_scantime_ +=
                (measures_.scan_->points.back().time / double(1000) - lidar_mean_scantime_) / scan_num_;

            if ((lidar_end_time_ - measures_.lidar_begin_time_) > 5 * lo::lidar_time_interval) {
                /// timestamp 有异常
                lidar_end_time_ = measures_.lidar_begin_time_ + lo::lidar_time_interval;
                lidar_mean_scantime_ = lo::lidar_time_interval;
            }
        }

        lo::lidar_time_interval = lidar_mean_scantime_;

        // LOG(INFO) << "recompute lidar end time: " << std::setprecision(14) << lidar_end_time_;
        measures_.lidar_end_time_ = lidar_end_time_;
        current_lidar_stats_.begin_time = measures_.lidar_begin_time_;
        current_lidar_stats_.end_time = measures_.lidar_end_time_;
        lidar_pushed_ = true;
    }

    if (last_timestamp_imu_ < lidar_end_time_) {
        return false;
    }

    /*** push imu_ data, and pop from imu_ buffer ***/
    double imu_time = imu_buffer_.front()->timestamp;
    measures_.imu_.clear();
    while ((!imu_buffer_.empty()) && (imu_time < lidar_end_time_)) {
        imu_time = imu_buffer_.front()->timestamp;
        if (imu_time > lidar_end_time_) {
            break;
        }

        measures_.imu_.push_back(imu_buffer_.front());

        imu_buffer_.pop_front();
    }

    current_max_imu_gap_ = 0.0;
    for (std::size_t i = 1; i < measures_.imu_.size(); ++i) {
        current_max_imu_gap_ =
            std::max(current_max_imu_gap_, measures_.imu_[i]->timestamp - measures_.imu_[i - 1]->timestamp);
    }

    lidar_buffer_.pop_front();
    time_buffer_.pop_front();
    if (!preprocess_time_buffer_ms_.empty()) preprocess_time_buffer_ms_.pop_front();
    if (!lidar_stats_buffer_.empty()) lidar_stats_buffer_.pop_front();
    lidar_pushed_ = false;

    // LOG(INFO) << "sync: " << std::setprecision(14) << measures_.lidar_begin_time_ << ", " <<
    // measures_.lidar_end_time_;

    return true;
}

/**
 * @brief 增量式地图构建：将当前帧点云添加到全局地图中
 *
 * 功能说明：
 * 1. 将当前帧点云从机体坐标系转换到世界坐标系
 * 2. 根据距离和邻近点信息决定哪些点需要添加到地图
 * 3. 实现自适应下采样，避免地图密度过大
 * 4. 更新IVox体素索引结构
 *
 * @note 该函数在每次激光雷达更新后调用，实现地图的实时增量构建
 * @note 使用体素栅格进行空间下采样，提高地图构建效率
 */
void LaserMapping::MapIncremental() {
    PointVector points_to_add;             // 需要添加到地图的点（经过智能下采样筛选）
    PointVector point_no_need_downsample;  // 无需下采样的点（稀疏区域的点，直接添加）

    size_t cur_pts = scan_down_lidar_->size();
    points_to_add.reserve(cur_pts);
    point_no_need_downsample.reserve(cur_pts);

    std::vector<size_t> index(cur_pts);
    for (size_t i = 0; i < cur_pts; ++i) {
        index[i] = i;
    }

    // TODO. 这里的逻辑不太对，导致体素边缘的点可能较为密集
    std::for_each(index.begin(), index.end(), [&](const size_t &i) {
        /* transform to world frame */
        PointBodyToWorld(scan_down_lidar_->points[i], scan_down_world_->points[i]);

        /* decide if need add to map */
        PointType &point_world = scan_down_world_->points[i];
        // 智能下采样：根据邻近点信息决定是否需要添加该点到地图
        if (!nearest_points_[i].empty() && flg_EKF_inited_) {
            const PointVector &points_near = nearest_points_[i];

            // 计算当前点所在体素的中心坐标
            Eigen::Vector3f center =
                ((point_world.getVector3fMap() / filter_size_map_min_).array().floor() + 0.5) * filter_size_map_min_;

            // 计算最近邻点到体素中心的距离
            Eigen::Vector3f dis_2_center = points_near[0].getVector3fMap() - center;

            // 如果最近邻点距离体素中心较远，说明该区域点云稀疏，直接添加点无需下采样
            if (fabs(dis_2_center.x()) > 0.5 * filter_size_map_min_ &&
                fabs(dis_2_center.y()) > 0.5 * filter_size_map_min_ &&
                fabs(dis_2_center.z()) > 0.5 * filter_size_map_min_) {
                point_no_need_downsample.emplace_back(point_world);
                return;
            }

            // 检查体素内是否已有更靠近中心的点（密集区域下采样）
            bool need_add = true;
            float dist = math::calc_dist(point_world.getVector3fMap(), center);  // 当前点到体素中心的距离
            if (points_near.size() >= fasterlio::NUM_MATCH_POINTS) {             // 只在邻近点足够多时进行下采样
                for (int readd_i = 0; readd_i < fasterlio::NUM_MATCH_POINTS; readd_i++) {
                    if (math::calc_dist(points_near[readd_i].getVector3fMap(), center) < dist + 1e-6) {
                        need_add = false;  // 已有更靠近中心的点，当前点无需添加（实现体素内点云去重）
                        break;
                    }
                }
            }

            // 只有需要时才添加到地图（实现自适应下采样）
            if (need_add) {
                points_to_add.emplace_back(point_world);  // 这并发可能有点问题
            }
        } else {
            // 初始化阶段或无邻近点时直接添加所有点
            points_to_add.emplace_back(point_world);
        }
    });

    Timer::Evaluate(
        [&, this]() {
            ivox_->AddPoints(points_to_add);
            ivox_->AddPoints(point_no_need_downsample);
        },
        "    IVox Add Points");
}

/**
 * @brief 激光雷达点云配准观测模型，用于传入ESKF内
 * @details 计算当前激光雷达点云与IVox局部地图之间的残差，并累加成ESKF需要的信息矩阵形式。
 *
 * 算法流程：
 * 1. 使用当前迭代状态s，将scan_down_lidar_从Lidar系变换到世界系；
 * 2. 在IVox地图中为每个点搜索最近邻，并尝试拟合局部平面；
 * 3. 用点到平面的距离构造点面残差和位姿雅可比；
 * 4. 可选地加入点到点ICP残差；
 * 5. 将所有残差累加为 H^T H 和 H^T r，交给ESKF::Update()求解位姿增量。
 *
 * @param s[in] 当前ESKF状态，包含位姿、速度、零偏等
 * @param obs[out] 观测模型结构体，填充valid_、HTH_、HTr_和残差统计信息
 */
double LaserMapping::PointInformationScale(const PointType &point, const Vec3d &plane_normal_world,
                                           const NavState &state) const {
    if (!point_noise_enabled_) return 1.0;
    Vec3d sensor_origin = Vec3d::Zero();
    if (const auto *sensor = multi_lidar_config_.FindLidar(point.lidar_id)) {
        sensor_origin = sensor->t_lidar_to_primary;
    }
    const Vec3d beam_primary = point.getVector3fMap().cast<double>() - sensor_origin;
    const double range = beam_primary.norm();
    if (range < 1e-6) return min_information_scale_;
    const Mat3d R_world_primary = state.rot_.matrix() * offset_R_lidar_fixed_;
    const Vec3d normal_primary = R_world_primary.transpose() * plane_normal_world.normalized();
    const double cos_incidence = std::clamp(normal_primary.dot(beam_primary / range), -1.0, 1.0);
    const double tangential_sigma = range * angular_noise_sigma_rad_;
    const double variance = range_noise_sigma_ * range_noise_sigma_ * cos_incidence * cos_incidence +
                            tangential_sigma * tangential_sigma * (1.0 - cos_incidence * cos_incidence);
    const double reference_variance = noise_reference_sigma_ * noise_reference_sigma_;
    return std::clamp(reference_variance / std::max(variance, 1e-9), min_information_scale_,
                      max_information_scale_);
}

void LaserMapping::ObsModel(NavState &s, ESKF::CustomObservationModel &obs) {
    int cnt_pts = scan_down_lidar_->size();

    // 并行处理当前帧点云时使用的索引数组，避免在lambda里依赖递增变量。
    std::vector<size_t> index(cnt_pts);
    for (size_t i = 0; i < index.size(); ++i) {
        index[i] = i;
    }

    // LOG(INFO) << "obs from state: " << s.pos_.transpose() << ", " << s.rot_.unit_quaternion().coeffs().transpose();

    Timer::Evaluate(
        [&, this]() {
            // 当前迭代状态s给出IMU到世界的位姿，offset_*是Lidar到IMU的固定外参。
            // 因此 R_wl/t_wl 表示当前Lidar帧到世界系的变换。
            Mat3f R_wl = (s.rot_.matrix() * offset_R_lidar_fixed_).cast<float>();
            Vec3f t_wl = (s.rot_ * offset_t_lidar_fixed_ + s.pos_).cast<float>();

            std::for_each(std::execution::par_unseq, index.begin(), index.end(), [&](const size_t &i) {
                PointType &point_lidar = scan_down_lidar_->points[i];
                PointType &point_world = scan_down_world_->points[i];

                /// 将点从当前帧Lidar坐标系变换到世界坐标系，用于在局部地图中查最近邻。
                Vec3f p_lidar = point_lidar.getVector3fMap();
                point_world.getVector3fMap() = R_wl * p_lidar + t_wl;
                point_world.intensity = point_lidar.intensity;

                auto &points_near = nearest_points_[i];
                points_near.clear();

                /// 在IVox地图中搜索当前世界系点附近的地图点，后续用这些点拟合局部平面。
                ivox_->GetClosestPoint(point_world, points_near, fasterlio::NUM_MATCH_POINTS);
                point_selected_surf_[i] = points_near.size() >= fasterlio::MIN_NUM_MATCH_POINTS;

                // 点到点ICP复用同一批最近邻。只有能找到足够最近邻的点才有资格参与点到点约束。
                point_selected_icp_[i] = point_selected_surf_[i];

                /// 最近邻数量足够时，尝试拟合局部平面，plane_coef_[i] = [nx, ny, nz, d]。
                if (point_selected_surf_[i]) {
                    point_selected_surf_[i] =
                        math::esti_plane(plane_coef_[i], points_near, fasterlio::ESTI_PLANE_THRESHOLD);
                }
                /// 平面拟合和有效性验证
                if (point_selected_surf_[i]) {
                    auto temp = point_world.getVector4fMap();
                    temp[3] = 1.0;
                    float pd2 = plane_coef_[i].dot(temp);  ///< 计算点到平面的距离（有符号）
                    // 根据点面残差和量测距离做经验筛选，剔除几何关系不可靠的匹配。
                    // 该条件等价于 FAST-LIO 中的 s = 1 - 0.9 * fabs(pd2) / sqrt(range), s > 0.9。
                    // 化简后得到 range > 81 * pd2^2，因此这里应使用 p_lidar.norm() 而不是 squaredNorm()。
                    // 含义是远处点允许稍大的点面残差，近处点需要更严格的平面一致性。
                    // 但是要踢掉点到面距离太离谱的匹配点
                    bool valid_corr = p_lidar.norm() > 81 * pd2 * pd2;
                    if (valid_corr) {
                        point_selected_surf_[i] = true;
                        residuals_[i] = pd2;
                    } else {
                        point_selected_surf_[i] = false;
                    }
                }
            });
        },
        "    ObsModel (Lidar Match)");

    effect_feat_surf_ = 0;
    effect_feat_icp_ = 0;

    // 将并行阶段筛选出的有效点压缩到corr_pts_/corr_norm_前部，便于后面只遍历有效点。
    // corr_pts_前三维保存Lidar系点坐标，第4维保存点面有符号残差。
    corr_pts_.resize(cnt_pts);
    corr_lidar_ids_.resize(cnt_pts);
    corr_norm_.resize(cnt_pts);
    for (int i = 0; i < cnt_pts; i++) {
        if (point_selected_surf_[i]) {
            corr_norm_[effect_feat_surf_] = plane_coef_[i];
            corr_pts_[effect_feat_surf_] = scan_down_lidar_->points[i].getVector4fMap();
            corr_pts_[effect_feat_surf_][3] = residuals_[i];
            corr_lidar_ids_[effect_feat_surf_] = scan_down_lidar_->points[i].lidar_id;

            effect_feat_surf_++;
        }

        if (point_selected_icp_[i]) {
            effect_feat_icp_++;
        }
    }

    corr_pts_.resize(effect_feat_surf_);
    corr_lidar_ids_.resize(effect_feat_surf_);
    corr_norm_.resize(effect_feat_surf_);

    // 有效点面约束太少时，观测模型不可用；ESKF::Update()会放弃本次更新并回退。
    if (effect_feat_surf_ < 20) {
        obs.valid_ = false;
        LOG(WARNING) << "No enough effective surface points: " << effect_feat_surf_ << ", icp: " << effect_feat_icp_
                     << ", required: " << 20;
        return;
    }

    index.resize(effect_feat_surf_);
    const Mat3f off_R = offset_R_lidar_fixed_.cast<float>();
    const Vec3f off_t = offset_t_lidar_fixed_.cast<float>();
    const Mat3f Rt = s.rot_.matrix().transpose().cast<float>();

    /// 点面ICP部分：每个有效点贡献一个标量残差 r = - point_to_plane_distance。
    /// 观测只约束6维位姿，因此最终累加到6x6的HTH_和6x1的HTr_。
    obs.HTH_.setZero();
    obs.HTr_.setZero();

    std::vector<Mat6d> JTJ(effect_feat_surf_);
    std::vector<Vec6d> JTr(effect_feat_surf_);

    std::vector<double> res_sq(index.size());

    std::for_each(std::execution::par_unseq, index.begin(), index.end(), [&](const size_t &i) {
        Vec3f point_this_be = corr_pts_[i].head<3>();      //< lidar坐标系下的点
        Vec3f point_this = off_R * point_this_be + off_t;  ///< IMU坐标系下的点
        // 旋转扰动对点坐标的影响由叉乘矩阵表达，后面用于构造姿态雅可比。
        Mat3f point_crossmat = math::SKEW_SYM_MATRIX(point_this);

        /*** get the normal vector of closest surface/corner ***/
        Vec3f norm_vec = corr_norm_[i].head<3>();

        /*** calculate the Measurement Jacobian matrix H ***/
        // 残差是世界系点到世界系平面的距离。平移雅可比就是平面法向量；
        // 姿态雅可比需要把世界系法向量转回IMU切空间，再与Lidar点在IMU系下的位置叉乘。
        Vec3f C(Rt * norm_vec);
        Vec3f A(point_crossmat * C);

        Eigen::Matrix<double, 1, ESKF::pose_obs_dim_> J;
        J.setZero();
        J << norm_vec[0], norm_vec[1], norm_vec[2], A[0], A[1], A[2];

        // corr_pts_[i][3]里保存的是点到平面的有符号距离pd2。
        // 这里取负号，是为了让后续求解的dx沿着减小残差的方向更新。
        float res = -corr_pts_[i][3];

        // double w = huber_weight(res);
        PointType weighted_point;
        weighted_point.x = point_this_be.x();
        weighted_point.y = point_this_be.y();
        weighted_point.z = point_this_be.z();
        weighted_point.lidar_id = corr_lidar_ids_[i];
        double w = PointInformationScale(weighted_point, norm_vec.cast<double>(), s);

        JTJ[i] = (J.transpose() * J).eval() * w;
        JTr[i] = J.transpose() * res * w;

        res_sq[i] = res * res;
    });

    // 并行阶段先把每个点的J^T J和J^T r保存下来，串行阶段再累加，避免并发写obs。
    for (int i = 0; i < index.size(); ++i) {
        obs.HTH_ += JTJ[i] * options_.plane_icp_weight_;
        obs.HTr_ += JTr[i] * options_.plane_icp_weight_;
    }

    // 残差统计用于ESKF迭代过程中的收敛判断和AA回退判断。
    if (!res_sq.empty()) {
        std::sort(res_sq.begin(), res_sq.end());
        obs.lidar_residual_mean_ = res_sq[res_sq.size() / 2];
        obs.lidar_residual_max_ = res_sq[res_sq.size() - 1];
        // LOG(INFO) << "residual mean: " << obs.lidar_residual_mean_ << ", max: " << obs.lidar_residual_max_
        //           << ", 85%: " << res_sq[res_sq.size() * 0.85];
    }

    /// 点到点ICP部分

    if (options_.enable_icp_part_) {
        // 点到点ICP为每个有效点贡献3维残差：当前点世界坐标 - 最近地图点世界坐标。
        // 它是点面约束的补充，权重由options_.icp_weight_控制。
        JTJ.resize(cnt_pts);
        JTr.resize(cnt_pts);

        std::vector<size_t> index(cnt_pts);
        for (size_t i = 0; i < index.size(); ++i) {
            index[i] = i;
        }

        std::for_each(std::execution::par_unseq, index.begin(), index.end(), [&](const size_t &i) {
            if (point_selected_icp_[i] == false) {
                return;
            }

            /// q是当前点在Lidar系下的坐标，qs是使用当前状态投影后的世界系坐标。
            Vec3d q = scan_down_lidar_->points[i].getVector3fMap().cast<double>();
            Vec3d qs = scan_down_world_->points[i].getVector3fMap().cast<double>();

            Eigen::Matrix<double, 3, ESKF::pose_obs_dim_> J;
            J.setZero();

            /// translation 部分
            J.block<3, 3>(0, 0) = Mat3d::Identity();

            /// rotation 部分
            J.block<3, 3>(0, 3) = -(s.rot_.matrix() * offset_R_lidar_fixed_) * SO3::hat(q);

            // 点到点残差：当前点世界坐标和最近地图点世界坐标之差。
            Vec3d e = qs - nearest_points_[i][0].getVector3fMap().cast<double>();

            // 过大的点到点残差通常对应错误最近邻，直接剔除。
            if (e.norm() > 0.5) {
                point_selected_icp_[i] = false;
                return;
            }

            JTJ[i] = J.transpose() * J;
            JTr[i] = -J.transpose() * e;
        });

        // 将点到点ICP的信息量累加到同一个6维位姿观测中。
        for (int i = 0; i < cnt_pts; ++i) {
            if (point_selected_icp_[i] == false) {
                continue;
            }
            obs.HTH_ += JTJ[i] * options_.icp_weight_;
            obs.HTr_ += JTr[i] * options_.icp_weight_;
        }
    }
}

///////////////////////////  private method /////////////////////////////////////////////////////////////////////

/**
 * @brief 根据全部关键帧点云拼接全局地图。
 * @param use_lio_pose true时使用前端LIO位姿，false时使用关键帧优化位姿。
 * @param use_voxel 是否对单帧关键帧点云和最终全局点云执行体素滤波。
 * @param res 体素滤波叶子尺寸，单位m。
 * @return 拼接并设置好PCD元信息的全局点云。
 */
CloudPtr LaserMapping::PrepareMapExportCloud(const CloudPtr& cloud) const {
    if (!cloud || !map_export_self_point_filter_config_.enabled) return cloud;
    CloudPtr filtered(new PointCloudType(*cloud));
    FilterSelfPoints(*filtered, map_export_self_point_filter_config_,
                     GetInitialLidarRotation().matrix());
    return filtered;
}

CloudPtr LaserMapping::GetGlobalMap(bool use_lio_pose, bool use_voxel, float res,
                                    bool apply_map_export_filter) {
    CloudPtr global_map(new PointCloudType);
    std::size_t map_export_removed = 0;

    /// 体素滤波器在关键帧级和全局地图级复用，分辨率由调用方指定。
    pcl::VoxelGrid<PointType> voxel;
    voxel.setLeafSize(res, res, res);

    /// 关键帧点云保存在Lidar坐标系下，拼接前需要先变换到IMU坐标系。
    SE3 T_imu_lidar(Eigen::Quaterniond(offset_R_lidar_fixed_).normalized(), offset_t_lidar_fixed_);

    for (auto &kf : all_keyframes_) {
        const CloudPtr source_cloud = kf->GetCloud();
        CloudPtr cloud = apply_map_export_filter
                             ? PrepareMapExportCloud(source_cloud)
                             : source_cloud;
        if (!cloud || cloud->empty()) continue;
        if (apply_map_export_filter && source_cloud) {
            map_export_removed += source_cloud->size() - cloud->size();
        }

        CloudPtr cloud_filter(new PointCloudType);

        /// 可选地先对单个关键帧点云降采样，降低全局拼接的点数和内存占用。
        if (use_voxel) {
            if (multi_lidar_config_.enabled) {
                cloud_filter = DownsamplePreservingSource(cloud, res);
            } else {
                voxel.setInputCloud(cloud);
                voxel.filter(*cloud_filter);
            }

        } else {
            cloud_filter = cloud;
        }
        if (!cloud_filter || cloud_filter->empty()) continue;

        CloudPtr cloud_trans(new PointCloudType);

        /// 保存调试地图时可使用原始LIO位姿；保存最终地图时通常使用回环优化后的位姿。
        if (use_lio_pose) {
            pcl::transformPointCloud(*cloud_filter, *cloud_trans, (kf->GetLIOPose() * T_imu_lidar).matrix());
        } else {
            pcl::transformPointCloud(*cloud_filter, *cloud_trans, (kf->GetOptPose() * T_imu_lidar).matrix());
        }

        /// 将当前关键帧点云累加到世界坐标系下的全局点云。
        *global_map += *cloud_trans;

        LOG(INFO) << "kf " << kf->GetID() << ", pose: " << kf->GetOptPose().translation().transpose();
    }

    CloudPtr global_map_filtered(new PointCloudType);
    /// 拼接完成后再做一次全局体素滤波，合并关键帧重叠区域中的冗余点。
    if (use_voxel) {
        if (multi_lidar_config_.enabled) {
            global_map_filtered = DownsamplePreservingSource(global_map, res);
        } else {
            voxel.setInputCloud(global_map);
            voxel.filter(*global_map_filtered);
        }
    } else {
        global_map_filtered = global_map;
    }

    /// 补齐PCL保存PCD时需要的点云组织信息。
    global_map_filtered->is_dense = false;
    global_map_filtered->height = 1;
    global_map_filtered->width = global_map_filtered->size();

    LOG(INFO) << "global map: " << global_map_filtered->size();
    if (apply_map_export_filter && map_export_self_point_filter_config_.enabled) {
        LOG(INFO) << "map-export self-point filter removed "
                  << map_export_removed << " keyframe points";
    }

    return global_map_filtered;
}

void LaserMapping::SaveMap() {
    /// 保存地图
    auto global_map = GetGlobalMap(true);

    pcl::io::savePCDFileBinaryCompressed("./data/lio.pcd", *global_map);

    LOG(INFO) << "lio map is saved to ./data/lio.pcd";
}

CloudPtr LaserMapping::GetRecentCloud() {
    if (lidar_buffer_.empty()) {
        return nullptr;
    }

    return lidar_buffer_.front();
}

CloudPtr LaserMapping::GetProjCloud() {
    auto cloud = scan_undistort_;
    ProjectKFs(cloud);
    return cloud;
}

}  // namespace lightning
