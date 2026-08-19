#include <pcl/common/transforms.h>
#include <pcl_conversions/pcl_conversions.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>

#include "core/localization/lidar_loc/lidar_loc.h"
#include "core/localization/localization.h"

#include <opencv2/highgui.hpp>

#include "core/localization/pose_graph/pgo.h"
#include "io/yaml_io.h"
#include "ui/pangolin_window.h"

namespace lightning::loc {

// ！ 构造函数
Localization::Localization(Options options) { options_ = options; }

// ！初始化函数
bool Localization::Init(const std::string& yaml_path, const std::string& global_map_path) {
    std::unique_lock<std::shared_mutex> lock(lifecycle_mutex_);
    if (lidar_loc_ != nullptr) {
        // 若已经启动，则变为初始化
        // Finish() joins the sensor worker, which may itself be waiting for
        // lifecycle access. Do not hold that lock while waiting for the worker.
        lock.unlock();
        Finish();
        lock.lock();
    }

    YAML_IO yaml(yaml_path);
    options_.with_ui_ = yaml.GetValue<bool>("system", "with_ui");

    /// lidar odom前端
    LaserMapping::Options opt_lio;
    opt_lio.is_in_slam_mode_ = false;

    lio_ = std::make_shared<LaserMapping>(opt_lio);
    if (!lio_->Init(yaml_path)) {
        LOG(ERROR) << "failed to init lio";
        return false;
    }

    /// 激光定位
    LidarLoc::Options lidar_loc_options;
    lidar_loc_options.update_dynamic_cloud_ = yaml.GetValue<bool>("lidar_loc", "update_dynamic_cloud");
    lidar_loc_options.force_2d_ = yaml.GetValue<bool>("lidar_loc", "force_2d");
    lidar_loc_options.map_option_.enable_dynamic_polygon_ = false;
    lidar_loc_options.map_option_.map_path_ = global_map_path;
    lidar_loc_ = std::make_shared<LidarLoc>(lidar_loc_options);

    if (options_.with_ui_) {
        ui_ = std::make_shared<ui::PangolinWindow>();
        ui_->SetCurrentScanSize(1);
        if (ui_->Init()) {
            lidar_loc_->SetUI(ui_);
        } else {
            LOG(ERROR) << "failed to init 3D UI, continue without Pangolin";
            ui_.reset();
            options_.with_ui_ = false;
        }

        // lio_->SetUI(ui_);
    }

    if (!lidar_loc_->Init(yaml_path)) {
        LOG(ERROR) << "failed to initialize lidar localization";
        return false;
    }

    /// pose graph
    pgo_ = std::make_shared<PGO>();
    pgo_->SetDebug(false);
    const YAML::Node localization_pgo = YAML::LoadFile(yaml_path)["localization_pgo"];
    pgo_->SetDrSmoothingEnabled(localization_pgo
                                    ? localization_pgo["enable_dr_smoothing"].as<bool>(true)
                                    : true);
    pgo_->SetDrExtrapolationEnabled(localization_pgo
                                        ? localization_pgo["enable_dr_extrapolation"].as<bool>(true)
                                        : true);

    ///  各模块的异步调用
    const YAML::Node system = YAML::LoadFile(yaml_path)["system"];
    options_.enable_lidar_loc_skip_ =
        system && system["enable_lidar_loc_skip"] ? system["enable_lidar_loc_skip"].as<bool>() : false;
    options_.enable_lidar_loc_rviz_ =
        system && system["enable_lidar_loc_rviz"] ? system["enable_lidar_loc_rviz"].as<bool>() : false;
    options_.lidar_loc_skip_num_ =
        system && system["lidar_loc_skip_num"] ? system["lidar_loc_skip_num"].as<int>() : 1;
    options_.enable_lidar_odom_skip_ =
        system && system["enable_lidar_odom_skip"] ? system["enable_lidar_odom_skip"].as<bool>() : false;
    options_.lidar_odom_skip_num_ =
        system && system["lidar_odom_skip_num"] ? system["lidar_odom_skip_num"].as<int>() : 1;
    options_.loc_on_kf_ = yaml.GetValue<bool>("lidar_loc", "loc_on_kf");
    imu_static_hold_enabled_ =
        system && system["enable_imu_static_hold"]
            ? system["enable_imu_static_hold"].as<bool>()
            : false;
    {
        std::lock_guard<std::mutex> static_lock(static_detector_mutex_);
        static_imu_window_.clear();
        static_gyro_sum_ = 0.0;
        static_gyro_sq_sum_ = 0.0;
        static_accel_sum_ = 0.0;
        static_accel_sq_sum_ = 0.0;
        last_static_lio_stamp_ = 0.0;
        last_static_lio_speed_ = 0.0;
        last_static_lio_reliable_ = false;
        last_valid_lidar_loc_stamp_ = 0.0;
        last_wheel_speed_stamp_ = 0.0;
        last_wheel_speed_arrival_ = {};
        last_wheel_speed_mps_ = 0.0;
        wheel_speed_observed_ = false;
        static_exit_count_ = 0;
        imu_static_hold_active_ = false;
    }

    const size_t sensor_queue_size =
        system && system["online_sensor_queue_size"]
            ? system["online_sensor_queue_size"].as<size_t>()
            : 10000;
    if (sensor_queue_size == 0) {
        LOG(ERROR) << "system.online_sensor_queue_size must be positive";
        return false;
    }
    sensor_proc_.SetMaxSize(sensor_queue_size);
    online_sensor_max_lag_sec_ =
        system && system["online_sensor_max_lag_sec"]
            ? system["online_sensor_max_lag_sec"].as<double>()
            : 0.0;
    online_sensor_resume_lag_sec_ =
        system && system["online_sensor_resume_lag_sec"]
            ? system["online_sensor_resume_lag_sec"].as<double>()
            : online_sensor_max_lag_sec_ * 0.5;
    if (online_sensor_max_lag_sec_ < 0.0 || online_sensor_resume_lag_sec_ < 0.0 ||
        (online_sensor_max_lag_sec_ > 0.0 &&
         online_sensor_resume_lag_sec_ >= online_sensor_max_lag_sec_)) {
        LOG(ERROR) << "online sensor lag thresholds must satisfy 0 <= resume < max";
        return false;
    }
    lidar_overload_throttled_ = false;
    LOG(INFO) << "online sensor queue size=" << sensor_queue_size
              << ", lidar overload max_lag_sec=" << online_sensor_max_lag_sec_
              << ", resume_lag_sec=" << online_sensor_resume_lag_sec_;
    // Online localization must operate on the freshest projected scan. A
    // backlog is harmful here: global relocalization can be expensive, and
    // replaying stale scans afterwards prevents timely confirmation.
    const size_t default_loc_queue_size = 1;
    const size_t loc_queue_size =
        system && system["online_lidar_loc_queue_size"]
            ? system["online_lidar_loc_queue_size"].as<size_t>()
            : default_loc_queue_size;
    if (loc_queue_size == 0) {
        LOG(ERROR) << "system.online_lidar_loc_queue_size must be positive";
        return false;
    }
    lidar_loc_proc_cloud_.SetMaxSize(loc_queue_size);

    sensor_proc_.SetName("传感器顺序队列");
    lidar_loc_proc_cloud_.SetName("激光定位");
    high_frequency_output_proc_.SetName("高频定位输出");
    high_frequency_output_proc_.SetMaxSize(1);

    // 允许跳帧
    lidar_loc_proc_cloud_.SetSkipParam(options_.enable_lidar_loc_skip_, options_.lidar_loc_skip_num_);

    sensor_proc_.SetProcFunc([this](const SensorInput& input) { ProcessSensorInput(input); });
    lidar_loc_proc_cloud_.SetProcFunc([this](CloudPtr cloud) { LidarLocProcCloud(cloud); });

    if (options_.online_mode_) {
        sensor_proc_.Start();
        lidar_loc_proc_cloud_.Start();
    }

    /// TODO: 发布
    auto publish_localization_result = [this](const LocalizationResult& res) {
        // if (loc_result_.timestamp_ > 0) {
        //             double loc_fps = 1.0 / (res.timestamp_ - loc_result_.timestamp_);
        //             // LOG_EVERY_N(INFO, 10) << "loc fps: " << loc_fps;
        //         }

        loc_result_ = res;

        if (tf_callback_ && loc_result_.valid_) {
            tf_callback_(loc_result_.ToGeoMsg());
        }

        if (ui_) {
            ui_->UpdateNavState(loc_result_.ToNavState());
            ui_->UpdateRecentPose(loc_result_.pose_);
        }

        if (localization_result_callback_) {
            localization_result_callback_(loc_result_);
        }
    };
    // Global optimization output can arrive behind the current extrapolated
    // timestamp. Keep it on a separate evaluation path so live ROS output
    // never rolls time backwards.
    pgo_->SetGlobalOutputHandleFunction([this](const LocalizationResult& res) {
        if (global_localization_result_callback_) global_localization_result_callback_(res);
    });
    high_frequency_output_proc_.SetProcFunc(publish_localization_result);
    pgo_->SetHighFrequencyGlobalOutputHandleFunction(
        [this, publish_localization_result](const LocalizationResult& result) {
            if (options_.online_mode_) {
                high_frequency_output_proc_.AddMessage(result);
            } else {
                publish_localization_result(result);
            }
        });
    if (options_.online_mode_) high_frequency_output_proc_.Start();

    /// 预处理器
    preprocess_.reset(new PointCloudPreprocess());
    preprocess_->Blind() = yaml.GetValue<double>("fasterlio", "blind");
    preprocess_->TimeScale() = yaml.GetValue<double>("fasterlio", "time_scale");
    int lidar_type = yaml.GetValue<int>("fasterlio", "lidar_type");
    preprocess_->NumScans() = yaml.GetValue<int>("fasterlio", "scan_line");
    preprocess_->PointFilterNum() = yaml.GetValue<int>("fasterlio", "point_filter_num");
    float height_max = yaml.GetValue<float>("roi", "height_max");
    float height_min = yaml.GetValue<float>("roi", "height_min");

    preprocess_->SetHeightROI(height_max, height_min);

    LOG(INFO) << "lidar_type " << lidar_type;
    if (lidar_type == 1) {
        preprocess_->SetLidarType(LidarType::AVIA);
        LOG(INFO) << "Using AVIA Lidar";
    } else if (lidar_type == 2) {
        preprocess_->SetLidarType(LidarType::VELO32);
        LOG(INFO) << "Using Velodyne 32 Lidar";
    } else if (lidar_type == 3) {
        preprocess_->SetLidarType(LidarType::OUST64);
        LOG(INFO) << "Using OUST 64 Lidar";
    } else if (lidar_type == 4) {
        preprocess_->SetLidarType(LidarType::ROBOSENSE);
        LOG(INFO) << "Using OUST 64 Lidar";
    } else {
        LOG(WARNING) << "unknown lidar_type";
    }

    return true;
}

void Localization::ProcessLidarMsg(const sensor_msgs::msg::PointCloud2::SharedPtr cloud) {
    std::shared_lock<std::shared_mutex> lifecycle_lock(lifecycle_mutex_);
    const int lidar_id = lio_ && lio_->IsMultiLidarEnabled() ? lio_->GetMultiLidarConfig().primary_lidar_id : 0;
    lifecycle_lock.unlock();
    ProcessLidarMsg(cloud, lidar_id);
}

void Localization::ProcessLidarMsg(const sensor_msgs::msg::PointCloud2::SharedPtr cloud, int lidar_id) {
    std::shared_lock<std::shared_mutex> lifecycle_lock(lifecycle_mutex_);
    const double timestamp = cloud
                                 ? static_cast<double>(cloud->header.stamp.sec) +
                                       static_cast<double>(cloud->header.stamp.nanosec) * 1e-9
                                 : 0.0;
    if (ShouldThrottleLidarInput(timestamp)) return;
    UL input_lock(input_mutex_);
    if (preprocess_ == nullptr || lidar_loc_ == nullptr || lio_ == nullptr || pgo_ == nullptr) {
        return;
    }

    // 串行模式
    CloudPtr laser_cloud(new PointCloudType);
    preprocess_->Process(cloud, laser_cloud);
    laser_cloud->header.stamp = cloud->header.stamp.sec * 1e9 + cloud->header.stamp.nanosec;

    bool process_lidar_odom = true;
    if (options_.enable_lidar_odom_skip_) {
        const int skip_num = options_.lidar_odom_skip_num_ > 0 ? options_.lidar_odom_skip_num_ : 1;
        process_lidar_odom = lidar_odom_skip_cnt_ == 0;
        lidar_odom_skip_cnt_ = (lidar_odom_skip_cnt_ + 1) % skip_num;
    }
    if (!process_lidar_odom) return;
    if (options_.online_mode_) {
        ObserveSensorEnqueued(static_cast<double>(laser_cloud->header.stamp) * 1e-9);
        sensor_proc_.AddMessage({nullptr, laser_cloud, lidar_id, false});
        return;
    }
    input_lock.unlock();
    lifecycle_lock.unlock();
    LidarOdomProcCloud(laser_cloud, lidar_id);
}

void Localization::ProcessLivoxLidarMsg(const livox_ros_driver2::msg::CustomMsg::SharedPtr cloud) {
    std::shared_lock<std::shared_mutex> lifecycle_lock(lifecycle_mutex_);
    const double timestamp = cloud
                                 ? static_cast<double>(cloud->header.stamp.sec) +
                                       static_cast<double>(cloud->header.stamp.nanosec) * 1e-9
                                 : 0.0;
    if (ShouldThrottleLidarInput(timestamp)) return;
    UL input_lock(input_mutex_);
    if (preprocess_ == nullptr || lidar_loc_ == nullptr || lio_ == nullptr || pgo_ == nullptr) {
        return;
    }

    // 串行模式
    CloudPtr laser_cloud(new PointCloudType);
    preprocess_->Process(cloud, laser_cloud);
    laser_cloud->header.stamp = cloud->header.stamp.sec * 1e9 + cloud->header.stamp.nanosec;

    const int lidar_id = lio_->IsMultiLidarEnabled() ? lio_->GetMultiLidarConfig().primary_lidar_id : 0;
    bool process_lidar_odom = true;
    if (options_.enable_lidar_odom_skip_) {
        const int skip_num = options_.lidar_odom_skip_num_ > 0 ? options_.lidar_odom_skip_num_ : 1;
        process_lidar_odom = lidar_odom_skip_cnt_ == 0;
        lidar_odom_skip_cnt_ = (lidar_odom_skip_cnt_ + 1) % skip_num;
    }
    if (!process_lidar_odom) return;
    if (options_.online_mode_) {
        ObserveSensorEnqueued(static_cast<double>(laser_cloud->header.stamp) * 1e-9);
        sensor_proc_.AddMessage({nullptr, laser_cloud, lidar_id, false});
        return;
    }
    input_lock.unlock();
    lifecycle_lock.unlock();
    LidarOdomProcCloud(laser_cloud, lidar_id);
}

void Localization::ProcessSensorInput(const SensorInput& input) {
    const double timestamp = input.is_imu && input.imu
                                 ? input.imu->timestamp
                                 : (input.cloud ? static_cast<double>(input.cloud->header.stamp) * 1e-9 : 0.0);
    if (input.is_imu) {
        ProcessIMUData(input.imu);
    } else {
        LidarOdomProcCloud(input.cloud, input.lidar_id);
    }
    // This frontier means callback completion, not merely dequeue/start.
    ObserveSensorProcessed(timestamp);
}

void Localization::LidarOdomProcCloud(CloudPtr cloud, int lidar_id) {
    // LIO must observe IMU and lidar in callback order. Running this update in
    // a second worker lets later IMU callbacks overtake the cloud and changes
    // the filter result relative to offline processing.
    std::shared_lock<std::shared_mutex> lifecycle_lock(lifecycle_mutex_);
    UL processing_lock(processing_mutex_);

    if (lio_ == nullptr) {
        return;
    }

    /// NOTE: 在NCLT这种数据集中，lio内部是有缓存的，它拿到的点云不一定是最新时刻的点云
    lio_->ProcessPointCloud2(cloud, lidar_id);
    DrainLioOutputs();
}

void Localization::DrainLioOutputs() {
    while (true) {
        const auto status = lio_->RunDetailed();
        if (status == LaserMapping::RunStatus::kNoData) break;
        if (status != LaserMapping::RunStatus::kOutput) continue;

        auto lo_state = lio_->GetState();
        ObserveLioForStaticDetector(lo_state);

        lidar_loc_->ProcessLO(lo_state);
        pgo_->ProcessLidarOdom(lo_state);

        // LOG(INFO) << "LO pose: " << std::setprecision(12) << lo_state.timestamp_ << " "
        //           << lo_state.GetPose().translation().transpose();

        /// 获得lio的关键帧
        auto scan = lio_->GetProjCloud();

        if (options_.loc_on_kf_) {
            auto kf = lio_->GetKeyframe();
            if (kf == lio_kf_) {
                /// 关键帧未更新，那就只更新IMU状态
                continue;
            }

            // if (ui_) {
            //     ui_->UpdateKF(kf);
            // }

            lio_kf_ = kf;

            // auto scan = lio_->GetScanUndist();

            if (options_.online_mode_) {
                lidar_loc_proc_cloud_.AddMessage(scan);
            } else {
                LidarLocProcCloud(scan);
            }
        } else {
            // auto scan = cloud;   // 这个cloud应该差一个外参

            if (options_.online_mode_) {
                lidar_loc_proc_cloud_.AddMessage(scan);
            } else {
                LidarLocProcCloud(scan);
            }
        }
    }
}

bool Localization::IsMultiLidarEnabled() const { return lio_ && lio_->IsMultiLidarEnabled(); }

const MultiLidarConfig& Localization::GetMultiLidarConfig() const {
    CHECK(lio_ != nullptr);
    return lio_->GetMultiLidarConfig();
}

SO3 Localization::GetInitialLidarRotation() const {
    CHECK(lio_ != nullptr);
    return lio_->GetInitialLidarRotation();
}

Localization::RuntimeStats Localization::GetRuntimeStats() const {
    std::lock_guard<std::mutex> lock(runtime_stats_mutex_);
    RuntimeStats stats = runtime_stats_;
    stats.sensor_queue_pending = sensor_proc_.PendingCount();
    stats.sensor_queue_dropped = sensor_proc_.DroppedCount();
    stats.sensor_queue_processed = sensor_proc_.ProcessedCount();
    stats.localization_queue_pending = lidar_loc_proc_cloud_.PendingCount();
    stats.localization_queue_dropped = lidar_loc_proc_cloud_.DroppedCount();
    stats.localization_queue_processed = lidar_loc_proc_cloud_.ProcessedCount();
    return stats;
}

void Localization::ObserveSensorEnqueued(double timestamp) {
    if (timestamp <= 0.0) return;
    std::lock_guard<std::mutex> lock(runtime_stats_mutex_);
    runtime_stats_.latest_enqueued_sensor_stamp =
        std::max(runtime_stats_.latest_enqueued_sensor_stamp, timestamp);
    runtime_stats_.current_sensor_lag_sec =
        runtime_stats_.latest_processed_sensor_stamp > 0.0
            ? std::max(0.0, runtime_stats_.latest_enqueued_sensor_stamp -
                                runtime_stats_.latest_processed_sensor_stamp)
            : 0.0;
    runtime_stats_.max_sensor_lag_sec =
        std::max(runtime_stats_.max_sensor_lag_sec, runtime_stats_.current_sensor_lag_sec);
}

bool Localization::ShouldThrottleLidarInput(double timestamp) {
    if (!options_.online_mode_ || online_sensor_max_lag_sec_ <= 0.0) return false;

    double lag_sec = 0.0;
    {
        std::lock_guard<std::mutex> lock(runtime_stats_mutex_);
        if (runtime_stats_.latest_enqueued_sensor_stamp > 0.0 &&
            runtime_stats_.latest_processed_sensor_stamp > 0.0) {
            lag_sec = std::max(0.0, runtime_stats_.latest_enqueued_sensor_stamp -
                                        runtime_stats_.latest_processed_sensor_stamp);
        }
    }

    bool throttled = lidar_overload_throttled_.load();
    if (!throttled && lag_sec >= online_sensor_max_lag_sec_) {
        lidar_overload_throttled_ = true;
        throttled = true;
        LOG(WARNING) << "sensor queue lag reached " << lag_sec
                     << " sec; keeping fresh lidar admissible while the bounded queue evicts stale input";
    } else if (throttled && lag_sec <= online_sensor_resume_lag_sec_) {
        lidar_overload_throttled_ = false;
        throttled = false;
        LOG(INFO) << "sensor queue lag recovered to " << lag_sec
                  << " sec; resuming lidar input";
    }

    if (throttled) {
        LOG_EVERY_N(WARNING, 100) << "sensor overload at lidar " << std::setprecision(14)
                                  << timestamp << "; lag_sec=" << lag_sec
                                  << ", queue_pending=" << sensor_proc_.PendingCount()
                                  << ", lidar remains admitted";
    }
    // Never latch into a state that rejects every future lidar frame. The
    // bounded FIFO already removes the oldest stale input under overload.
    return false;
}

void Localization::ObserveSensorProcessed(double timestamp) {
    if (timestamp <= 0.0) return;
    std::lock_guard<std::mutex> lock(runtime_stats_mutex_);
    constexpr double kSevereRollbackSec = 1.0;
    const double rollback = runtime_stats_.latest_processed_sensor_stamp - timestamp;
    if (runtime_stats_.latest_processed_sensor_stamp > 0.0 && rollback > kSevereRollbackSec) {
        ++runtime_stats_.severe_timestamp_rollback_count;
        runtime_stats_.worst_timestamp_rollback_sec =
            std::max(runtime_stats_.worst_timestamp_rollback_sec, rollback);
    }
    // ROS callbacks from multiple lidar topics may arrive out of timestamp
    // order. Keep the processed frontier monotonic so a late stale message is
    // still diagnosed as a rollback without manufacturing queue lag or
    // retriggering overload admission control.
    runtime_stats_.latest_processed_sensor_stamp =
        std::max(runtime_stats_.latest_processed_sensor_stamp, timestamp);
    runtime_stats_.current_sensor_lag_sec =
        std::max(0.0, runtime_stats_.latest_enqueued_sensor_stamp -
                          runtime_stats_.latest_processed_sensor_stamp);
    runtime_stats_.max_sensor_lag_sec =
        std::max(runtime_stats_.max_sensor_lag_sec, runtime_stats_.current_sensor_lag_sec);
}

void Localization::ObserveLioForStaticDetector(const NavState& state) {
    if (!imu_static_hold_enabled_ || state.timestamp_ <= 0.0) return;
    std::lock_guard<std::mutex> lock(static_detector_mutex_);
    last_static_lio_stamp_ = state.timestamp_;
    last_static_lio_speed_ = state.GetVel().norm();
    last_static_lio_reliable_ = state.lidar_odom_reliable_;
}

void Localization::ObserveLidarLocForStaticDetector(const LocalizationResult& result) {
    // LidarLoc owns lidar_loc_valid_, but valid_ is only populated later by
    // PGO. Requiring valid_ here permanently prevented online static entry.
    if (!imu_static_hold_enabled_ || !result.lidar_loc_valid_ ||
        result.timestamp_ <= 0.0) {
        return;
    }
    std::lock_guard<std::mutex> lock(static_detector_mutex_);
    last_valid_lidar_loc_stamp_ = result.timestamp_;
}

bool Localization::UpdateImuStaticState(const IMUPtr& imu) {
    if (!imu_static_hold_enabled_ || !imu || imu->timestamp <= 0.0) return false;

    constexpr double kWindowSec = 1.0;
    constexpr double kMinWindowSec = 0.8;
    constexpr double kEnterGyroMean = 0.025;
    constexpr double kEnterGyroStd = 0.010;
    constexpr double kEnterAccelCv = 0.025;
    constexpr double kEnterLioSpeed = 0.05;
    constexpr double kExitLioSpeed = 0.08;
    constexpr double kExitLioSpeedWithZeroCan = 0.30;
    constexpr double kMaxLioAgeSec = 0.4;
    constexpr double kMaxValidLocAgeSec = 0.8;
    constexpr double kExitGyro = 0.05;
    constexpr double kExitAccelRatio = 0.05;
    constexpr double kMaxWheelSpeedAgeSec = 0.25;
    constexpr double kEnterWheelSpeed = 0.03;
    constexpr double kExitWheelSpeed = 0.05;
    constexpr int kExitSamples = 3;

    const StaticImuSample sample{imu->timestamp, imu->angular_velocity.norm(),
                                 imu->linear_acceleration.norm()};
    std::lock_guard<std::mutex> lock(static_detector_mutex_);
    static_imu_window_.push_back(sample);
    static_gyro_sum_ += sample.gyro_norm;
    static_gyro_sq_sum_ += sample.gyro_norm * sample.gyro_norm;
    static_accel_sum_ += sample.accel_norm;
    static_accel_sq_sum_ += sample.accel_norm * sample.accel_norm;
    while (!static_imu_window_.empty() &&
           sample.timestamp - static_imu_window_.front().timestamp > kWindowSec) {
        const auto& old = static_imu_window_.front();
        static_gyro_sum_ -= old.gyro_norm;
        static_gyro_sq_sum_ -= old.gyro_norm * old.gyro_norm;
        static_accel_sum_ -= old.accel_norm;
        static_accel_sq_sum_ -= old.accel_norm * old.accel_norm;
        static_imu_window_.pop_front();
    }

    const double count = static_cast<double>(static_imu_window_.size());
    if (count < 2.0) return imu_static_hold_active_;
    const double gyro_mean = static_gyro_sum_ / count;
    const double accel_mean = static_accel_sum_ / count;
    const double gyro_std = std::sqrt(std::max(0.0, static_gyro_sq_sum_ / count - gyro_mean * gyro_mean));
    const double accel_std = std::sqrt(std::max(0.0, static_accel_sq_sum_ / count - accel_mean * accel_mean));
    const double accel_scale = std::max(1e-6, accel_mean);
    const double accel_cv = accel_std / accel_scale;
    const double accel_delta_ratio = std::abs(sample.accel_norm - accel_mean) / accel_scale;
    const double window_span = sample.timestamp - static_imu_window_.front().timestamp;
    const double lio_age = sample.timestamp - last_static_lio_stamp_;
    const double loc_age = sample.timestamp - last_valid_lidar_loc_stamp_;
    const bool fresh_lio = lio_age >= 0.0 && lio_age <= kMaxLioAgeSec;
    const bool fresh_valid_loc = loc_age >= 0.0 && loc_age <= kMaxValidLocAgeSec;
    // CAN is published by the other Orin. Its header clock can have a fixed
    // offset from the LiDAR/IMU clock, so cross-sensor header subtraction is
    // not a valid freshness test. Steady callback-arrival time detects an
    // actual CAN silence without depending on clock synchronization.
    const double wheel_arrival_age = wheel_speed_observed_
        ? std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                        last_wheel_speed_arrival_).count()
        : std::numeric_limits<double>::infinity();
    const bool fresh_wheel_speed = wheel_speed_observed_ &&
                                   wheel_arrival_age <= kMaxWheelSpeedAgeSec;
    const bool wheel_reports_stationary =
        fresh_wheel_speed && std::abs(last_wheel_speed_mps_) < kEnterWheelSpeed;

    if (!imu_static_hold_active_) {
        const bool stable_imu_window = gyro_mean < kEnterGyroMean &&
                                       gyro_std < kEnterGyroStd &&
                                       accel_cv < kEnterAccelCv;
        // A running loader has substantially more stationary IMU vibration
        // than the no-CAN bags. Once the motor-speed signal is present, use a
        // debounced zero wheel speed as the stationary observation and retain
        // low-speed LIO plus a fresh map match as independent safeguards. Old
        // bags without CAN continue to require the strict IMU window.
        const bool stationary_observation = fresh_wheel_speed
                                                ? wheel_reports_stationary
                                                : stable_imu_window;
        if (window_span >= kMinWindowSec && stationary_observation && fresh_lio &&
            fresh_valid_loc &&
            last_static_lio_reliable_ &&
            last_static_lio_speed_ < kEnterLioSpeed) {
            imu_static_hold_active_ = true;
            static_exit_count_ = 0;
            LOG(WARNING) << "enter conservative IMU static hold at " << std::setprecision(14)
                         << sample.timestamp << ", lio_speed=" << last_static_lio_speed_
                         << ", wheel_speed="
                         << (fresh_wheel_speed ? last_wheel_speed_mps_ : 0.0)
                         << ", wheel_observed=" << wheel_speed_observed_
                         << ", gyro_mean=" << gyro_mean << ", accel_cv=" << accel_cv;
        }
        return imu_static_hold_active_;
    }

    const bool inertial_motion = sample.gyro_norm > kExitGyro ||
                                 accel_delta_ratio > kExitAccelRatio;
    // Fresh zero CAN suppresses engine-vibration false exits. A moving CAN
    // sample or LIO motion still releases the hold within three IMU samples;
    // if CAN becomes stale, the IMU fallback is active again.
    const double lio_exit_speed = wheel_reports_stationary
                                      ? kExitLioSpeedWithZeroCan
                                      : kExitLioSpeed;
    const bool moving = (fresh_wheel_speed &&
                         std::abs(last_wheel_speed_mps_) > kExitWheelSpeed) ||
                        (fresh_lio && last_static_lio_speed_ > lio_exit_speed) ||
                        (!wheel_reports_stationary && inertial_motion);
    static_exit_count_ = moving ? static_exit_count_ + 1 : 0;
    if (static_exit_count_ >= kExitSamples) {
        imu_static_hold_active_ = false;
        static_exit_count_ = 0;
        LOG(WARNING) << "exit conservative IMU static hold at " << std::setprecision(14)
                     << sample.timestamp << ", gyro=" << sample.gyro_norm
                     << ", accel_delta_ratio=" << accel_delta_ratio
                     << ", lio_speed=" << last_static_lio_speed_
                     << ", wheel_speed=" << last_wheel_speed_mps_;
    }
    return imu_static_hold_active_;
}

void Localization::ProcessWheelSpeed(double timestamp, double longitudinal_speed_mps) {
    if (!imu_static_hold_enabled_ || !std::isfinite(longitudinal_speed_mps)) {
        return;
    }
    std::lock_guard<std::mutex> lock(static_detector_mutex_);
    last_wheel_speed_stamp_ = timestamp;
    last_wheel_speed_arrival_ = std::chrono::steady_clock::now();
    last_wheel_speed_mps_ = longitudinal_speed_mps;
    wheel_speed_observed_ = true;
}

void Localization::LidarLocProcCloud(CloudPtr scan_undist) {
    std::shared_lock<std::shared_mutex> lifecycle_lock(lifecycle_mutex_);
    if (lidar_loc_ == nullptr || pgo_ == nullptr) return;

    lidar_loc_->ProcessCloud(scan_undist);

    auto res = lidar_loc_->GetLocalizationResult();
    ObserveLidarLocForStaticDetector(res);
    const auto match_stats = lidar_loc_->GetLastMatchStats();
    {
        std::lock_guard<std::mutex> lock(runtime_stats_mutex_);
        if (match_stats.relocalization_attempted) ++runtime_stats_.relocalization_attempt_count;
        if (match_stats.relocalization_accepted) ++runtime_stats_.relocalization_accept_count;
        if (match_stats.relocalization_attempted) {
            runtime_stats_.last_relocalization_candidate_found = match_stats.relocalization_candidate_found;
            runtime_stats_.last_relocalization_accepted = match_stats.relocalization_accepted;
            runtime_stats_.last_relocalization_candidate_id = match_stats.relocalization_candidate_id;
            runtime_stats_.last_relocalization_score = match_stats.relocalization_score;
            runtime_stats_.last_relocalization_search_time_ms = match_stats.relocalization_search_time_ms;
            runtime_stats_.last_relocalization_reason = match_stats.relocalization_reason;
        }
    }
    if (match_stats.relocalization_accepted) {
        pgo_->Reset();
        LOG(WARNING) << "reset localization PGO after accepted global relocalization";
    }
    if (processed_cloud_callback_) {
        processed_cloud_callback_(scan_undist, res);
    }
    pgo_->ProcessLidarLoc(res);

    if (ui_) {
        // Twi with Til, here pose means Twl, thus Til=I
        ui_->UpdateScan(scan_undist, res.pose_);
    }

    if (loc_state_callback_) {
        auto loc_state = std::make_shared<std_msgs::msg::Int32>();
        loc_state->data = static_cast<int>(res.status_);
        LOG(INFO) << "loc_state: " << loc_state->data;
        loc_state_callback_(*loc_state);
    }

    // cv::Mat img(100, 100, CV_8UC3, cv::Scalar(255, 255, 255));
    // cv::imshow("img", img);
    // cv::waitKey(0);
}

void Localization::ProcessIMUMsg(IMUPtr imu) {
    std::shared_lock<std::shared_mutex> lifecycle_lock(lifecycle_mutex_);
    if (options_.online_mode_) {
        ObserveSensorEnqueued(imu ? imu->timestamp : 0.0);
        sensor_proc_.AddMessage({imu, nullptr, 0, true});
        return;
    }
    lifecycle_lock.unlock();
    ProcessIMUData(std::move(imu));
}

void Localization::ProcessIMUData(IMUPtr imu) {
    std::shared_lock<std::shared_mutex> lifecycle_lock(lifecycle_mutex_);
    UL lock(processing_mutex_);

    if (lidar_loc_ == nullptr || lio_ == nullptr || pgo_ == nullptr) {
        return;
    }

    double this_imu_time = imu->timestamp;
    if (last_imu_time_ > 0 && this_imu_time < last_imu_time_) {
        LOG(WARNING) << "IMU 时间异常：" << this_imu_time << ", last: " << last_imu_time_;
    }
    last_imu_time_ = this_imu_time;

    /// 里程计处理IMU
    lio_->ProcessIMU(imu);

    // A scan waiting for this IMU has an earlier timestamp than the current
    // IMU prediction. Publish the completed LIO output first so PGO observes
    // relative poses in sensor-time order.
    DrainLioOutputs();

    /// 这里需要 IMU predict，否则没法process DR了
    auto dr_state = lio_->GetIMUState();

    if (!dr_state.pose_is_ok_) {
        return;
    }

    if (UpdateImuStaticState(imu)) {
        dr_state.is_parking_ = true;
        dr_state.SetVel(Vec3d::Zero());
        // Prevent the high-frequency ESKF prediction from integrating a
        // stationary bias into unbounded velocity and position drift.
        lio_->SetIMUVelocity(Vec3d::Zero());
    }

    /// 如果没有odm, 用lio替代DR

    // LOG(INFO) << "dr state: " << std::setprecision(12) << dr_state.timestamp_ << " "
    //           << dr_state.GetPose().translation().transpose()
    //           << ", q=" << dr_state.GetPose().unit_quaternion().coeffs().transpose();

    // Lidar localization must continue matching while parked so freshness and
    // relocalization are based on real scans, not on a synthetic parking flag.
    auto lidar_loc_dr_state = dr_state;
    lidar_loc_dr_state.is_parking_ = false;
    lidar_loc_->ProcessDR(lidar_loc_dr_state);
    pgo_->ProcessDR(dr_state);
}

// void Localization::ProcessOdomMsg(const nav_msgs::msg::Odometry::SharedPtr odom_msg) {
//     std::shared_lock<std::shared_mutex> lifecycle_lock(lifecycle_mutex_);
//     UL lock(processing_mutex_);
//
//     if (lidar_loc_ == nullptr || lio_ == nullptr || pgo_ == nullptr) {
//         return;
//     }
//     double this_odom_time = ToSec(odom_msg->header.stamp);
//     if (last_odom_time_ > 0 && this_odom_time < last_odom_time_) {
//         LOG(WARNING) << "Odom Time Abnormal:" << this_odom_time << ", last: " << last_odom_time_;
//     }
//     last_odom_time_ = this_odom_time;
//
//     lio_->ProcessOdometry(odom_msg);
//
//     if (!lio_->GetbOdomHF()) {
//         return;
//     }
//
//     auto dr_state = lio_->GetStateHF(mapping::FasterLioMapping::kHFStateOdomFiltered);
//
//     constexpr auto kThVbrbStill = 0.03;  // 0.08;
//     constexpr auto kThOmegaStill = 0.03;
//     if (dr_state.Getvwi().norm() < kThVbrbStill && dr_state.Getwii().norm() < kThOmegaStill) {
//         dr_state.is_parking_ = true;
//         dr_state.Setvwi(Vec3d::Zero());
//         dr_state.Setwii(Vec3d::Zero());
//     }
//
//     lidar_loc_->ProcessDR(dr_state);
//     pgo_->ProcessDR(dr_state);
// }

void Localization::Finish() {
    sensor_proc_.Quit();
    lidar_loc_proc_cloud_.Quit();
    high_frequency_output_proc_.Quit();
    std::unique_lock<std::shared_mutex> lifecycle_lock(lifecycle_mutex_);
    if (lidar_loc_) lidar_loc_->Finish();
    if (ui_) {
        ui_->Quit();
    }
}

void Localization::SetExternalPose(const Eigen::Quaterniond& q, const Eigen::Vector3d& t) {
    std::shared_lock<std::shared_mutex> lifecycle_lock(lifecycle_mutex_);
    UL lock(processing_mutex_);
    /// 设置外部重定位的pose
    if (lidar_loc_) {
        lidar_loc_->SetInitialPose(SE3(q, t));
        if (pgo_) pgo_->Reset();
    }
}

void Localization::SetTFCallback(Localization::TFCallback&& callback) { tf_callback_ = callback; }

void Localization::SetLocalizationResultCallback(Localization::LocalizationResultCallback&& callback) {
    localization_result_callback_ = std::move(callback);
}

void Localization::SetGlobalLocalizationResultCallback(Localization::LocalizationResultCallback&& callback) {
    global_localization_result_callback_ = std::move(callback);
}

void Localization::SetProcessedCloudCallback(Localization::ProcessedCloudCallback&& callback) {
    processed_cloud_callback_ = std::move(callback);
}

void Localization::SetLocStateCallback(Localization::LocStateCallback&& callback) {
    loc_state_callback_ = std::move(callback);
}

}  // namespace lightning::loc
