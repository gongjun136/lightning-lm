#include <pcl/common/transforms.h>
#include <pcl_conversions/pcl_conversions.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iterator>
#include <limits>
#include <sstream>

#include "common/debug_event.h"
#include "core/localization/lidar_loc/lidar_loc.h"
#include "core/localization/localization.h"

#include <opencv2/highgui.hpp>

#include "core/localization/pose_graph/pgo.h"
#include "io/yaml_io.h"
#include "ui/pangolin_window.h"
#include "utils/causal_trace.h"

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
    live_output_timestamp_gate_.Reset();
    last_imu_time_ = 0.0;

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
        static_wheel_speed_history_.clear();
        last_wheel_speed_stamp_ = 0.0;
        wheel_speed_observed_ = false;
        wheel_timestamp_mismatch_reported_ = false;
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
    lidar_loc_proc_cloud_.SetProcFunc([this](const LidarLocInput& input) { LidarLocProcCloud(input); });

    if (options_.online_mode_) {
        {
            std::lock_guard<std::mutex> arrival_lock(arrival_mutex_);
            primary_arrivals_.clear();
        }
        sensor_proc_.Start();
        lidar_loc_proc_cloud_.Start();
    }

    /// TODO: 发布
    auto publish_localization_result = [this](const LocalizationResult& res) {
        std::lock_guard<std::mutex> dispatch_lock(live_output_dispatch_mutex_);
        const auto timestamp_decision = live_output_timestamp_gate_.Observe(res.timestamp_);
        if (!timestamp_decision.accepted) {
            std::ostringstream message;
            message << "Dropped non-monotonic live localization output: timestamp="
                    << std::setprecision(16) << res.timestamp_
                    << ", last_accepted=" << timestamp_decision.reference_timestamp
                    << ", rollback_sec=" << timestamp_decision.lag_sec;
            LOG(WARNING) << message.str();
            debug_event::EmitThrottled("live_output_timestamp_rollback", message.str(),
                                       std::chrono::seconds(1));
            return;
        }
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
    const auto received_at = std::chrono::steady_clock::now();
    const bool profiling_enabled = profiling::ComputeProfilingEnabled();
    profiling::Stopwatch outer_profile_timer(profiling_enabled);
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

    RememberPrimaryArrival(lidar_id, static_cast<std::uint64_t>(cloud->header.stamp.sec) * 1000000000ULL +
                               cloud->header.stamp.nanosec, received_at);
    // 串行模式
    CloudPtr laser_cloud(new PointCloudType);
    profiling::Stopwatch preprocess_profile_timer(profiling_enabled);
    preprocess_->Process(cloud, laser_cloud);
    const profiling::TimingSample preprocess_timing = preprocess_profile_timer.Stop();
    laser_cloud->header.stamp = cloud->header.stamp.sec * 1e9 + cloud->header.stamp.nanosec;

    bool process_lidar_odom = true;
    if (options_.enable_lidar_odom_skip_) {
        const int skip_num = options_.lidar_odom_skip_num_ > 0 ? options_.lidar_odom_skip_num_ : 1;
        process_lidar_odom = lidar_odom_skip_cnt_ == 0;
        lidar_odom_skip_cnt_ = (lidar_odom_skip_cnt_ + 1) % skip_num;
    }
    if (!process_lidar_odom) return;
    profiling::Stopwatch dispatch_profile_timer(profiling_enabled);
    if (options_.online_mode_) {
        ObserveSensorEnqueued(static_cast<double>(laser_cloud->header.stamp) * 1e-9);
        sensor_proc_.AddMessage({nullptr, laser_cloud, lidar_id, false,
                                profiling_enabled ? std::chrono::steady_clock::now()
                                                  : std::chrono::steady_clock::time_point{}});
        input_lock.unlock();
        lifecycle_lock.unlock();
    } else {
        input_lock.unlock();
        lifecycle_lock.unlock();
        LidarOdomProcCloud(laser_cloud, lidar_id);
    }
    const profiling::TimingSample dispatch_timing = dispatch_profile_timer.Stop();
    const profiling::TimingSample outer_timing = outer_profile_timer.Stop();
    if (profiling_enabled) {
        const auto summary =
            lidar_input_timing_window_.Add({preprocess_timing, dispatch_timing, outer_timing});
        if (summary) {
            LOG(INFO) << "COMPUTE_BENCH_SUMMARY module=lidar_input source=pointcloud2"
                      << " timestamp_s=" << timestamp << " lidar_id=" << lidar_id
                      << " last_output_points=" << laser_cloud->size()
                      << ' ' << profiling::FormatTimingWindow(*summary)
                      << " sensor_queue_pending=" << sensor_proc_.PendingCount()
                      << " sensor_queue_dropped_total=" << sensor_proc_.DroppedCount();
        }
    }
}

void Localization::ProcessLivoxLidarMsg(const livox_ros_driver2::msg::CustomMsg::SharedPtr cloud) {
    const auto received_at = std::chrono::steady_clock::now();
    const bool profiling_enabled = profiling::ComputeProfilingEnabled();
    profiling::Stopwatch outer_profile_timer(profiling_enabled);
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
    profiling::Stopwatch preprocess_profile_timer(profiling_enabled);
    preprocess_->Process(cloud, laser_cloud);
    const profiling::TimingSample preprocess_timing = preprocess_profile_timer.Stop();
    laser_cloud->header.stamp = cloud->header.stamp.sec * 1e9 + cloud->header.stamp.nanosec;

    const int lidar_id = lio_->IsMultiLidarEnabled() ? lio_->GetMultiLidarConfig().primary_lidar_id : 0;
    RememberPrimaryArrival(lidar_id, laser_cloud->header.stamp, received_at);
    bool process_lidar_odom = true;
    if (options_.enable_lidar_odom_skip_) {
        const int skip_num = options_.lidar_odom_skip_num_ > 0 ? options_.lidar_odom_skip_num_ : 1;
        process_lidar_odom = lidar_odom_skip_cnt_ == 0;
        lidar_odom_skip_cnt_ = (lidar_odom_skip_cnt_ + 1) % skip_num;
    }
    if (!process_lidar_odom) return;
    profiling::Stopwatch dispatch_profile_timer(profiling_enabled);
    if (options_.online_mode_) {
        ObserveSensorEnqueued(static_cast<double>(laser_cloud->header.stamp) * 1e-9);
        sensor_proc_.AddMessage({nullptr, laser_cloud, lidar_id, false,
                                profiling_enabled ? std::chrono::steady_clock::now()
                                                  : std::chrono::steady_clock::time_point{}});
        input_lock.unlock();
        lifecycle_lock.unlock();
    } else {
        input_lock.unlock();
        lifecycle_lock.unlock();
        LidarOdomProcCloud(laser_cloud, lidar_id);
    }
    const profiling::TimingSample dispatch_timing = dispatch_profile_timer.Stop();
    const profiling::TimingSample outer_timing = outer_profile_timer.Stop();
    if (profiling_enabled) {
        const auto summary =
            lidar_input_timing_window_.Add({preprocess_timing, dispatch_timing, outer_timing});
        if (summary) {
            LOG(INFO) << "COMPUTE_BENCH_SUMMARY module=lidar_input source=livox_custom"
                      << " timestamp_s=" << timestamp << " lidar_id=" << lidar_id
                      << " last_output_points=" << laser_cloud->size()
                      << ' ' << profiling::FormatTimingWindow(*summary)
                      << " sensor_queue_pending=" << sensor_proc_.PendingCount()
                      << " sensor_queue_dropped_total=" << sensor_proc_.DroppedCount();
        }
    }
}

void Localization::RememberPrimaryArrival(int lidar_id, std::uint64_t stamp,
                                         std::chrono::steady_clock::time_point received_at) {
    const int primary = lio_->IsMultiLidarEnabled() ? lio_->GetMultiLidarConfig().primary_lidar_id : 0;
    if (lidar_id != primary) return;
    std::lock_guard<std::mutex> lock(arrival_mutex_);
    // Keep the first arrival for duplicate timestamps. At 10 Hz the bound
    // covers 25.6 seconds; drops never turn this into an unbounded ledger.
    primary_arrivals_.emplace(stamp, received_at);
    while (primary_arrivals_.size() > 256) primary_arrivals_.erase(primary_arrivals_.begin());
}

void Localization::ProcessSensorInput(const SensorInput& input) {
    const bool trace = profiling::CausalTraceEnabled();
    profiling::Stopwatch trace_timer(trace);
    const double trace_queue_ms = trace && input.enqueued_at != std::chrono::steady_clock::time_point{}
        ? std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - input.enqueued_at).count()
        : -1.0;
    if (profiling::ComputeProfilingEnabled() &&
        input.enqueued_at != std::chrono::steady_clock::time_point{}) {
        const auto dequeued_at = std::chrono::steady_clock::now();
        profiling::TimingSample age;
        age.valid = true;
        age.wall_ms = std::chrono::duration<double, std::milli>(dequeued_at - input.enqueued_at).count();
        auto& window = input.is_imu ? imu_queue_timing_window_ : lidar_queue_timing_window_;
        if (const auto summary = window.Add({age}, dequeued_at)) {
            const auto& wall = summary->stages.front().wall_ms;
            LOG(INFO) << std::fixed << std::setprecision(6)
                      << "COMPUTE_BENCH_SUMMARY module=sensor_queue source="
                      << (input.is_imu ? "imu" : "lidar")
                      << " window_s=" << summary->window_s << " samples=" << wall.count
                      << " queue_wait_mean_ms=" << wall.mean
                      << " queue_wait_p95_ms=" << wall.p95
                      << " queue_wait_p99_ms=" << wall.p99
                      << " queue_wait_max_ms=" << wall.max;
        }
    }
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
    if (trace) {
        const auto elapsed = trace_timer.Stop();
        // Normal IMU callbacks need no detailed trace; retain every lidar and
        // slow IMU callback, without altering existing all-message summaries.
        if (!input.is_imu || trace_queue_ms >= 5.0 || elapsed.wall_ms >= 20.0) {
            profiling::CausalEvent event;
            event.kind = "sensor_process"; event.source = input.is_imu ? "imu" : "lidar";
            event.stamp = timestamp;
            event.frame_id = input.cloud ? input.cloud->header.stamp : 0;
            event.duration_ms = elapsed.wall_ms; event.thread_cpu_ms = elapsed.thread_cpu_ms;
            event.queue_wait_ms = trace_queue_ms;
            profiling::RecordCausalEvent(event);
        }
    }
}

void Localization::LidarOdomProcCloud(CloudPtr cloud, int lidar_id) {
    // LIO must observe IMU and lidar in callback order. Running this update in
    // a second worker lets later IMU callbacks overtake the cloud and changes
    // the filter result relative to offline processing.
    const auto wait_started = profiling::CausalTraceEnabled() ? std::chrono::steady_clock::now()
                                                            : std::chrono::steady_clock::time_point{};
    const double trace_stamp = cloud ? cloud->header.stamp * 1e-9 : 0.0;
    std::shared_lock<std::shared_mutex> lifecycle_lock(lifecycle_mutex_);
    profiling::TraceLockWait("lidar_lifecycle", trace_stamp, wait_started);
    const auto processing_started = wait_started == std::chrono::steady_clock::time_point{}
        ? wait_started : std::chrono::steady_clock::now();
    UL processing_lock(processing_mutex_);
    profiling::TraceLockWait("lidar_processing", trace_stamp, processing_started);

    if (lio_ == nullptr) {
        return;
    }
    double latest_input_stamp = 0.0;
    {
        std::lock_guard<std::mutex> stats_lock(runtime_stats_mutex_);
        latest_input_stamp = runtime_stats_.latest_enqueued_sensor_stamp;
    }
    lio_->SetLatestInputTimestamp(latest_input_stamp);

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
        LidarLocInput lidar_loc_input;
        lidar_loc_input.localization_cloud = scan;
        lidar_loc_input.publication_cloud = lio_->GetPublicationCloud();
        lidar_loc_input.frame_stats = lio_->GetCurrentFrameStats();
        {
            std::lock_guard<std::mutex> lock(arrival_mutex_);
            const auto stamp = scan->header.stamp;
            // The fused header is the earliest source, up to match_tolerance
            // before the primary header. Half a microsecond covers double conversion.
            auto arrival = primary_arrivals_.lower_bound(stamp > 500 ? stamp - 500 : 0);
            const auto tolerance = static_cast<std::uint64_t>(
                lio_->GetMultiLidarConfig().match_tolerance * 1e9) + 500;
            if (arrival != primary_arrivals_.end() && arrival->first <= stamp + tolerance) {
                lidar_loc_input.primary_received_at = arrival->second;
                primary_arrivals_.erase(primary_arrivals_.begin(), std::next(arrival));
            }
        }
        lidar_loc_input.publication_eligible = lio_->CanPublishCurrentCloud();
        const auto& selection = lio_->GetCurrentLidarSelection();
        const auto& multi_config = lio_->GetMultiLidarConfig();
        const auto wheel_stats = lio_->GetWheelSpeedDrStats();
        {
            std::lock_guard<std::mutex> stats_lock(runtime_stats_mutex_);
            runtime_stats_.adaptive_lidar_load_enabled =
                multi_config.adaptive_load.enabled;
            runtime_stats_.adaptive_lidar_degradation_step =
                lio_->GetAdaptiveLidarDegradationStep();
            runtime_stats_.selected_lidar_point_stride = selection.point_stride;
            runtime_stats_.selected_lidar_ids = selection.lidar_ids;
            runtime_stats_.current_frame_lidar_ids =
                lidar_loc_input.frame_stats.present_lidar_ids;
            runtime_stats_.current_frame_missing_lidar_ids =
                lidar_loc_input.frame_stats.missing_lidar_ids;
            runtime_stats_.lidar_correction_age_sec = lio_->GetLastLidarLatencySec();
            runtime_stats_.last_lio_processing_ms = lio_->GetLastFrameProcessingMs();
            runtime_stats_.adaptive_lidar_stale_drop_count =
                lio_->GetAdaptiveStaleDropCount();
            runtime_stats_.cloud_publish_min_lidars = static_cast<std::uint32_t>(
                multi_config.adaptive_load.cloud_publish_min_lidars);
            runtime_stats_.cloud_publish_eligible = lidar_loc_input.publication_eligible;
            runtime_stats_.wheel_speed_dr_enabled = lio_->GetWheelSpeedDrConfig().enabled;
            runtime_stats_.wheel_speed_dr_stats = wheel_stats;
        }

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
                if (profiling::ComputeProfilingEnabled()) {
                    lidar_loc_input.enqueued_at = std::chrono::steady_clock::now();
                }
                lidar_loc_proc_cloud_.AddMessage(lidar_loc_input);
            } else {
                LidarLocProcCloud(lidar_loc_input);
            }
        } else {
            // auto scan = cloud;   // 这个cloud应该差一个外参

            if (options_.online_mode_) {
                if (profiling::ComputeProfilingEnabled()) {
                    lidar_loc_input.enqueued_at = std::chrono::steady_clock::now();
                }
                lidar_loc_proc_cloud_.AddMessage(lidar_loc_input);
            } else {
                LidarLocProcCloud(lidar_loc_input);
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
    RuntimeStats stats;
    {
        std::lock_guard<std::mutex> lock(runtime_stats_mutex_);
        stats = runtime_stats_;
    }
    {
        std::shared_lock<std::shared_mutex> lifecycle_lock(lifecycle_mutex_);
        if (lio_ != nullptr) {
            stats.lidar_correction_age_sec = lio_->GetLastLidarLatencySec();
            stats.last_lio_processing_ms = lio_->GetLastFrameProcessingMs();
            stats.adaptive_lidar_stale_drop_count =
                lio_->GetAdaptiveStaleDropCount();
            stats.wheel_speed_dr_enabled = lio_->GetWheelSpeedDrConfig().enabled;
            stats.wheel_speed_dr_stats = lio_->GetWheelSpeedDrStats();
        }
    }
    stats.sensor_queue_pending = sensor_proc_.PendingCount();
    stats.sensor_queue_dropped = sensor_proc_.DroppedCount();
    stats.sensor_queue_processed = sensor_proc_.ProcessedCount();
    stats.localization_queue_pending = lidar_loc_proc_cloud_.PendingCount();
    stats.localization_queue_dropped = lidar_loc_proc_cloud_.DroppedCount();
    stats.localization_queue_processed = lidar_loc_proc_cloud_.ProcessedCount();
    stats.high_frequency_queue_pending = high_frequency_output_proc_.PendingCount();
    stats.high_frequency_queue_dropped = high_frequency_output_proc_.DroppedCount();
    stats.high_frequency_queue_processed = high_frequency_output_proc_.ProcessedCount();
    stats.live_output_non_monotonic_drop_count =
        live_output_timestamp_gate_.RejectedCount();
    stats.worst_live_output_timestamp_rollback_sec =
        live_output_timestamp_gate_.WorstRollbackSec();
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
    debug_event::ReportState("sensor_queue_overload", throttled,
                             "Sensor queue overloaded",
                             "Sensor queue recovered from overload",
                             std::chrono::seconds(2));

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
    bool severe_rollback = false;
    double rollback = 0.0;
    {
        std::lock_guard<std::mutex> lock(runtime_stats_mutex_);
        constexpr double kSevereRollbackSec = 1.0;
        rollback = runtime_stats_.latest_processed_sensor_stamp - timestamp;
        if (runtime_stats_.latest_processed_sensor_stamp > 0.0 &&
            rollback > kSevereRollbackSec) {
            severe_rollback = true;
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
            std::max(runtime_stats_.max_sensor_lag_sec,
                     runtime_stats_.current_sensor_lag_sec);
    }
    if (severe_rollback) {
        std::ostringstream message;
        message << "Processed sensor timestamp rolled back severely: rollback_sec="
                << rollback << ", timestamp=" << std::setprecision(16) << timestamp;
        debug_event::EmitThrottled("processed_sensor_timestamp_rollback", message.str(),
                                   std::chrono::seconds(1));
    }
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
    // CAN callbacks may be ahead of the serialized IMU worker. Select the
    // latest causal sample, not the latest arrival (and never a future zero).
    const StaticWheelSpeedSample* wheel_sample = nullptr;
    for (auto it = static_wheel_speed_history_.rbegin(); it != static_wheel_speed_history_.rend(); ++it) {
        if (it->timestamp <= sample.timestamp) {
            wheel_sample = &*it;
            break;
        }
    }
    const double wheel_age = wheel_sample ? sample.timestamp - wheel_sample->timestamp
                                         : std::numeric_limits<double>::infinity();
    const double wheel_speed_mps = wheel_sample ? wheel_sample->speed_mps : 0.0;
    const double wheel_imu_timestamp_delta = last_wheel_speed_stamp_ - sample.timestamp;
    const bool fresh_wheel_speed = wheel_sample && wheel_age <= kMaxWheelSpeedAgeSec;
    if (wheel_speed_observed_ && !fresh_wheel_speed &&
        !wheel_timestamp_mismatch_reported_) {
        wheel_timestamp_mismatch_reported_ = true;
        LOG(WARNING) << "ignore CAN wheel speed: latest CAN-IMU Header delta="
                     << std::setprecision(14) << wheel_imu_timestamp_delta
                     << " sec, selected_history_age=" << wheel_age
                     << " sec; no causal sample within " << kMaxWheelSpeedAgeSec
                     << " sec; using IMU/LIO fallback";
    } else if (fresh_wheel_speed && wheel_timestamp_mismatch_reported_) {
        wheel_timestamp_mismatch_reported_ = false;
        LOG(INFO) << "CAN wheel-speed history recovered: selected_history_age="
                  << std::setprecision(14) << wheel_age
                  << " sec, latest CAN-IMU delta=" << wheel_imu_timestamp_delta << " sec";
    }
    debug_event::ReportState("can_imu_timestamp_mismatch",
                             wheel_speed_observed_ && !fresh_wheel_speed,
                             "No fresh causal CAN sample at IMU epoch; using IMU/LIO fallback",
                             "Fresh causal CAN sample recovered",
                             std::chrono::seconds(1));
    const bool wheel_reports_stationary =
        fresh_wheel_speed && std::abs(wheel_speed_mps) < kEnterWheelSpeed;
    const auto trace_decision = [&](bool hold, const char* reason) {
        if (!profiling::CausalTraceEnabled()) return;
        profiling::CausalEvent event;
        event.kind = "static_decision"; event.source = "static_detector"; event.reason = reason;
        event.stamp = sample.timestamp;
        if (wheel_sample) {
            event.input_stamp = wheel_sample->timestamp;
            event.measured_v = wheel_speed_mps;
        }
        event.hold = hold;
        profiling::RecordCausalEvent(event);
    };

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
            trace_decision(true, fresh_wheel_speed ? "can_stationary" : "imu_stationary");
            LOG(WARNING) << "enter conservative IMU static hold at " << std::setprecision(14)
                         << sample.timestamp << ", lio_speed=" << last_static_lio_speed_
                         << ", wheel_speed="
                         << (fresh_wheel_speed ? wheel_speed_mps : 0.0)
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
                         std::abs(wheel_speed_mps) > kExitWheelSpeed) ||
                        (fresh_lio && last_static_lio_speed_ > lio_exit_speed) ||
                        (!wheel_reports_stationary && inertial_motion);
    static_exit_count_ = moving ? static_exit_count_ + 1 : 0;
    if (static_exit_count_ >= kExitSamples) {
        imu_static_hold_active_ = false;
        static_exit_count_ = 0;
        // These flags describe the final decision sample, not all three debounce samples.
        const unsigned flags = (fresh_wheel_speed && std::abs(wheel_speed_mps) > kExitWheelSpeed ? 1u : 0u) |
            (fresh_lio && last_static_lio_speed_ > lio_exit_speed ? 2u : 0u) |
            (!wheel_reports_stationary && inertial_motion ? 4u : 0u);
        static const char* reasons[] = {"none", "can", "lio", "can+lio", "imu", "can+imu", "lio+imu", "can+lio+imu"};
        trace_decision(false, reasons[flags]);
        LOG(WARNING) << "exit conservative IMU static hold at " << std::setprecision(14)
                     << sample.timestamp << ", gyro=" << sample.gyro_norm
                     << ", accel_delta_ratio=" << accel_delta_ratio
                     << ", lio_speed=" << last_static_lio_speed_
                     << ", wheel_speed=" << wheel_speed_mps;
    }
    return imu_static_hold_active_;
}

void Localization::ProcessWheelSpeed(double timestamp, double longitudinal_speed_mps,
                                     double motor_torque_nm) {
    if (!std::isfinite(timestamp) || timestamp <= 0.0 || !std::isfinite(longitudinal_speed_mps) ||
        !std::isfinite(motor_torque_nm)) return;
    if (profiling::CausalTraceEnabled()) {
        profiling::CausalEvent event;
        event.kind = "can_input"; event.source = "callback";
        event.stamp = timestamp; event.measured_v = longitudinal_speed_mps;
        event.torque = motor_torque_nm;
        profiling::RecordCausalEvent(event);
    }
    std::shared_lock<std::shared_mutex> lifecycle_lock(lifecycle_mutex_);
    if (lio_ != nullptr) {
        lio_->ProcessWheelSpeed(timestamp, longitudinal_speed_mps, motor_torque_nm);
    }
    if (!imu_static_hold_enabled_) return;
    std::lock_guard<std::mutex> lock(static_detector_mutex_);
    if (timestamp <= last_wheel_speed_stamp_) return;
    last_wheel_speed_stamp_ = timestamp;
    static_wheel_speed_history_.push_back({timestamp, longitudinal_speed_mps});
    constexpr double kHistorySec = 2.0;
    constexpr std::size_t kMaxHistorySamples = 512;
    while (!static_wheel_speed_history_.empty() &&
           (timestamp - static_wheel_speed_history_.front().timestamp > kHistorySec ||
            static_wheel_speed_history_.size() > kMaxHistorySamples)) {
        static_wheel_speed_history_.pop_front();
    }
    wheel_speed_observed_ = true;
}

void Localization::LidarLocProcCloud(const LidarLocInput& input) {
    const bool profiling_enabled = profiling::ComputeProfilingEnabled();
    const double queue_wait_ms = profiling_enabled &&
                                        input.enqueued_at != std::chrono::steady_clock::time_point{}
                                    ? std::chrono::duration<double, std::milli>(
                                          std::chrono::steady_clock::now() - input.enqueued_at).count()
                                    : -1.0;
    profiling::Stopwatch outer_profile_timer(profiling_enabled);
    std::shared_lock<std::shared_mutex> lifecycle_lock(lifecycle_mutex_);
    if (lidar_loc_ == nullptr || pgo_ == nullptr) return;
    const CloudPtr& scan_undist = input.localization_cloud;

    profiling::Stopwatch process_cloud_profile_timer(profiling_enabled);
    lidar_loc_->ProcessCloud(scan_undist);
    const profiling::TimingSample process_cloud_timing = process_cloud_profile_timer.Stop();

    profiling::Stopwatch result_profile_timer(profiling_enabled);
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
    // LidarLoc owns lidar_loc_valid_; valid_ belongs to fused/DR outputs.
    lio_->SetLocalizationGood(res.lidar_loc_valid_ && res.status_ == LocalizationStatus::GOOD);
    const profiling::TimingSample result_timing = result_profile_timer.Stop();
    profiling::Stopwatch cloud_callback_profile_timer(profiling_enabled);
    if (processed_cloud_callback_) {
        processed_cloud_callback_(input.publication_cloud, res, input.frame_stats,
                                  input.publication_eligible);
    }
    const profiling::TimingSample cloud_callback_timing = cloud_callback_profile_timer.Stop();
    profiling::Stopwatch pgo_profile_timer(profiling_enabled);
    pgo_->ProcessLidarLoc(res);
    const profiling::TimingSample pgo_timing = pgo_profile_timer.Stop();

    double primary_to_pgo_ms = -1.0;
    if (input.primary_received_at != std::chrono::steady_clock::time_point{}) {
        const auto completed_at = std::chrono::steady_clock::now();
        primary_to_pgo_ms = std::chrono::duration<double, std::milli>(
            completed_at - input.primary_received_at).count();
        if (res.lidar_loc_valid_ && res.status_ == LocalizationStatus::GOOD && !match_stats.relocalization_accepted) {
            if (primary_to_pgo_ms >= 100.0) ++lidar_deadline_misses_;
            profiling::TimingSample sample;
            sample.valid = true;
            sample.wall_ms = primary_to_pgo_ms;
            const auto summary = lidar_end_to_end_window_.Add({sample}, completed_at);
            if (summary) {
                const auto& wall = summary->stages.front().wall_ms;
                // Low-rate production telemetry; CPU is not measured by this wall-clock ledger.
                LOG(INFO) << std::fixed << std::setprecision(6)
                          << "COMPUTE_BENCH_SUMMARY module=lidar_end_to_end phase=tracking"
                          << " window_s=" << summary->window_s << " samples=" << wall.count
                          << " mean_ms=" << wall.mean << " p95_ms=" << wall.p95
                          << " p99_ms=" << wall.p99 << " max_ms=" << wall.max
                          << " healthy_observation_rate_hz=" << wall.count / summary->window_s
                          << " deadline_misses_total=" << lidar_deadline_misses_
                          << " sensor_queue_dropped_total=" << sensor_proc_.DroppedCount()
                          << " lidar_queue_dropped_total=" << lidar_loc_proc_cloud_.DroppedCount();
            }
        }
    }

    if (ui_) {
        // Twi with Til, here pose means Twl, thus Til=I
        ui_->UpdateScan(scan_undist, res.pose_);
    }

    if (loc_state_callback_) {
        auto loc_state = std::make_shared<std_msgs::msg::Int32>();
        loc_state->data = static_cast<int>(res.status_);
        if (!profiling::ReduceNonessentialOverhead()) {
            LOG(INFO) << "loc_state: " << loc_state->data;
        }
        loc_state_callback_(*loc_state);
    }
    lifecycle_lock.unlock();

    if (profiling_enabled) {
        const profiling::TimingSample outer_timing = outer_profile_timer.Stop();
        LOG(INFO) << std::fixed << std::setprecision(6) << "COMPUTE_BENCH_FRAME module=lidar_loc_pipeline"
                  << " frame_id=" << scan_undist->header.stamp
                  << " primary_to_pgo_ms=" << primary_to_pgo_ms
                  << " queue_wait_ms=" << queue_wait_ms
                  << " timestamp_s=" << res.timestamp_
                  << " localization_points=" << (scan_undist ? scan_undist->size() : 0)
                  << " lidar_valid=" << res.lidar_loc_valid_
                  << " publication_points="
                  << (input.publication_cloud ? input.publication_cloud->size() : 0)
                  << " status=" << static_cast<int>(res.status_)
                  << " valid=" << res.valid_
                  << " confidence=" << res.confidence_
                  << " relocalization_attempted=" << match_stats.relocalization_attempted
                  << " relocalization_accepted=" << match_stats.relocalization_accepted
                  << ' ' << profiling::FormatTimingSample("process_cloud", process_cloud_timing)
                  << ' ' << profiling::FormatTimingSample("result_bookkeeping", result_timing)
                  << ' ' << profiling::FormatTimingSample("cloud_callback", cloud_callback_timing)
                  << ' ' << profiling::FormatTimingSample("pgo", pgo_timing)
                  << ' ' << profiling::FormatTimingSample("outer", outer_timing)
                  << " queue_pending=" << lidar_loc_proc_cloud_.PendingCount()
                  << " queue_dropped_total=" << lidar_loc_proc_cloud_.DroppedCount()
                  << " queue_processed_total=" << lidar_loc_proc_cloud_.ProcessedCount();
    }

    // cv::Mat img(100, 100, CV_8UC3, cv::Scalar(255, 255, 255));
    // cv::imshow("img", img);
    // cv::waitKey(0);
}

void Localization::ProcessIMUMsg(IMUPtr imu) {
    std::shared_lock<std::shared_mutex> lifecycle_lock(lifecycle_mutex_);
    if (options_.online_mode_) {
        ObserveSensorEnqueued(imu ? imu->timestamp : 0.0);
        sensor_proc_.AddMessage({imu, nullptr, 0, true,
                                profiling::ComputeProfilingEnabled() ? std::chrono::steady_clock::now()
                                                                     : std::chrono::steady_clock::time_point{}});
        return;
    }
    lifecycle_lock.unlock();
    ProcessIMUData(std::move(imu));
}

void Localization::ProcessIMUData(IMUPtr imu) {
    const bool profiling_enabled = profiling::ComputeProfilingEnabled();
    profiling::Stopwatch outer_profile_timer(profiling_enabled);
    const auto wait_started = profiling::CausalTraceEnabled() ? std::chrono::steady_clock::now()
                                                            : std::chrono::steady_clock::time_point{};
    const double trace_stamp = imu ? imu->timestamp : 0.0;
    std::shared_lock<std::shared_mutex> lifecycle_lock(lifecycle_mutex_);
    profiling::TraceLockWait("imu_lifecycle", trace_stamp, wait_started);
    const auto processing_started = wait_started == std::chrono::steady_clock::time_point{}
        ? wait_started : std::chrono::steady_clock::now();
    UL lock(processing_mutex_);
    profiling::TraceLockWait("imu_processing", trace_stamp, processing_started);

    if (lidar_loc_ == nullptr || lio_ == nullptr || pgo_ == nullptr) {
        return;
    }

    if (!imu || !std::isfinite(imu->timestamp) || imu->timestamp <= 0.0 ||
        !imu->angular_velocity.allFinite() || !imu->linear_acceleration.allFinite()) return;
    double this_imu_time = imu->timestamp;
    if (last_imu_time_ > 0 && this_imu_time < last_imu_time_) {
        LOG(WARNING) << "IMU 时间异常：" << this_imu_time << ", last: " << last_imu_time_;
        std::ostringstream message;
        message << "Rejected IMU timestamp rollback: rollback_sec="
                << last_imu_time_ - this_imu_time;
        debug_event::EmitThrottled("imu_timestamp_rollback", message.str(),
                                   std::chrono::seconds(1));
        return;
    }
    if (this_imu_time == last_imu_time_) return;
    last_imu_time_ = this_imu_time;

    /// 里程计处理IMU
    double latest_input_stamp = 0.0;
    {
        std::lock_guard<std::mutex> stats_lock(runtime_stats_mutex_);
        latest_input_stamp = runtime_stats_.latest_enqueued_sensor_stamp;
    }
    lio_->SetLatestInputTimestamp(latest_input_stamp);
    profiling::Stopwatch lio_imu_profile_timer(profiling_enabled);
    lio_->ProcessIMU(imu);
    const profiling::TimingSample lio_imu_timing = lio_imu_profile_timer.Stop();

    // A scan waiting for this IMU has an earlier timestamp than the current
    // IMU prediction. Publish the completed LIO output first so PGO observes
    // relative poses in sensor-time order.
    profiling::Stopwatch drain_lio_profile_timer(profiling_enabled);
    DrainLioOutputs();
    const profiling::TimingSample drain_lio_timing = drain_lio_profile_timer.Stop();

    /// 这里需要 IMU predict，否则没法process DR了
    profiling::Stopwatch state_profile_timer(profiling_enabled);
    auto dr_state = lio_->GetIMUState();

    profiling::TimingSample state_timing;
    profiling::TimingSample lidar_loc_dr_timing;
    profiling::TimingSample pgo_dr_timing;
    const auto emit_imu_profile = [&](const char* phase) {
        const profiling::TimingSample outer_timing = outer_profile_timer.Stop();
        lock.unlock();
        lifecycle_lock.unlock();
        if (!profiling_enabled) return;
        const auto summary = imu_dr_timing_window_.Add(
            {lio_imu_timing, drain_lio_timing, state_timing, lidar_loc_dr_timing,
             pgo_dr_timing, outer_timing});
        if (summary) {
            LOG(INFO) << "COMPUTE_BENCH_SUMMARY module=imu_dr_pipeline"
                      << " phase=" << phase
                      << " timestamp_s=" << (imu ? imu->timestamp : 0.0)
                      << ' ' << profiling::FormatTimingWindow(*summary)
                      << " sensor_queue_pending=" << sensor_proc_.PendingCount()
                      << " sensor_queue_dropped_total=" << sensor_proc_.DroppedCount()
                      << " lidar_loc_queue_pending=" << lidar_loc_proc_cloud_.PendingCount()
                      << " lidar_loc_queue_dropped_total=" << lidar_loc_proc_cloud_.DroppedCount()
                      << " output_queue_pending=" << high_frequency_output_proc_.PendingCount()
                      << " output_queue_dropped_total=" << high_frequency_output_proc_.DroppedCount()
                      << " output_queue_processed_total=" << high_frequency_output_proc_.ProcessedCount();
        }
    };

    if (!dr_state.pose_is_ok_) {
        state_timing = state_profile_timer.Stop();
        emit_imu_profile("waiting_for_pose");
        return;
    }

    const bool static_hold_active = UpdateImuStaticState(imu);
    debug_event::ReportState(
        "imu_static_hold", static_hold_active,
        "Entered conservative IMU static hold",
        "Exited conservative IMU static hold",
        std::chrono::seconds(1));
    lio_->SetIMUStaticHold(static_hold_active);
    dr_state = lio_->GetIMUState();
    if (static_hold_active) {
        dr_state.is_parking_ = true;
    }
    state_timing = state_profile_timer.Stop();

    /// 如果没有odm, 用lio替代DR

    // LOG(INFO) << "dr state: " << std::setprecision(12) << dr_state.timestamp_ << " "
    //           << dr_state.GetPose().translation().transpose()
    //           << ", q=" << dr_state.GetPose().unit_quaternion().coeffs().transpose();

    // Lidar localization must continue matching while parked so freshness and
    // relocalization are based on real scans, not on a synthetic parking flag.
    auto lidar_loc_dr_state = dr_state;
    lidar_loc_dr_state.is_parking_ = false;
    profiling::Stopwatch lidar_loc_dr_profile_timer(profiling_enabled);
    lidar_loc_->ProcessDR(lidar_loc_dr_state);
    lidar_loc_dr_timing = lidar_loc_dr_profile_timer.Stop();
    profiling::Stopwatch pgo_dr_profile_timer(profiling_enabled);
    pgo_->ProcessDR(dr_state);
    pgo_dr_timing = pgo_dr_profile_timer.Stop();
    emit_imu_profile("tracking");
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
