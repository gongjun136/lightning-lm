#include <pcl/common/transforms.h>
#include <pcl_conversions/pcl_conversions.h>

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
    UL lock(global_mutex_);
    if (lidar_loc_ != nullptr) {
        // 若已经启动，则变为初始化
        // Finish() joins the sensor worker, which may itself be waiting for
        // global_mutex_. Do not hold that mutex while waiting for the worker.
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

    sensor_proc_.SetMaxSize(10000);
    // Multi-lidar relocalization can spend several seconds building the first
    // BTC query. Keep that bounded startup backlog instead of dropping the
    // scans that offline mode naturally retains while processing is paused.
    const size_t default_loc_queue_size = lio_->IsMultiLidarEnabled() ? 100 : 10;
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
    pgo_->SetHighFrequencyGlobalOutputHandleFunction(publish_localization_result);

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
    const int lidar_id = lio_ && lio_->IsMultiLidarEnabled() ? lio_->GetMultiLidarConfig().primary_lidar_id : 0;
    ProcessLidarMsg(cloud, lidar_id);
}

void Localization::ProcessLidarMsg(const sensor_msgs::msg::PointCloud2::SharedPtr cloud, int lidar_id) {
    UL lock(global_mutex_);
    if (lidar_loc_ == nullptr || lio_ == nullptr || pgo_ == nullptr) {
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
        sensor_proc_.AddMessage({nullptr, laser_cloud, lidar_id, false});
        return;
    }
    lock.unlock();
    LidarOdomProcCloud(laser_cloud, lidar_id);
}

void Localization::ProcessLivoxLidarMsg(const livox_ros_driver2::msg::CustomMsg::SharedPtr cloud) {
    UL lock(global_mutex_);
    if (lidar_loc_ == nullptr || lio_ == nullptr || pgo_ == nullptr) {
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
        sensor_proc_.AddMessage({nullptr, laser_cloud, lidar_id, false});
        return;
    }
    lock.unlock();
    LidarOdomProcCloud(laser_cloud, lidar_id);
}

void Localization::ProcessSensorInput(const SensorInput& input) {
    if (input.is_imu) {
        ProcessIMUData(input.imu);
    } else {
        LidarOdomProcCloud(input.cloud, input.lidar_id);
    }
}

void Localization::LidarOdomProcCloud(CloudPtr cloud, int lidar_id) {
    // LIO must observe IMU and lidar in callback order. Running this update in
    // a second worker lets later IMU callbacks overtake the cloud and changes
    // the filter result relative to offline processing.
    UL processing_lock(global_mutex_);

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

void Localization::LidarLocProcCloud(CloudPtr scan_undist) {
    lidar_loc_->ProcessCloud(scan_undist);

    auto res = lidar_loc_->GetLocalizationResult();
    if (lidar_loc_->GetLastMatchStats().relocalization_accepted) {
        pgo_->Reset();
        LOG(WARNING) << "reset localization PGO after accepted BTC relocalization";
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
    if (options_.online_mode_) {
        sensor_proc_.AddMessage({imu, nullptr, 0, true});
        return;
    }
    ProcessIMUData(std::move(imu));
}

void Localization::ProcessIMUData(IMUPtr imu) {
    UL lock(global_mutex_);

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

    // /// 停车判定
    // constexpr auto kThVbrbStill = 0.05;  // 0.08;
    // constexpr auto kThOmegaStill = 0.05;

    // if (dr_state.GetVel().norm() < kThVbrbStill && imu->angular_velocity.norm() < kThOmegaStill) {
    //     dr_state.is_parking_ = true;
    //     dr_state.SetVel(Vec3d::Zero());
    // }

    /// 如果没有odm, 用lio替代DR

    // LOG(INFO) << "dr state: " << std::setprecision(12) << dr_state.timestamp_ << " "
    //           << dr_state.GetPose().translation().transpose()
    //           << ", q=" << dr_state.GetPose().unit_quaternion().coeffs().transpose();

    lidar_loc_->ProcessDR(dr_state);
    pgo_->ProcessDR(dr_state);
}

// void Localization::ProcessOdomMsg(const nav_msgs::msg::Odometry::SharedPtr odom_msg) {
//     UL lock(global_mutex_);
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
    if (lidar_loc_) lidar_loc_->Finish();
    if (ui_) {
        ui_->Quit();
    }
}

void Localization::SetExternalPose(const Eigen::Quaterniond& q, const Eigen::Vector3d& t) {
    UL lock(global_mutex_);
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
