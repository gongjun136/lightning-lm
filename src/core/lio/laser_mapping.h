#ifndef FASTER_LIO_LASER_MAPPING_H
#define FASTER_LIO_LASER_MAPPING_H

#include <pcl/filters/voxel_grid.h>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <limits>
#include <mutex>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <thread>

#include "common/eigen_types.h"
#include "common/imu.h"
#include "common/imu_body_velocity.h"
#include "common/keyframe.h"
#include "common/options.h"
#include "core/ivox3d/ivox3d.h"
#include "core/lio/eskf.hpp"
#include "core/lio/imu_processing.hpp"
#include "core/lio/multi_lidar_fusion.h"
#include "pointcloud_preprocess.h"

#include "livox_ros_driver2/msg/custom_msg.hpp"

namespace lightning {

namespace ui {
class PangolinWindow;
}

struct WheelSpeedDrConfig {
    bool enabled = false;  ///< 是否启用轮速观测更新。
    bool use_body_axis = true;  ///< 是否沿车体纵轴构建速度观测。
    bool lever_arm_enabled = false;  ///< 是否进行主雷达到后轴的速度杆臂补偿。
    double lever_gyro_std_radps = 0.02;  ///< 杆臂速度补偿使用的角速度噪声标准差（rad/s）。
    Vec3d rear_to_lidar_body = Vec3d(2.199, 0.0, 2.740);  ///< 后轴到主雷达的车体系杆臂向量（m）。
    double base_std_mps = 0.35;  ///< 正常行驶时的轮速观测标准差（m/s）。
    double stationary_std_mps = 0.08;  ///< 静止时的轮速观测标准差（m/s）。
    double stationary_speed_threshold_mps = 0.03;  ///< 判定近静止轮速的阈值（m/s）。
    double max_age_sec = 0.25;  ///< 轮速观测相对更新时刻允许的最大年龄（s）。
    double future_tolerance_sec = 0.0;  ///< 允许轮速时间戳超前更新时刻的容差（s）。
    double max_acceleration_mps2 = 4.0;  ///< 相邻轮速观测允许的最大加速度（m/s²）。
    double max_abs_innovation_mps = 1.5;  ///< 轮速观测允许的最大绝对新息（m/s）。
    double normalized_innovation_squared_gate = 9.0;  ///< 归一化新息平方的拒绝阈值。
    double torque_reference_nm = 1000.0;  ///< 根据扭矩调整轮速噪声的参考值（N·m）。
    double torque_std_scale = 3.0;  ///< 扭矩较大时对观测标准差的缩放系数。
    double max_velocity_step_mps = 0.35;  ///< 轮速更新允许的最大速度修正（m/s）。
};

struct WheelSpeedDrStats {
    std::uint64_t input_count = 0;  ///< 接收到的轮速观测数。
    std::uint64_t invalid_input_count = 0;  ///< 因数值无效而拒绝的轮速观测数。
    std::uint64_t timestamp_reject_count = 0;  ///< 因时间戳无效或不递增而拒绝的轮速观测数。
    std::uint64_t acceleration_reject_count = 0;  ///< 因加速度超限而拒绝的轮速观测数。
    std::uint64_t lidar_filter_accepted_count = 0;  ///< 点云时刻 ESKF 接受的轮速更新数。
    std::uint64_t lidar_filter_rejected_count = 0;  ///< 点云时刻 ESKF 拒绝的轮速更新数。
    std::uint64_t lidar_filter_stale_count = 0;  ///< 点云时刻 ESKF 遇到的过期轮速数。
    std::uint64_t imu_filter_accepted_count = 0;  ///< 高频 IMU ESKF 接受的轮速更新数。
    std::uint64_t imu_filter_rejected_count = 0;  ///< 高频 IMU ESKF 拒绝的轮速更新数。
    std::uint64_t imu_filter_stale_count = 0;  ///< 高频 IMU ESKF 遇到的过期轮速数。
    double last_measurement_mps = 0.0;  ///< 最近参与门控的轮速量测（m/s）。
    double last_predicted_mps = 0.0;  ///< 最近轮速更新前的预测速度（m/s）。
    double last_innovation_mps = 0.0;  ///< 最近轮速量测与预测之差（m/s）。
    double last_standard_deviation_mps = 0.0;  ///< 最近轮速更新使用的观测标准差（m/s）。
    double last_normalized_innovation_squared = 0.0;  ///< 最近轮速更新的归一化新息平方。
};

/**
 * @brief LIO前端主流程：点云预处理、IMU同步、ESKF更新和局部地图维护。
 *
 * LaserMapping串起整个前端定位流程：
 * 1. 接收ROS/Livox/已预处理点云和IMU数据，并放入缓存队列；
 * 2. 将一帧Lidar与覆盖其扫描时间的IMU同步成MeasureGroup；
 * 3. 调用ImuProcess完成IMU预测和点云运动畸变补偿；
 * 4. 用当前点云和IVox局部地图构造Lidar观测，迭代更新ESKF；
 * 5. 根据运动阈值创建关键帧，并增量更新局部地图。
 *
 * 同步上有一个容易误解的点：bag中的点云时间戳通常是扫描开始时间，
 * 但处理一帧点云需要等到扫描结束时间之后的IMU。因此缓存中最新点云往往还不能处理，
 * 实际处理的通常是队列中已经等到足够IMU覆盖的那一帧。
 */
/** @see @ref laser_mapping_module */
class LaserMapping {
   public:
    enum class RunStatus { kNoData, kConsumed, kOutput };

    /**
     * @brief LaserMapping运行选项。
     *
     * 这些参数主要控制前端配准权重、关键帧创建阈值和关键帧投影策略。
     * 传感器、体素分辨率、噪声等运行参数由Init()从yaml中读取。
     */
    struct Options {
        Options() {}

        bool is_in_slam_mode_ = true;  ///< SLAM模式会保存关键帧；定位模式下可用时间阈值强制创建关键帧

        bool enable_icp_part_ = true;    ///< 是否在点面约束外额外加入点到点ICP约束
        double plane_icp_weight_ = 1.0;  ///< 点面ICP残差在ESKF观测信息矩阵中的权重
        double icp_weight_ = 100;        ///< 点到点ICP残差在ESKF观测信息矩阵中的权重

        int min_pts = 300;  ///< 触发有效前端处理所需的最少点数。

        /// 关键帧创建阈值：当前帧相对上一关键帧的平移或旋转超过阈值时创建新关键帧。
        double kf_dis_th_ = 2.0;                    // 平移阈值，单位m
        double kf_angle_th_ = 15 * M_PI / 180.0;    ///< 旋转阈值，单位rad

        bool proj_kfs_ = false;  ///< 是否将附近关键帧点云投影到当前帧，用于显示或辅助匹配
        int max_proj_kfs_ = 5;   ///< 保留用于投影的关键帧数量上限
    };

    /// 类内有Eigen固定大小矩阵成员，使用对齐new避免内存对齐问题。
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /// IVox局部地图类型，保存世界系点云并支持快速最近邻查询。
    using IVoxType = IVox<3, IVoxNodeType::DEFAULT, PointType>;

    LaserMapping(Options options = Options());
    ~LaserMapping() {
        scan_down_lidar_ = nullptr;
        scan_undistort_full_ = nullptr;
        scan_undistort_ = nullptr;
        scan_down_world_ = nullptr;
        LOG(INFO) << "laser mapping deconstruct";
    }

    /// 从yaml加载参数并初始化预处理器、IMU处理器、IVox局部地图和ESKF观测函数。
    bool Init(const std::string &config_yaml);

    /// 处理缓存中的一帧同步数据，完成去畸变、ESKF更新、关键帧判断和地图维护。
    bool Run();
    RunStatus RunDetailed();

    // 三个ProcessPointCloud2函数处理逻辑：
    // 1.时间戳检查和回环检测：检查时间戳是否倒退，如果是则清空缓冲区
    // 2.点云预处理：调用 preprocess_->Process() 进行统一格式 (CloudPtr)转换和预处理
    // 3.数据缓存：将处理后的点云和时间戳存入缓冲区
    // 4.性能监控：使用 Timer::Evaluate 记录预处理耗时
    // callbacks of lidar and imu
    /// 处理标准ROS2 PointCloud2点云消息，完成预处理后写入Lidar缓存队列。
    bool ProcessPointCloud2(const sensor_msgs::msg::PointCloud2::SharedPtr &msg);
    bool ProcessPointCloud2(const sensor_msgs::msg::PointCloud2::SharedPtr &msg, int lidar_id);

    /// 处理Livox自定义点云消息，完成预处理后写入Lidar缓存队列。
    bool ProcessPointCloud2(const livox_ros_driver2::msg::CustomMsg::SharedPtr &msg);
    bool ProcessPointCloud2(const livox_ros_driver2::msg::CustomMsg::SharedPtr &msg, int lidar_id);

    /// 处理已经转换好的点云，直接写入Lidar缓存队列。
    bool ProcessPointCloud2(CloudPtr cloud);
    bool ProcessPointCloud2(CloudPtr cloud, int lidar_id);

    void FlushMultiLidar();

    /// 处理一条IMU消息：写入IMU缓存，并在IMU初始化后维护一份高频kf_imu_状态供UI显示。
    void ProcessIMU(const lightning::IMUPtr &msg_in);

    /// Buffer a signed body-forward speed converted from motor rpm.
    void ProcessWheelSpeed(double timestamp, double longitudinal_speed_mps, double motor_torque_nm);
    const WheelSpeedDrConfig& GetWheelSpeedDrConfig() const { return wheel_speed_dr_config_; }
    WheelSpeedDrStats GetWheelSpeedDrStats() const;
    /// 保存前端的地图
    void SaveMap();

    /// 绑定UI窗口，用于在Run()/ProcessIMU()中更新状态和点云显示。
    void SetUI(std::shared_ptr<ui::PangolinWindow> ui) { ui_ = ui; }

    /// 获取关键帧
    Keyframe::Ptr GetKeyframe() const { return last_kf_; }

    /// 获取激光的状态
    NavState GetState() const { return state_point_; }
    bool IsInitialized() const { return !flg_first_scan_ && p_imu_->IsIMUInited(); }
    bool IsTrackingHealthy() const { return last_tracking_healthy_; }
    bool IsMultiLidarEnabled() const { return multi_lidar_config_.enabled; }
    const MultiLidarConfig &GetMultiLidarConfig() const { return multi_lidar_config_; }
    const MultiLidarFrameStats &GetCurrentFrameStats() const { return current_lidar_stats_; }
    CloudPtr GetPublicationCloud() const { return scan_undistort_full_; }
    bool CanPublishCurrentCloud() const {
        return adaptive_lidar_load_controller_.CanPublishCloud(current_lidar_stats_);
    }
    void SetLocalizationGood(bool good) { adaptive_lidar_load_controller_.SetLocalizationGood(good); }
    void SetLatestInputTimestamp(double timestamp) {
        double current = latest_input_sensor_timestamp_.load(std::memory_order_relaxed);
        while (timestamp > current &&
               !latest_input_sensor_timestamp_.compare_exchange_weak(
                   current, timestamp, std::memory_order_relaxed)) {
        }
    }
    int GetAdaptiveLidarDegradationStep() const {
        return adaptive_lidar_load_controller_.DegradationStep();
    }
    const AdaptiveLidarSelection& GetCurrentLidarSelection() const {
        return current_lidar_selection_;
    }
    double GetLastLidarLatencySec() const {
        return last_lidar_latency_sec_.load(std::memory_order_relaxed);
    }
    double GetLastFrameProcessingMs() const {
        return last_frame_processing_ms_.load(std::memory_order_relaxed);
    }
    std::size_t GetAdaptiveStaleDropCount() const {
        return adaptive_stale_drop_count_.load(std::memory_order_relaxed);
    }
    std::size_t GetMultiLidarLateDropCount() const { return multi_lidar_assembler_.LateDropCount(); }
    std::size_t GetMultiLidarDuplicateDropCount() const { return multi_lidar_assembler_.DuplicateDropCount(); }
    std::size_t GetMultiLidarToleranceDropCount() const { return multi_lidar_assembler_.ToleranceDropCount(); }
    std::size_t GetMultiLidarInvalidDropCount() const { return multi_lidar_assembler_.InvalidDropCount(); }
    std::size_t GetMultiLidarAssembledFrameCount() const { return multi_lidar_assembler_.EmittedFrameCount(); }
    std::size_t GetMultiLidarInsufficientDropCount() const {
        return multi_lidar_assembler_.InsufficientLidarDropCount();
    }
    std::size_t GetPendingLidarFrameCount() const { return lidar_buffer_.size(); }
    std::size_t GetPreImuDropCount() const { return pre_imu_drop_count_; }
    double GetLastFrameBeginTime() const { return measures_.lidar_begin_time_; }
    double GetLastFrameEndTime() const { return measures_.lidar_end_time_; }
    const Mat3d &GetLidarToImuRotation() const { return offset_R_lidar_fixed_; }
    const Vec3d &GetLidarToImuTranslation() const { return offset_t_lidar_fixed_; }
    // Initial level vehicle defines body axes, matching publication convention.
    // This is not an extrinsic calibration for an arbitrary sloped start.
    SO3 GetImuToBodyRotation() const { return p_imu_->GetInitialRotation(); }
    SO3 GetInitialLidarRotation() const {
        return p_imu_->GetInitialRotation() *
               SO3(Eigen::Quaterniond(offset_R_lidar_fixed_).normalized());
    }

    /// 获取IMU最新时刻状态；IMU未初始化时返回pose_is_ok_=false的无效状态。
    NavState GetIMUState() const {
        if (p_imu_->IsIMUInited()) {
            auto state = kf_imu_.GetX();
            state.rear_axle_speed_offset_ = high_frequency_static_hold_ ? 0.0 : GetRearAxleSpeedOffset(state);
            return state;
        } else {
            NavState s;
            s.pose_is_ok_ = false;
            return s;
        }
    }

    /// One-shot velocity assignment; use SetIMUStaticHold for a persistent hold.
    void SetIMUVelocity(const Vec3d& velocity) {
        auto state = kf_imu_.GetX();
        state.SetVel(velocity);
        kf_imu_.ChangeX(state);
    }

    /// Persist the external detector's hold across prediction and LIO rebuilds.
    /// The main LIO filter remains unconstrained so motion can release the hold.
    void SetIMUStaticHold(bool active);

    /// 获取最近一次去畸变后的点云，点仍在当前Lidar坐标系下。
    CloudPtr GetScanUndist() const { return scan_undistort_; }
    /// 获取当前去畸变点云叠加附近投影关键帧后的点云。
    CloudPtr GetProjCloud();

    /// 获取最新的点云
    CloudPtr GetRecentCloud();

    /// 获取全部关键帧，用于后端、保存地图或调试显示。
    std::vector<Keyframe::Ptr> GetAllKeyframes() { return all_keyframes_; }

    /// Copy a keyframe cloud and remove map-export-only articulated vehicle points.
    CloudPtr PrepareMapExportCloud(const CloudPtr& cloud) const;

    /**
     * @brief 根据关键帧点云拼接全局地图。
     *
     * @param use_lio_pose true时使用前端LIO位姿，false时使用关键帧优化位姿。
     * @param use_voxel 是否对每个关键帧点云和最终地图做体素滤波。
     * @param res 体素滤波分辨率，单位m。
     * @return 拼接后的全局点云。
     */
    CloudPtr GetGlobalMap(bool use_lio_pose, bool use_voxel = true,
                          float res = 0.1,
                          bool apply_map_export_filter = false);

   private:
    /// 从Lidar/IMU缓存队列中取出一帧时间覆盖完整的同步数据，写入measures_。
    bool SyncPackages();

    /// Lidar观测模型：根据当前状态匹配局部地图，构造ESKF需要的HTH/HTr。
    void ObsModel(NavState &s, ESKF::CustomObservationModel &obs);

    /// 将当前帧Lidar坐标系下的点变换到世界系，使用state_point_和固定外参。
    inline void PointBodyToWorld(const PointType &pi, PointType &po) {
        Vec3d p_global(state_point_.rot_ *
                           (offset_R_lidar_fixed_ * pi.getVector3fMap().cast<double>() + offset_t_lidar_fixed_) +
                       state_point_.pos_);

        po.x = p_global(0);
        po.y = p_global(1);
        po.z = p_global(2);
        po.intensity = pi.intensity;
        po.time = pi.time;
        po.lidar_id = pi.lidar_id;
    }

    bool EnqueueCloud(double timestamp, CloudPtr cloud, const MultiLidarFrameStats *stats = nullptr,
                      double preprocess_ms = 0.0);
    friend class LaserMappingTimingTestPeer;
    bool DrainAssembledFrames();
    void RebuildHighFrequencyState();
    void ApplyHighFrequencyStaticHold();
    double GetRearAxleSpeedOffset(const NavState& state) const;
    bool ApplyWheelSpeedObservation(ESKF& filter, double state_timestamp,
                                    double& last_applied_observation_timestamp,
                                    bool high_frequency_filter,
                                    const char* trace_source);
    void ResetWheelSpeedIntegrationBridge(double timestamp);
    int gyro_prefilter_samples_ = 1;  ///< 陀螺仪预滤波的累计样本数。
    std::deque<Vec3d> raw_gyro_window_;  ///< 原始角速度样本窗口，用于预滤波统计。
    struct WheelSpeedSample {
        double timestamp = 0.0;  ///< 轮速观测的传感器时间戳（s）。
        double speed_mps = 0.0;  ///< 带符号的车体纵向轮速（m/s）。
        double torque_nm = 0.0;  ///< 与轮速同帧的电机扭矩（N·m）。
    };
    double PointInformationScale(const PointType &point, const Vec3d &plane_normal_world,
                                 const NavState &state) const;

    /// 将当前帧降采样点云增量加入IVox局部地图，并根据邻近点做自适应下采样。
    void MapIncremental();

    /// 从yaml读取传感器、预处理、IVox、ESKF噪声和关键帧等参数。
    bool LoadParamsFromYAML(const std::string &yaml);

    /// 创建关键帧
    void MakeKF();

    /// 将附近的关键帧投影至cloud中
    void ProjectKFs(CloudPtr cloud, int size_limit = 1000);

   private:
    Options options_;  ///< 前端工作模式、约束权重与关键帧配置。

    /// 核心模块
    IVoxType::Options ivox_options_;
    std::shared_ptr<IVoxType> ivox_ = nullptr;                    ///< IVox局部地图，支持地图点插入和最近邻查询
    std::shared_ptr<PointCloudPreprocess> preprocess_ = nullptr;  ///< 点云预处理模块，统一不同雷达消息格式
    std::shared_ptr<ImuProcess> p_imu_ = nullptr;                 ///< IMU初始化、预测和点云去畸变模块

    /// 局部地图相关
    double filter_size_map_min_ = 0;  // 地图体素滤波分辨率（m），控制地图点的密度

    /// Lidar-IMU固定外参和地图保存路径
    std::vector<double> extrinT_{3, 0.0};  // yaml读取的Lidar到IMU平移外参
    std::vector<double> extrinR_{9, 0.0};  ///< yaml读取的Lidar到IMU旋转外参
    Mat3d offset_R_lidar_fixed_ = Mat3d::Identity();  ///< Lidar到IMU的旋转矩阵形式
    Vec3d offset_t_lidar_fixed_ = Vec3d::Zero();      ///< Lidar到IMU的平移向量形式
    std::string map_file_path_;  ///< 地图文件的保存路径。
    double filter_size_scan_ = 0.0;  ///< 当前扫描的体素降采样分辨率（m）。
    MultiLidarConfig multi_lidar_config_;  ///< 雷达外参、组帧与自适应降载配置。
    MultiLidarFrameAssembler multi_lidar_assembler_;  ///< 按时间窗口合并多雷达输入的组帧器。
    SelfPointFilterConfig self_point_filter_config_;  ///< 前端点云的车体自反射过滤区域。
    SelfPointFilterConfig map_export_self_point_filter_config_;  ///< 地图导出时使用的车体自反射过滤区域。

    bool point_noise_enabled_ = false;  ///< 是否根据量测距离与角度噪声调整点面信息权重。
    double range_noise_sigma_ = 0.02;  ///< 雷达距离量测噪声标准差（m）。
    double angular_noise_sigma_rad_ = 0.05 * M_PI / 180.0;  ///< 雷达角度量测噪声标准差（rad）。
    double noise_reference_sigma_ = 0.05;  ///< 点面信息权重归一化使用的参考噪声标准差（m）。
    double min_information_scale_ = 0.25;  ///< 单点信息权重缩放的下限。
    double max_information_scale_ = 4.0;  ///< 单点信息权重缩放的上限。

    std::vector<Keyframe::Ptr> all_keyframes_;  ///< 所有关键帧的存储列表
    Keyframe::Ptr last_kf_ = nullptr;           ///< 最近的关键帧指针（用于快速访问）
    int kf_id_ = 0;                             ///< 关键帧ID计数器（唯一标识每个关键帧）

    /// 当前帧点云缓存
    CloudPtr scan_undistort_full_{new PointCloudType()};  // 所有有效雷达的去畸变点云，用于发布
    CloudPtr scan_undistort_{new PointCloudType()};       ///< 动态选择后的定位点云
    CloudPtr scan_down_lidar_{new PointCloudType()};  ///< 当前帧降采样点云，Lidar坐标系
    CloudPtr scan_down_world_{new PointCloudType()};  ///< 当前帧降采样点云，世界坐标系
    pcl::VoxelGrid<PointType> voxel_scan_;            ///< 当前帧点云体素滤波器

    /// 点面相关
    std::vector<PointVector> nearest_points_;  // 当前帧每个点在IVox地图中的最近邻点
    std::vector<Vec4f> corr_pts_;              ///< 点面内点：[x,y,z,残差]，点坐标在Lidar系
    std::vector<std::uint8_t> corr_lidar_ids_;  ///< 与点面内点对应的原始雷达 ID。
    std::vector<Vec4f> corr_norm_;             ///< 点面内点对应平面：[nx,ny,nz,d]，平面在世界系
    std::vector<float> residuals_;             ///< 点到平面有符号残差
    // Fixed blocks reduce allocation and keep summation independent of worker count.
    std::vector<Mat6d> observation_hessian_blocks_;  ///< 按固定点块累计的观测 Hessian，便于并行归并。
    std::vector<Vec6d> observation_gradient_blocks_;  ///< 与观测 Hessian 点块对应的梯度累加值。
    std::vector<double> observation_residual_squared_;  ///< 各点块的残差平方和。
    int observation_evaluations_ = 0;  ///< 当前帧观测模型的计算次数。
    std::vector<char> point_selected_surf_;    ///< 点面约束是否有效
    std::vector<Vec4f> plane_coef_;            ///< 每个点局部拟合出的平面参数

    /// 点到点相关
    std::vector<char> point_selected_icp_;  // 点到点ICP约束是否有效

    std::mutex mtx_buffer_;          ///< 保护Lidar/IMU缓存队列的互斥锁
    std::deque<double> time_buffer_; ///< 与lidar_buffer_一一对应的点云起始时间

    std::deque<PointCloudType::Ptr> lidar_buffer_;  ///< 激光雷达数据缓冲队列（用于与IMU时间同步）
    std::deque<lightning::IMUPtr> imu_buffer_;      ///< IMU数据缓冲队列（高频数据，用于状态预测）
    IMUPtr high_frequency_previous_imu_;  ///< Raw left endpoint, independent of the filter epoch.
    IMUPtr last_synchronized_imu_;        ///< Retain the sample removed at the scan boundary.
    std::deque<MultiLidarFrameStats> lidar_stats_buffer_;  ///< 与点云队列对应的多雷达组帧统计。
    std::deque<double> preprocess_time_buffer_ms_;  ///< 与点云队列对应的预处理耗时（ms）。
    double current_preprocess_ms_ = 0.0;  ///< 当前同步帧的预处理耗时（ms）。

    /// options
    bool keep_first_imu_estimation_ = false;  // 在没有建立地图前，是否要使用前几帧的IMU状态
    double timediff_lidar_wrt_imu_ = 0.0;     ///< 激光雷达与IMU之间的时间偏移量，用于时间同步
    double last_timestamp_lidar_ = 0;         ///< 上一帧激光雷达数据的结束时间戳，用于检查数据间隙
    double lidar_end_time_ = 0;               ///< 当前帧激光雷达的结束时间戳（本地计算，用于IMU同步）
    double last_timestamp_imu_ = -1.0;        ///< 上一帧处理过的IMU数据时间戳，用于IMU数据筛选
    double first_lidar_time_ = 0.0;           ///< 系统启动后第一帧激光雷达的时间戳，作为参考基准
    bool lidar_pushed_ = false;               ///< 标记当前帧的激光雷达数据是否已添加到处理队列中

    bool enable_skip_lidar_ = true;  ///< 雷达是否需要跳帧
    int skip_lidar_num_ = 5;         ///< 每隔多少帧跳一个雷达
    int skip_lidar_cnt_ = 0;         ///< 跳帧计数器

    /// statistics and flags ///
    int scan_count_ = 0;                // 总扫描帧数统计
    int publish_count_ = 0;             ///< 发布次数统计（用于控制发布频率）
    bool flg_first_scan_ = true;        ///< 是否为第一帧扫描（用于初始化判断）
    bool flg_EKF_inited_ = false;       ///< ESKF滤波器是否已初始化（影响是否进行观测更新）
    double lidar_mean_scantime_ = 0.0;  ///< 激光雷达平均扫描时间（用于时间统计和性能监控）
    int scan_num_ = 0;                  ///< 当前扫描序列号
    int effect_feat_surf_ = 0 /**< 本帧有效点面观测数量。 */, frame_num_ = 0 /**< 前端累计处理的帧数。 */, effect_feat_icp_ = 0 /**< 本帧有效点到点 ICP 观测数量。 */;
    MultiLidarFrameStats current_lidar_stats_;  ///< 当前已同步点云的雷达完整性与点数统计。
    AdaptiveLidarLoadController adaptive_lidar_load_controller_;  ///< 根据处理耗时与跟踪健康度调整雷达和点数预算。
    AdaptiveLidarSelection current_lidar_selection_;  ///< 当前帧选中的雷达、采样步长和点数上限。
    std::atomic<double> latest_input_sensor_timestamp_{0.0};  ///< 最近输入的传感器时间戳（s），用于处理滞后估计。
    std::atomic<double> last_lidar_latency_sec_{0.0};  ///< 最近点云输出相对最新输入的滞后（s）。
    std::atomic<double> last_frame_processing_ms_{0.0};  ///< 最近一帧前端处理耗时（ms）。
    std::atomic<std::size_t> adaptive_stale_drop_count_{0};  ///< 因自适应时效门槛丢弃的点云帧数。
    bool last_tracking_healthy_ = false;  ///< 最近一帧前端跟踪是否健康。
    double current_max_imu_gap_ = 0.0;  ///< 当前同步帧内最大的 IMU 时间间隙（s）。
    std::size_t pre_imu_drop_count_ = 0;  ///< 因尚无可用 IMU 而丢弃的点云帧数。

    WheelSpeedDrConfig wheel_speed_dr_config_;  ///< 轮速观测的噪声、时效、杆臂及新息门控配置。
    mutable std::mutex wheel_speed_mutex_;  ///< 保护轮速观测队列与更新统计。
    std::deque<WheelSpeedSample> wheel_speed_buffer_;  ///< 带传感器时间戳的轮速观测队列。
    WheelSpeedDrStats wheel_speed_dr_stats_;  ///< 轮速输入校验与两路 ESKF 更新统计。
    double last_raw_wheel_speed_timestamp_ = std::numeric_limits<double>::lowest();  ///< 最近原始轮速观测的时间戳（s）。
    double last_raw_wheel_speed_mps_ = 0.0;  ///< 最近原始轮速观测值（m/s），用于加速度检查。
    double last_lidar_filter_wheel_timestamp_ = std::numeric_limits<double>::lowest();  ///< 点云时刻 ESKF 最近使用的轮速观测时间（s）。
    double last_imu_filter_wheel_timestamp_ = std::numeric_limits<double>::lowest();  ///< 高频 IMU ESKF 最近使用的轮速观测时间（s）。

    double last_lidar_time_ = 0;  ///< 上一帧激光雷达时间戳（用于时间同步和断流检测）

    ///////////////////////// EKF inputs and output ///////////////////////////////////////////////////////
    MeasureGroup measures_;  ///< SyncPackages输出的一帧同步Lidar和IMU数据

    ESKF kf_;      ///< 点云时刻的IMU状态，用于畸变矫正+雷达里程计观测更新
    ESKF kf_imu_;  ///< imu 最新时刻的eskf状态，提供UI的高频位姿输出
    bool high_frequency_static_hold_ = false;  ///< 高频 ESKF 是否保持静止状态。

    NavState state_point_;  ///< 当前Lidar帧结束时刻的前端状态

    bool use_aa_ = false;  ///< ESKF观测更新是否使用Anderson Acceleration
    bool propagate_velocity_ = false;  ///< 是否在ESKF名义状态中传播速度
    bool lidar_update_pose_only_ = false;  ///< Lidar观测是否只修正位姿
    bool lidar_update_inertial_states_ = true;  ///< Lidar观测是否间接修正bg/ba/gravity
    int lio_parallel_threads_ = 1;  ///< LIO匹配及Livox预处理共享的OpenMP线程上限
    double max_update_velocity_step_ = 2.0;  ///< 一次观测更新允许的速度改变量上限（m/s）。
    double max_update_gyro_bias_step_ = 0.05;  ///< 一次观测更新允许的陀螺仪偏置改变量上限（rad/s）。
    double max_update_acc_bias_step_ = 0.5;  ///< 一次观测更新允许的加速度计偏置改变量上限（m/s²）。
    double max_update_gravity_step_ = 0.05;  ///< 一次观测更新允许的重力状态改变量上限。
    bool adaptive_velocity_propagation_ = false;  ///< 是否根据速度新息自适应开启速度传播。
    double velocity_innovation_ema_alpha_ = 0.05;  ///< 速度新息指数滑动均值的更新系数。
    double velocity_innovation_enable_threshold_ = 0.16;  ///< 触发速度传播的速度新息均值阈值。
    double velocity_innovation_disable_threshold_ = 0.10;  ///< 关闭速度传播的速度新息均值阈值。
    int velocity_propagation_max_active_updates_ = 200;  ///< 速度传播连续启用的最大更新次数。
    int velocity_propagation_cooldown_updates_ = 100;  ///< 关闭速度传播后的冷却更新次数。
    bool velocity_innovation_initialized_ = false;  ///< 速度新息滑动统计是否已初始化。
    double velocity_innovation_ema_ = 0.0;  ///< 当前速度新息的指数滑动均值。
    bool velocity_propagation_active_ = false;  ///< 当前是否开启速度传播。
    int velocity_propagation_active_updates_ = 0;  ///< 本轮速度传播已持续的更新次数。
    int velocity_propagation_cooldown_remaining_ = 0;  ///< 距离允许再次开启速度传播的剩余更新次数。
    bool velocity_propagation_safety_lockout_ = false;  ///< 因安全条件而锁定速度传播的状态。

    std::list<Keyframe::Ptr> proj_kfs_;  ///< 投影到当前帧的关键帧

    std::shared_ptr<ui::PangolinWindow> ui_ = nullptr;  ///< 可选UI，用于显示状态和点云
};

}  // namespace lightning

#endif  // FASTER_LIO_LASER_MAPPING_H
