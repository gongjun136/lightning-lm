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
    bool enabled = false;
    double base_std_mps = 0.35;
    double stationary_std_mps = 0.08;
    double stationary_speed_threshold_mps = 0.03;
    double max_age_sec = 0.25;
    double future_tolerance_sec = 0.03;
    double max_acceleration_mps2 = 4.0;
    double max_abs_innovation_mps = 1.5;
    double normalized_innovation_squared_gate = 9.0;
    double torque_reference_nm = 1000.0;
    double torque_std_scale = 3.0;
    double max_velocity_step_mps = 0.35;
};

struct WheelSpeedDrStats {
    std::uint64_t input_count = 0;
    std::uint64_t invalid_input_count = 0;
    std::uint64_t timestamp_reject_count = 0;
    std::uint64_t acceleration_reject_count = 0;
    std::uint64_t lidar_filter_accepted_count = 0;
    std::uint64_t lidar_filter_rejected_count = 0;
    std::uint64_t lidar_filter_stale_count = 0;
    std::uint64_t imu_filter_accepted_count = 0;
    std::uint64_t imu_filter_rejected_count = 0;
    std::uint64_t imu_filter_stale_count = 0;
    double last_measurement_mps = 0.0;
    double last_predicted_mps = 0.0;
    double last_innovation_mps = 0.0;
    double last_standard_deviation_mps = 0.0;
    double last_normalized_innovation_squared = 0.0;
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

        bool is_in_slam_mode_ = true;  // SLAM模式会保存关键帧；定位模式下可用时间阈值强制创建关键帧

        bool enable_icp_part_ = true;    // 是否在点面约束外额外加入点到点ICP约束
        double plane_icp_weight_ = 1.0;  // 点面ICP残差在ESKF观测信息矩阵中的权重
        double icp_weight_ = 100;        // 点到点ICP残差在ESKF观测信息矩阵中的权重

        int min_pts = 300;  // 降采样后进入配准流程所需的最少点数

        /// 关键帧创建阈值：当前帧相对上一关键帧的平移或旋转超过阈值时创建新关键帧。
        double kf_dis_th_ = 2.0;                    // 平移阈值，单位m
        double kf_angle_th_ = 15 * M_PI / 180.0;    // 旋转阈值，单位rad

        bool proj_kfs_ = false;  // 是否将附近关键帧点云投影到当前帧，用于显示或辅助匹配
        int max_proj_kfs_ = 5;   // 保留用于投影的关键帧数量上限
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
    SO3 GetInitialLidarRotation() const {
        return p_imu_->GetInitialRotation() *
               SO3(Eigen::Quaterniond(offset_R_lidar_fixed_).normalized());
    }

    /// 获取IMU最新时刻状态；IMU未初始化时返回pose_is_ok_=false的无效状态。
    NavState GetIMUState() const {
        if (p_imu_->IsIMUInited()) {
            return kf_imu_.GetX();
        } else {
            NavState s;
            s.pose_is_ok_ = false;
            return s;
        }
    }

    /// Apply a conservative zero-velocity constraint to the high-frequency
    /// prediction state after an external stationary detector has fired.
    void SetIMUVelocity(const Vec3d& velocity) {
        auto state = kf_imu_.GetX();
        state.SetVel(velocity);
        kf_imu_.ChangeX(state);
    }

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
    bool DrainAssembledFrames();
    bool ApplyWheelSpeedObservation(ESKF& filter, double state_timestamp,
                                    double& last_applied_observation_timestamp,
                                    bool high_frequency_filter);
    void ResetWheelSpeedIntegrationBridge(double timestamp);
    struct WheelSpeedSample {
        double timestamp = 0.0;
        double speed_mps = 0.0;
        double torque_nm = 0.0;
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
    Options options_;

    /// 核心模块
    IVoxType::Options ivox_options_;
    std::shared_ptr<IVoxType> ivox_ = nullptr;                    // IVox局部地图，支持地图点插入和最近邻查询
    std::shared_ptr<PointCloudPreprocess> preprocess_ = nullptr;  // 点云预处理模块，统一不同雷达消息格式
    std::shared_ptr<ImuProcess> p_imu_ = nullptr;                 // IMU初始化、预测和点云去畸变模块

    /// 局部地图相关
    double filter_size_map_min_ = 0;  // 地图体素滤波分辨率（m），控制地图点的密度

    /// Lidar-IMU固定外参和地图保存路径
    std::vector<double> extrinT_{3, 0.0};  // yaml读取的Lidar到IMU平移外参
    std::vector<double> extrinR_{9, 0.0};  // yaml读取的Lidar到IMU旋转外参
    Mat3d offset_R_lidar_fixed_ = Mat3d::Identity();  // Lidar到IMU的旋转矩阵形式
    Vec3d offset_t_lidar_fixed_ = Vec3d::Zero();      // Lidar到IMU的平移向量形式
    std::string map_file_path_;
    double filter_size_scan_ = 0.0;
    MultiLidarConfig multi_lidar_config_;
    MultiLidarFrameAssembler multi_lidar_assembler_;
    SelfPointFilterConfig self_point_filter_config_;
    SelfPointFilterConfig map_export_self_point_filter_config_;

    bool point_noise_enabled_ = false;
    double range_noise_sigma_ = 0.02;
    double angular_noise_sigma_rad_ = 0.05 * M_PI / 180.0;
    double noise_reference_sigma_ = 0.05;
    double min_information_scale_ = 0.25;
    double max_information_scale_ = 4.0;

    std::vector<Keyframe::Ptr> all_keyframes_;  // 所有关键帧的存储列表
    Keyframe::Ptr last_kf_ = nullptr;           // 最近的关键帧指针（用于快速访问）
    int kf_id_ = 0;                             // 关键帧ID计数器（唯一标识每个关键帧）

    /// 当前帧点云缓存
    CloudPtr scan_undistort_full_{new PointCloudType()};  // 所有有效雷达的去畸变点云，用于发布
    CloudPtr scan_undistort_{new PointCloudType()};       // 动态选择后的定位点云
    CloudPtr scan_down_lidar_{new PointCloudType()};  // 当前帧降采样点云，Lidar坐标系
    CloudPtr scan_down_world_{new PointCloudType()};  // 当前帧降采样点云，世界坐标系
    pcl::VoxelGrid<PointType> voxel_scan_;            // 当前帧点云体素滤波器

    /// 点面相关
    std::vector<PointVector> nearest_points_;  // 当前帧每个点在IVox地图中的最近邻点
    std::vector<Vec4f> corr_pts_;              // 点面内点：[x,y,z,残差]，点坐标在Lidar系
    std::vector<std::uint8_t> corr_lidar_ids_;
    std::vector<Vec4f> corr_norm_;             // 点面内点对应平面：[nx,ny,nz,d]，平面在世界系
    std::vector<float> residuals_;             // 点到平面有符号残差
    std::vector<char> point_selected_surf_;    // 点面约束是否有效
    std::vector<Vec4f> plane_coef_;            // 每个点局部拟合出的平面参数

    /// 点到点相关
    std::vector<char> point_selected_icp_;  // 点到点ICP约束是否有效

    std::mutex mtx_buffer_;          // 保护Lidar/IMU缓存队列的互斥锁
    std::deque<double> time_buffer_; // 与lidar_buffer_一一对应的点云起始时间

    std::deque<PointCloudType::Ptr> lidar_buffer_;  // 激光雷达数据缓冲队列（用于与IMU时间同步）
    std::deque<lightning::IMUPtr> imu_buffer_;      // IMU数据缓冲队列（高频数据，用于状态预测）
    std::deque<MultiLidarFrameStats> lidar_stats_buffer_;
    std::deque<double> preprocess_time_buffer_ms_;
    double current_preprocess_ms_ = 0.0;

    /// options
    bool keep_first_imu_estimation_ = false;  // 在没有建立地图前，是否要使用前几帧的IMU状态
    double timediff_lidar_wrt_imu_ = 0.0;     // 激光雷达与IMU之间的时间偏移量，用于时间同步
    double last_timestamp_lidar_ = 0;         // 上一帧激光雷达数据的结束时间戳，用于检查数据间隙
    double lidar_end_time_ = 0;               // 当前帧激光雷达的结束时间戳（本地计算，用于IMU同步）
    double last_timestamp_imu_ = -1.0;        // 上一帧处理过的IMU数据时间戳，用于IMU数据筛选
    double first_lidar_time_ = 0.0;           // 系统启动后第一帧激光雷达的时间戳，作为参考基准
    bool lidar_pushed_ = false;               // 标记当前帧的激光雷达数据是否已添加到处理队列中

    bool enable_skip_lidar_ = true;  // 雷达是否需要跳帧
    int skip_lidar_num_ = 5;         // 每隔多少帧跳一个雷达
    int skip_lidar_cnt_ = 0;         // 跳帧计数器

    /// statistics and flags ///
    int scan_count_ = 0;                // 总扫描帧数统计
    int publish_count_ = 0;             // 发布次数统计（用于控制发布频率）
    bool flg_first_scan_ = true;        // 是否为第一帧扫描（用于初始化判断）
    bool flg_EKF_inited_ = false;       // ESKF滤波器是否已初始化（影响是否进行观测更新）
    double lidar_mean_scantime_ = 0.0;  // 激光雷达平均扫描时间（用于时间统计和性能监控）
    int scan_num_ = 0;                  // 当前扫描序列号
    int effect_feat_surf_ = 0, frame_num_ = 0, effect_feat_icp_ = 0;
    MultiLidarFrameStats current_lidar_stats_;
    AdaptiveLidarLoadController adaptive_lidar_load_controller_;
    AdaptiveLidarSelection current_lidar_selection_;
    std::atomic<double> latest_input_sensor_timestamp_{0.0};
    std::atomic<double> last_lidar_latency_sec_{0.0};
    std::atomic<double> last_frame_processing_ms_{0.0};
    std::atomic<std::size_t> adaptive_stale_drop_count_{0};
    bool last_tracking_healthy_ = false;
    double current_max_imu_gap_ = 0.0;
    std::size_t pre_imu_drop_count_ = 0;

    WheelSpeedDrConfig wheel_speed_dr_config_;
    mutable std::mutex wheel_speed_mutex_;
    std::deque<WheelSpeedSample> wheel_speed_buffer_;
    WheelSpeedDrStats wheel_speed_dr_stats_;
    double last_raw_wheel_speed_timestamp_ = std::numeric_limits<double>::lowest();
    double last_raw_wheel_speed_mps_ = 0.0;
    double last_lidar_filter_wheel_timestamp_ = std::numeric_limits<double>::lowest();
    double last_imu_filter_wheel_timestamp_ = std::numeric_limits<double>::lowest();

    double last_lidar_time_ = 0;  // 上一帧激光雷达时间戳（用于时间同步和断流检测）

    ///////////////////////// EKF inputs and output ///////////////////////////////////////////////////////
    MeasureGroup measures_;  // SyncPackages输出的一帧同步Lidar和IMU数据

    ESKF kf_;      // 点云时刻的IMU状态，用于畸变矫正+雷达里程计观测更新
    ESKF kf_imu_;  // imu 最新时刻的eskf状态，提供UI的高频位姿输出

    NavState state_point_;  // 当前Lidar帧结束时刻的前端状态

    bool use_aa_ = false;  // ESKF观测更新是否使用Anderson Acceleration
    bool propagate_velocity_ = false;  // 是否在ESKF名义状态中传播速度
    bool lidar_update_pose_only_ = false;  // Lidar观测是否只修正位姿
    bool lidar_update_inertial_states_ = true;  // Lidar观测是否间接修正bg/ba/gravity
    double max_update_velocity_step_ = 2.0;
    double max_update_gyro_bias_step_ = 0.05;
    double max_update_acc_bias_step_ = 0.5;
    double max_update_gravity_step_ = 0.05;
    bool adaptive_velocity_propagation_ = false;
    double velocity_innovation_ema_alpha_ = 0.05;
    double velocity_innovation_enable_threshold_ = 0.16;
    double velocity_innovation_disable_threshold_ = 0.10;
    int velocity_propagation_max_active_updates_ = 200;
    int velocity_propagation_cooldown_updates_ = 100;
    bool velocity_innovation_initialized_ = false;
    double velocity_innovation_ema_ = 0.0;
    bool velocity_propagation_active_ = false;
    int velocity_propagation_active_updates_ = 0;
    int velocity_propagation_cooldown_remaining_ = 0;
    bool velocity_propagation_safety_lockout_ = false;

    std::list<Keyframe::Ptr> proj_kfs_;  // 投影到当前帧的关键帧

    std::shared_ptr<ui::PangolinWindow> ui_ = nullptr;  // 可选UI，用于显示状态和点云
};

}  // namespace lightning

#endif  // FASTER_LIO_LASER_MAPPING_H
