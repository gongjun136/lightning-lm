#ifndef FASTER_LIO_LASER_MAPPING_H
#define FASTER_LIO_LASER_MAPPING_H

#include <pcl/filters/voxel_grid.h>
#include <condition_variable>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <thread>

#include "common/eigen_types.h"
#include "common/imu.h"
#include "common/keyframe.h"
#include "common/options.h"
#include "core/ivox3d/ivox3d.h"
#include "core/lio/eskf.hpp"
#include "core/lio/imu_processing.hpp"
#include "pointcloud_preprocess.h"

#include "livox_ros_driver2/msg/custom_msg.hpp"

namespace lightning {

namespace ui {
class PangolinWindow;
}

/**
 * laser mapping
 * 目前有个问题：点云在缓存之后，实际处理的并不是最新的那个点云（通常是buffer里的前一个），这是因为bag里的点云用的开始时间戳，导致
 * 点云的结束时间要比IMU多0.1s左右。为了同步最近的IMU，就只能处理缓冲队列里的那个点云，而不是最新的点云
 */
class LaserMapping {
   public:
    struct Options {
        Options() {}

        bool is_in_slam_mode_ = true;  // 是否在slam模式下

        /// 关键帧阈值
        double kf_dis_th_ = 2.0;
        double kf_angle_th_ = 15 * M_PI / 180.0;
    };

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using IVoxType = IVox<3, IVoxNodeType::DEFAULT, PointType>;

    LaserMapping(Options options = Options());
    ~LaserMapping() {
        scan_down_lidar_ = nullptr;
        scan_undistort_ = nullptr;
        scan_down_world_ = nullptr;
        LOG(INFO) << "laser mapping deconstruct";
    }

    /// init without ros
    bool Init(const std::string &config_yaml);

    bool Run();

    // 三个ProcessPointCloud2函数处理逻辑：
    // 1.时间戳检查和回环检测：检查时间戳是否倒退，如果是则清空缓冲区
    // 2.点云预处理：调用 preprocess_->Process() 进行统一格式 (CloudPtr)转换和预处理
    // 3.数据缓存：将处理后的点云和时间戳存入缓冲区
    // 4.性能监控：使用 Timer::Evaluate 记录预处理耗时
    // callbacks of lidar and imu
    /// 处理ROS2的点云
    void ProcessPointCloud2(const sensor_msgs::msg::PointCloud2::SharedPtr &msg);

    /// 处理livox的点云
    void ProcessPointCloud2(const livox_ros_driver2::msg::CustomMsg::SharedPtr &msg);

    /// 如果已经做了预处理，也可以直接处理点云
    void ProcessPointCloud2(CloudPtr cloud);

    void ProcessIMU(const lightning::IMUPtr &msg_in);

    /// 保存前端的地图
    void SaveMap();

    void SetUI(std::shared_ptr<ui::PangolinWindow> ui) { ui_ = ui; }

    /// 获取关键帧
    Keyframe::Ptr GetKeyframe() const { return last_kf_; }

    /// 获取激光的状态
    NavState GetState() const { return state_point_; }

    /// 获取IMU状态
    NavState GetIMUState() const {
        if (p_imu_->IsIMUInited()) {
            return kf_imu_.GetX();
        } else {
            NavState s;
            s.pose_is_ok_ = false;
            return s;
        }
    }

    CloudPtr GetScanUndist() const { return scan_undistort_; }

    /// 获取最新的点云
    CloudPtr GetRecentCloud();

    std::vector<Keyframe::Ptr> GetAllKeyframes() { return all_keyframes_; }

    /**
     * 计算全局地图
     * @param use_lio_pose
     * @return
     */
    CloudPtr GetGlobalMap(bool use_lio_pose, bool use_voxel = true, float res = 0.1);

   private:
    // sync lidar with imu
    bool SyncPackages();

    void ObsModel(NavState &s, ESKF::CustomObservationModel &obs);

    inline void PointLidarToWorld(const PointType &pi, PointType &po) {
        Vec3d p_global(state_point_.rot_ * (state_point_.offset_R_lidar_ * pi.getVector3fMap().cast<double>() +
                                            state_point_.offset_t_lidar_) +
                       state_point_.pos_);

        po.x = p_global(0);
        po.y = p_global(1);
        po.z = p_global(2);
        po.intensity = pi.intensity;
    }

    void MapIncremental();

    bool LoadParamsFromYAML(const std::string &yaml);

    /// 创建关键帧
    void MakeKF();

   private:
    Options options_;

    /// modules
    IVoxType::Options ivox_options_;
    std::shared_ptr<IVoxType> ivox_ = nullptr;                    // localmap in ivox
    std::shared_ptr<PointCloudPreprocess> preprocess_ = nullptr;  // point cloud preprocess
    std::shared_ptr<ImuProcess> p_imu_ = nullptr;                 // imu process

    /// local map related
    double filter_size_map_min_ = 0;  // 地图体素滤波分辨率（m），控制地图点的密度

    /// params
    std::vector<double> extrinT_{3, 0.0};  // lidar-imu translation
    std::vector<double> extrinR_{9, 0.0};  // lidar-imu rotation
    std::string map_file_path_;

    std::vector<Keyframe::Ptr> all_keyframes_;  // 所有关键帧的存储列表
    Keyframe::Ptr last_kf_ = nullptr;           // 最近的关键帧指针（用于快速访问）
    int kf_id_ = 0;                             // 关键帧ID计数器（唯一标识每个关键帧）

    /// point clouds data
    CloudPtr scan_undistort_{new PointCloudType()};   // scan after undistortion in lidar
    CloudPtr scan_down_lidar_{new PointCloudType()};  // downsampled scan in lidar
    CloudPtr scan_down_world_{new PointCloudType()};  // downsampled scan in world
    std::vector<PointVector> nearest_points_;         // nearest points of current scan in world
    std::vector<Vec4f> corr_pts_;                     // 内点：有效匹配点 [x,y,z,残差]，lidar系
    std::vector<Vec4f> corr_norm_;                    // 内点：对应平面法向量 [nx,ny,nz,d]，world系
    pcl::VoxelGrid<PointType> voxel_scan_;            // voxel filter for current scan

    std::vector<float> residuals_;  // point-to-plane residuals
    // [gj-2025-11-28] bool -> uint8_t
    std::vector<uint8_t> point_selected_surf_;  // selected points (uint8_t for thread safety)
    std::vector<Vec4f> plane_coef_;             // plane coeffs

    std::mutex mtx_buffer_;
    std::deque<double> time_buffer_;

    std::deque<PointCloudType::Ptr> lidar_buffer_;  // 激光雷达数据缓冲队列（用于与IMU时间同步）
    std::deque<lightning::IMUPtr> imu_buffer_;      // IMU数据缓冲队列（高频数据，用于状态预测）

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
    int skip_lidar_cnt_ = 0;

    /// statistics and flags ///
    int scan_count_ = 0;                // 总扫描帧数统计
    int publish_count_ = 0;             // 发布次数统计（用于控制发布频率）
    bool flg_first_scan_ = true;        // 是否为第一帧扫描（用于初始化判断）
    bool flg_EKF_inited_ = false;       // ESKF滤波器是否已初始化（影响是否进行观测更新）
    double lidar_mean_scantime_ = 0.0;  // 激光雷达平均扫描时间（用于时间统计和性能监控）
    int scan_num_ = 0;                  // 当前扫描序列号
    int effect_feat_num_ = 0;           // 有效特征点数量（成功匹配的点云特征数）
    int frame_num_ = 0;                 // 总处理帧数

    double last_lidar_time_ = 0;  // 上一帧激光雷达时间戳（用于时间同步和断流检测）

    ///////////////////////// EKF inputs and output ///////////////////////////////////////////////////////
    MeasureGroup measures_;  // sync IMU and lidar scan

    ESKF kf_;      // 点云时刻的IMU状态，用于畸变矫正+雷达里程计观测更新
    ESKF kf_imu_;  // imu 最新时刻的eskf状态，提供UI的高频位姿输出

    NavState state_point_;  // ekf current state

    Vec3d pos_lidar_;               // lidar position after eskf update
    SO3 euler_cur_;                 // rotation in euler angles
    bool extrinsic_est_en_ = true;  // 是否估计lidar和IMU外参
    bool use_aa_ = false;           // use anderson acceleration?

    std::shared_ptr<ui::PangolinWindow> ui_ = nullptr;
};

}  // namespace lightning

#endif  // FASTER_LIO_LASER_MAPPING_H