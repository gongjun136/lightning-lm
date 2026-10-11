#ifndef FASTER_LIO_POINTCLOUD_PROCESSING_H
#define FASTER_LIO_POINTCLOUD_PROCESSING_H

#include <pcl_conversions/pcl_conversions.h>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include "common/measure_group.h"
#include "common/point_def.h"
#include "livox_ros_driver2/msg/custom_msg.hpp"

namespace lightning {

enum class LidarType { AVIA = 1, VELO32 = 2, OUST64 = 3, ROBOSENSE = 4 };

/**
 * point cloud preprocess
 * just unify the point format from livox/velodyne to PCL
 *
 * 预处理程序
 * 主要是对各种不同的雷达处理时间戳差异
 */
class PointCloudPreprocess {
   public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    PointCloudPreprocess() = default;
    ~PointCloudPreprocess() = default;

    /// processors
    void Process(const sensor_msgs::msg::PointCloud2 ::SharedPtr &msg, PointCloudType::Ptr &pcl_out);

    void Process(const livox_ros_driver2::msg::CustomMsg::SharedPtr &cloud, PointCloudType::Ptr &pcl_out);

    void Set(LidarType lid_type, double bld, int pfilt_num);

    // accessors
    double &Blind() { return blind_; }
    int &NumScans() { return num_scans_; }
    int &PointFilterNum() { return point_filter_num_; }
    float &TimeScale() { return time_scale_; }
    double &LivoxPointTimeScale() { return livox_point_time_scale_; }
    int &ParallelThreads() { return parallel_threads_; }
    LidarType GetLidarType() const { return lidar_type_; }
    void SetLidarType(LidarType lt) { lidar_type_ = lt; }

    void SetHeightROI(float height_max, float height_min) {
        height_max_ = height_max;
        height_min_ = height_min;
    }

   private:
    void LivoxPointCloud2Handler(const sensor_msgs::msg::PointCloud2 ::SharedPtr &msg);
    void Oust64Handler(const sensor_msgs::msg::PointCloud2 ::SharedPtr &msg);
    void RoboSenseHandler(const sensor_msgs::msg::PointCloud2 ::SharedPtr &msg);
    void VelodyneHandler(const sensor_msgs::msg::PointCloud2 ::SharedPtr &msg);

    PointCloudType cloud_full_ /**< 保留的完整预处理点云缓存。 */,
        cloud_out_ /**< 经过有效性、距离和点采样过滤的输出点云缓存。 */;

    LidarType lidar_type_ = LidarType::AVIA;  ///< 选择点格式与时间戳解释方式的雷达类型。
    int point_filter_num_ = 1;  ///< 点序列的采样步长。
    int num_scans_ = 6;  ///< 雷达扫描线数量。
    double blind_ = 0.01;  ///< 丢弃近距离盲区点的距离阈值（m）。
    float time_scale_ = 1e-3;  ///< 标准点云时间字段到内部时间单位的换算系数。
    double livox_point_time_scale_ = 1.0;  ///< Livox 点时间偏移到内部时间单位的换算系数。
    int parallel_threads_ = 1;  ///< 预处理允许使用的并行线程数。
    bool given_offset_time_ = false;  ///< 输入是否提供逐点时间偏移。

    float height_max_ = 1.0;  ///< 预处理高度 ROI 的上界（m）。
    float height_min_ = -1.0;  ///< 预处理高度 ROI 的下界（m）。
};
}  // namespace lightning

#endif
