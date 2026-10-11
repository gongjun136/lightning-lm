#pragma once

#include <cstddef>
#include <string>

#include <yaml-cpp/yaml.h>

#include "common/eigen_types.h"
#include "common/point_def.h"

namespace lightning::map_frame {

inline constexpr char kMetadataFilename[] = "map_frame.yaml";

struct ExportOptions {
    bool normalize_start_ground_z = false;  ///< 是否把起点地面归一化到导出地图参考平面。
    double start_ground_search_radius = 10.0;  ///< 起点附近地面估计的搜索半径（m）。
    double expected_sensor_height = 2.74;  ///< 用于预测地面高度的传感器安装高度（m）。
    double expected_height_tolerance = 0.8;  ///< 候选地面相对预期高度的允许偏差（m）。
    double histogram_bin_size = 0.05;  ///< 地面高度直方图的分箱宽度（m）。
    double inlier_tolerance = 0.20;  ///< 本次地面拟合的内点距离容差（m）。
    double max_median_absolute_deviation = 0.10;  ///< 地面估计允许的最大中位绝对偏差（m）。
    std::size_t min_ground_inliers = 100;  ///< 接受地面估计所需的最少内点数。
    double max_ground_tilt_deg = 10.0;  ///< 允许的最大地面倾斜角（deg）。

    bool export_pgm = false;  ///< 是否导出导航用 PGM 与 YAML。
    double pgm_resolution = 0.05;  ///< 导航栅格的像素分辨率（m）。
    double pgm_min_height = 2.5;  ///< 导出障碍点高度范围的下界（m）。
    double pgm_max_height = 3.0;  ///< 导出障碍点高度范围的上界（m）。
    double pgm_obstacle_dilation_radius = 0.0;  ///< 栅格障碍物膨胀半径（m）。
    int pgm_min_observation_frames = 3;  ///< 确认障碍物所需的最少独立观测帧数。
};

struct Metadata {
    bool normalized = false;  ///< 地图是否应用过起点地面归一化变换。
    std::string transform_id;  ///< 用于核对地图包与描述子坐标关系的变换标识。
    SE3 T_export_slam;  ///< 将原 SLAM 地图坐标变换到导出地图坐标的固定变换。
    double ground_z_slam = 0.0;  ///< 原 SLAM 坐标系中估计的地面高度（m）。
    double ground_z_export = 0.0;  ///< 导出坐标系中的地面参考高度（m）。
    Vec2d start_xy = Vec2d::Zero();  ///< 地面估计使用的原地图起点水平位置（m）。
    double expected_ground_z = 0.0;  ///< 由传感器高度推算的预期地面高度（m）。
    double search_radius = 0.0;  ///< 本次地面估计实际使用的搜索半径（m）。
    double histogram_bin_size = 0.0;  ///< 地面高度直方图的分箱宽度（m）。
    double inlier_tolerance = 0.0;  ///< 本次地面拟合的内点距离容差（m）。
    double median_absolute_deviation = 0.0;  ///< 地面高度估计的中位绝对偏差（m）。
    Vec3d ground_normal_slam = Vec3d::UnitZ();  ///< 原 SLAM 坐标系中的地面单位法向。
    double ground_plane_offset_slam = 0.0;  ///< 原 SLAM 坐标系地面平面的偏置项（m）。
    double ground_tilt_deg = 0.0;  ///< 本次估计的地面倾斜角（deg）。
    std::size_t candidate_count = 0;  ///< 参与地面估计的候选点数。
    std::size_t inlier_count = 0;  ///< 通过地面拟合的内点数。
};

bool ReadExportOptions(const YAML::Node& root, ExportOptions& options, std::string& error);

bool EstimateStartGroundFrame(const CloudPtr& map, const SE3& start_pose,
                              const ExportOptions& options, Metadata& metadata,
                              std::string& error);

void TransformCloudInPlace(const Metadata& metadata, const CloudPtr& cloud);

inline SE3 TransformPose(const Metadata& metadata, const SE3& pose) {
    return metadata.normalized ? metadata.T_export_slam * pose : pose;
}

bool SaveMetadata(const std::string& map_directory, const Metadata& metadata,
                  std::string& error);
bool LoadMetadata(const std::string& map_directory, Metadata& metadata,
                  std::string& error);

YAML::Node MakeTransformReference(const Metadata& metadata);
bool ReadTransformReference(const YAML::Node& node, Metadata& metadata,
                            std::string& error);
bool SameTransform(const Metadata& lhs, const Metadata& rhs, double tolerance = 1e-9);

}  // namespace lightning::map_frame
