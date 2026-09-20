#pragma once

#include <cstddef>
#include <string>

#include <yaml-cpp/yaml.h>

#include "common/eigen_types.h"
#include "common/point_def.h"

namespace lightning::map_frame {

inline constexpr char kMetadataFilename[] = "map_frame.yaml";

struct ExportOptions {
    bool normalize_start_ground_z = false;
    double start_ground_search_radius = 10.0;
    double expected_sensor_height = 2.74;
    double expected_height_tolerance = 0.8;
    double histogram_bin_size = 0.05;
    double inlier_tolerance = 0.20;
    double max_median_absolute_deviation = 0.10;
    std::size_t min_ground_inliers = 100;
    double max_ground_tilt_deg = 10.0;

    bool export_pgm = false;
    double pgm_resolution = 0.05;
    double pgm_min_height = 2.5;
    double pgm_max_height = 3.0;
    double pgm_obstacle_dilation_radius = 0.0;
    int pgm_min_observation_frames = 3;
};

struct Metadata {
    bool normalized = false;
    std::string transform_id;
    SE3 T_export_slam;
    double ground_z_slam = 0.0;
    double ground_z_export = 0.0;
    Vec2d start_xy = Vec2d::Zero();
    double expected_ground_z = 0.0;
    double search_radius = 0.0;
    double histogram_bin_size = 0.0;
    double inlier_tolerance = 0.0;
    double median_absolute_deviation = 0.0;
    Vec3d ground_normal_slam = Vec3d::UnitZ();
    double ground_plane_offset_slam = 0.0;
    double ground_tilt_deg = 0.0;
    std::size_t candidate_count = 0;
    std::size_t inlier_count = 0;
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
