#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <vector>

#include <opencv2/imgcodecs.hpp>
#include <yaml-cpp/yaml.h>

#include "core/maps/navigation_map_export.h"

namespace {

void AddPoint(const lightning::CloudPtr& cloud, float x, float y, float z) {
    lightning::PointType point;
    point.x = x;
    point.y = y;
    point.z = z;
    cloud->push_back(point);
}

}  // namespace

int main() {
    const auto cloud = std::make_shared<lightning::PointCloudType>();
    AddPoint(cloud, -0.02F, -0.02F, 0.0F);
    AddPoint(cloud, 0.11F, 0.11F, 0.0F);
    AddPoint(cloud, 0.01F, 0.06F, 2.75F);

    lightning::map_frame::ExportOptions options;
    options.export_pgm = true;
    options.pgm_resolution = 0.05;
    options.pgm_min_height = 2.5;
    options.pgm_max_height = 3.0;
    options.pgm_obstacle_dilation_radius = 0.0;
    options.pgm_min_observation_frames = 3;

    const auto unique =
        std::chrono::steady_clock::now().time_since_epoch().count();
    const std::filesystem::path directory =
        std::filesystem::temp_directory_path() /
        ("lightning_navigation_map_test_" + std::to_string(unique));
    std::filesystem::create_directories(directory);

    const auto common_ray_cloud = std::make_shared<lightning::PointCloudType>();
    AddPoint(common_ray_cloud, 0.01F, 0.06F, 2.75F);
    const auto extra_ray_cloud =
        std::make_shared<lightning::PointCloudType>(*common_ray_cloud);
    AddPoint(extra_ray_cloud, 0.11F, -0.02F, 0.0F);
    const lightning::SE3 identity(lightning::Quatd::Identity(),
                                  lightning::Vec3d::Zero());
    const std::vector<lightning::navigation_map::RaycastFrame> raycast_frames{
        {extra_ray_cloud, identity},
        {common_ray_cloud, identity},
        {common_ray_cloud, identity}};
    const lightning::navigation_map::SensorOrigins sensor_origins{
        {0, lightning::Vec3d::Zero()}};
    lightning::navigation_map::ExportResult result;
    std::string error;
    if (!lightning::navigation_map::ExportPgmAndYaml(
            cloud, raycast_frames, sensor_origins, directory.string(), options,
            result, error)) {
        std::cerr << error << std::endl;
        return 1;
    }
    if (result.width != 4 || result.height != 4 ||
        std::abs(result.origin_x + 0.05) > 1e-9 ||
        std::abs(result.origin_y + 0.05) > 1e-9 ||
        result.occupied_cells != 1 || result.free_cells == 0 ||
        result.unknown_cells == 0 || result.ray_count != 4) {
        std::cerr << "unexpected PGM geometry or occupied-cell count" << std::endl;
        return 2;
    }

    const cv::Mat image = cv::imread(
        (directory / "map.pgm").string(), cv::IMREAD_GRAYSCALE);
    if (image.empty() || image.rows != 4 || image.cols != 4 ||
        image.at<unsigned char>(1, 1) != 0 ||
        image.at<unsigned char>(2, 1) != 255 ||
        image.at<unsigned char>(3, 3) != 128 ||
        cv::countNonZero(image == 0) != 1 ||
        cv::countNonZero(image == 128) == 0 ||
        cv::countNonZero(image == 255) == 0) {
        std::cerr << "PGM pixels do not follow ROS lower-left origin convention"
                  << std::endl;
        return 3;
    }

    const YAML::Node manifest = YAML::LoadFile((directory / "map.yaml").string());
    const auto origin = manifest["origin"].as<std::vector<double>>();
    if (manifest["image"].as<std::string>() != "map.pgm" ||
        std::abs(manifest["resolution"].as<double>() - 0.05) > 1e-12 ||
        origin.size() != 3 || std::abs(origin[0] + 0.05) > 1e-9 ||
        std::abs(origin[1] + 0.05) > 1e-9 ||
        manifest["negate"].as<int>() != 0) {
        std::cerr << "map.yaml is not ROS navigation compatible" << std::endl;
        return 4;
    }

    std::error_code cleanup_error;
    std::filesystem::remove_all(directory, cleanup_error);
    std::cout << "navigation_map_export_test passed" << std::endl;
    return 0;
}
