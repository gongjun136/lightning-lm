#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>

#include "core/maps/map_frame.h"

int main() {
    auto cloud = std::make_shared<lightning::PointCloudType>();
    constexpr double slope_x = 0.05;
    constexpr double slope_y = -0.035;
    for (int x = -50; x <= 50; ++x) {
        for (int y = -50; y <= 50; ++y) {
            lightning::PointType point;
            point.x = 0.1F * x;
            point.y = 0.1F * y;
            point.z = static_cast<float>(
                1.25 + slope_x * point.x + slope_y * point.y +
                0.002 * ((x + y) % 5));
            cloud->push_back(point);
        }
    }
    for (int index = 0; index < 500; ++index) {
        lightning::PointType outlier;
        outlier.x = 0.01F * index;
        outlier.y = 0.0F;
        outlier.z = static_cast<float>(-4.0 + 0.02 * index);
        cloud->push_back(outlier);
    }

    lightning::map_frame::ExportOptions options;
    options.normalize_start_ground_z = true;
    options.expected_sensor_height = 2.74;
    lightning::map_frame::Metadata metadata;
    std::string error;
    const lightning::SE3 start_pose(
        lightning::Quatd::Identity(), lightning::Vec3d(0.0, 0.0, 3.99));
    if (!lightning::map_frame::EstimateStartGroundFrame(
            cloud, start_pose, options, metadata, error)) {
        std::cerr << error << std::endl;
        return 1;
    }
    if (std::abs(metadata.ground_z_slam - 1.25) > 0.03 ||
        metadata.ground_tilt_deg < 2.0 || metadata.ground_tilt_deg > 5.0) {
        std::cerr << "unexpected ground estimate: " << metadata.ground_z_slam
                  << ", tilt=" << metadata.ground_tilt_deg
                  << std::endl;
        return 2;
    }
    const lightning::Vec3d expected_normal =
        lightning::Vec3d(-slope_x, -slope_y, 1.0).normalized();
    if ((metadata.ground_normal_slam - expected_normal).norm() > 0.02 ||
        (metadata.T_export_slam.unit_quaternion() *
             metadata.ground_normal_slam -
         lightning::Vec3d::UnitZ())
                .norm() > 1e-6) {
        std::cerr << "ground normal was not aligned to +Z" << std::endl;
        return 2;
    }

    lightning::map_frame::TransformCloudInPlace(metadata, cloud);
    std::vector<double> ground;
    for (const auto& point : cloud->points) {
        if (std::abs(point.z) < 0.1) ground.push_back(point.z);
    }
    if (ground.size() < 10000) {
        std::cerr << "ground cloud was not translated near Z=0" << std::endl;
        return 3;
    }

    const auto unique =
        std::chrono::steady_clock::now().time_since_epoch().count();
    const std::filesystem::path directory =
        std::filesystem::temp_directory_path() /
        ("lightning_map_frame_test_" + std::to_string(unique));
    std::filesystem::create_directories(directory);
    if (!lightning::map_frame::SaveMetadata(directory.string(), metadata, error)) {
        std::cerr << error << std::endl;
        return 4;
    }
    lightning::map_frame::Metadata loaded;
    if (!lightning::map_frame::LoadMetadata(directory.string(), loaded, error) ||
        !lightning::map_frame::SameTransform(metadata, loaded)) {
        std::cerr << "metadata round trip failed: " << error << std::endl;
        return 5;
    }
    std::error_code cleanup_error;
    std::filesystem::remove_all(directory, cleanup_error);
    std::cout << "map_frame_test passed: ground_z=" << metadata.ground_z_slam
              << ", inliers=" << metadata.inlier_count << std::endl;
    return 0;
}
