#include <chrono>
#include <filesystem>
#include <iostream>

#include "common/pcd_io.h"
#include "common/point_def.h"

int main(int argc, char** argv) {
    const auto unique = std::chrono::steady_clock::now().time_since_epoch().count();
    const std::filesystem::path directory = argc > 1
        ? std::filesystem::path(argv[1])
        : std::filesystem::temp_directory_path() /
              ("lightning_world_pcd_test_" + std::to_string(unique));
    std::filesystem::create_directories(directory);

    lightning::PointCloudType source;
    for (int index = 0; index < 4; ++index) {
        lightning::PointType point;
        point.x = 10.0F + index;
        point.y = 20.0F - index;
        point.z = 3.5F + 0.25F * index;
        point.intensity = 12.5F + index;
        point.time = 1791517198.125 + index * 0.125;
        point.lidar_id = static_cast<std::uint8_t>(index % 3);
        source.push_back(point);
    }
    source.width = 2;
    source.height = 2;
    source.sensor_origin_ = Eigen::Vector4f(2.3F, -3.1F, 3.5F, 0.0F);
    source.sensor_orientation_ = Eigen::Quaternionf(
        Eigen::AngleAxisf(0.05F, Eigen::Vector3f::UnitZ()));
    const auto source_origin = source.sensor_origin_;
    const auto source_orientation = source.sensor_orientation_;

    const auto path = directory / "world.pcd";
    if (lightning::pcd_io::SaveWorldCloudBinaryCompressed(path.string(), source) != 0) {
        std::cerr << "world PCD export failed" << std::endl;
        return 1;
    }
    lightning::PointCloudType loaded;
    if (pcl::io::loadPCDFile(path.string(), loaded) != 0 ||
        loaded.size() != source.size() || loaded.width != source.width ||
        loaded.height != source.height || !loaded.sensor_origin_.isZero() ||
        !loaded.sensor_orientation_.isApprox(Eigen::Quaternionf::Identity())) {
        std::cerr << "world PCD must retain dimensions and use an identity VIEWPOINT" << std::endl;
        return 2;
    }
    for (std::size_t index = 0; index < source.size(); ++index) {
        const auto& expected = source[index];
        const auto& actual = loaded[index];
        if (actual.x != expected.x || actual.y != expected.y || actual.z != expected.z ||
            actual.intensity != expected.intensity || actual.time != expected.time ||
            actual.lidar_id != expected.lidar_id) {
            std::cerr << "world PCD changed a point coordinate or attribute" << std::endl;
            return 3;
        }
    }
    if (!source.sensor_origin_.isApprox(source_origin) ||
        !source.sensor_orientation_.isApprox(source_orientation)) {
        std::cerr << "export changed in-memory sensor metadata" << std::endl;
        return 4;
    }

    // Alignment overlays use RGB clouds and may inherit a loaded PCD viewpoint.
    pcl::PointCloud<pcl::PointXYZRGB> rgb;
    pcl::PointXYZRGB point;
    point.x = 1.0F; point.y = -2.0F; point.z = 3.0F;
    point.r = 123; point.g = 45; point.b = 67;
    rgb.push_back(point);
    rgb.sensor_origin_ = source_origin;
    rgb.sensor_orientation_ = source_orientation;
    const auto rgb_path = directory / "world_rgb.pcd";
    pcl::PointCloud<pcl::PointXYZRGB> loaded_rgb;
    if (lightning::pcd_io::SaveWorldCloudBinaryCompressed(rgb_path.string(), rgb) != 0 ||
        pcl::io::loadPCDFile(rgb_path.string(), loaded_rgb) != 0 || loaded_rgb.size() != 1 ||
        !loaded_rgb.sensor_origin_.isZero() ||
        !loaded_rgb.sensor_orientation_.isApprox(Eigen::Quaternionf::Identity()) ||
        loaded_rgb[0].x != point.x || loaded_rgb[0].y != point.y || loaded_rgb[0].z != point.z ||
        loaded_rgb[0].r != point.r || loaded_rgb[0].g != point.g || loaded_rgb[0].b != point.b) {
        std::cerr << "RGB world PCD round trip failed" << std::endl;
        return 5;
    }
    if (argc == 1) {
        std::error_code cleanup_error;
        std::filesystem::remove_all(directory, cleanup_error);
    }
    std::cout << "world_pcd_export_test passed: identity VIEWPOINT; XYZ, fields and sensor metadata preserved"
              << std::endl;
    return 0;
}
