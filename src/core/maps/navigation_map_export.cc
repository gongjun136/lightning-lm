#include "core/maps/navigation_map_export.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <vector>

#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <yaml-cpp/yaml.h>

namespace lightning::navigation_map {

namespace {

constexpr int kRayBinCount = 3600;
constexpr int kMaxRayGapBins = 5;

struct RayEndpoint {
    double range_squared = -1.0;
    cv::Point cell;
};

struct SensorRays {
    Vec3d origin_map = Vec3d::Zero();
    std::vector<RayEndpoint> bins = std::vector<RayEndpoint>(kRayBinCount);
};

bool ToImageCell(const Vec3d& point, double origin_x, double origin_y,
                 double resolution, int height, cv::Point& cell) {
    const double grid_x = std::floor((point.x() - origin_x) / resolution);
    const double grid_y = std::floor((point.y() - origin_y) / resolution);
    if (!std::isfinite(grid_x) || !std::isfinite(grid_y) ||
        grid_x < std::numeric_limits<int>::lowest() ||
        grid_x > std::numeric_limits<int>::max() ||
        grid_y < std::numeric_limits<int>::lowest() ||
        grid_y > std::numeric_limits<int>::max()) {
        return false;
    }
    const auto grid_x_integer = static_cast<long long>(grid_x);
    const auto grid_y_integer = static_cast<long long>(grid_y);
    const long long image_y = static_cast<long long>(height) - 1 - grid_y_integer;
    if (image_y < std::numeric_limits<int>::lowest() ||
        image_y > std::numeric_limits<int>::max()) {
        return false;
    }
    cell.x = static_cast<int>(grid_x_integer);
    cell.y = static_cast<int>(image_y);
    return true;
}

void FillRaySector(const cv::Point& origin,
                   const std::vector<cv::Point>& endpoints,
                   cv::Mat& known_mask) {
    if (endpoints.empty()) return;
    if (endpoints.size() == 1) {
        cv::line(known_mask, origin, endpoints.front(), cv::Scalar(255));
        return;
    }
    std::vector<cv::Point> polygon;
    polygon.reserve(endpoints.size() + 1);
    polygon.push_back(origin);
    polygon.insert(polygon.end(), endpoints.begin(), endpoints.end());
    const cv::Point* polygon_data = polygon.data();
    const int polygon_size = static_cast<int>(polygon.size());
    cv::fillPoly(known_mask, &polygon_data, &polygon_size, 1, cv::Scalar(255));
}

std::size_t FillSensorRays(const cv::Point& origin,
                           const std::vector<RayEndpoint>& bins,
                           cv::Mat& known_mask) {
    std::vector<int> occupied_bins;
    occupied_bins.reserve(bins.size());
    for (int index = 0; index < static_cast<int>(bins.size()); ++index) {
        if (bins[index].range_squared >= 0.0) occupied_bins.push_back(index);
    }
    if (occupied_bins.empty()) return 0;
    if (occupied_bins.size() == 1) {
        FillRaySector(origin, {bins[occupied_bins.front()].cell}, known_mask);
        return 1;
    }

    std::size_t largest_gap_after = 0;
    int largest_gap = -1;
    for (std::size_t index = 0; index < occupied_bins.size(); ++index) {
        const int current = occupied_bins[index];
        const int next = occupied_bins[(index + 1) % occupied_bins.size()];
        const int gap = (next - current + kRayBinCount) % kRayBinCount;
        if (gap > largest_gap) {
            largest_gap = gap;
            largest_gap_after = index;
        }
    }

    const std::size_t start = (largest_gap_after + 1) % occupied_bins.size();
    std::vector<cv::Point> sector;
    int previous_bin = -1;
    for (std::size_t offset = 0; offset < occupied_bins.size(); ++offset) {
        const int bin = occupied_bins[(start + offset) % occupied_bins.size()];
        if (previous_bin >= 0) {
            const int gap = (bin - previous_bin + kRayBinCount) % kRayBinCount;
            if (gap > kMaxRayGapBins) {
                FillRaySector(origin, sector, known_mask);
                sector.clear();
            }
        }
        sector.push_back(bins[bin].cell);
        previous_bin = bin;
    }
    FillRaySector(origin, sector, known_mask);
    return occupied_bins.size();
}

std::size_t ClearObservedRays(
    const std::vector<RaycastFrame>& frames,
    const SensorOrigins& sensor_origins_primary, double origin_x,
    double origin_y, double resolution, int height,
    int min_observation_frames, cv::Mat& known_mask) {
    const double pi = std::acos(-1.0);
    std::size_t ray_count = 0;
    cv::Mat observation_count(known_mask.size(), CV_8UC1, cv::Scalar(0));
    cv::Mat frame_mask(known_mask.size(), CV_8UC1, cv::Scalar(0));
    for (const auto& frame : frames) {
        if (!frame.cloud || frame.cloud->empty()) continue;
        frame_mask.setTo(0);
        std::map<std::uint8_t, SensorRays> rays_by_sensor;
        for (const auto& point : frame.cloud->points) {
            if (!std::isfinite(point.x) || !std::isfinite(point.y) ||
                !std::isfinite(point.z)) {
                continue;
            }
            auto rays = rays_by_sensor.find(point.lidar_id);
            if (rays == rays_by_sensor.end()) {
                const auto origin = sensor_origins_primary.find(point.lidar_id);
                const Vec3d origin_primary =
                    origin == sensor_origins_primary.end() ? Vec3d::Zero()
                                                           : origin->second;
                SensorRays sensor_rays;
                sensor_rays.origin_map = frame.T_map_primary * origin_primary;
                rays = rays_by_sensor
                           .emplace(point.lidar_id, std::move(sensor_rays))
                           .first;
            }

            const Vec3d endpoint_map =
                frame.T_map_primary *
                Vec3d(point.x, point.y, point.z);
            const Vec2d delta =
                endpoint_map.head<2>() - rays->second.origin_map.head<2>();
            const double range_squared = delta.squaredNorm();
            if (!endpoint_map.allFinite() || !std::isfinite(range_squared) ||
                range_squared < resolution * resolution) {
                continue;
            }
            cv::Point endpoint_cell;
            if (!ToImageCell(endpoint_map, origin_x, origin_y, resolution,
                             height, endpoint_cell)) {
                continue;
            }
            const double angle = std::atan2(delta.y(), delta.x());
            const int bin = std::clamp(
                static_cast<int>(std::floor((angle + pi) * kRayBinCount /
                                            (2.0 * pi))),
                0, kRayBinCount - 1);
            if (range_squared > rays->second.bins[bin].range_squared) {
                rays->second.bins[bin].range_squared = range_squared;
                rays->second.bins[bin].cell = endpoint_cell;
            }
        }

        std::size_t frame_ray_count = 0;
        for (const auto& [lidar_id, sensor_rays] : rays_by_sensor) {
            (void)lidar_id;
            cv::Point origin_cell;
            if (!ToImageCell(sensor_rays.origin_map, origin_x, origin_y,
                             resolution, height, origin_cell)) {
                continue;
            }
            frame_ray_count +=
                FillSensorRays(origin_cell, sensor_rays.bins, frame_mask);
        }
        if (frame_ray_count > 0) {
            cv::add(observation_count, cv::Scalar(1), observation_count,
                    frame_mask);
            ray_count += frame_ray_count;
        }
    }
    cv::compare(observation_count, min_observation_frames - 1, known_mask,
                cv::CMP_GT);
    return ray_count;
}

}  // namespace

bool ExportPgmAndYaml(const CloudPtr& map,
                      const std::vector<RaycastFrame>& raycast_frames,
                      const SensorOrigins& sensor_origins_primary,
                      const std::string& map_directory,
                      const map_frame::ExportOptions& options,
                      ExportResult& result, std::string& error) {
    result = ExportResult{};
    if (!map || map->empty()) {
        error = "cannot export a PGM from an empty global map";
        return false;
    }

    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    for (const auto& point : map->points) {
        if (!std::isfinite(point.x) || !std::isfinite(point.y) ||
            !std::isfinite(point.z)) {
            continue;
        }
        min_x = std::min(min_x, static_cast<double>(point.x));
        min_y = std::min(min_y, static_cast<double>(point.y));
        max_x = std::max(max_x, static_cast<double>(point.x));
        max_y = std::max(max_y, static_cast<double>(point.y));
    }
    if (!std::isfinite(min_x) || !std::isfinite(min_y) ||
        !std::isfinite(max_x) || !std::isfinite(max_y)) {
        error = "the global map has no finite points";
        return false;
    }

    const double resolution = options.pgm_resolution;
    const double origin_x = std::floor(min_x / resolution) * resolution;
    const double origin_y = std::floor(min_y / resolution) * resolution;
    const double width_cells = std::floor((max_x - origin_x) / resolution) + 1.0;
    const double height_cells = std::floor((max_y - origin_y) / resolution) + 1.0;
    if (!std::isfinite(width_cells) || !std::isfinite(height_cells) ||
        width_cells <= 0.0 || height_cells <= 0.0 ||
        width_cells > std::numeric_limits<int>::max() ||
        height_cells > std::numeric_limits<int>::max() ||
        width_cells * height_cells > 1.0e9) {
        error = "the requested PGM dimensions are invalid or too large";
        return false;
    }

    const int width = static_cast<int>(width_cells);
    const int height = static_cast<int>(height_cells);
    cv::Mat obstacle_mask(height, width, CV_8UC1, cv::Scalar(0));
    for (const auto& point : map->points) {
        if (!std::isfinite(point.x) || !std::isfinite(point.y) ||
            !std::isfinite(point.z)) {
            continue;
        }
        const int x = static_cast<int>(
            std::floor((static_cast<double>(point.x) - origin_x) / resolution));
        const int y = static_cast<int>(
            std::floor((static_cast<double>(point.y) - origin_y) / resolution));
        if (x >= 0 && x < width && y >= 0 && y < height) {
            const int image_y = height - 1 - y;
            if (point.z >= options.pgm_min_height &&
                point.z <= options.pgm_max_height) {
                obstacle_mask.at<unsigned char>(image_y, x) = 255;
            }
        }
    }

    cv::Mat known_mask(height, width, CV_8UC1, cv::Scalar(0));
    const std::size_t ray_count = ClearObservedRays(
        raycast_frames, sensor_origins_primary, origin_x, origin_y, resolution,
        height, options.pgm_min_observation_frames, known_mask);
    if (ray_count == 0) {
        error = "trajectory ray clearing produced no valid rays";
        return false;
    }

    if (options.pgm_obstacle_dilation_radius > 0.0) {
        const int radius_pixels = static_cast<int>(
            std::ceil(options.pgm_obstacle_dilation_radius / resolution));
        const int diameter = 2 * radius_pixels + 1;
        const cv::Mat kernel = cv::getStructuringElement(
            cv::MORPH_ELLIPSE, cv::Size(diameter, diameter));
        cv::dilate(obstacle_mask, obstacle_mask, kernel);
    }

    constexpr unsigned char kUnknown = 128;
    cv::Mat image(height, width, CV_8UC1, cv::Scalar(kUnknown));
    image.setTo(255, known_mask);
    image.setTo(0, obstacle_mask);
    const std::filesystem::path directory(map_directory);
    const std::filesystem::path pgm_path = directory / "map.pgm";
    const std::filesystem::path yaml_path = directory / "map.yaml";
    if (!cv::imwrite(pgm_path.string(), image)) {
        error = "failed to write navigation PGM: " + pgm_path.string();
        return false;
    }

    try {
        YAML::Node manifest;
        manifest["image"] = "map.pgm";
        manifest["mode"] = "trinary";
        manifest["resolution"] = resolution;
        manifest["origin"] = std::vector<double>{origin_x, origin_y, 0.0};
        manifest["negate"] = 0;
        manifest["occupied_thresh"] = 0.65;
        manifest["free_thresh"] = 0.25;
        manifest["width"] = width;
        manifest["height"] = height;
        std::ofstream output(yaml_path);
        if (!output.is_open()) {
            error = "failed to open navigation map manifest: " + yaml_path.string();
            return false;
        }
        output << manifest;
        output.close();
        if (!output) {
            error = "failed to write navigation map manifest: " + yaml_path.string();
            return false;
        }
    } catch (const std::exception& exception) {
        error = std::string("failed to write navigation map manifest: ") +
                exception.what();
        return false;
    }

    result.width = width;
    result.height = height;
    result.origin_x = origin_x;
    result.origin_y = origin_y;
    result.occupied_cells = static_cast<std::size_t>(cv::countNonZero(obstacle_mask));
    result.free_cells = static_cast<std::size_t>(cv::countNonZero(image == 255));
    result.unknown_cells = static_cast<std::size_t>(cv::countNonZero(image == kUnknown));
    result.ray_count = ray_count;
    return true;
}

}  // namespace lightning::navigation_map
