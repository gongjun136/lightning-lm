#include "core/maps/map_frame.h"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

namespace lightning::map_frame {
namespace {

template <typename T>
T ReadOr(const YAML::Node& node, const char* key, const T& fallback) {
    return node && node[key] ? node[key].as<T>() : fallback;
}

double Median(std::vector<double> values) {
    const std::size_t middle = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    const double upper = values[middle];
    if (values.size() % 2 != 0) return upper;
    std::nth_element(values.begin(), values.begin() + middle - 1, values.begin() + middle);
    return 0.5 * (values[middle - 1] + upper);
}

bool FitPlane(const std::vector<Vec3d>& points, Vec3d& normal, double& offset) {
    if (points.size() < 3) return false;
    Vec3d centroid = Vec3d::Zero();
    for (const auto& point : points) centroid += point;
    centroid /= static_cast<double>(points.size());

    Mat3d covariance = Mat3d::Zero();
    for (const auto& point : points) {
        const Vec3d centered = point - centroid;
        covariance += centered * centered.transpose();
    }
    covariance /= static_cast<double>(points.size());
    const Eigen::SelfAdjointEigenSolver<Mat3d> solver(covariance);
    if (solver.info() != Eigen::Success) return false;
    normal = solver.eigenvectors().col(0).normalized();
    if (!normal.allFinite()) return false;
    if (normal.z() < 0.0) normal = -normal;
    offset = -normal.dot(centroid);
    return std::isfinite(offset);
}

std::string MakeTransformId(const SE3& transform) {
    std::ostringstream id;
    const Vec3d translation = transform.translation();
    const Quatd quaternion = transform.unit_quaternion();
    id << "start_ground_plane:" << std::fixed << std::setprecision(12)
       << translation.x() << "," << translation.y() << "," << translation.z()
       << "," << quaternion.x() << "," << quaternion.y() << ","
       << quaternion.z() << "," << quaternion.w();
    return id.str();
}

bool ReadFiniteVector(const YAML::Node& node, const char* key, std::size_t size,
                      std::vector<double>& values, std::string& error) {
    if (!node || !node[key]) {
        error = std::string("map-frame metadata is missing ") + key;
        return false;
    }
    values = node[key].as<std::vector<double>>();
    if (values.size() != size) {
        error = std::string("map-frame metadata ") + key + " has invalid size";
        return false;
    }
    for (const double value : values) {
        if (!std::isfinite(value)) {
            error = std::string("map-frame metadata ") + key + " contains a non-finite value";
            return false;
        }
    }
    return true;
}

}  // namespace

bool ReadExportOptions(const YAML::Node& root, ExportOptions& options, std::string& error) {
    try {
        const YAML::Node node = root["map_export"];
        if (!node) return true;
        options.normalize_start_ground_z =
            ReadOr<bool>(node, "normalize_start_ground_z", options.normalize_start_ground_z);
        options.start_ground_search_radius =
            ReadOr<double>(node, "start_ground_search_radius", options.start_ground_search_radius);
        options.expected_sensor_height =
            ReadOr<double>(node, "expected_sensor_height", options.expected_sensor_height);
        options.expected_height_tolerance =
            ReadOr<double>(node, "expected_height_tolerance", options.expected_height_tolerance);
        options.histogram_bin_size =
            ReadOr<double>(node, "histogram_bin_size", options.histogram_bin_size);
        options.inlier_tolerance =
            ReadOr<double>(node, "inlier_tolerance", options.inlier_tolerance);
        options.max_median_absolute_deviation = ReadOr<double>(
            node, "max_median_absolute_deviation", options.max_median_absolute_deviation);
        options.min_ground_inliers =
            ReadOr<std::size_t>(node, "min_ground_inliers", options.min_ground_inliers);
        options.max_ground_tilt_deg =
            ReadOr<double>(node, "max_ground_tilt_deg", options.max_ground_tilt_deg);
        options.export_pgm = ReadOr<bool>(node, "export_pgm", options.export_pgm);
        options.pgm_resolution =
            ReadOr<double>(node, "pgm_resolution", options.pgm_resolution);
        options.pgm_min_height =
            ReadOr<double>(node, "pgm_min_height", options.pgm_min_height);
        options.pgm_max_height =
            ReadOr<double>(node, "pgm_max_height", options.pgm_max_height);
        options.pgm_obstacle_dilation_radius = ReadOr<double>(
            node, "pgm_obstacle_dilation_radius", options.pgm_obstacle_dilation_radius);
        options.pgm_min_observation_frames = ReadOr<int>(
            node, "pgm_min_observation_frames",
            options.pgm_min_observation_frames);
    } catch (const YAML::Exception& exception) {
        error = std::string("invalid map_export configuration: ") + exception.what();
        return false;
    }

    if (options.normalize_start_ground_z &&
        (!std::isfinite(options.start_ground_search_radius) ||
         options.start_ground_search_radius <= 0.0 ||
         !std::isfinite(options.expected_sensor_height) ||
         options.expected_sensor_height <= 0.0 ||
         !std::isfinite(options.expected_height_tolerance) ||
         options.expected_height_tolerance <= 0.0 ||
         !std::isfinite(options.histogram_bin_size) ||
         options.histogram_bin_size <= 0.0 ||
         !std::isfinite(options.inlier_tolerance) ||
         options.inlier_tolerance < options.histogram_bin_size ||
         !std::isfinite(options.max_median_absolute_deviation) ||
         options.max_median_absolute_deviation <= 0.0 ||
         options.min_ground_inliers == 0 ||
         !std::isfinite(options.max_ground_tilt_deg) ||
         options.max_ground_tilt_deg <= 0.0 ||
         options.max_ground_tilt_deg >= 90.0)) {
        error = "invalid start-ground map_export parameters";
        return false;
    }
    if (options.export_pgm &&
        (!std::isfinite(options.pgm_resolution) || options.pgm_resolution <= 0.0 ||
         !std::isfinite(options.pgm_min_height) ||
         !std::isfinite(options.pgm_max_height) ||
         options.pgm_min_height >= options.pgm_max_height ||
         !std::isfinite(options.pgm_obstacle_dilation_radius) ||
         options.pgm_obstacle_dilation_radius < 0.0 ||
         options.pgm_min_observation_frames <= 0 ||
         options.pgm_min_observation_frames > 255)) {
        error = "invalid PGM map_export parameters";
        return false;
    }
    return true;
}

bool EstimateStartGroundFrame(const CloudPtr& map, const SE3& start_pose,
                              const ExportOptions& options, Metadata& metadata,
                              std::string& error) {
    metadata = Metadata{};
    if (!options.normalize_start_ground_z) return true;
    if (!map || map->empty()) {
        error = "cannot estimate the start ground from an empty global map";
        return false;
    }

    const Vec3d start = start_pose.translation();
    if (!start.allFinite()) {
        error = "the start pose contains a non-finite translation";
        return false;
    }

    const double expected_ground_z = start.z() - options.expected_sensor_height;
    const double minimum_z = expected_ground_z - options.expected_height_tolerance;
    const double maximum_z = expected_ground_z + options.expected_height_tolerance;
    const double radius_squared =
        options.start_ground_search_radius * options.start_ground_search_radius;
    const std::size_t bin_count = static_cast<std::size_t>(
                                      std::ceil((maximum_z - minimum_z) /
                                                options.histogram_bin_size)) +
                                  1;
    std::vector<std::size_t> histogram(bin_count, 0);
    std::vector<Vec3d> candidates;
    candidates.reserve(map->size() / 20);

    for (const auto& point : map->points) {
        if (!std::isfinite(point.x) || !std::isfinite(point.y) ||
            !std::isfinite(point.z)) {
            continue;
        }
        const double dx = static_cast<double>(point.x) - start.x();
        const double dy = static_cast<double>(point.y) - start.y();
        const double z = point.z;
        if (dx * dx + dy * dy > radius_squared || z < minimum_z || z > maximum_z) {
            continue;
        }
        candidates.emplace_back(point.x, point.y, point.z);
        const auto bin = std::min(
            bin_count - 1,
            static_cast<std::size_t>(std::floor((z - minimum_z) /
                                                options.histogram_bin_size)));
        ++histogram[bin];
    }

    if (candidates.size() < options.min_ground_inliers) {
        error = "too few start-area points in the expected ground-height band";
        return false;
    }

    std::size_t best_bin = 0;
    for (std::size_t bin = 1; bin < histogram.size(); ++bin) {
        if (histogram[bin] > histogram[best_bin]) {
            best_bin = bin;
        } else if (histogram[bin] == histogram[best_bin]) {
            const double current_center =
                minimum_z + (static_cast<double>(bin) + 0.5) * options.histogram_bin_size;
            const double best_center =
                minimum_z + (static_cast<double>(best_bin) + 0.5) * options.histogram_bin_size;
            if (std::abs(current_center - expected_ground_z) <
                std::abs(best_center - expected_ground_z)) {
                best_bin = bin;
            }
        }
    }

    const double mode_center =
        minimum_z + (static_cast<double>(best_bin) + 0.5) * options.histogram_bin_size;
    std::vector<Vec3d> plane_inliers;
    plane_inliers.reserve(candidates.size());
    for (const auto& point : candidates) {
        if (std::abs(point.z() - mode_center) <= options.inlier_tolerance) {
            plane_inliers.push_back(point);
        }
    }
    if (plane_inliers.size() < options.min_ground_inliers) {
        error = "the dominant start-ground height has too few inliers";
        return false;
    }

    Vec3d ground_normal;
    double ground_offset = 0.0;
    if (!FitPlane(plane_inliers, ground_normal, ground_offset)) {
        error = "failed to fit the start-ground plane";
        return false;
    }

    for (int iteration = 0; iteration < 2; ++iteration) {
        plane_inliers.clear();
        for (const auto& point : candidates) {
            if (std::abs(ground_normal.dot(point) + ground_offset) <=
                options.inlier_tolerance) {
                plane_inliers.push_back(point);
            }
        }
        if (plane_inliers.size() < options.min_ground_inliers ||
            !FitPlane(plane_inliers, ground_normal, ground_offset)) {
            error = "the start-ground plane has too few inliers";
            return false;
        }
    }

    const double tilt_deg =
        std::acos(std::clamp(ground_normal.z(), -1.0, 1.0)) * 180.0 /
        std::acos(-1.0);
    if (!std::isfinite(tilt_deg) || tilt_deg > options.max_ground_tilt_deg) {
        error = "the estimated start-ground tilt exceeds max_ground_tilt_deg";
        return false;
    }
    if (std::abs(ground_normal.z()) < 1e-6) {
        error = "the estimated start-ground plane is vertical";
        return false;
    }

    std::vector<double> plane_residuals;
    plane_residuals.reserve(plane_inliers.size());
    for (const auto& point : plane_inliers) {
        plane_residuals.push_back(std::abs(ground_normal.dot(point) + ground_offset));
    }
    const double residual_median = Median(plane_residuals);
    if (!std::isfinite(residual_median) ||
        residual_median > options.max_median_absolute_deviation) {
        error = "the start-ground plane residual is not sufficiently concentrated";
        return false;
    }

    const double ground_z =
        -(ground_normal.x() * start.x() + ground_normal.y() * start.y() +
          ground_offset) /
        ground_normal.z();
    Quatd rotation = Quatd::FromTwoVectors(ground_normal, Vec3d::UnitZ());
    rotation.normalize();
    const SE3 transform(rotation, Vec3d(0.0, 0.0, ground_offset));
    metadata.normalized = true;
    metadata.transform_id = MakeTransformId(transform);
    metadata.T_export_slam = transform;
    metadata.ground_z_slam = ground_z;
    metadata.ground_z_export = 0.0;
    metadata.start_xy = start.head<2>();
    metadata.expected_ground_z = expected_ground_z;
    metadata.search_radius = options.start_ground_search_radius;
    metadata.histogram_bin_size = options.histogram_bin_size;
    metadata.inlier_tolerance = options.inlier_tolerance;
    metadata.median_absolute_deviation = residual_median;
    metadata.ground_normal_slam = ground_normal;
    metadata.ground_plane_offset_slam = ground_offset;
    metadata.ground_tilt_deg = tilt_deg;
    metadata.candidate_count = candidates.size();
    metadata.inlier_count = plane_inliers.size();
    return true;
}

void TransformCloudInPlace(const Metadata& metadata, const CloudPtr& cloud) {
    if (!metadata.normalized || !cloud) return;
    for (auto& point : cloud->points) {
        const Vec3d transformed = metadata.T_export_slam * Vec3d(point.x, point.y, point.z);
        point.x = static_cast<float>(transformed.x());
        point.y = static_cast<float>(transformed.y());
        point.z = static_cast<float>(transformed.z());
    }
    const Vec3d sensor_origin = metadata.T_export_slam *
                                cloud->sensor_origin_.head<3>().cast<double>();
    cloud->sensor_origin_.head<3>() = sensor_origin.cast<float>();
    cloud->sensor_orientation_ =
        (metadata.T_export_slam.unit_quaternion().cast<float>() *
         cloud->sensor_orientation_)
            .normalized();
}

YAML::Node MakeTransformReference(const Metadata& metadata) {
    YAML::Node node;
    node["transform_id"] = metadata.transform_id;
    node["normalization"] = "start_ground_plane";
    const Vec3d translation = metadata.T_export_slam.translation();
    const Quatd quaternion = metadata.T_export_slam.unit_quaternion();
    node["translation_xyz"] =
        std::vector<double>{translation.x(), translation.y(), translation.z()};
    node["quaternion_xyzw"] =
        std::vector<double>{quaternion.x(), quaternion.y(), quaternion.z(), quaternion.w()};
    return node;
}

bool ReadTransformReference(const YAML::Node& node, Metadata& metadata,
                            std::string& error) {
    try {
        if (!node || !node["transform_id"] || !node["normalization"]) {
            error = "map-frame transform reference is incomplete or unsupported";
            return false;
        }
        const std::string normalization = node["normalization"].as<std::string>();
        if (normalization != "start_ground_z" &&
            normalization != "start_ground_plane") {
            error = "map-frame transform reference is incomplete or unsupported";
            return false;
        }
        std::vector<double> translation;
        std::vector<double> quaternion_xyzw;
        if (!ReadFiniteVector(node, "translation_xyz", 3, translation, error) ||
            !ReadFiniteVector(node, "quaternion_xyzw", 4, quaternion_xyzw, error)) {
            return false;
        }
        Quatd quaternion(quaternion_xyzw[3], quaternion_xyzw[0],
                         quaternion_xyzw[1], quaternion_xyzw[2]);
        if (quaternion.norm() < 1e-9) {
            error = "map-frame transform quaternion has zero norm";
            return false;
        }
        quaternion.normalize();
        const double tilt_deg = quaternion.angularDistance(Quatd::Identity()) *
                                180.0 / std::acos(-1.0);
        if (std::abs(translation[0]) > 1e-9 || std::abs(translation[1]) > 1e-9 ||
            !std::isfinite(tilt_deg) || tilt_deg >= 90.0 ||
            (normalization == "start_ground_z" && tilt_deg > 1e-9)) {
            error = "map-frame transform is not a supported ground normalization";
            return false;
        }
        metadata = Metadata{};
        metadata.normalized = true;
        metadata.transform_id = node["transform_id"].as<std::string>();
        metadata.T_export_slam =
            SE3(quaternion, Vec3d(translation[0], translation[1], translation[2]));
        return true;
    } catch (const YAML::Exception& exception) {
        error = std::string("invalid map-frame transform reference: ") + exception.what();
        return false;
    }
}

bool SaveMetadata(const std::string& map_directory, const Metadata& metadata,
                  std::string& error) {
    if (!metadata.normalized) return true;
    try {
        YAML::Node root;
        root["schema_version"] = 1;
        root["complete"] = true;
        root["map_frame"] = MakeTransformReference(metadata);
        root["ground_z_slam"] = metadata.ground_z_slam;
        root["ground_z_export"] = metadata.ground_z_export;
        YAML::Node estimator;
        estimator["start_xy"] =
            std::vector<double>{metadata.start_xy.x(), metadata.start_xy.y()};
        estimator["expected_ground_z"] = metadata.expected_ground_z;
        estimator["search_radius"] = metadata.search_radius;
        estimator["histogram_bin_size"] = metadata.histogram_bin_size;
        estimator["inlier_tolerance"] = metadata.inlier_tolerance;
        estimator["median_absolute_deviation"] = metadata.median_absolute_deviation;
        estimator["ground_normal_slam"] =
            std::vector<double>{metadata.ground_normal_slam.x(),
                                metadata.ground_normal_slam.y(),
                                metadata.ground_normal_slam.z()};
        estimator["ground_plane_offset_slam"] = metadata.ground_plane_offset_slam;
        estimator["ground_tilt_deg"] = metadata.ground_tilt_deg;
        estimator["candidate_count"] = metadata.candidate_count;
        estimator["inlier_count"] = metadata.inlier_count;
        root["estimator"] = estimator;

        const std::filesystem::path path =
            std::filesystem::path(map_directory) / kMetadataFilename;
        std::ofstream output(path);
        if (!output.is_open()) {
            error = "failed to open map-frame metadata: " + path.string();
            return false;
        }
        output << root;
        output.close();
        if (!output) {
            error = "failed to write map-frame metadata: " + path.string();
            return false;
        }
        return true;
    } catch (const std::exception& exception) {
        error = std::string("failed to save map-frame metadata: ") + exception.what();
        return false;
    }
}

bool LoadMetadata(const std::string& map_directory, Metadata& metadata,
                  std::string& error) {
    const std::filesystem::path path =
        std::filesystem::path(map_directory) / kMetadataFilename;
    if (!std::filesystem::exists(path)) {
        error = "required map-frame metadata is missing: " + path.string();
        return false;
    }
    try {
        const YAML::Node root = YAML::LoadFile(path.string());
        if (!root["schema_version"] || root["schema_version"].as<int>() != 1 ||
            !root["complete"] || !root["complete"].as<bool>()) {
            error = "map-frame metadata is incomplete or uses an unsupported schema";
            return false;
        }
        if (!ReadTransformReference(root["map_frame"], metadata, error)) return false;
        if (!root["ground_z_slam"] || !root["ground_z_export"]) {
            error = "map-frame metadata is missing ground-height fields";
            return false;
        }
        metadata.ground_z_slam = root["ground_z_slam"].as<double>();
        metadata.ground_z_export = root["ground_z_export"].as<double>();
        if (!std::isfinite(metadata.ground_z_slam) ||
            !std::isfinite(metadata.ground_z_export) ||
            std::abs(metadata.ground_z_export) > 1e-9) {
            error = "map-frame ground height and transform are inconsistent";
            return false;
        }
        const YAML::Node estimator = root["estimator"];
        if (estimator) {
            const auto start_xy =
                ReadOr<std::vector<double>>(estimator, "start_xy", {});
            if (start_xy.size() == 2) metadata.start_xy = Vec2d(start_xy[0], start_xy[1]);
            metadata.expected_ground_z =
                ReadOr<double>(estimator, "expected_ground_z", 0.0);
            metadata.search_radius = ReadOr<double>(estimator, "search_radius", 0.0);
            metadata.histogram_bin_size =
                ReadOr<double>(estimator, "histogram_bin_size", 0.0);
            metadata.inlier_tolerance =
                ReadOr<double>(estimator, "inlier_tolerance", 0.0);
            metadata.median_absolute_deviation =
                ReadOr<double>(estimator, "median_absolute_deviation", 0.0);
            const auto ground_normal =
                ReadOr<std::vector<double>>(estimator, "ground_normal_slam", {});
            if (!ground_normal.empty()) {
                if (ground_normal.size() != 3) {
                    error = "map-frame metadata has an invalid ground normal";
                    return false;
                }
                metadata.ground_normal_slam =
                    Vec3d(ground_normal[0], ground_normal[1], ground_normal[2]);
            }
            metadata.ground_plane_offset_slam =
                ReadOr<double>(estimator, "ground_plane_offset_slam",
                               -metadata.ground_z_slam);
            metadata.ground_tilt_deg =
                ReadOr<double>(estimator, "ground_tilt_deg", 0.0);
            metadata.candidate_count =
                ReadOr<std::size_t>(estimator, "candidate_count", 0);
            metadata.inlier_count =
                ReadOr<std::size_t>(estimator, "inlier_count", 0);
        }
        const Vec3d start_ground(metadata.start_xy.x(), metadata.start_xy.y(),
                                 metadata.ground_z_slam);
        const Vec3d transformed_ground = metadata.T_export_slam * start_ground;
        if (!metadata.ground_normal_slam.allFinite() ||
            std::abs(metadata.ground_normal_slam.norm() - 1.0) > 1e-6 ||
            !std::isfinite(metadata.ground_plane_offset_slam) ||
            !std::isfinite(metadata.ground_tilt_deg) ||
            std::abs(transformed_ground.z() - metadata.ground_z_export) > 1e-6) {
            error = "map-frame ground plane and transform are inconsistent";
            return false;
        }
        return true;
    } catch (const std::exception& exception) {
        error = std::string("failed to load map-frame metadata: ") + exception.what();
        return false;
    }
}

bool SameTransform(const Metadata& lhs, const Metadata& rhs, double tolerance) {
    if (!lhs.normalized || !rhs.normalized ||
        lhs.transform_id != rhs.transform_id) {
        return false;
    }
    return (lhs.T_export_slam.translation() -
            rhs.T_export_slam.translation()).norm() <= tolerance &&
           lhs.T_export_slam.unit_quaternion().angularDistance(
               rhs.T_export_slam.unit_quaternion()) <= tolerance;
}

}  // namespace lightning::map_frame
