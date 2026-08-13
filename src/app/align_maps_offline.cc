#include <gflags/gflags.h>
#include <glog/logging.h>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <pcl/common/centroid.h>
#include <pcl/common/transforms.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/registration/icp.h>
#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "common/eigen_types.h"
#include "common/point_def.h"
#include "core/localization/btc_relocalizer.h"
#include "core/localization/solid_relocalizer.h"

DEFINE_string(old_map, "", "old annotated Lightning-LM map package");
DEFINE_string(new_map, "", "new Lightning-LM map package used by localization");
DEFINE_string(input_bag, "", "mapping ROS 2 bag directory (used to verify the data time range)");
DEFINE_string(config, "", "online localization YAML to update after successful alignment");
DEFINE_string(output_dir, "", "directory for metrics, overlays, and review artifacts");
DEFINE_bool(update_config, true, "write output.fixed_map_transform after all gates pass");
DEFINE_int32(max_static_submaps, 8, "maximum evenly spaced static BTC submaps used for validation");
DEFINE_double(static_translation_threshold, 0.15, "maximum start-pose displacement considered static (m)");
DEFINE_double(static_rotation_threshold_deg, 2.0, "maximum start-pose rotation considered static (deg)");
DEFINE_double(coarse_voxel, 1.5, "global coarse-map voxel size (m)");
DEFINE_double(fine_voxel, 0.45, "global validation-map voxel size (m)");

namespace fs = std::filesystem;
using Cloud = pcl::PointCloud<pcl::PointXYZI>;
using CloudPtr = Cloud::Ptr;

namespace {

constexpr double kRadToDeg = 180.0 / M_PI;

struct DatabaseEntry {
    int descriptor_id = -1;
    double timestamp = 0.0;
    lightning::SE3 T_map_lidar;
    fs::path cloud_path;
};

struct DirectionalMetrics {
    double overlap_05 = 0.0;
    double overlap_10 = 0.0;
    double overlap_20 = 0.0;
    double rmse_10 = std::numeric_limits<double>::infinity();
    double median_10 = std::numeric_limits<double>::infinity();
    double p95_10 = std::numeric_limits<double>::infinity();
    std::size_t points = 0;
};

struct SymmetricMetrics {
    DirectionalMetrics new_to_old;
    DirectionalMetrics old_to_new;
    double overlap_05 = 0.0;
    double overlap_10 = 0.0;
    double overlap_20 = 0.0;
    double rmse_10 = std::numeric_limits<double>::infinity();
};

struct StaticWindowMetric {
    int descriptor_id = -1;
    double timestamp = 0.0;
    DirectionalMetrics metrics;
    double correction_translation = std::numeric_limits<double>::infinity();
    double correction_rotation_deg = std::numeric_limits<double>::infinity();
    lightning::SE3 refined_transform;
};

struct StaticSummary {
    double mean_overlap_05 = 0.0;
    double mean_overlap_10 = 0.0;
    double min_overlap_10 = 0.0;
    double mean_rmse_10 = std::numeric_limits<double>::infinity();
    double median_correction_translation = std::numeric_limits<double>::infinity();
    double median_correction_rotation_deg = std::numeric_limits<double>::infinity();
    double max_transform_spread_translation = std::numeric_limits<double>::infinity();
    double max_transform_spread_rotation_deg = std::numeric_limits<double>::infinity();
    std::vector<StaticWindowMetric> windows;
};

struct Candidate {
    std::string origin;
    int query_descriptor_id = -1;
    double retrieval_score = 0.0;
    lightning::SE3 seed;
    lightning::SE3 transform;
    double rough_overlap_20 = 0.0;
    bool converged = false;
    SymmetricMetrics global;
    StaticSummary static_summary;
    double regional_p10_overlap_10 = 0.0;
    double quality = -std::numeric_limits<double>::infinity();
};

double RotationDistanceDeg(const lightning::SE3& left, const lightning::SE3& right) {
    const lightning::Mat3d relative = left.rotationMatrix().transpose() * right.rotationMatrix();
    const double cosine = std::clamp((relative.trace() - 1.0) * 0.5, -1.0, 1.0);
    return std::acos(cosine) * kRadToDeg;
}

double TranslationDistance(const lightning::SE3& left, const lightning::SE3& right) {
    return (left.translation() - right.translation()).norm();
}

lightning::SE3 ReadPose(const YAML::Node& entry) {
    const auto pose = entry["pose_xyzw"].as<std::vector<double>>();
    if (pose.size() != 7) throw std::runtime_error("database pose_xyzw must contain 7 values");
    lightning::Quatd quaternion(pose[6], pose[3], pose[4], pose[5]);
    if (quaternion.norm() < 1e-9) throw std::runtime_error("database pose has zero quaternion");
    quaternion.normalize();
    return lightning::SE3(quaternion, lightning::Vec3d(pose[0], pose[1], pose[2]));
}

std::vector<DatabaseEntry> LoadDatabaseEntries(const fs::path& map_path) {
    const fs::path directory = map_path / "btc_relocalization";
    const fs::path manifest_path = directory / "database.yaml";
    if (!fs::exists(manifest_path)) {
        throw std::runtime_error("BTC database not found: " + manifest_path.string());
    }
    const YAML::Node manifest = YAML::LoadFile(manifest_path.string());
    const YAML::Node entries = manifest["entries"];
    if (!entries || !entries.IsSequence() || entries.size() == 0) {
        throw std::runtime_error("BTC database has no entries: " + manifest_path.string());
    }
    std::vector<DatabaseEntry> result;
    result.reserve(entries.size());
    for (std::size_t index = 0; index < entries.size(); ++index) {
        const YAML::Node entry = entries[index];
        DatabaseEntry value;
        value.descriptor_id = entry["descriptor_id"].as<int>();
        value.timestamp = entry["timestamp"].as<double>();
        value.T_map_lidar = ReadPose(entry);
        value.cloud_path = directory / entry["cloud"].as<std::string>();
        if (!fs::exists(value.cloud_path)) {
            throw std::runtime_error("submap cloud not found: " + value.cloud_path.string());
        }
        result.push_back(std::move(value));
    }
    return result;
}

std::vector<std::size_t> DetectStaticEntryIndices(const std::vector<DatabaseEntry>& entries) {
    if (entries.empty()) return {};
    std::size_t last_static = 0;
    int consecutive_motion = 0;
    for (std::size_t index = 1; index < entries.size(); ++index) {
        const bool is_static =
            TranslationDistance(entries.front().T_map_lidar, entries[index].T_map_lidar) <=
                FLAGS_static_translation_threshold &&
            RotationDistanceDeg(entries.front().T_map_lidar, entries[index].T_map_lidar) <=
                FLAGS_static_rotation_threshold_deg;
        if (is_static) {
            last_static = index;
            consecutive_motion = 0;
        } else if (++consecutive_motion >= 3) {
            break;
        }
    }
    const std::size_t count = last_static + 1;
    const std::size_t requested = static_cast<std::size_t>(std::max(1, FLAGS_max_static_submaps));
    const std::size_t selected_count = std::min(count, requested);
    std::vector<std::size_t> selected;
    selected.reserve(selected_count);
    if (selected_count == 1) return {0};
    for (std::size_t i = 0; i < selected_count; ++i) {
        selected.push_back(static_cast<std::size_t>(std::llround(
            static_cast<double>(i) * static_cast<double>(count - 1) /
            static_cast<double>(selected_count - 1))));
    }
    selected.erase(std::unique(selected.begin(), selected.end()), selected.end());
    return selected;
}

std::pair<double, double> ReadBagTimeRange(const fs::path& bag_path) {
    const fs::path metadata_path = fs::is_directory(bag_path) ? bag_path / "metadata.yaml" : fs::path();
    if (metadata_path.empty() || !fs::exists(metadata_path)) {
        throw std::runtime_error("ROS 2 bag metadata.yaml not found: " + bag_path.string());
    }
    const YAML::Node info = YAML::LoadFile(metadata_path.string())["rosbag2_bagfile_information"];
    const double start = info["starting_time"]["nanoseconds_since_epoch"].as<double>() * 1e-9;
    const double duration = info["duration"]["nanoseconds"].as<double>() * 1e-9;
    return {start, start + duration};
}

CloudPtr LoadCloud(const fs::path& path) {
    CloudPtr cloud(new Cloud);
    if (pcl::io::loadPCDFile(path.string(), *cloud) != 0 || cloud->empty()) {
        throw std::runtime_error("failed to load point cloud: " + path.string());
    }
    return cloud;
}

CloudPtr LimitCloud(const CloudPtr& input, std::size_t maximum) {
    if (!input || input->size() <= maximum) return input;
    CloudPtr output(new Cloud);
    const std::size_t stride = (input->size() + maximum - 1) / maximum;
    output->reserve(maximum);
    for (std::size_t index = 0; index < input->size(); index += stride) {
        output->push_back((*input)[index]);
    }
    return output;
}

CloudPtr Downsample(const CloudPtr& input, double leaf, std::size_t maximum = 250000) {
    CloudPtr output(new Cloud);
    pcl::VoxelGrid<pcl::PointXYZI> voxel;
    const float size = static_cast<float>(leaf);
    voxel.setLeafSize(size, size, size);
    voxel.setInputCloud(input);
    voxel.filter(*output);
    return LimitCloud(output, maximum);
}

lightning::CloudPtr ToLightningCloud(const CloudPtr& input) {
    lightning::CloudPtr output(new lightning::PointCloudType);
    output->reserve(input->size());
    for (const auto& source : input->points) {
        lightning::PointType point;
        point.x = source.x;
        point.y = source.y;
        point.z = source.z;
        point.intensity = source.intensity;
        output->push_back(point);
    }
    return output;
}

lightning::SE3 MatrixToSE3(const Eigen::Matrix4f& matrix) {
    lightning::Mat3d rotation = matrix.block<3, 3>(0, 0).cast<double>();
    Eigen::JacobiSVD<lightning::Mat3d> svd(rotation, Eigen::ComputeFullU | Eigen::ComputeFullV);
    rotation = svd.matrixU() * svd.matrixV().transpose();
    if (rotation.determinant() < 0.0) {
        lightning::Mat3d u = svd.matrixU();
        u.col(2) *= -1.0;
        rotation = u * svd.matrixV().transpose();
    }
    return lightning::SE3(lightning::Quatd(rotation).normalized(),
                          matrix.block<3, 1>(0, 3).cast<double>());
}

double Percentile(std::vector<double> values, double fraction) {
    if (values.empty()) return std::numeric_limits<double>::infinity();
    const std::size_t index = std::min(values.size() - 1,
        static_cast<std::size_t>(std::floor(fraction * static_cast<double>(values.size() - 1))));
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(index), values.end());
    return values[index];
}

DirectionalMetrics EvaluateDirectional(const CloudPtr& source, const CloudPtr& target,
                                       const lightning::SE3& T_target_source) {
    DirectionalMetrics result;
    if (!source || !target || source->empty() || target->empty()) return result;
    pcl::KdTreeFLANN<pcl::PointXYZI> tree;
    tree.setInputCloud(target);
    std::vector<double> inliers;
    inliers.reserve(source->size());
    double squared_sum = 0.0;
    std::size_t count05 = 0, count10 = 0, count20 = 0;
    std::vector<int> indices(1);
    std::vector<float> squared_distances(1);
    for (const auto& point : source->points) {
        const lightning::Vec3d transformed =
            T_target_source * lightning::Vec3d(point.x, point.y, point.z);
        pcl::PointXYZI query;
        query.x = static_cast<float>(transformed.x());
        query.y = static_cast<float>(transformed.y());
        query.z = static_cast<float>(transformed.z());
        if (tree.nearestKSearch(query, 1, indices, squared_distances) <= 0) continue;
        const double distance = std::sqrt(squared_distances[0]);
        if (distance <= 0.5) ++count05;
        if (distance <= 1.0) {
            ++count10;
            squared_sum += distance * distance;
            inliers.push_back(distance);
        }
        if (distance <= 2.0) ++count20;
    }
    result.points = source->size();
    const double denominator = static_cast<double>(std::max<std::size_t>(1, source->size()));
    result.overlap_05 = static_cast<double>(count05) / denominator;
    result.overlap_10 = static_cast<double>(count10) / denominator;
    result.overlap_20 = static_cast<double>(count20) / denominator;
    if (!inliers.empty()) {
        result.rmse_10 = std::sqrt(squared_sum / static_cast<double>(inliers.size()));
        result.median_10 = Percentile(inliers, 0.50);
        result.p95_10 = Percentile(inliers, 0.95);
    }
    return result;
}

SymmetricMetrics EvaluateSymmetric(const CloudPtr& new_cloud, const CloudPtr& old_cloud,
                                   const lightning::SE3& T_old_new) {
    SymmetricMetrics result;
    result.new_to_old = EvaluateDirectional(new_cloud, old_cloud, T_old_new);
    result.old_to_new = EvaluateDirectional(old_cloud, new_cloud, T_old_new.inverse());
    result.overlap_05 = std::min(result.new_to_old.overlap_05, result.old_to_new.overlap_05);
    result.overlap_10 = std::min(result.new_to_old.overlap_10, result.old_to_new.overlap_10);
    result.overlap_20 = std::min(result.new_to_old.overlap_20, result.old_to_new.overlap_20);
    result.rmse_10 = std::max(result.new_to_old.rmse_10, result.old_to_new.rmse_10);
    return result;
}

bool RefineIcp(const CloudPtr& source, const CloudPtr& target, const lightning::SE3& initial,
               double maximum_correspondence, int iterations, lightning::SE3& refined) {
    pcl::IterativeClosestPoint<pcl::PointXYZI, pcl::PointXYZI> icp;
    icp.setInputSource(source);
    icp.setInputTarget(target);
    icp.setMaximumIterations(iterations);
    icp.setMaxCorrespondenceDistance(maximum_correspondence);
    icp.setTransformationEpsilon(1e-8);
    icp.setEuclideanFitnessEpsilon(1e-7);
    Cloud aligned;
    icp.align(aligned, initial.matrix().cast<float>());
    if (!icp.hasConverged() || !icp.getFinalTransformation().allFinite()) return false;
    refined = MatrixToSE3(icp.getFinalTransformation());
    return true;
}

std::vector<Candidate> BuildPcaSeeds(const CloudPtr& new_cloud, const CloudPtr& old_cloud) {
    Eigen::Vector4f new_centroid4, old_centroid4;
    pcl::compute3DCentroid(*new_cloud, new_centroid4);
    pcl::compute3DCentroid(*old_cloud, old_centroid4);
    Eigen::Matrix3f new_covariance, old_covariance;
    pcl::computeCovarianceMatrixNormalized(*new_cloud, new_centroid4, new_covariance);
    pcl::computeCovarianceMatrixNormalized(*old_cloud, old_centroid4, old_covariance);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> new_solver(new_covariance);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> old_solver(old_covariance);
    Eigen::Matrix3d new_basis, old_basis;
    for (int column = 0; column < 3; ++column) {
        new_basis.col(column) = new_solver.eigenvectors().col(2 - column).cast<double>();
        old_basis.col(column) = old_solver.eigenvectors().col(2 - column).cast<double>();
    }
    if (new_basis.determinant() < 0.0) new_basis.col(2) *= -1.0;
    if (old_basis.determinant() < 0.0) old_basis.col(2) *= -1.0;

    std::array<int, 3> permutation{0, 1, 2};
    std::vector<Candidate> result;
    int seed_id = 0;
    do {
        for (int signs = 0; signs < 8; ++signs) {
            lightning::Mat3d mapping = lightning::Mat3d::Zero();
            for (int row = 0; row < 3; ++row) {
                mapping(row, permutation[row]) = (signs & (1 << row)) ? -1.0 : 1.0;
            }
            const lightning::Mat3d rotation = old_basis * mapping * new_basis.transpose();
            if (rotation.determinant() < 0.999) continue;
            const lightning::Vec3d new_centroid = new_centroid4.head<3>().cast<double>();
            const lightning::Vec3d old_centroid = old_centroid4.head<3>().cast<double>();
            Candidate candidate;
            candidate.origin = "global_pca_" + std::to_string(seed_id++);
            candidate.seed = lightning::SE3(lightning::Quatd(rotation).normalized(),
                                             old_centroid - rotation * new_centroid);
            candidate.transform = candidate.seed;
            result.push_back(std::move(candidate));
        }
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return result;
}

fs::path WriteSearchConfig(const fs::path& config_path, const fs::path& output_dir) {
    YAML::Node root = YAML::LoadFile(config_path.string());
    YAML::Node relocalization = root["relocalization"];
    relocalization["enabled"] = true;
    relocalization["query_submap_sizes"] = std::vector<int>{1};
    relocalization["query_stride"] = 1;
    relocalization["top_k"] = 20;
    relocalization["min_btc_score"] = 0.10;
    YAML::Node solid = relocalization["solid"];
    solid["query_submap_sizes"] = std::vector<int>{1};
    solid["query_stride"] = 1;
    solid["top_k"] = 20;
    solid["retrieval_pool_size"] = 200;
    solid["min_similarity"] = 0.35;
    solid["refine_with_icp"] = true;
    solid["icp_batch_size"] = 20;
    solid["icp_workers"] = 8;
    solid["icp_yaw_hypothesis_offsets_deg"] = std::vector<double>{0.0, 90.0, -90.0, 180.0};
    const fs::path path = output_dir / "alignment_search_config.yaml";
    std::ofstream stream(path);
    stream << root;
    if (!stream.good()) throw std::runtime_error("failed to write search config: " + path.string());
    return path;
}

template <typename Relocalizer>
void CollectRelocalizationSeeds(const std::string& name, Relocalizer& relocalizer,
                                const std::vector<DatabaseEntry>& entries,
                                const std::vector<std::size_t>& selected,
                                std::vector<Candidate>& candidates) {
    for (const std::size_t index : selected) {
        const auto cloud = ToLightningCloud(LoadCloud(entries[index].cloud_path));
        relocalizer.ResetQuery();
        const auto result = relocalizer.AddFrame(
            cloud, entries[index].T_map_lidar, entries[index].timestamp);
        if (!result) continue;
        for (const auto& match : result->candidates) {
            Candidate candidate;
            candidate.origin = name;
            candidate.query_descriptor_id = entries[index].descriptor_id;
            candidate.retrieval_score = match.score;
            candidate.seed = match.T_world_imu * entries[index].T_map_lidar.inverse();
            candidate.transform = candidate.seed;
            candidates.push_back(std::move(candidate));
        }
    }
}

CloudPtr BuildStaticCloudInNewMap(const std::vector<DatabaseEntry>& entries,
                                  const std::vector<std::size_t>& selected) {
    CloudPtr combined(new Cloud);
    for (const std::size_t index : selected) {
        const CloudPtr local = LoadCloud(entries[index].cloud_path);
        Cloud transformed;
        pcl::transformPointCloud(*local, transformed, entries[index].T_map_lidar.matrix().cast<float>());
        *combined += transformed;
    }
    return Downsample(combined, 0.30, 180000);
}

StaticSummary EvaluateStaticWindows(const std::vector<DatabaseEntry>& entries,
                                    const std::vector<std::size_t>& selected,
                                    const CloudPtr& old_map, const lightning::SE3& T_old_new) {
    StaticSummary summary;
    std::vector<double> corrections_t, corrections_r, rmse;
    summary.min_overlap_10 = 1.0;
    for (const std::size_t index : selected) {
        const CloudPtr local = Downsample(LoadCloud(entries[index].cloud_path), 0.30, 60000);
        StaticWindowMetric window;
        window.descriptor_id = entries[index].descriptor_id;
        window.timestamp = entries[index].timestamp;
        const lightning::SE3 initial_old_lidar = T_old_new * entries[index].T_map_lidar;
        window.metrics = EvaluateDirectional(local, old_map, initial_old_lidar);
        lightning::SE3 refined_old_lidar = initial_old_lidar;
        if (RefineIcp(local, old_map, initial_old_lidar, 1.25, 50, refined_old_lidar)) {
            window.refined_transform = refined_old_lidar * entries[index].T_map_lidar.inverse();
            window.correction_translation = TranslationDistance(T_old_new, window.refined_transform);
            window.correction_rotation_deg = RotationDistanceDeg(T_old_new, window.refined_transform);
        }
        summary.mean_overlap_05 += window.metrics.overlap_05;
        summary.mean_overlap_10 += window.metrics.overlap_10;
        summary.min_overlap_10 = std::min(summary.min_overlap_10, window.metrics.overlap_10);
        if (std::isfinite(window.metrics.rmse_10)) rmse.push_back(window.metrics.rmse_10);
        if (std::isfinite(window.correction_translation)) corrections_t.push_back(window.correction_translation);
        if (std::isfinite(window.correction_rotation_deg)) corrections_r.push_back(window.correction_rotation_deg);
        summary.windows.push_back(std::move(window));
    }
    const double denominator = static_cast<double>(std::max<std::size_t>(1, summary.windows.size()));
    summary.mean_overlap_05 /= denominator;
    summary.mean_overlap_10 /= denominator;
    if (!rmse.empty()) summary.mean_rmse_10 = std::accumulate(rmse.begin(), rmse.end(), 0.0) / rmse.size();
    summary.median_correction_translation = Percentile(corrections_t, 0.50);
    summary.median_correction_rotation_deg = Percentile(corrections_r, 0.50);
    summary.max_transform_spread_translation = 0.0;
    summary.max_transform_spread_rotation_deg = 0.0;
    for (std::size_t i = 0; i < summary.windows.size(); ++i) {
        for (std::size_t j = i + 1; j < summary.windows.size(); ++j) {
            if (!std::isfinite(summary.windows[i].correction_translation) ||
                !std::isfinite(summary.windows[j].correction_translation)) continue;
            summary.max_transform_spread_translation = std::max(
                summary.max_transform_spread_translation,
                TranslationDistance(summary.windows[i].refined_transform,
                                    summary.windows[j].refined_transform));
            summary.max_transform_spread_rotation_deg = std::max(
                summary.max_transform_spread_rotation_deg,
                RotationDistanceDeg(summary.windows[i].refined_transform,
                                    summary.windows[j].refined_transform));
        }
    }
    return summary;
}

double RegionalP10Overlap(const CloudPtr& new_map, const CloudPtr& old_map,
                          const lightning::SE3& T_old_new) {
    if (!new_map || new_map->empty()) return 0.0;
    float min_x = std::numeric_limits<float>::infinity(), min_y = min_x;
    float max_x = -min_x, max_y = -min_x;
    for (const auto& point : new_map->points) {
        min_x = std::min(min_x, point.x); max_x = std::max(max_x, point.x);
        min_y = std::min(min_y, point.y); max_y = std::max(max_y, point.y);
    }
    constexpr int kGrid = 4;
    std::array<std::vector<std::size_t>, kGrid * kGrid> cells;
    for (std::size_t index = 0; index < new_map->size(); ++index) {
        const auto& point = (*new_map)[index];
        const int x = std::clamp(static_cast<int>(kGrid * (point.x - min_x) /
                          std::max(1e-6F, max_x - min_x)), 0, kGrid - 1);
        const int y = std::clamp(static_cast<int>(kGrid * (point.y - min_y) /
                          std::max(1e-6F, max_y - min_y)), 0, kGrid - 1);
        cells[static_cast<std::size_t>(y * kGrid + x)].push_back(index);
    }
    pcl::KdTreeFLANN<pcl::PointXYZI> tree;
    tree.setInputCloud(old_map);
    std::vector<double> overlaps;
    for (const auto& cell : cells) {
        if (cell.size() < 100) continue;
        std::size_t inside = 0;
        std::vector<int> indices(1); std::vector<float> distances(1);
        for (const std::size_t index : cell) {
            const auto& point = (*new_map)[index];
            const lightning::Vec3d transformed =
                T_old_new * lightning::Vec3d(point.x, point.y, point.z);
            pcl::PointXYZI query;
            query.x = transformed.x(); query.y = transformed.y(); query.z = transformed.z();
            if (tree.nearestKSearch(query, 1, indices, distances) > 0 && distances[0] <= 1.0F) ++inside;
        }
        overlaps.push_back(static_cast<double>(inside) / static_cast<double>(cell.size()));
    }
    return overlaps.empty() ? 0.0 : Percentile(overlaps, 0.10);
}

YAML::Node PoseNode(const lightning::SE3& pose) {
    YAML::Node node;
    const auto q = pose.unit_quaternion();
    node["translation_xyz"] = std::vector<double>{pose.translation().x(), pose.translation().y(), pose.translation().z()};
    node["quaternion_xyzw"] = std::vector<double>{q.x(), q.y(), q.z(), q.w()};
    return node;
}

YAML::Node MetricsNode(const DirectionalMetrics& metrics) {
    YAML::Node node;
    node["points"] = metrics.points;
    node["overlap_0_5m"] = metrics.overlap_05;
    node["overlap_1_0m"] = metrics.overlap_10;
    node["overlap_2_0m"] = metrics.overlap_20;
    node["rmse_within_1_0m"] = metrics.rmse_10;
    node["median_within_1_0m"] = metrics.median_10;
    node["p95_within_1_0m"] = metrics.p95_10;
    return node;
}

void SaveReviewClouds(const fs::path& output_dir, const CloudPtr& old_map,
                      const CloudPtr& new_map, const CloudPtr& static_new,
                      const lightning::SE3& transform) {
    CloudPtr aligned_new(new Cloud), aligned_static(new Cloud);
    pcl::transformPointCloud(*new_map, *aligned_new, transform.matrix().cast<float>());
    pcl::transformPointCloud(*static_new, *aligned_static, transform.matrix().cast<float>());
    pcl::io::savePCDFileBinaryCompressed((output_dir / "old_map_review.pcd").string(), *old_map);
    pcl::io::savePCDFileBinaryCompressed((output_dir / "new_map_aligned_review.pcd").string(), *aligned_new);
    pcl::io::savePCDFileBinaryCompressed((output_dir / "static_window_aligned.pcd").string(), *aligned_static);

    pcl::PointCloud<pcl::PointXYZRGB> overlay;
    overlay.reserve(old_map->size() + aligned_new->size() + aligned_static->size());
    auto append = [&](const CloudPtr& cloud, std::uint8_t r, std::uint8_t g, std::uint8_t b) {
        for (const auto& point : cloud->points) {
            pcl::PointXYZRGB output;
            output.x = point.x; output.y = point.y; output.z = point.z;
            output.r = r; output.g = g; output.b = b;
            overlay.push_back(output);
        }
    };
    append(old_map, 180, 180, 180);
    append(aligned_new, 230, 70, 55);
    append(aligned_static, 30, 160, 255);
    pcl::io::savePCDFileBinaryCompressed((output_dir / "alignment_overlay_rgb.pcd").string(), overlay);

    float min_x = std::numeric_limits<float>::infinity(), min_y = min_x;
    float max_x = -min_x, max_y = -min_x;
    auto bounds = [&](const CloudPtr& cloud) {
        for (const auto& point : cloud->points) {
            min_x = std::min(min_x, point.x); max_x = std::max(max_x, point.x);
            min_y = std::min(min_y, point.y); max_y = std::max(max_y, point.y);
        }
    };
    bounds(old_map); bounds(aligned_new);
    constexpr int kSize = 1800, kPadding = 40;
    const double scale = std::min((kSize - 2.0 * kPadding) / std::max(1.0F, max_x - min_x),
                                  (kSize - 2.0 * kPadding) / std::max(1.0F, max_y - min_y));
    cv::Mat image(kSize, kSize, CV_8UC3, cv::Scalar(18, 18, 18));
    auto draw = [&](const CloudPtr& cloud, const cv::Scalar& color, int radius) {
        for (const auto& point : cloud->points) {
            const int x = static_cast<int>(kPadding + (point.x - min_x) * scale);
            const int y = static_cast<int>(kSize - kPadding - (point.y - min_y) * scale);
            if (x >= 0 && x < kSize && y >= 0 && y < kSize) cv::circle(image, {x, y}, radius, color, -1);
        }
    };
    draw(old_map, cv::Scalar(155, 155, 155), 1);
    draw(aligned_new, cv::Scalar(50, 65, 235), 1);
    draw(aligned_static, cv::Scalar(255, 165, 25), 2);
    cv::putText(image, "old=gray  aligned-new=red  static=blue", {45, 32},
                cv::FONT_HERSHEY_SIMPLEX, 0.75, cv::Scalar(245, 245, 245), 2);
    cv::imwrite((output_dir / "alignment_birdseye.png").string(), image);
}

std::string QuoteYaml(const std::string& value) {
    std::string escaped;
    for (const char character : value) {
        escaped += character == '\'' ? "''" : std::string(1, character);
    }
    return "'" + escaped + "'";
}

bool PatchConfig(const fs::path& config_path, const fs::path& output_dir,
                 const lightning::SE3& transform, double quality, std::string& error) {
    std::ifstream input(config_path);
    if (!input.is_open()) { error = "failed to read config"; return false; }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        lines.push_back(line);
    }
    fs::copy_file(config_path, output_dir / "localization_config_before.yaml",
                  fs::copy_options::overwrite_existing);
    const auto q = transform.unit_quaternion();
    std::ostringstream block;
    block << std::setprecision(15)
          << "  fixed_map_transform:\n"
          << "    enabled: true\n"
          << "    convention: target_from_localization\n"
          << "    source_frame: localization_map\n"
          << "    target_frame: map\n"
          << "    translation_xyz: [" << transform.translation().x() << ", "
          << transform.translation().y() << ", " << transform.translation().z() << "]\n"
          << "    quaternion_xyzw: [" << q.x() << ", " << q.y() << ", " << q.z() << ", " << q.w() << "]\n"
          << "    provenance:\n"
          << "      generated_by: align_maps_offline\n"
          << "      report: " << QuoteYaml((output_dir / "alignment_report.yaml").string()) << "\n"
          << "      quality: " << quality;
    std::vector<std::string> block_lines;
    std::istringstream block_stream(block.str());
    while (std::getline(block_stream, line)) block_lines.push_back(line);

    auto top_level = [](const std::string& value) {
        return !value.empty() && value[0] != ' ' && value[0] != '\t' && value[0] != '#';
    };
    std::size_t output_start = lines.size();
    for (std::size_t index = 0; index < lines.size(); ++index) {
        if (lines[index].rfind("output:", 0) == 0) { output_start = index; break; }
    }
    if (output_start == lines.size()) {
        if (!lines.empty() && !lines.back().empty()) lines.push_back("");
        lines.push_back("output:");
        lines.push_back("  map_frame: map");
        lines.insert(lines.end(), block_lines.begin(), block_lines.end());
    } else {
        std::size_t output_end = lines.size();
        for (std::size_t index = output_start + 1; index < lines.size(); ++index) {
            if (top_level(lines[index])) { output_end = index; break; }
        }
        std::size_t fixed_start = output_end;
        bool has_map_frame = false;
        for (std::size_t index = output_start + 1; index < output_end; ++index) {
            if (lines[index].rfind("  map_frame:", 0) == 0) has_map_frame = true;
            if (lines[index].rfind("  fixed_map_transform:", 0) == 0) { fixed_start = index; break; }
        }
        if (fixed_start < output_end) {
            std::size_t fixed_end = fixed_start + 1;
            while (fixed_end < output_end &&
                   (lines[fixed_end].empty() || lines[fixed_end][0] == ' ' || lines[fixed_end][0] == '\t')) {
                if (!lines[fixed_end].empty() && lines[fixed_end].rfind("  ", 0) == 0 &&
                    lines[fixed_end].rfind("    ", 0) != 0) break;
                ++fixed_end;
            }
            lines.erase(lines.begin() + static_cast<std::ptrdiff_t>(fixed_start),
                        lines.begin() + static_cast<std::ptrdiff_t>(fixed_end));
            output_end -= fixed_end - fixed_start;
        }
        std::size_t insertion = output_start + 1;
        if (!has_map_frame) {
            lines.insert(lines.begin() + static_cast<std::ptrdiff_t>(insertion++), "  map_frame: map");
        }
        lines.insert(lines.begin() + static_cast<std::ptrdiff_t>(insertion), block_lines.begin(), block_lines.end());
    }
    const fs::path temporary = config_path.string() + ".map_alignment.tmp";
    std::ofstream output(temporary);
    for (const auto& value : lines) output << value << '\n';
    output.close();
    try {
        YAML::LoadFile(temporary.string());
        fs::rename(temporary, config_path);
        fs::copy_file(config_path, output_dir / "localization_config_after.yaml",
                      fs::copy_options::overwrite_existing);
        return true;
    } catch (const std::exception& exception) {
        std::error_code ignored;
        fs::remove(temporary, ignored);
        error = exception.what();
        return false;
    }
}

}  // namespace

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_colorlogtostderr = true;
    FLAGS_stderrthreshold = google::INFO;
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    if (FLAGS_old_map.empty() || FLAGS_new_map.empty() || FLAGS_input_bag.empty() ||
        FLAGS_config.empty() || FLAGS_output_dir.empty()) {
        LOG(ERROR) << "--old_map, --new_map, --input_bag, --config, and --output_dir are required";
        return 2;
    }

    const fs::path old_map = fs::absolute(FLAGS_old_map);
    const fs::path new_map = fs::absolute(FLAGS_new_map);
    const fs::path bag = fs::absolute(FLAGS_input_bag);
    const fs::path config = fs::absolute(FLAGS_config);
    const fs::path output_dir = fs::absolute(FLAGS_output_dir);
    try {
        fs::create_directories(output_dir);
        const auto bag_range = ReadBagTimeRange(bag);
        const auto new_entries = LoadDatabaseEntries(new_map);
        const auto static_indices = DetectStaticEntryIndices(new_entries);
        if (static_indices.size() < 4) throw std::runtime_error("fewer than four stable start submaps detected");
        const double static_start = new_entries[static_indices.front()].timestamp;
        const double static_end = new_entries[static_indices.back()].timestamp;
        if (static_start < bag_range.first - 1.0 || static_end > bag_range.second + 1.0 ||
            static_start - bag_range.first > 60.0 || static_end - static_start < 3.0) {
            throw std::runtime_error("detected static submaps are not a valid beginning-of-bag window");
        }
        LOG(INFO) << "detected static map window " << std::fixed << std::setprecision(3)
                  << static_start << " .. " << static_end << " using " << static_indices.size()
                  << " submaps";

        const fs::path search_config = WriteSearchConfig(config, output_dir);
        std::vector<Candidate> seeds;

        const fs::path old_solid_manifest = old_map / "solid_relocalization" / "database.yaml";
        if (!fs::exists(old_solid_manifest)) {
            LOG(INFO) << "old map has no SOLiD database; generating it from BTC submaps";
            if (!lightning::loc::SolidRelocalizer::BuildDatabase(search_config.string(), old_map.string())) {
                LOG(WARNING) << "failed to generate old-map SOLiD database; continuing with fallbacks";
            }
        }
        if (fs::exists(old_solid_manifest)) {
            lightning::loc::SolidRelocalizer solid;
            if (solid.Init(search_config.string(), old_map.string(), lightning::SE3())) {
                CollectRelocalizationSeeds("solid", solid, new_entries, static_indices, seeds);
            }
        }
        lightning::loc::BtcRelocalizer btc;
        if (btc.Init(search_config.string(), old_map.string(), lightning::SE3())) {
            CollectRelocalizationSeeds("btc", btc, new_entries, static_indices, seeds);
        }
        LOG(INFO) << "descriptor candidate seeds: " << seeds.size();

        const CloudPtr old_full = LoadCloud(old_map / "global.pcd");
        const CloudPtr new_full = LoadCloud(new_map / "global.pcd");
        const CloudPtr old_coarse = Downsample(old_full, FLAGS_coarse_voxel, 120000);
        const CloudPtr new_coarse = Downsample(new_full, FLAGS_coarse_voxel, 120000);
        const CloudPtr old_fine = Downsample(old_full, FLAGS_fine_voxel, 250000);
        const CloudPtr new_fine = Downsample(new_full, FLAGS_fine_voxel, 250000);
        auto pca_seeds = BuildPcaSeeds(new_coarse, old_coarse);
        seeds.insert(seeds.end(), pca_seeds.begin(), pca_seeds.end());

        for (auto& candidate : seeds) {
            candidate.rough_overlap_20 =
                EvaluateDirectional(new_coarse, old_coarse, candidate.seed).overlap_20;
        }
        std::stable_sort(seeds.begin(), seeds.end(), [](const Candidate& left, const Candidate& right) {
            if (left.rough_overlap_20 != right.rough_overlap_20) return left.rough_overlap_20 > right.rough_overlap_20;
            return left.retrieval_score > right.retrieval_score;
        });
        std::vector<Candidate> selected_seeds;
        for (const auto& candidate : seeds) {
            const bool duplicate = std::any_of(selected_seeds.begin(), selected_seeds.end(), [&](const Candidate& selected) {
                return TranslationDistance(candidate.seed, selected.seed) < 0.75 &&
                       RotationDistanceDeg(candidate.seed, selected.seed) < 4.0;
            });
            const bool weak_global_pca = candidate.origin.rfind("global_pca_", 0) == 0 &&
                                         candidate.rough_overlap_20 < 0.03;
            if (!duplicate && !weak_global_pca) {
                selected_seeds.push_back(candidate);
            }
            if (selected_seeds.size() >= 16) break;
        }
        if (selected_seeds.empty()) throw std::runtime_error("no plausible descriptor or global geometry seed");

        const CloudPtr static_new = BuildStaticCloudInNewMap(new_entries, static_indices);
        std::vector<Candidate> evaluated;
        for (auto candidate : selected_seeds) {
            lightning::SE3 coarse = candidate.seed;
            if (!RefineIcp(new_coarse, old_coarse, candidate.seed, 6.0, 60, coarse)) continue;
            lightning::SE3 fine = coarse;
            if (!RefineIcp(new_fine, old_fine, coarse, 1.5, 70, fine)) continue;
            lightning::SE3 static_refined = fine;
            RefineIcp(static_new, old_fine, fine, 1.25, 60, static_refined);
            candidate.transform = static_refined;
            candidate.converged = true;
            candidate.global = EvaluateSymmetric(new_fine, old_fine, candidate.transform);
            candidate.static_summary = EvaluateStaticWindows(
                new_entries, static_indices, old_fine, candidate.transform);
            candidate.regional_p10_overlap_10 = RegionalP10Overlap(new_fine, old_fine, candidate.transform);
            candidate.quality = 0.30 * candidate.global.overlap_05 +
                                0.50 * candidate.global.overlap_10 +
                                0.35 * candidate.static_summary.mean_overlap_05 -
                                0.20 * candidate.global.rmse_10 -
                                0.15 * candidate.static_summary.mean_rmse_10 -
                                0.02 * candidate.static_summary.median_correction_translation;
            evaluated.push_back(std::move(candidate));
        }
        if (evaluated.empty()) throw std::runtime_error("all candidate refinements failed");
        std::stable_sort(evaluated.begin(), evaluated.end(), [](const Candidate& left, const Candidate& right) {
            return left.quality > right.quality;
        });
        const Candidate& best = evaluated.front();

        int solid_support = 0, btc_support = 0;
        for (const auto& seed : seeds) {
            if (TranslationDistance(seed.seed, best.transform) <= 2.0 &&
                RotationDistanceDeg(seed.seed, best.transform) <= 8.0) {
                if (seed.origin == "solid") ++solid_support;
                if (seed.origin == "btc") ++btc_support;
            }
        }
        bool ambiguous = false;
        for (std::size_t index = 1; index < evaluated.size(); ++index) {
            const auto& alternative = evaluated[index];
            if (TranslationDistance(best.transform, alternative.transform) > 1.0 ||
                RotationDistanceDeg(best.transform, alternative.transform) > 5.0) {
                ambiguous = alternative.quality >= best.quality - 0.05;
                break;
            }
        }

        std::vector<std::string> failed_gates;
        auto gate = [&](bool condition, const std::string& name) { if (!condition) failed_gates.push_back(name); };
        gate(best.global.overlap_05 >= 0.35, "global_symmetric_overlap_0_5m_below_0_35");
        gate(best.global.overlap_10 >= 0.55, "global_symmetric_overlap_1_0m_below_0_55");
        gate(best.global.rmse_10 <= 0.45, "global_symmetric_rmse_above_0_45m");
        gate(best.regional_p10_overlap_10 >= 0.25, "regional_p10_overlap_1_0m_below_0_25");
        gate(best.static_summary.mean_overlap_05 >= 0.55, "static_mean_overlap_0_5m_below_0_55");
        gate(best.static_summary.mean_overlap_10 >= 0.75, "static_mean_overlap_1_0m_below_0_75");
        gate(best.static_summary.min_overlap_10 >= 0.65, "static_min_overlap_1_0m_below_0_65");
        gate(best.static_summary.mean_rmse_10 <= 0.35, "static_rmse_above_0_35m");
        gate(best.static_summary.median_correction_translation <= 0.35,
             "static_icp_median_correction_above_0_35m");
        gate(best.static_summary.median_correction_rotation_deg <= 2.0,
             "static_icp_median_rotation_above_2deg");
        gate(best.static_summary.max_transform_spread_translation <= 0.50,
             "static_window_transform_spread_above_0_50m");
        gate(best.static_summary.max_transform_spread_rotation_deg <= 3.0,
             "static_window_rotation_spread_above_3deg");
        const bool descriptor_support = solid_support >= 2 || btc_support >= 2;
        const bool strong_geometry_fallback = best.global.overlap_10 >= 0.70 &&
            best.static_summary.mean_overlap_05 >= 0.70 && best.regional_p10_overlap_10 >= 0.40;
        gate(descriptor_support || strong_geometry_fallback, "insufficient_independent_candidate_support");
        gate(!ambiguous, "ambiguous_distinct_transform_with_similar_quality");
        bool success = failed_gates.empty();

        SaveReviewClouds(output_dir, old_fine, new_fine, static_new, best.transform);
        std::ofstream candidates_csv(output_dir / "alignment_candidates.csv");
        candidates_csv << "rank,origin,query_descriptor_id,retrieval_score,rough_overlap_2m,quality,"
                          "global_overlap_0_5m,global_overlap_1m,global_rmse_1m,static_overlap_0_5m,"
                          "static_overlap_1m,regional_p10_overlap_1m,tx,ty,tz,qx,qy,qz,qw\n";
        for (std::size_t index = 0; index < evaluated.size(); ++index) {
            const auto& candidate = evaluated[index];
            const auto q = candidate.transform.unit_quaternion();
            candidates_csv << index + 1 << ',' << candidate.origin << ',' << candidate.query_descriptor_id << ','
                           << candidate.retrieval_score << ',' << candidate.rough_overlap_20 << ','
                           << candidate.quality << ',' << candidate.global.overlap_05 << ','
                           << candidate.global.overlap_10 << ',' << candidate.global.rmse_10 << ','
                           << candidate.static_summary.mean_overlap_05 << ','
                           << candidate.static_summary.mean_overlap_10 << ','
                           << candidate.regional_p10_overlap_10 << ',' << candidate.transform.translation().x() << ','
                           << candidate.transform.translation().y() << ',' << candidate.transform.translation().z() << ','
                           << q.x() << ',' << q.y() << ',' << q.z() << ',' << q.w() << '\n';
        }

        YAML::Node report;
        report["schema_version"] = 1;
        report["success"] = success;
        report["old_map"] = old_map.string();
        report["new_map"] = new_map.string();
        report["input_bag"] = bag.string();
        report["config"] = config.string();
        report["transform_convention"] = "T_old_from_new; left-multiply online output poses only";
        report["transform"] = PoseNode(best.transform);
        report["quality"] = best.quality;
        report["candidate_support"]["solid"] = solid_support;
        report["candidate_support"]["btc"] = btc_support;
        report["candidate_support"]["strong_geometry_fallback"] = strong_geometry_fallback;
        report["ambiguous"] = ambiguous;
        report["bag_time_range"] = std::vector<double>{bag_range.first, bag_range.second};
        report["static_window"]["start"] = static_start;
        report["static_window"]["end"] = static_end;
        report["static_window"]["duration"] = static_end - static_start;
        report["static_window"]["selected_submaps"] = static_indices.size();
        report["global_metrics"]["new_to_old"] = MetricsNode(best.global.new_to_old);
        report["global_metrics"]["old_to_new"] = MetricsNode(best.global.old_to_new);
        report["global_metrics"]["symmetric_overlap_0_5m"] = best.global.overlap_05;
        report["global_metrics"]["symmetric_overlap_1_0m"] = best.global.overlap_10;
        report["global_metrics"]["symmetric_overlap_2_0m"] = best.global.overlap_20;
        report["global_metrics"]["symmetric_rmse_within_1_0m"] = best.global.rmse_10;
        report["global_metrics"]["regional_p10_overlap_1_0m"] = best.regional_p10_overlap_10;
        report["static_metrics"]["mean_overlap_0_5m"] = best.static_summary.mean_overlap_05;
        report["static_metrics"]["mean_overlap_1_0m"] = best.static_summary.mean_overlap_10;
        report["static_metrics"]["min_overlap_1_0m"] = best.static_summary.min_overlap_10;
        report["static_metrics"]["mean_rmse_within_1_0m"] = best.static_summary.mean_rmse_10;
        report["static_metrics"]["median_icp_correction_m"] = best.static_summary.median_correction_translation;
        report["static_metrics"]["median_icp_correction_deg"] = best.static_summary.median_correction_rotation_deg;
        report["static_metrics"]["max_transform_spread_m"] = best.static_summary.max_transform_spread_translation;
        report["static_metrics"]["max_transform_spread_deg"] = best.static_summary.max_transform_spread_rotation_deg;
        for (const auto& window : best.static_summary.windows) {
            YAML::Node node;
            node["descriptor_id"] = window.descriptor_id;
            node["timestamp"] = window.timestamp;
            node["metrics"] = MetricsNode(window.metrics);
            node["icp_correction_m"] = window.correction_translation;
            node["icp_correction_deg"] = window.correction_rotation_deg;
            report["static_windows"].push_back(node);
        }
        for (const auto& failure : failed_gates) report["failed_gates"].push_back(failure);
        report["thresholds"]["global_symmetric_overlap_0_5m"] = 0.35;
        report["thresholds"]["global_symmetric_overlap_1_0m"] = 0.55;
        report["thresholds"]["global_rmse_1_0m"] = 0.45;
        report["thresholds"]["static_mean_overlap_0_5m"] = 0.55;
        report["thresholds"]["static_mean_overlap_1_0m"] = 0.75;
        report["artifacts"] = std::vector<std::string>{"alignment_birdseye.png", "alignment_overlay_rgb.pcd",
            "old_map_review.pcd", "new_map_aligned_review.pcd", "static_window_aligned.pcd",
            "alignment_candidates.csv"};

        bool config_written = false;
        std::string config_error;
        if (success && FLAGS_update_config) {
            config_written = PatchConfig(config, output_dir, best.transform, best.quality, config_error);
            if (!config_written) {
                success = false;
                report["success"] = false;
                report["failed_gates"].push_back("config_update_failed: " + config_error);
            }
        }
        report["config_updated"] = config_written;
        std::ofstream report_stream(output_dir / "alignment_report.yaml");
        report_stream << std::setprecision(15) << report;

        std::ofstream review(output_dir / "MANUAL_REVIEW.md");
        review << "# Map alignment manual review\n\n"
               << "Automatic result: **" << (success ? "PASS" : "FAIL") << "**\n\n"
               << "- Inspect `alignment_birdseye.png` for global red/gray agreement.\n"
               << "- Open `alignment_overlay_rgb.pcd`; old is gray, aligned new is red, static data is blue.\n"
               << "- Review `alignment_report.yaml` and per-candidate `alignment_candidates.csv`.\n"
               << "- The transform convention is `T_old_from_new`; online localization remains in the new map internally.\n";

        LOG(INFO) << "MAP_ALIGNMENT_RESULT success=" << success
                  << " config_updated=" << config_written
                  << " quality=" << best.quality
                  << " global_overlap_1m=" << best.global.overlap_10
                  << " static_overlap_0.5m=" << best.static_summary.mean_overlap_05
                  << " output=" << output_dir;
        if (!success) {
            for (const auto& failure : failed_gates) LOG(ERROR) << "alignment gate failed: " << failure;
        }
        return success ? 0 : 1;
    } catch (const std::exception& exception) {
        std::error_code ignored;
        fs::create_directories(output_dir, ignored);
        std::ofstream failure(output_dir / "fatal_error.txt");
        failure << exception.what() << '\n';
        LOG(ERROR) << "map alignment failed: " << exception.what();
        return 2;
    }
}
