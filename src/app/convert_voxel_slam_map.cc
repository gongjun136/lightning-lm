#include <gflags/gflags.h>
#include <glog/logging.h>

#include <pcl/io/pcd_io.h>

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "common/point_def.h"
#include "core/lightning_math.hpp"
#include "core/maps/tiled_map.h"

DEFINE_string(input_dir, "", "Voxel-SLAM output directory containing numbered PCDs and alidarState.txt");
DEFINE_string(pose_file, "", "optimized Voxel-SLAM pose file; defaults to INPUT_DIR/alidarState.txt");
DEFINE_string(output_map_dir, "", "new Lightning tiled-map directory");
DEFINE_double(voxel_size, 0.1, "global and per-chunk voxel size in metres");
DEFINE_int32(frame_stride, 1, "use every Nth optimized scan");
DEFINE_int32(max_frames, 0, "optional maximum number of selected scans; disabled when <= 0");

namespace {

struct TimedPose {
    double timestamp = 0.0;
    lightning::SE3 pose;
};

struct VoxelKey {
    std::int64_t x = 0;
    std::int64_t y = 0;
    std::int64_t z = 0;

    bool operator==(const VoxelKey& other) const { return x == other.x && y == other.y && z == other.z; }
};

struct VoxelKeyHash {
    std::size_t operator()(const VoxelKey& key) const {
        std::size_t seed = std::hash<std::int64_t>{}(key.x);
        seed ^= std::hash<std::int64_t>{}(key.y) + 0x9e3779b9 + (seed << 6U) + (seed >> 2U);
        seed ^= std::hash<std::int64_t>{}(key.z) + 0x9e3779b9 + (seed << 6U) + (seed >> 2U);
        return seed;
    }
};

struct VoxelAccumulator {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double intensity = 0.0;
    double timestamp = 0.0;
    std::size_t count = 0;
};

bool ReadVoxelPoses(const std::filesystem::path& path, std::vector<TimedPose>& poses) {
    std::ifstream input(path);
    if (!input.is_open()) return false;

    std::string line;
    double previous_timestamp = -std::numeric_limits<double>::infinity();
    int line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty() || line.front() == '#') continue;
        std::istringstream stream(line);
        double timestamp, x, y, z, qx, qy, qz, qw;
        if (!(stream >> timestamp >> x >> y >> z >> qx >> qy >> qz >> qw)) {
            LOG(ERROR) << path << ':' << line_number << ": invalid pose row";
            return false;
        }
        if (!std::isfinite(timestamp) || !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z) ||
            !std::isfinite(qx) || !std::isfinite(qy) || !std::isfinite(qz) || !std::isfinite(qw) ||
            timestamp <= previous_timestamp) {
            LOG(ERROR) << path << ':' << line_number << ": non-finite or non-monotonic pose";
            return false;
        }
        lightning::Quatd quaternion(qw, qx, qy, qz);
        if (std::abs(quaternion.norm() - 1.0) > 1e-3) {
            LOG(ERROR) << path << ':' << line_number << ": invalid quaternion norm " << quaternion.norm();
            return false;
        }
        poses.push_back({timestamp, lightning::SE3(quaternion.normalized(), lightning::Vec3d(x, y, z))});
        previous_timestamp = timestamp;
    }
    return !poses.empty();
}

}  // namespace

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_colorlogtostderr = true;
    FLAGS_stderrthreshold = google::INFO;
    google::ParseCommandLineFlags(&argc, &argv, true);

    if (FLAGS_input_dir.empty() || FLAGS_output_map_dir.empty() || FLAGS_voxel_size <= 0.0 ||
        FLAGS_frame_stride <= 0) {
        LOG(ERROR) << "input_dir, output_map_dir, positive voxel_size and positive frame_stride are required";
        return 2;
    }

    const std::filesystem::path input_dir = std::filesystem::absolute(FLAGS_input_dir);
    const std::filesystem::path pose_path = FLAGS_pose_file.empty()
                                                ? input_dir / "alidarState.txt"
                                                : std::filesystem::absolute(FLAGS_pose_file);
    const std::filesystem::path output_dir = std::filesystem::absolute(FLAGS_output_map_dir);
    if (!std::filesystem::is_directory(input_dir) || !std::filesystem::is_regular_file(pose_path)) {
        LOG(ERROR) << "missing Voxel-SLAM input directory or pose file";
        return 2;
    }
    if (std::filesystem::exists(output_dir) && !std::filesystem::is_empty(output_dir)) {
        LOG(ERROR) << "refusing to overwrite non-empty output directory: " << output_dir;
        return 2;
    }

    std::vector<TimedPose> poses;
    if (!ReadVoxelPoses(pose_path, poses)) {
        LOG(ERROR) << "failed to read optimized poses from " << pose_path;
        return 3;
    }

    std::unordered_map<VoxelKey, VoxelAccumulator, VoxelKeyHash> voxels;
    voxels.reserve(1U << 20U);
    std::size_t loaded_frames = 0;
    std::size_t input_points = 0;
    for (std::size_t index = 0; index < poses.size(); index += FLAGS_frame_stride) {
        if (FLAGS_max_frames > 0 && loaded_frames >= static_cast<std::size_t>(FLAGS_max_frames)) break;
        const std::filesystem::path cloud_path = input_dir / (std::to_string(index) + ".pcd");
        pcl::PointCloud<pcl::PointXYZI> local_cloud;
        if (pcl::io::loadPCDFile(cloud_path.string(), local_cloud) != 0 || local_cloud.empty()) {
            LOG(ERROR) << "failed to load required scan " << cloud_path;
            return 3;
        }
        input_points += local_cloud.size();
        for (const auto& local_point : local_cloud) {
            const lightning::Vec3d world = poses[index].pose * local_point.getVector3fMap().cast<double>();
            const VoxelKey key{
                static_cast<std::int64_t>(std::floor(world.x() / FLAGS_voxel_size)),
                static_cast<std::int64_t>(std::floor(world.y() / FLAGS_voxel_size)),
                static_cast<std::int64_t>(std::floor(world.z() / FLAGS_voxel_size)),
            };
            auto& accumulator = voxels[key];
            accumulator.x += world.x();
            accumulator.y += world.y();
            accumulator.z += world.z();
            accumulator.intensity += local_point.intensity;
            accumulator.timestamp += poses[index].timestamp;
            ++accumulator.count;
        }
        ++loaded_frames;
    }
    if (voxels.empty()) {
        LOG(ERROR) << "no map points were reconstructed";
        return 3;
    }

    auto filtered_map = std::make_shared<lightning::PointCloudType>();
    filtered_map->reserve(voxels.size());
    for (const auto& [key, accumulator] : voxels) {
        (void)key;
        const double inverse_count = 1.0 / static_cast<double>(accumulator.count);
        lightning::PointType point;
        point.x = static_cast<float>(accumulator.x * inverse_count);
        point.y = static_cast<float>(accumulator.y * inverse_count);
        point.z = static_cast<float>(accumulator.z * inverse_count);
        point.data[3] = 1.0F;
        point.intensity = static_cast<float>(accumulator.intensity * inverse_count);
        point.time = accumulator.timestamp * inverse_count;
        filtered_map->push_back(point);
    }
    voxels.clear();
    voxels.rehash(0);

    std::filesystem::create_directories(output_dir);
    lightning::TiledMap::Options options;
    options.map_path_ = output_dir.string();
    options.voxel_size_in_chunk_ = static_cast<float>(FLAGS_voxel_size);
    lightning::TiledMap tiled_map(options);
    if (!tiled_map.ConvertFromFullPCD(filtered_map, poses.front().pose, output_dir.string())) {
        LOG(ERROR) << "failed to export tiled map";
        return 4;
    }
    const std::filesystem::path global_path = output_dir / "global.pcd";
    if (pcl::io::savePCDFileBinaryCompressed(global_path.string(), *filtered_map) != 0) {
        LOG(ERROR) << "failed to save reconstructed global map";
        return 4;
    }

    std::ofstream metadata(output_dir / "conversion_metadata.txt");
    metadata << std::setprecision(17)
             << "source=voxel_slam_optimized_scans\n"
             << "input_dir=" << input_dir.string() << "\n"
             << "pose_file=" << pose_path.string() << "\n"
             << "pose_count=" << poses.size() << "\n"
             << "loaded_frames=" << loaded_frames << "\n"
             << "frame_stride=" << FLAGS_frame_stride << "\n"
             << "input_points=" << input_points << "\n"
             << "filtered_points=" << filtered_map->size() << "\n"
             << "voxel_size_m=" << FLAGS_voxel_size << "\n"
             << "voxelization=incremental_global_centroid_hash\n";
    metadata.close();
    LOG(INFO) << "converted " << loaded_frames << " Voxel-SLAM scans (" << input_points << " points) into "
              << filtered_map->size() << " tiled-map points at " << output_dir;
    return 0;
}
