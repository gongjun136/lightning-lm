#include "core/localization/lidar_loc/lidar_loc.h"
#include "core/localization/lidar_loc/pclomp/voxel_grid_covariance_omp_impl.hpp"

#include <cstdlib>
#include <iostream>
#include <random>

namespace lightning::loc {
struct LidarLocRegressionAccess {
    static auto Matcher(LidarLoc& loc) { return loc.pcl_ndt_; }
    static bool RefreshMap(LidarLoc& loc) {
        loc.map_ = std::make_shared<TiledMap>();
        return loc.UpdateGlobalMap();
    }
    static void SetMap(LidarLoc& loc, const CloudPtr& cloud) {
        loc.current_lo_pose_set_ = true;
        loc.current_lo_pose_ = SE3();
        loc.last_lo_pose_set_ = true;
        loc.last_lo_pose_ = SE3();
        loc.last_abs_pose_ = SE3();
        loc.relocalization_static_map_ = cloud;
        loc.relocalization_map_min_ = Vec3d::Constant(-2.0);
        loc.relocalization_map_max_ = Vec3d::Constant(5.0);
        auto tree = std::make_unique<pcl::KdTreeFLANN<PointType>>();
        tree->setInputCloud(cloud);
        loc.relocalization_kdtrees_.push_back(std::move(tree));
    }
    static bool MapAgrees(LidarLoc& loc, const CloudPtr& scan, const SE3& pose) {
        return loc.EvaluateRelocalizationMapConsistency(scan, pose, 0, 0.75, 2000).passed;
    }
    static bool RawOdomAgrees(LidarLoc& loc, const SE3& pose) {
        double delta = 0.0;
        return loc.CheckLidarOdomValid(pose, delta);
    }
    static SE3 LastPose(const LidarLoc& loc) { return loc.last_abs_pose_; }
};
}  // namespace lightning::loc

namespace {
void Require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        std::exit(1);
    }
}
}

int main() {
    using namespace lightning;
    using namespace lightning::loc;
    using Access = LidarLocRegressionAccess;
    LidarLoc loc;
    Require(Access::RefreshMap(loc), "refresh tracking map matcher");
    auto ndt = Access::Matcher(loc);
    Require(ndt->getMaximumIterations() == 20, "map refresh retains the 20-iteration tracking budget");
    Require(std::abs(ndt->getTransformationEpsilon() - 0.01) < 1e-9,
            "map refresh retains the tracking convergence tolerance");
    Require(std::abs(ndt->getOulierRatio() - 0.45) < 1e-9,
            "map refresh retains the NDT score normalization");
    ndt->setNumThreads(1);
    CloudPtr map(new PointCloudType);
    std::mt19937 rng(20260925);
    std::normal_distribution<float> noise(0.0f, 0.09f);
    for (int x = 0; x < 4; ++x) {
        for (int y = 0; y < 4; ++y) {
            for (int z = 0; z < 3; ++z) {
                for (int n = 0; n < 40; ++n) {
                    PointType point{};
                    // This implementation rounds voxel keys: integer locations
                    // are cell centers, while half-integers are boundaries.
                    point.x = x + noise(rng);
                    point.y = y + noise(rng);
                    point.z = z + noise(rng);
                    point.data[3] = 1.0f;
                    map->push_back(point);
                }
            }
        }
    }
    ndt->AddTarget(map);
    ndt->ComputeTargetGrids();
    ndt->setInputSource(map);
    ndt->setMaximumIterations(1);
    PointCloudType output;
    Eigen::Matrix4f guess = Eigen::Matrix4f::Identity();
    // Start within this narrow synthetic Gaussian's convergence basin, while
    // still requiring more than one 0.1 m step to recover the map pose.
    guess(0, 3) = 0.12f;
    ndt->align(output, guess);
    Require(ndt->getFinalNumIteration() == 1, "NDT respects its actual iteration budget");
    Require(!ndt->hasConverged(), "iteration exhaustion is not convergence");
    ndt->setMaximumIterations(20);
    ndt->align(output, guess);
    Require(ndt->hasConverged(), "a sufficient iteration budget converges on the same scan");
    std::cout << "NDT regression: iterations=" << ndt->getFinalNumIteration()
              << " translation=" << ndt->getFinalTransformation().block<3, 1>(0, 3).transpose()
              << " score=" << ndt->getTransformationProbability() << std::endl;
    Require(ndt->getFinalTransformation().block<3, 1>(0, 3).norm() < 0.05,
            "converged registration recovers the map pose");

    Access::SetMap(loc, map);
    Require(Access::MapAgrees(loc, map, SE3()), "whole-scan map validation accepts correct alignment");
    Require(!Access::MapAgrees(loc, map, SE3(SO3(), Vec3d(40.0, 0.0, 0.0))),
            "whole-scan map validation rejects a 40 m displaced pose independently of NDT score");
    Require(!Access::RawOdomAgrees(loc, SE3(SO3(), Vec3d(2.0, 0.0, 0.0))),
            "raw 2 m correction is rejected before smoothing can turn it into 0.2 m");
    Require(Access::LastPose(loc).translation().norm() == 0.0,
            "rejected raw registration does not contaminate the next pose prior");
    std::cout << "lidar tracking health regression tests passed\n";
}
