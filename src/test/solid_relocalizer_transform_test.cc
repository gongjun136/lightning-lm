#include <algorithm>
#include <cmath>
#include <iostream>

#include "core/localization/solid_relocalizer.h"

namespace {

constexpr double kTolerance = 1e-8;

lightning::Quatd Rotation(double roll, double pitch, double yaw) {
    return lightning::Quatd(
        lightning::AngAxisd(yaw, lightning::Vec3d::UnitZ()) *
        lightning::AngAxisd(pitch, lightning::Vec3d::UnitY()) *
        lightning::AngAxisd(roll, lightning::Vec3d::UnitX()));
}

bool PreservesMatchedPosition(const lightning::Vec3d& odom_translation) {
    const lightning::SE3 T_odom_current_lidar(
        Rotation(-0.09, 0.13, -1.1), odom_translation);
    const lightning::SE3 T_map_current_lidar(
        Rotation(0.02, -0.03, 0.7),
        lightning::Vec3d(25.6, 18.4, 3.5));

    const lightning::SE3 T_map_odom = lightning::loc::MakePlanarMapOdom(
        T_map_current_lidar, T_odom_current_lidar);
    const lightning::SE3 projected_lidar =
        T_map_odom * T_odom_current_lidar;
    if ((projected_lidar.translation() -
         T_map_current_lidar.translation()).norm() > kTolerance) {
        return false;
    }

    const lightning::Mat3d rotation = T_map_odom.rotationMatrix();
    return std::abs(rotation(2, 0)) <= kTolerance &&
           std::abs(rotation(2, 1)) <= kTolerance &&
           std::abs(rotation(0, 2)) <= kTolerance &&
           std::abs(rotation(1, 2)) <= kTolerance &&
           std::abs(rotation(2, 2) - 1.0) <= kTolerance;
}

bool CandidateSortingIsSafe() {
    lightning::loc::RelocalizationCandidateVector candidates;
    for (int index = 0; index < 513; ++index) {
        lightning::loc::RelocalizationCandidate candidate;
        candidate.candidate_id = index;
        candidate.score = static_cast<double>(index % 17);
        candidate.T_world_imu = lightning::SE3(
            Rotation(0.01 * (index % 3), -0.02 * (index % 5), 0.03 * index),
            lightning::Vec3d(index * 1e6, -index * 2e6, index));
        candidates.push_back(std::move(candidate));
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const auto& left, const auto& right) {
                  if (left.score != right.score) {
                      return left.score > right.score;
                  }
                  if (left.candidate_id != right.candidate_id) {
                      return left.candidate_id < right.candidate_id;
                  }
                  return left.query_submap_size < right.query_submap_size;
              });
    return candidates.size() == 513 && candidates.front().score == 16.0 &&
           candidates.back().score == 0.0;
}

}  // namespace

int main() {
    if (!PreservesMatchedPosition(lightning::Vec3d(4.0, -8.0, 1.5)) ||
        !PreservesMatchedPosition(
            lightning::Vec3d(1.4e7, -3.1e6, 1.9e7))) {
        std::cerr << "planar map<-odom projection moved the matched lidar position\n";
        return 1;
    }
    if (!CandidateSortingIsSafe()) {
        std::cerr << "aligned relocalization candidate sorting failed\n";
        return 2;
    }
    return 0;
}
