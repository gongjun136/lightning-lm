#include "core/localization/pose_graph/pgo.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

void Require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << std::endl;
        std::exit(1);
    }
}

}  // namespace

int main() {
    using namespace lightning;
    using namespace lightning::loc;

    PGO pgo;
    pgo.SetDebug(false);
    pgo.SetDrSmoothingEnabled(false);
    pgo.SetDrExtrapolationEnabled(false);

    LocalizationResult output;
    int output_count = 0;
    LocalizationResult high_frequency_output;
    pgo.SetGlobalOutputHandleFunction([&](const LocalizationResult& result) {
        output = result;
        ++output_count;
    });
    pgo.SetHighFrequencyGlobalOutputHandleFunction(
        [&](const LocalizationResult& result) { high_frequency_output = result; });

    const SE3 T_world_odom(SO3::exp(Vec3d(0.15, -0.08, 1.1)), Vec3d(12.0, -18.0, 2.5));
    NavState seed_lidar_odom;
    seed_lidar_odom.timestamp_ = 999.9;
    seed_lidar_odom.confidence_ = 1.0;
    seed_lidar_odom.pose_is_ok_ = true;
    seed_lidar_odom.lidar_odom_reliable_ = true;
    seed_lidar_odom.SetPose(SE3());
    Require(pgo.ProcessLidarOdom(seed_lidar_odom), "accept seed lidar odometry");

    for (int index = 0; index < 20; ++index) {
        const double timestamp = 1000.0 + 0.1 * index;
        const SE3 T_odom_imu(SO3::exp(Vec3d(0.001 * index, -0.0005 * index, 0.002 * index)),
                             Vec3d(0.12 * index, 0.02 * std::sin(0.2 * index), 0.01 * index));

        NavState lidar_odom;
        lidar_odom.timestamp_ = timestamp;
        lidar_odom.confidence_ = 1.0;
        lidar_odom.pose_is_ok_ = true;
        lidar_odom.lidar_odom_reliable_ = true;
        lidar_odom.SetPose(T_odom_imu);
        Require(pgo.ProcessLidarOdom(lidar_odom), "accept synthetic lidar odometry");
        Require(pgo.ProcessDR(lidar_odom), "accept synthetic DR pose");

        LocalizationResult lidar_loc;
        lidar_loc.timestamp_ = timestamp;
        lidar_loc.pose_ = T_world_odom * T_odom_imu;
        lidar_loc.valid_ = true;
        lidar_loc.lidar_loc_valid_ = true;
        lidar_loc.lidar_loc_odom_error_normal_ = true;
        lidar_loc.lidar_loc_smooth_flag_ = true;
        lidar_loc.confidence_ = 1.0;
        lidar_loc.status_ = LocalizationStatus::GOOD;
        Require(pgo.ProcessLidarLoc(lidar_loc), "accept synthetic map localization");
        Require(output_count == index + 1, "produce one PGO result per map-localization frame");
        Require(output.valid_, "PGO result is valid");
        Require((output.pose_.translation() - lidar_loc.pose_.translation()).norm() < 1e-3,
                "non-identity global translation remains stable across incremental window replacement");
        Require((output.pose_.so3().inverse() * lidar_loc.pose_.so3()).log().norm() < 1e-3,
                "non-identity global rotation remains stable across incremental window replacement");
        Require(high_frequency_output.valid_, "high-frequency PGO result is valid");
        Require((high_frequency_output.pose_.translation() - lidar_loc.pose_.translation()).norm() < 1e-3,
                "disabled DR extrapolation preserves scan-time map translation");
    }

    std::cout << "localization_pgo_test passed" << std::endl;
    return 0;
}
