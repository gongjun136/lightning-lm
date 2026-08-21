#include "core/localization/pose_graph/pgo.h"
#include "core/localization/pose_graph/smoother.h"

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

    PoseSmoother smoother;
    Require(smoother.PushDRPose(SE3(Quatd::Identity(), Vec3d::Zero()), 100.0, 3.0),
            "smoother accepts the first timed DR pose");
    Require(smoother.PushDRPose(
                SE3(Quatd::Identity(), Vec3d(0.63, 0.0, 0.0)), 100.2, 3.0),
            "time-aware smoother accepts plausible 0.63 m motion over 0.2 s");
    Require(smoother.PushDRPose(
                SE3(Quatd::Identity(), Vec3d(0.63, 0.0, 0.0)), 100.2, 3.0),
            "smoother treats a repeated PubResult DR sample as idempotent");
    Require(!smoother.PushDRPose(
                SE3(Quatd::Identity(), Vec3d(10.0, 0.0, 0.0)), 100.3, 3.0),
            "time-aware smoother rejects an implausible pose jump");

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

    PGO velocity_pgo;
    velocity_pgo.SetDebug(false);
    velocity_pgo.SetDrSmoothingEnabled(false);
    velocity_pgo.SetDrExtrapolationEnabled(true);
    LocalizationResult velocity_output;
    int velocity_output_count = 0;
    velocity_pgo.SetHighFrequencyGlobalOutputHandleFunction(
        [&](const LocalizationResult& result) {
            velocity_output = result;
            ++velocity_output_count;
        });

    const SO3 body_rotation = SO3::exp(Vec3d(0.0, 0.0, 0.4));
    NavState velocity_seed_before;
    velocity_seed_before.timestamp_ = 1999.99;
    velocity_seed_before.pose_is_ok_ = true;
    velocity_seed_before.lidar_odom_reliable_ = true;
    velocity_seed_before.rot_ = body_rotation;
    velocity_seed_before.SetVel(body_rotation * Vec3d(0.5, 0.0, 0.0));
    NavState velocity_seed_after = velocity_seed_before;
    velocity_seed_after.timestamp_ = 2000.01;
    Require(velocity_pgo.ProcessLidarOdom(velocity_seed_before), "accept first velocity-test lidar odometry");
    Require(velocity_pgo.ProcessLidarOdom(velocity_seed_after), "accept second velocity-test lidar odometry");
    Require(velocity_pgo.ProcessDR(velocity_seed_before), "accept first velocity-test DR");
    Require(velocity_pgo.ProcessDR(velocity_seed_after), "accept second velocity-test DR");

    LocalizationResult velocity_loc;
    velocity_loc.timestamp_ = 2000.0;
    velocity_loc.pose_ = velocity_seed_before.GetPose();
    velocity_loc.valid_ = true;
    velocity_loc.lidar_loc_valid_ = true;
    velocity_loc.lidar_loc_odom_error_normal_ = true;
    velocity_loc.lidar_loc_smooth_flag_ = true;
    velocity_loc.confidence_ = 1.0;
    velocity_loc.status_ = LocalizationStatus::GOOD;
    Require(velocity_pgo.ProcessLidarLoc(velocity_loc), "initialize velocity-test PGO result");

    NavState latest_dr = velocity_seed_after;
    latest_dr.timestamp_ = 2000.02;
    const Vec3d latest_body_velocity(2.25, -0.4, 0.1);
    latest_dr.SetVel(body_rotation * latest_body_velocity);
    Require(velocity_pgo.ProcessDR(latest_dr), "accept high-frequency DR velocity update");
    Require(std::abs(velocity_output.timestamp_ - latest_dr.timestamp_) < 1e-9,
            "high-frequency output advances to the latest DR timestamp");
    Require((velocity_output.vel_b_ - latest_body_velocity).norm() < 1e-9,
            "high-frequency output refreshes body velocity from the latest DR state");

    NavState parked_dr = latest_dr;
    parked_dr.timestamp_ = 2000.03;
    parked_dr.is_parking_ = true;
    parked_dr.SetVel(Vec3d::Zero());
    Require(velocity_pgo.ProcessDR(parked_dr), "accept first stationary DR state");
    const SE3 parked_output_pose = velocity_output.pose_;
    parked_dr.timestamp_ = 2000.04;
    parked_dr.pos_ += Vec3d(10.0, 0.0, 0.0);
    Require(velocity_pgo.ProcessDR(parked_dr), "accept drifting DR while stationary hold is active");
    Require((velocity_output.pose_.translation() - parked_output_pose.translation()).norm() < 1e-9,
            "stationary hold rejects subsequent DR position drift");
    Require(velocity_output.is_parking_ && velocity_output.vel_b_.norm() < 1e-9,
            "stationary output is explicitly marked and has zero velocity");

    LocalizationResult parked_map_update;
    parked_map_update.timestamp_ = 2000.035;
    parked_map_update.pose_ = SE3(body_rotation, Vec3d(1.0, 2.0, 3.0));
    parked_map_update.valid_ = false;
    parked_map_update.lidar_loc_valid_ = true;
    parked_map_update.status_ = LocalizationStatus::GOOD;
    const int output_count_before_map_update = velocity_output_count;
    Require(velocity_pgo.ProcessLidarLoc(parked_map_update),
            "continue processing valid map matches while stationary");
    Require(velocity_output_count == output_count_before_map_update,
            "parked lidar match does not publish on the lidar timestamp");
    parked_dr.timestamp_ = 2000.045;
    Require(velocity_pgo.ProcessDR(parked_dr), "publish the held pose on the next DR sample");
    Require(std::abs(velocity_output.timestamp_ - parked_dr.timestamp_) < 1e-9,
            "parked output timestamps remain monotonic");
    Require((velocity_output.pose_.translation() - parked_output_pose.translation()).norm() < 1e-9,
            "stationary output rejects scan-to-map pose jitter");

    NavState moving_dr = parked_dr;
    moving_dr.timestamp_ = 2000.05;
    moving_dr.is_parking_ = false;
    moving_dr.SetVel(body_rotation * latest_body_velocity);
    Require(velocity_pgo.ProcessDR(moving_dr), "resume PGO extrapolation after stationary hold");
    Require(!velocity_output.is_parking_, "motion releases the PGO stationary output");

    // Reproduce a hold long enough to evict the pre-hold fused timestamp from
    // the bounded DR queue.  The first moving result must continue from the
    // held sensor-time epoch rather than publishing that evicted old result.
    NavState long_parked_dr = moving_dr;
    long_parked_dr.is_parking_ = true;
    long_parked_dr.SetVel(Vec3d::Zero());
    for (int index = 1; index <= 10010; ++index) {
        long_parked_dr.timestamp_ = moving_dr.timestamp_ + 0.01 * index;
        Require(velocity_pgo.ProcessDR(long_parked_dr),
                "accept DR throughout long stationary hold");
    }
    const double last_long_parked_stamp = velocity_output.timestamp_;
    NavState long_hold_exit = long_parked_dr;
    long_hold_exit.timestamp_ += 0.01;
    long_hold_exit.is_parking_ = false;
    long_hold_exit.SetVel(body_rotation * latest_body_velocity);
    Require(velocity_pgo.ProcessDR(long_hold_exit),
            "resume PGO after the pre-hold DR epoch was evicted");
    Require(velocity_output.timestamp_ > last_long_parked_stamp &&
                std::abs(velocity_output.timestamp_ - long_hold_exit.timestamp_) < 1e-9,
            "long stationary hold exit does not release the frozen old timestamp");
    Require(!velocity_output.is_parking_,
            "long stationary hold exit clears the parking marker");

    PGO parked_initialization_pgo;
    parked_initialization_pgo.SetDebug(false);
    parked_initialization_pgo.SetDrSmoothingEnabled(false);
    parked_initialization_pgo.SetDrExtrapolationEnabled(false);
    LocalizationResult parked_initialization_output;
    parked_initialization_pgo.SetHighFrequencyGlobalOutputHandleFunction(
        [&](const LocalizationResult& result) { parked_initialization_output = result; });
    NavState parked_initialization_state;
    parked_initialization_state.timestamp_ = 2500.0;
    parked_initialization_state.pose_is_ok_ = true;
    parked_initialization_state.lidar_odom_reliable_ = true;
    parked_initialization_state.is_parking_ = true;
    Require(parked_initialization_pgo.ProcessLidarOdom(parked_initialization_state),
            "accept parked initialization lidar odometry");
    Require(parked_initialization_pgo.ProcessDR(parked_initialization_state),
            "accept parked initialization DR");
    parked_initialization_state.timestamp_ = 2500.2;
    Require(parked_initialization_pgo.ProcessLidarOdom(parked_initialization_state),
            "accept second parked initialization lidar odometry");
    Require(parked_initialization_pgo.ProcessDR(parked_initialization_state),
            "accept second parked initialization DR");
    LocalizationResult parked_initialization_loc;
    parked_initialization_loc.timestamp_ = 2500.2;
    parked_initialization_loc.pose_ = T_world_odom;
    parked_initialization_loc.valid_ = false;
    parked_initialization_loc.lidar_loc_valid_ = true;
    parked_initialization_loc.lidar_loc_odom_error_normal_ = true;
    parked_initialization_loc.lidar_loc_smooth_flag_ = true;
    parked_initialization_loc.confidence_ = 1.0;
    parked_initialization_loc.status_ = LocalizationStatus::GOOD;
    Require(parked_initialization_pgo.ProcessLidarLoc(parked_initialization_loc),
            "initialize PGO while stationary hold is requested");
    Require(parked_initialization_output.valid_,
            "stationary request does not suppress the first fused map pose");
    parked_initialization_state.timestamp_ = 2500.3;
    Require(parked_initialization_pgo.ProcessDR(parked_initialization_state),
            "enter stationary output after parked initialization");
    Require(parked_initialization_output.valid_ && parked_initialization_output.is_parking_,
            "parked initialization transitions to a valid stationary output");

    PGO recovery_pgo;
    recovery_pgo.SetDebug(false);
    recovery_pgo.SetDrSmoothingEnabled(false);
    recovery_pgo.SetDrExtrapolationEnabled(false);
    NavState recovery_relative_pose;
    recovery_relative_pose.timestamp_ = 3000.0;
    recovery_relative_pose.pose_is_ok_ = true;
    recovery_relative_pose.lidar_odom_reliable_ = true;
    Require(recovery_pgo.ProcessLidarOdom(recovery_relative_pose), "accept recovery-test lidar odometry");
    NavState rolled_back_pose = recovery_relative_pose;
    rolled_back_pose.timestamp_ = 2999.0;
    Require(!recovery_pgo.ProcessLidarOdom(rolled_back_pose), "reject lidar odometry timestamp rollback");

    LocalizationResult recovery_loc;
    recovery_loc.timestamp_ = 3000.0;
    recovery_loc.pose_ = SE3();
    recovery_loc.valid_ = true;
    recovery_loc.lidar_loc_valid_ = true;
    recovery_loc.lidar_loc_odom_error_normal_ = true;
    recovery_loc.lidar_loc_smooth_flag_ = true;
    recovery_loc.confidence_ = 1.0;
    recovery_loc.status_ = LocalizationStatus::GOOD;
    recovery_relative_pose.timestamp_ = 3000.5;
    Require(recovery_pgo.ProcessLidarOdom(recovery_relative_pose), "advance lidar odometry before reset");
    Require(recovery_pgo.ProcessDR(recovery_relative_pose), "advance DR before reset");
    Require(recovery_pgo.Reset(), "reset PGO for global relocalization");
    recovery_loc.timestamp_ = 3000.4;
    Require(!recovery_pgo.ProcessLidarLoc(recovery_loc), "drop queued localization from before reset watermark");
    recovery_relative_pose.timestamp_ = 3000.6;
    Require(recovery_pgo.ProcessLidarOdom(recovery_relative_pose), "accept first fresh lidar odometry after reset");
    Require(recovery_pgo.ProcessDR(recovery_relative_pose), "accept first fresh DR after reset");
    recovery_loc.timestamp_ = 3000.65;
    Require(!recovery_pgo.ProcessLidarLoc(recovery_loc), "wait for enough fresh relative poses to interpolate");
    recovery_relative_pose.timestamp_ = 3000.8;
    Require(recovery_pgo.ProcessLidarOdom(recovery_relative_pose), "accept second fresh lidar odometry after reset");
    Require(recovery_pgo.ProcessDR(recovery_relative_pose), "accept second fresh DR after reset");
    recovery_loc.timestamp_ = 3000.7;
    Require(recovery_pgo.ProcessLidarLoc(recovery_loc), "resume PGO with only fresh relative poses");

    std::cout << "localization_pgo_test passed" << std::endl;
    return 0;
}
