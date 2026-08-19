#include "core/localization/localization.h"

#include <chrono>
#include <future>
#include <iostream>
#include <mutex>

namespace lightning::loc {

class LocalizationLockingTestPeer {
   public:
    static std::unique_lock<std::mutex> LockBackend(Localization& localization) {
        return std::unique_lock<std::mutex>(localization.processing_mutex_);
    }

    static void ConfigureOverloadGuard(Localization& localization, double max_lag, double resume_lag) {
        localization.options_.online_mode_ = true;
        localization.online_sensor_max_lag_sec_ = max_lag;
        localization.online_sensor_resume_lag_sec_ = resume_lag;
    }

    static void SetSensorLag(Localization& localization, double enqueued, double processed) {
        std::lock_guard<std::mutex> lock(localization.runtime_stats_mutex_);
        localization.runtime_stats_.latest_enqueued_sensor_stamp = enqueued;
        localization.runtime_stats_.latest_processed_sensor_stamp = processed;
    }

    static bool ShouldThrottleLidar(Localization& localization, double timestamp) {
        return localization.ShouldThrottleLidarInput(timestamp);
    }

    static void ObserveSensorEnqueued(Localization& localization, double timestamp) {
        localization.ObserveSensorEnqueued(timestamp);
    }

    static void ObserveSensorProcessed(Localization& localization, double timestamp) {
        localization.ObserveSensorProcessed(timestamp);
    }

    static void EnableImuStaticHold(Localization& localization) {
        localization.imu_static_hold_enabled_ = true;
    }

    static void ObserveLio(Localization& localization, const NavState& state) {
        localization.ObserveLioForStaticDetector(state);
    }

    static void ObserveLidarLoc(Localization& localization, const LocalizationResult& result) {
        localization.ObserveLidarLocForStaticDetector(result);
    }

    static bool UpdateImuStaticState(Localization& localization, const IMUPtr& imu) {
        return localization.UpdateImuStaticState(imu);
    }
};

}  // namespace lightning::loc

int main() {
    using namespace std::chrono_literals;
    using lightning::loc::Localization;
    using lightning::loc::LocalizationLockingTestPeer;

    Localization localization;
    auto backend_lock = LocalizationLockingTestPeer::LockBackend(localization);
    auto cloud = std::make_shared<sensor_msgs::msg::PointCloud2>();
    std::promise<void> callback_started;
    auto callback_started_future = callback_started.get_future();

    auto callback = std::async(std::launch::async, [&]() {
        callback_started.set_value();
        localization.ProcessLidarMsg(cloud, 0);
    });
    if (callback_started_future.wait_for(2s) != std::future_status::ready) {
        backend_lock.unlock();
        callback.wait();
        std::cerr << "Point-cloud callback thread did not start" << std::endl;
        return 1;
    }
    const bool callback_blocked = callback.wait_for(1s) != std::future_status::ready;

    // Always release the lock before waiting for the future so a regression
    // reports a failure instead of hanging the test process.
    backend_lock.unlock();
    callback.wait();

    if (callback_blocked) {
        std::cerr << "Point-cloud input waited for the backend processing lock" << std::endl;
        return 1;
    }

    Localization guarded_localization;
    LocalizationLockingTestPeer::ConfigureOverloadGuard(guarded_localization, 0.5, 0.1);
    LocalizationLockingTestPeer::SetSensorLag(guarded_localization, 10.0, 9.4);
    if (LocalizationLockingTestPeer::ShouldThrottleLidar(guarded_localization, 10.0)) {
        std::cerr << "overload guard rejected a fresh lidar above the maximum lag" << std::endl;
        return 1;
    }
    LocalizationLockingTestPeer::SetSensorLag(guarded_localization, 10.0, 9.7);
    if (LocalizationLockingTestPeer::ShouldThrottleLidar(guarded_localization, 10.1)) {
        std::cerr << "overload warning state rejected a fresh lidar" << std::endl;
        return 1;
    }
    LocalizationLockingTestPeer::SetSensorLag(guarded_localization, 10.0, 9.95);
    if (LocalizationLockingTestPeer::ShouldThrottleLidar(guarded_localization, 10.2)) {
        std::cerr << "overload guard did not resume below the recovery lag" << std::endl;
        return 1;
    }
    if (guarded_localization.GetRuntimeStats().sensor_queue_dropped != 0) {
        std::cerr << "overload warning incorrectly counted an admission drop" << std::endl;
        return 1;
    }

    Localization reordered_localization;
    LocalizationLockingTestPeer::ObserveSensorEnqueued(reordered_localization, 20.0);
    LocalizationLockingTestPeer::ObserveSensorProcessed(reordered_localization, 20.0);
    LocalizationLockingTestPeer::ObserveSensorProcessed(reordered_localization, 18.0);
    const auto reordered_stats = reordered_localization.GetRuntimeStats();
    if (reordered_stats.latest_processed_sensor_stamp != 20.0 ||
        reordered_stats.current_sensor_lag_sec != 0.0 ||
        reordered_stats.severe_timestamp_rollback_count != 1 ||
        reordered_stats.worst_timestamp_rollback_sec != 2.0) {
        std::cerr << "out-of-order processing moved the monotonic sensor frontier" << std::endl;
        return 1;
    }

    Localization static_localization;
    LocalizationLockingTestPeer::EnableImuStaticHold(static_localization);
    bool static_active = false;
    for (int index = 0; index <= 120; ++index) {
        const double timestamp = 100.0 + 0.01 * index;
        if (index % 10 == 0) {
            lightning::NavState lio_state;
            lio_state.timestamp_ = timestamp;
            lio_state.pose_is_ok_ = true;
            lio_state.lidar_odom_reliable_ = true;
            lio_state.SetVel(lightning::Vec3d::Zero());
            LocalizationLockingTestPeer::ObserveLio(static_localization, lio_state);

            lightning::loc::LocalizationResult lidar_loc;
            lidar_loc.timestamp_ = timestamp;
            // Raw LidarLoc results set lidar_loc_valid_; PGO sets valid_ only
            // after this observation point.
            lidar_loc.valid_ = false;
            lidar_loc.lidar_loc_valid_ = true;
            LocalizationLockingTestPeer::ObserveLidarLoc(static_localization, lidar_loc);
        }
        auto imu = std::make_shared<lightning::IMU>();
        imu->timestamp = timestamp;
        imu->angular_velocity = lightning::Vec3d(0.002, -0.001, 0.003);
        imu->linear_acceleration = lightning::Vec3d(0.0, 0.0, 1.0);
        static_active = LocalizationLockingTestPeer::UpdateImuStaticState(static_localization, imu);
    }
    if (!static_active) {
        std::cerr << "stable IMU with fresh zero-speed lidar evidence did not enter static hold" << std::endl;
        return 1;
    }
    static_localization.ProcessWheelSpeed(101.21, 0.2);
    for (int index = 1; index <= 3; ++index) {
        auto imu = std::make_shared<lightning::IMU>();
        imu->timestamp = 101.2 + 0.01 * index;
        imu->angular_velocity = lightning::Vec3d(0.002, -0.001, 0.003);
        imu->linear_acceleration = lightning::Vec3d(0.0, 0.0, 1.0);
        static_active = LocalizationLockingTestPeer::UpdateImuStaticState(static_localization, imu);
    }
    if (static_active) {
        std::cerr << "fresh moving CAN speed did not release static hold" << std::endl;
        return 1;
    }

    Localization can_static_localization;
    LocalizationLockingTestPeer::EnableImuStaticHold(can_static_localization);
    for (int index = 0; index <= 120; ++index) {
        const double timestamp = 200.0 + 0.01 * index;
        if (index % 10 == 0) {
            lightning::NavState lio_state;
            lio_state.timestamp_ = timestamp;
            lio_state.pose_is_ok_ = true;
            lio_state.lidar_odom_reliable_ = true;
            lio_state.SetVel(lightning::Vec3d::Zero());
            LocalizationLockingTestPeer::ObserveLio(can_static_localization, lio_state);

            lightning::loc::LocalizationResult lidar_loc;
            lidar_loc.timestamp_ = timestamp;
            lidar_loc.lidar_loc_valid_ = true;
            LocalizationLockingTestPeer::ObserveLidarLoc(can_static_localization, lidar_loc);
        }
        if (index % 2 == 0) {
            // CAN is sourced from the other Orin and its header clock is not
            // required to share the LiDAR/IMU epoch.
            can_static_localization.ProcessWheelSpeed(timestamp - 35.0, 0.0);
        }
        auto imu = std::make_shared<lightning::IMU>();
        imu->timestamp = timestamp;
        const double vibration = index % 2 == 0 ? 0.09 : 0.04;
        imu->angular_velocity = lightning::Vec3d(vibration, 0.01, -0.02);
        imu->linear_acceleration = lightning::Vec3d(0.0, 0.0,
                                                     index % 2 == 0 ? 1.08 : 0.94);
        static_active = LocalizationLockingTestPeer::UpdateImuStaticState(
            can_static_localization, imu);
    }
    if (!static_active) {
        std::cerr << "fresh zero CAN with low-speed lidar evidence did not override engine vibration"
                  << std::endl;
        return 1;
    }
    lightning::NavState vibrating_lio;
    vibrating_lio.timestamp_ = 201.21;
    vibrating_lio.pose_is_ok_ = true;
    vibrating_lio.lidar_odom_reliable_ = true;
    vibrating_lio.SetVel(lightning::Vec3d(0.18, 0.0, 0.0));
    LocalizationLockingTestPeer::ObserveLio(can_static_localization, vibrating_lio);
    for (int index = 1; index <= 3; ++index) {
        const double timestamp = 201.2 + 0.01 * index;
        can_static_localization.ProcessWheelSpeed(timestamp - 35.0, 0.0);
        auto imu = std::make_shared<lightning::IMU>();
        imu->timestamp = timestamp;
        imu->angular_velocity = lightning::Vec3d(0.09, 0.01, -0.02);
        imu->linear_acceleration = lightning::Vec3d(0.0, 0.0, 1.08);
        static_active = LocalizationLockingTestPeer::UpdateImuStaticState(
            can_static_localization, imu);
    }
    if (!static_active) {
        std::cerr << "fresh zero CAN did not suppress a stationary LIO vibration spike"
                  << std::endl;
        return 1;
    }
    return 0;
}
