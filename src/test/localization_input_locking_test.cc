#include "core/localization/localization.h"

#include <chrono>
#include <future>
#include <iostream>
#include <mutex>
#include <limits>

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
    static std::size_t WheelHistorySize(const Localization& localization) {
        return localization.static_wheel_speed_history_.size();
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
            can_static_localization.ProcessWheelSpeed(timestamp, 0.0);
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
        can_static_localization.ProcessWheelSpeed(timestamp, 0.0);
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

    Localization offset_can_localization;
    LocalizationLockingTestPeer::EnableImuStaticHold(offset_can_localization);
    for (int index = 0; index <= 120; ++index) {
        const double timestamp = 300.0 + 0.01 * index;
        if (index % 10 == 0) {
            lightning::NavState lio_state;
            lio_state.timestamp_ = timestamp;
            lio_state.pose_is_ok_ = true;
            lio_state.lidar_odom_reliable_ = true;
            lio_state.SetVel(lightning::Vec3d::Zero());
            LocalizationLockingTestPeer::ObserveLio(offset_can_localization, lio_state);

            lightning::loc::LocalizationResult lidar_loc;
            lidar_loc.timestamp_ = timestamp;
            lidar_loc.lidar_loc_valid_ = true;
            LocalizationLockingTestPeer::ObserveLidarLoc(offset_can_localization, lidar_loc);
        }
        if (index % 2 == 0) {
            offset_can_localization.ProcessWheelSpeed(timestamp - 37.0, 0.0);
        }
        auto imu = std::make_shared<lightning::IMU>();
        imu->timestamp = timestamp;
        const double vibration = index % 2 == 0 ? 0.09 : 0.04;
        imu->angular_velocity = lightning::Vec3d(vibration, 0.01, -0.02);
        imu->linear_acceleration = lightning::Vec3d(
            0.0, 0.0, index % 2 == 0 ? 1.08 : 0.94);
        static_active = LocalizationLockingTestPeer::UpdateImuStaticState(
            offset_can_localization, imu);
    }
    if (static_active) {
        std::cerr << "37-second-offset CAN incorrectly overrode the IMU fallback"
                  << std::endl;
        return 1;
    }
    offset_can_localization.ProcessWheelSpeed(301.21, 0.0);
    auto aligned_imu = std::make_shared<lightning::IMU>();
    aligned_imu->timestamp = 301.21;
    aligned_imu->angular_velocity = lightning::Vec3d(0.09, 0.01, -0.02);
    aligned_imu->linear_acceleration = lightning::Vec3d(0.0, 0.0, 1.08);
    static_active = LocalizationLockingTestPeer::UpdateImuStaticState(
        offset_can_localization, aligned_imu);
    if (!static_active) {
        std::cerr << "timestamp-aligned zero CAN did not recover the stationary observation"
                  << std::endl;
        return 1;
    }

    // CAN callbacks can run ahead of the serialized IMU worker. A future
    // sample must neither hide a valid historical zero nor release the hold.
    Localization queued_can_localization;
    LocalizationLockingTestPeer::EnableImuStaticHold(queued_can_localization);
    for (int index = 0; index <= 150; ++index) {
        queued_can_localization.ProcessWheelSpeed(400.0 + 0.01 * index,
                                                  index < 125 ? 0.0 : 0.5);
    }
    const auto observe_queued_imu = [&](double stamp) {
        lightning::NavState lio;
        lio.timestamp_ = stamp;
        lio.lidar_odom_reliable_ = true;
        LocalizationLockingTestPeer::ObserveLio(queued_can_localization, lio);
        lightning::loc::LocalizationResult loc;
        loc.timestamp_ = stamp;
        loc.lidar_loc_valid_ = true;
        LocalizationLockingTestPeer::ObserveLidarLoc(queued_can_localization, loc);
        auto imu = std::make_shared<lightning::IMU>();
        imu->timestamp = stamp;
        imu->angular_velocity = lightning::Vec3d(0.09, 0.01, -0.02);
        imu->linear_acceleration = lightning::Vec3d(0.0, 0.0, 9.81);
        return LocalizationLockingTestPeer::UpdateImuStaticState(queued_can_localization, imu);
    };
    for (int index = 0; index <= 120; ++index) {
        static_active = observe_queued_imu(400.0 + 0.01 * index);
    }
    if (!static_active) {
        std::cerr << "future CAN hid the valid stationary history" << std::endl;
        return 1;
    }
    for (int index = 121; index <= 124; ++index) {
        if (!observe_queued_imu(400.0 + 0.01 * index)) {
            std::cerr << "future moving CAN released the hold before its epoch" << std::endl;
            return 1;
        }
    }
    for (int index = 125; index <= 127; ++index) {
        static_active = observe_queued_imu(400.0 + 0.01 * index);
    }
    if (static_active) {
        std::cerr << "historical moving CAN did not release hold in three samples" << std::endl;
        return 1;
    }
    Localization future_zero_localization;
    LocalizationLockingTestPeer::EnableImuStaticHold(future_zero_localization);
    future_zero_localization.ProcessWheelSpeed(500.0, 0.5);
    future_zero_localization.ProcessWheelSpeed(501.1, 0.0);
    for (int index = 0; index <= 95; ++index) {
        const double stamp = 500.0 + 0.01 * index;
        lightning::NavState lio;
        lio.timestamp_ = stamp;
        lio.lidar_odom_reliable_ = true;
        LocalizationLockingTestPeer::ObserveLio(future_zero_localization, lio);
        lightning::loc::LocalizationResult loc;
        loc.timestamp_ = stamp;
        loc.lidar_loc_valid_ = true;
        LocalizationLockingTestPeer::ObserveLidarLoc(future_zero_localization, loc);
        auto imu = std::make_shared<lightning::IMU>();
        imu->timestamp = stamp;
        imu->angular_velocity = lightning::Vec3d(0.09, 0.0, 0.0);
        imu->linear_acceleration = lightning::Vec3d(0.0, 0.0, 9.81);
        if (LocalizationLockingTestPeer::UpdateImuStaticState(future_zero_localization, imu)) {
            std::cerr << "future zero CAN incorrectly engaged static hold" << std::endl;
            return 1;
        }
    }
    Localization bounded_can;
    LocalizationLockingTestPeer::EnableImuStaticHold(bounded_can);
    bounded_can.ProcessWheelSpeed(600.0, 0.0);
    bounded_can.ProcessWheelSpeed(600.0, 0.5);
    bounded_can.ProcessWheelSpeed(599.9, 0.5);
    bounded_can.ProcessWheelSpeed(std::numeric_limits<double>::quiet_NaN(), 0.0);
    bounded_can.ProcessWheelSpeed(600.1, std::numeric_limits<double>::infinity());
    if (LocalizationLockingTestPeer::WheelHistorySize(bounded_can) != 1) {
        std::cerr << "invalid or nonmonotonic CAN changed history" << std::endl;
        return 1;
    }
    for (int i = 1; i <= 1000; ++i) bounded_can.ProcessWheelSpeed(600.0 + i * 0.001, 0.0);
    if (LocalizationLockingTestPeer::WheelHistorySize(bounded_can) > 512) {
        std::cerr << "CAN history sample bound exceeded" << std::endl;
        return 1;
    }
    bounded_can.ProcessWheelSpeed(604.0, 0.0);
    if (LocalizationLockingTestPeer::WheelHistorySize(bounded_can) != 1) {
        std::cerr << "CAN history age bound exceeded" << std::endl;
        return 1;
    }
    return 0;
}
