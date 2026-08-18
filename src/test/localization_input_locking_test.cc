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
    if (!LocalizationLockingTestPeer::ShouldThrottleLidar(guarded_localization, 10.0)) {
        std::cerr << "overload guard did not throttle above the maximum lag" << std::endl;
        return 1;
    }
    LocalizationLockingTestPeer::SetSensorLag(guarded_localization, 10.0, 9.7);
    if (!LocalizationLockingTestPeer::ShouldThrottleLidar(guarded_localization, 10.1)) {
        std::cerr << "overload guard did not preserve hysteresis" << std::endl;
        return 1;
    }
    LocalizationLockingTestPeer::SetSensorLag(guarded_localization, 10.0, 9.95);
    if (LocalizationLockingTestPeer::ShouldThrottleLidar(guarded_localization, 10.2)) {
        std::cerr << "overload guard did not resume below the recovery lag" << std::endl;
        return 1;
    }
    if (guarded_localization.GetRuntimeStats().sensor_queue_dropped != 2) {
        std::cerr << "overload guard drops were not exposed in runtime diagnostics" << std::endl;
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
    return 0;
}
