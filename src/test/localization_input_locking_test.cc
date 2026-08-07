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
    return 0;
}
