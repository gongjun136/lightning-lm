#pragma once

namespace lightning {
class AsyncOnlineLogger;
// Install after InitGoogleLogging/flag parsing, before any ROS/algorithm
// workers. Destroy after those workers stop. Fatal/crash logging remains glog's
// emergency synchronous path. Normal shutdown drains accepted logs and can wait
// for storage; live localization never performs the writer's I/O.
class OnlineLogging {
 public:
    OnlineLogging();
    ~OnlineLogging();
    OnlineLogging(const OnlineLogging&) = delete;
    OnlineLogging& operator=(const OnlineLogging&) = delete;
 private:
    AsyncOnlineLogger* logger_;  // Owned by glog; explicitly drained before shutdown.
};
}  // namespace lightning
