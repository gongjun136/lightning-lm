#pragma once

namespace lightning {
class AsyncOnlineLogger;
// Install after InitGoogleLogging/flag parsing, before any ROS/algorithm
// workers. Destroy after those workers stop. Fatal/crash logging remains glog's
// emergency synchronous path. Normal shutdown drains accepted logs and can wait
// for storage; live localization never performs the writer's I/O.
// LIGHTNING_LM_COMPUTE_PROFILE=1 routes detailed timing into an independent
// bounded writer. Its new file is LIGHTNING_LM_COMPUTE_PROFILE_PATH, defaulting
// to ./compute_profile.log for direct binary invocations; existing files fail
// startup rather than being overwritten. End-to-end health stays on stderr.
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
