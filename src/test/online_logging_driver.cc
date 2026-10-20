#include "utils/online_logging.h"
#include <glog/logging.h>
#include <chrono>
#include <cstdio>
#include <memory>
#include <string>
#include <thread>
#include <vector>

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    const bool synchronous = argc > 1 && std::string(argv[1]) == "sync";
    std::unique_ptr<lightning::OnlineLogging> logging;
    if (synchronous) FLAGS_logtostderr = true;
    else logging = std::make_unique<lightning::OnlineLogging>();
    const bool flood = argc > 1;
    const auto start = std::chrono::steady_clock::now();
    std::vector<std::thread> producers;
    for (int t = 0; t < 4; ++t) {
        producers.emplace_back([t, flood] {
            for (int i = 0; i < (flood ? 10000 : 100); ++i) {
                LOG(INFO) << "TEST_ROW producer=" << t << " item=" << i << " " << std::string(256, 'x');
            }
        });
    }
    for (auto& thread : producers) thread.join();
    LOG(WARNING) << "TEST_WARNING";
    LOG(ERROR) << "TEST_ERROR";
    google::FlushLogFiles(google::INFO);  // Must not turn into a storage wait.
    const auto ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start).count();
    // Independent pipe: proves producer/Flush completion while stderr is stalled.
    std::printf("PRODUCERS_DONE elapsed_ms=%.3f\n", ms);
    std::fflush(stdout);
}
