#include "utils/online_logging.h"
#include <glog/logging.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
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
    const std::string mode = argc > 1 ? argv[1] : "";
    if (mode.rfind("profile", 0) == 0) {
        const bool many = mode == "profile_flood";
        for (int i = 0; i < (many ? 10000 : 25); ++i) {
            for (const char* module : {"multi_lidar_assemble", "lio_pipeline", "lidar_loc",
                                       "lidar_loc_pipeline", "cloud_publish"}) {
                LOG(INFO) << "COMPUTE_BENCH_FRAME module=" << module << " frame_id=" << i
                          << " outer_wall_ms=1.25 outer_thread_cpu_ms=0.75 outer_process_cpu_ms_concurrent=2.5";
            }
            LOG(INFO) << "LIO_BENCH_FRAME frame_id=" << i << " scan_match_ms=0.75";
            LOG(INFO) << "COMPUTE_BENCH_EVENT module=pgo_optimization frame_id=" << i << " outer_wall_ms=0.5";
            LOG(INFO) << "COMPUTE_BENCH_SUMMARY module=imu_dr_pipeline samples=100 outer_p99_ms=0.5";
            if (!many) std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        if (mode == "profile_long")
            LOG(INFO) << "COMPUTE_BENCH_FRAME module=long_record " << std::string(10000, 'x');
        LOG(INFO) << "COMPUTE_BENCH_SUMMARY module=lidar_end_to_end phase=tracking samples=100 p99_ms=5";
        LOG(WARNING) << "TEST_WARNING";
        LOG(ERROR) << "TEST_ERROR";
        std::printf("PRODUCERS_DONE\n");
        std::fflush(stdout);
        if (mode == "profile_crash") {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            std::_Exit(77);
        }
        return 0;
    }
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
