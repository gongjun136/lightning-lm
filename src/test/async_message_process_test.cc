#include <atomic>
#include <chrono>
#include <condition_variable>
#include <iostream>
#include <mutex>
#include <thread>

#include "core/system/async_message_process.h"

int main() {
    using namespace std::chrono_literals;

    lightning::sys::AsyncMessageProcess<int> processor;
    std::atomic<int> processed{0};
    std::mutex started_mutex;
    std::condition_variable started_cv;
    bool first_started = false;

    processor.SetMaxSize(100);
    processor.SetProcFunc([&](int value) {
        if (value == 0) {
            {
                std::lock_guard<std::mutex> lock(started_mutex);
                first_started = true;
            }
            started_cv.notify_one();
            std::this_thread::sleep_for(50ms);
        }
        ++processed;
    });
    processor.Start();
    processor.AddMessage(0);

    {
        std::unique_lock<std::mutex> lock(started_mutex);
        if (!started_cv.wait_for(lock, 2s, [&] { return first_started; })) {
            std::cerr << "worker did not start" << std::endl;
            processor.Quit();
            return 1;
        }
    }

    constexpr int kMessages = 50;
    for (int value = 1; value < kMessages; ++value) processor.AddMessage(value);
    processor.Quit();

    if (processed != kMessages) {
        std::cerr << "Quit dropped pending messages: processed=" << processed << ", expected=" << kMessages
                  << std::endl;
        return 1;
    }
    return 0;
}
