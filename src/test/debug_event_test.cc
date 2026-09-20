#include <cstdlib>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "common/debug_event.h"

namespace {

void Require(bool condition, const char* message) {
    if (condition) return;
    std::cerr << "debug_event_test failed: " << message << std::endl;
    std::exit(1);
}

}  // namespace

int main() {
    std::mutex messages_mutex;
    std::vector<std::string> messages;
    lightning::debug_event::SetSink([&](const std::string& message) {
        std::lock_guard<std::mutex> lock(messages_mutex);
        messages.push_back(message);
    });

    lightning::debug_event::Emit("ready");
    lightning::debug_event::EmitThrottled(
        "noisy", "first", std::chrono::hours(1));
    lightning::debug_event::EmitThrottled(
        "noisy", "suppressed", std::chrono::hours(1));

    lightning::debug_event::ReportState(
        "health", false, "failed", "recovered", std::chrono::milliseconds(0));
    lightning::debug_event::ReportState(
        "health", true, "failed", "recovered", std::chrono::milliseconds(0));
    lightning::debug_event::ReportState(
        "health", true, "failed again", "recovered", std::chrono::milliseconds(0));
    lightning::debug_event::ReportState(
        "health", false, "failed", "recovered", std::chrono::milliseconds(0));

    std::vector<std::thread> workers;
    for (int i = 0; i < 8; ++i) {
        workers.emplace_back([]() {
            lightning::debug_event::EmitThrottled(
                "concurrent", "one concurrent event", std::chrono::hours(1));
        });
    }
    for (auto& worker : workers) worker.join();

    {
        std::lock_guard<std::mutex> lock(messages_mutex);
        Require(messages.size() == 5, "unexpected event count");
        Require(messages[0] == "ready", "immediate event was not delivered");
        Require(messages[1] == "first", "throttled event did not preserve the first value");
        Require(messages[2] == "failed", "active transition was not delivered");
        Require(messages[3] == "recovered", "recovery transition was not delivered");
        Require(messages[4] == "one concurrent event", "concurrent throttling emitted more than once");
    }

    lightning::debug_event::ClearSink();
    lightning::debug_event::Emit("ignored");
    {
        std::lock_guard<std::mutex> lock(messages_mutex);
        Require(messages.size() == 5, "ClearSink did not detach the sink");
    }

    std::cout << "debug_event_test passed" << std::endl;
    return 0;
}
