#include <cstdlib>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include <glog/logging.h>

#include "common/debug_event.h"
#include "core/localization/pose_graph/smoother.h"

namespace {

void Require(bool condition, const char* message) {
    if (condition) return;
    std::cerr << "debug_event_test failed: " << message << std::endl;
    std::exit(1);
}

}  // namespace

int main() {
    google::InitGoogleLogging("debug_event_test");
    std::mutex messages_mutex;
    std::vector<std::string> messages;
    lightning::debug_event::SetSink([&](const std::string& message) {
        std::lock_guard<std::mutex> lock(messages_mutex);
        messages.push_back(message);
    });

    lightning::debug_event::Emit("在线定位节点已就绪");
    lightning::debug_event::EmitThrottled(
        "noisy", "首次异常 count=1", std::chrono::hours(1));
    lightning::debug_event::EmitThrottled(
        "noisy", "限频抑制的异常", std::chrono::hours(1));

    lightning::debug_event::ReportState(
        "health", false, "定位已丢失", "定位已恢复", std::chrono::milliseconds(0));
    lightning::debug_event::ReportState(
        "health", true, "定位已丢失", "定位已恢复", std::chrono::milliseconds(0));
    lightning::debug_event::ReportState(
        "health", true, "定位仍然丢失", "定位已恢复", std::chrono::milliseconds(0));
    lightning::debug_event::ReportState(
        "health", false, "定位已丢失", "定位已恢复", std::chrono::milliseconds(0));

    std::vector<std::thread> workers;
    for (int i = 0; i < 8; ++i) {
        workers.emplace_back([]() {
            lightning::debug_event::EmitThrottled(
                "concurrent", "并发异常仅发送一次", std::chrono::hours(1));
        });
    }
    for (auto& worker : workers) worker.join();

    {
        std::lock_guard<std::mutex> lock(messages_mutex);
        Require(messages.size() == 5, "unexpected event count");
        Require(messages[0] == "在线定位节点已就绪", "Chinese immediate event was not preserved");
        Require(messages[1] == "首次异常 count=1", "Chinese throttled event or numeric field changed");
        Require(messages[2] == "定位已丢失", "Chinese active transition was not preserved");
        Require(messages[3] == "定位已恢复", "Chinese recovery transition was not preserved");
        Require(messages[4] == "并发异常仅发送一次", "Chinese concurrent event changed or was duplicated");
    }

    lightning::debug_event::ClearSink();
    lightning::debug_event::Emit("ignored");
    {
        std::lock_guard<std::mutex> lock(messages_mutex);
        Require(messages.size() == 5, "ClearSink did not detach the sink");
    }

    // Exercise actual production descriptions, not only dispatcher fixtures.
    messages.clear();
    lightning::debug_event::SetSink([&](const std::string& message) {
        messages.push_back(message);
    });
    lightning::loc::PoseSmoother invalid_metadata;
    Require(!invalid_metadata.PushDRPose(lightning::SE3(), -1.0, 0.0),
            "invalid DR metadata unexpectedly accepted");
    Require(!invalid_metadata.PushDRPose(lightning::SE3(), -2.0, 0.0),
            "repeated invalid DR metadata unexpectedly accepted");
    Require(messages.size() == 1, "production event throttling changed");
    Require(messages.back() == "位姿平滑器拒绝无效的 DR 元数据",
            "production invalid-metadata description is not Chinese");

    lightning::loc::PoseSmoother timestamp_rollback;
    Require(timestamp_rollback.PushDRPose(lightning::SE3(), 1.0, 0.0), "initial DR pose rejected");
    Require(!timestamp_rollback.PushDRPose(lightning::SE3(), 0.5, 0.0), "DR rollback accepted");
    Require(messages.size() == 2 && messages.back() == "位姿平滑器拒绝非递增的 DR 时间戳",
            "production timestamp-rollback description is not Chinese");

    lightning::loc::PoseSmoother motion_jump;
    Require(motion_jump.PushDRPose(lightning::SE3(), 1.0, 0.0), "initial motion pose rejected");
    const lightning::SE3 jumped_pose(lightning::SO3(), lightning::Vec3d(10.0, 0.0, 0.0));
    Require(!motion_jump.PushDRPose(jumped_pose, 1.1, 0.0), "DR motion jump accepted");
    Require(messages.size() == 3 &&
                messages.back().find("位姿平滑器拒绝运动突变: distance_m=10") == 0 &&
                messages.back().find(", allowed_distance_m=") != std::string::npos,
            "production dynamic description or numeric field names changed");
    lightning::debug_event::ClearSink();

    std::cout << "debug_event_test passed" << std::endl;
    return 0;
}
