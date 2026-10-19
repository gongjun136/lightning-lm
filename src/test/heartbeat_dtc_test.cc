#include "utils/functional_safety_heartbeat.h"
#include <stdexcept>
#include <iostream>

using lightning::functional_safety::HeartbeatConfig;
using lightning::functional_safety::HeartbeatPublisher;
using Message = diagnostic_monitor_interfaces::msg::NodeHeartbeat;

int main(int argc, char** argv) {
    rclcpp::init(argc, argv);
    auto node = std::make_shared<rclcpp::Node>("heartbeat_dtc_test");
    Message latest;
    size_t received = 0;
    auto sub = node->create_subscription<Message>("/test/heartbeat_dtc", rclcpp::QoS(10).best_effort(),
        [&](Message::ConstSharedPtr msg) { latest = *msg; ++received; });
    HeartbeatConfig config;
    config.node_id = 6;
    config.topic = "/test/heartbeat_dtc";
    HeartbeatPublisher publisher(*node, config);
    publisher.Stop(); // Deterministic manual publications through the production serializer.
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
    while (sub->get_publisher_count() == 0 && std::chrono::steady_clock::now() < deadline)
        rclcpp::spin_some(node);
    auto check = [&](std::vector<uint32_t> expected, uint32_t work) {
        rclcpp::spin_some(node);
        const auto before = received;
        const auto end = std::chrono::steady_clock::now() + std::chrono::seconds(5);
        // Graph discovery can precede DDS matching; retry the same complete snapshot.
        do {
            publisher.PublishNow();
            rclcpp::spin_some(node);
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        } while ((received == before || latest.active_fault != expected || latest.work_seq != work) &&
                 std::chrono::steady_clock::now() < end);
        if (received == before || latest.active_fault != expected || latest.work_seq != work ||
            latest.node_id != 6 || latest.boot_id != publisher.BootId()) {
            std::cerr << "received=" << received << " before=" << before
                      << " node=" << latest.node_id << " work=" << latest.work_seq
                      << " active_count=" << latest.active_fault.size()
                      << " expected_count=" << expected.size() << '\n';
            throw std::runtime_error("heartbeat fault snapshot/work identity mismatch");
        }
    };
    check({}, 0); // Startup UNKNOWN does not invent faults.
    publisher.UpdateFaults({0x40450198, 0x30811813, 0x30811813, 0}, {});
    check({0x30811813, 0x40450198}, 0); // Sorted, unique, no zero.
    const auto first_seq = latest.heartbeat_seq;
    publisher.UpdateFaults({}, {}); // UNKNOWN / missing observations.
    check({0x30811813, 0x40450198}, 0);
    if (latest.heartbeat_seq <= first_seq) throw std::runtime_error("heartbeat sequence stalled");
    publisher.UpdateFaults({}, {0x30811813});
    check({0x40450198}, 0); // Only the explicitly passing code recovers.
    publisher.RecordWork();
    publisher.UpdateFaults({}, {0x40450198});
    check({}, 1);
    publisher.UpdateFaults({0x30811813}, {0x30811813});
    check({0x30811813}, 1); // Conflicting batch cannot clear a failure.
    std::atomic_bool updating{true};
    std::thread updater([&]() {
        while (updating.load()) {
            publisher.UpdateFaults({0x40450198}, {0x30811813});
            publisher.UpdateFaults({0x30811813}, {0x40450198});
        }
    });
    bool coherent = true;
    for (int i = 0; i < 30; ++i) {
        const auto before = received;
        publisher.PublishNow();
        const auto end = std::chrono::steady_clock::now() + std::chrono::seconds(2);
        while (received == before && std::chrono::steady_clock::now() < end) rclcpp::spin_some(node);
        coherent &= received > before && latest.active_fault.size() == 1 &&
            (latest.active_fault[0] == 0x30811813 || latest.active_fault[0] == 0x40450198);
    }
    updating.store(false);
    updater.join();
    if (!coherent) throw std::runtime_error("concurrent update exposed a partial fault batch");
    HeartbeatConfig restart_config = config;
    restart_config.topic = "/test/restarted_heartbeat_dtc";
    HeartbeatPublisher restarted(*node, restart_config);
    restarted.Stop();
    if (restarted.BootId() == publisher.BootId() || restarted.BootId() == 0)
        throw std::runtime_error("restart must have distinct boot identity");
    rclcpp::shutdown();
    std::cout << "heartbeat_dtc_test passed\n";
}
