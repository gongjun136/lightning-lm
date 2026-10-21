#include "utils/functional_safety_heartbeat.h"
#include "utils/localization_dtc_reporting.h"
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
    // Exercise the actual localization detector -> production adapter -> ROS
    // heartbeat path, rather than supplying fault codes only to the serializer.
    using namespace lightning;
    using namespace lightning::loc::dtc;
    Monitor monitor(Config{}, Clock::time_point{});
    auto feed = [&](int ms, double x, bool large_covariance = false) {
        const auto now = Clock::time_point{} + std::chrono::milliseconds(ms);
        const double stamp = 100.0 + ms * 0.001;
        const SE3 pose(SO3(), Vec3d(x, 0, 0));
        Mat6d covariance = 0.0001 * Mat6d::Identity();
        if (large_covariance) covariance(1, 1) = 0.31 * 0.31;
        monitor.ObserveSpeed(stamp, 1.0, now);
        monitor.ObserveGyro(stamp, 0.0, now);
        monitor.ObservePose(stamp, pose, pose, covariance, true, false, 1, now);
        functional_safety::UpdateLocalizationFaults(publisher, monitor.Evaluate(false, now));
    };
    functional_safety::UpdateLocalizationFaults(publisher, monitor.Evaluate(false, Clock::time_point{}));
    check({}, 0);
    for (int ms = 0; ms <= 1000; ms += 20) feed(ms, ms * 0.001);
    check({}, 0);
    feed(1020, 3.02, true);
    check({0x30811813, 0x40450198}, 0);
    functional_safety::UpdateLocalizationFaults(publisher,
        monitor.Evaluate(false, Clock::time_point{} + std::chrono::milliseconds(1400)));
    check({0x30811813, 0x40450198}, 0); // UNKNOWN/stale cannot clear an established fault.
    feed(1420, 3.42);
    check({0x30811813}, 0); // Fresh uncertainty PASS clears only its own code.
    for (int ms = 1440; ms <= 2020; ms += 20) feed(ms, 2.0 + ms * 0.001);
    check({}, 0);
    // A yaw fault can recover without optional CAN, but a translation fault
    // must remain active until translation itself is verified again.
    Monitor optional_monitor(Config{}, Clock::time_point{});
    auto optional_feed = [&](int ms, double x, double yaw, bool can) {
        const auto now = Clock::time_point{} + std::chrono::milliseconds(ms);
        const double stamp = 100.0 + ms * 0.001;
        const SE3 pose(SO3::exp(Vec3d(0, 0, yaw)), Vec3d(x, 0, 0));
        if (can) optional_monitor.ObserveSpeed(stamp, 0.0, now);
        optional_monitor.ObserveGyro(stamp, 0.0, now);
        optional_monitor.ObservePose(stamp, pose, pose, 0.0001 * Mat6d::Identity(), true, false, 1, now);
        functional_safety::UpdateLocalizationFaults(publisher, optional_monitor.Evaluate(false, now));
    };
    for (int ms = 0; ms <= 1000; ms += 20) optional_feed(ms, 0, 0, false);
    check({}, 0);
    optional_feed(1020, 0, .2, false);
    check({0x30811813}, 0);
    for (int ms = 1040; ms <= 1800; ms += 20) optional_feed(ms, 0, 0, false);
    check({}, 0);
    for (int ms = 1820; ms <= 2500; ms += 20) optional_feed(ms, 0, 0, true);
    optional_feed(2520, 2, 0, true);
    check({0x30811813}, 0);
    for (int ms = 2540; ms <= 3300; ms += 20) optional_feed(ms, 2, 0, false);
    check({0x30811813}, 0);
    for (int ms = 3320; ms <= 4000; ms += 20) optional_feed(ms, 2, 0, true);
    check({}, 0);
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
