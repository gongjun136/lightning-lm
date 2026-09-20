#ifndef LIGHTNING_FUNCTIONAL_SAFETY_HEARTBEAT_H
#define LIGHTNING_FUNCTIONAL_SAFETY_HEARTBEAT_H

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <mutex>
#include <string>
#include <thread>

#include <diagnostic_monitor_interfaces/msg/node_heartbeat.hpp>
#include <rclcpp/rclcpp.hpp>

namespace lightning::functional_safety {

enum class NodeState : std::uint8_t {
    kNotReady = 0,  // Localization initialization is in progress.
    kIdle = 1,      // Initialized, but no valid pose has been published.
    kRunning = 2,   // Localization is healthy and publishing poses.
    kDegraded = 3,  // P1: localization is degraded and DR remains available.
    kFault = 4,     // P0: localization is lost.
};

struct HeartbeatConfig {
    std::uint16_t node_id = 0;
    std::string topic;
    std::chrono::milliseconds period{100};
    std::size_t qos_depth = 1;
    bool reliable = false;
};

class HeartbeatPublisher {
   public:
    HeartbeatPublisher(rclcpp::Node& node, HeartbeatConfig config);
    ~HeartbeatPublisher();

    HeartbeatPublisher(const HeartbeatPublisher&) = delete;
    HeartbeatPublisher& operator=(const HeartbeatPublisher&) = delete;

    void SetState(NodeState state);
    void RecordWork();
    void PublishNow();
    void Stop();

    std::uint16_t NodeId() const { return config_.node_id; }
    std::uint64_t BootId() const { return boot_id_; }
    std::uint32_t WorkSequence() const { return work_seq_.load(std::memory_order_relaxed); }
    NodeState State() const {
        return static_cast<NodeState>(state_.load(std::memory_order_relaxed));
    }

    static std::uint64_t GenerateBootId();
    static std::uint64_t UnixMicrosecondsNow();
    static std::uint32_t Advance(std::atomic<std::uint32_t>& sequence);

   private:
    void Run();

    HeartbeatConfig config_;
    rclcpp::Publisher<diagnostic_monitor_interfaces::msg::NodeHeartbeat>::SharedPtr publisher_;
    const std::uint64_t boot_id_;
    std::atomic<std::uint32_t> heartbeat_seq_{0};
    std::atomic<std::uint32_t> work_seq_{0};
    std::atomic<std::uint8_t> state_{static_cast<std::uint8_t>(NodeState::kNotReady)};
    std::atomic_bool stop_requested_{false};
    std::mutex wait_mutex_;
    std::condition_variable wait_condition_;
    std::thread worker_;
};

}  // namespace lightning::functional_safety

#endif  // LIGHTNING_FUNCTIONAL_SAFETY_HEARTBEAT_H
