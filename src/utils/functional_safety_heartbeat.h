#ifndef LIGHTNING_FUNCTIONAL_SAFETY_HEARTBEAT_H
#define LIGHTNING_FUNCTIONAL_SAFETY_HEARTBEAT_H

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <mutex>
#include <string>
#include <thread>
#include <set>
#include <vector>

#include <diagnostic_monitor_interfaces/msg/node_heartbeat.hpp>
#include <rclcpp/rclcpp.hpp>

namespace lightning::functional_safety {

enum class NodeState : std::uint8_t {
    kNotReady = 0,  // Localization initialization is in progress.
    kIdle = 1,      // Initialized, but no valid pose has been published.
    kRunning = 2,   // Caller-reported progress state, not a fresh diagnostic PASS.
    kDegraded = 3,  // P1: localization is degraded and DR remains available.
    kFault = 4,     // P0: localization is lost.
};

struct HeartbeatConfig {
    std::uint16_t node_id = 0;  ///< 监测端区分节点的标识。
    std::string topic;  ///< 心跳消息发布话题。
    std::chrono::milliseconds period{100};  ///< 心跳发布周期，使用毫秒表示。
    std::size_t qos_depth = 1;  ///< 心跳发布器的 QoS 队列深度。
    bool reliable = false;  ///< 是否采用可靠传输 QoS。
};

class HeartbeatPublisher {
   public:
    HeartbeatPublisher(rclcpp::Node& node, HeartbeatConfig config);
    ~HeartbeatPublisher();

    HeartbeatPublisher(const HeartbeatPublisher&) = delete;
    HeartbeatPublisher& operator=(const HeartbeatPublisher&) = delete;

    void SetState(NodeState state);
    // Apply one evaluated batch atomically. UNKNOWN/omitted codes retain state.
    // Callers supply only owned, non-historical codes with fresh evidence.
    void UpdateFaults(const std::vector<std::uint32_t>& failed,
                      const std::vector<std::uint32_t>& passed);
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

    HeartbeatConfig config_;  ///< 心跳标识、发布周期与 QoS 配置。
    rclcpp::Publisher<diagnostic_monitor_interfaces::msg::NodeHeartbeat>::SharedPtr publisher_;  ///< 发布 NodeHeartbeat 消息的 ROS 发布器。
    const std::uint64_t boot_id_;  ///< 本次进程启动的唯一标识。
    std::atomic<std::uint32_t> heartbeat_seq_{0};  ///< 心跳发布序列号。
    std::atomic<std::uint32_t> work_seq_{0};  ///< 有效定位工作完成后推进的工作序列号。
    std::atomic<std::uint8_t> state_{static_cast<std::uint8_t>(NodeState::kNotReady)};  ///< 当前节点运行状态，供监测端识别未就绪、运行或故障。
    std::atomic_bool stop_requested_{false};  ///< 请求心跳线程退出的标志。
    std::mutex fault_mutex_;
    std::set<std::uint32_t> active_faults_;
    std::mutex publish_mutex_;
    std::mutex wait_mutex_;  ///< 配合条件变量等待发布周期或停止请求。
    std::condition_variable wait_condition_;  ///< 心跳周期等待及停止时唤醒使用的条件变量。
    std::thread worker_;  ///< 独立于算法处理线程的心跳发布线程。
};

}  // namespace lightning::functional_safety

#endif  // LIGHTNING_FUNCTIONAL_SAFETY_HEARTBEAT_H
