#pragma once

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>
#include <utility>

namespace lightning {

/**
 * @brief 多生产线程、单消费线程的有界异步批处理队列，供因果追踪和在线日志共用。
 * @details 生产者通过 TryPush() 将记录加入 pending_；后台线程周期性交换 pending_ 与
 * batch_，解锁后将整个批次交给 consumer_。生产者不等待队列锁、消费者或磁盘写入，
 * 入队锁竞争、达到容量阈值或正在停止时直接丢弃本次记录。
 * @tparam T 必须具有可赋值的 std::uint64_t sequence 成员；序号按成功入队顺序分配，
 * 不表示不同线程中业务事件发生的先后顺序。
 * @note 每个实例仅支持一轮 Start()/Stop()。先 Start() 再启动生产者；所有生产者结束后，
 * 由消费线程之外的所有者调用 Stop()，生命周期操作由所有者串行执行。
 */
template <typename T>
class BoundedAsyncQueue {
 public:
    // 在后台线程同步执行的批处理回调，可能收到空批次；批次内容仅在回调期间有效。
    // 消费者负责实际处理记录，并保证异常不会逸出回调。
    using Consumer = std::function<void(const std::vector<T>&)>;

    /**
     * @brief 配置容量和消费者，并为两个缓冲区预留空间；此时不启动线程。
     * @param capacity pending_ 的容量上限；后台正在处理的 batch_ 单独占用一组缓冲。
     * @param consumer 后台线程调用的批处理函数。
     * @param critical_reserve 为关键记录预留的槽位数，普通记录只能使用剩余容量。
     * @pre capacity > 0、critical_reserve <= capacity 且 consumer 有效，由调用者保证。
     */
    BoundedAsyncQueue(std::size_t capacity, Consumer consumer,
                      std::size_t critical_reserve = 0)
        : capacity_(capacity), reserve_(critical_reserve), consumer_(std::move(consumer)) {
        pending_.reserve(capacity_);
        batch_.reserve(capacity_);
    }
    // 析构时等待消费线程退出；调用者须事先结束所有生产者。
    ~BoundedAsyncQueue() { Stop(); }

    // 启动唯一的消费线程；同一实例只能调用一次。
    void Start() { worker_ = std::thread([this] { Run(); }); }

    /**
     * @brief 尝试接收一条记录，成功返回 true，失败则累计丢弃数并返回 false。
     * @param row 待接收记录；成功入队时会覆盖其 sequence 字段。
     * @param critical 关键记录可使用全部容量；普通记录在 capacity_ - reserve_ 处限流。
     * @note 预留容量只影响接收阈值，关键记录仍按入队顺序消费，也可能因锁竞争而丢弃。
     */
    bool TryPush(T row, bool critical = false) {
        // received_ 统计所有尝试，包括随后因锁竞争、停止或容量限制而丢弃的记录。
        ++received_;
        std::unique_lock<std::mutex> lock(mutex_, std::try_to_lock);
        if (!lock.owns_lock() || stopping_ ||
            pending_.size() >= (critical ? capacity_ : capacity_ - reserve_)) {
            ++dropped_;
            if (critical) ++critical_dropped_;
            return false;
        }
        // 序号分配与插入共用同一把锁，使序号顺序与成功入队顺序一致。
        row.sequence = next_sequence_++;
        pending_.push_back(std::move(row));
        return true;
    }

    // 标记停止接收，唤醒消费者交付剩余批次，再等待线程退出。
    // 已结束后重复调用会直接返回；此函数不得由消费者回调自身调用。
    void Stop() {
        if (!worker_.joinable()) return;
        { std::lock_guard<std::mutex> lock(mutex_); stopping_ = true; }
        wake_.notify_one();
        worker_.join();
    }
    /// @brief 累计入队尝试次数，包含成功和失败。
    std::uint64_t Received() const { return received_.load(); }
    /// @brief 累计入队失败次数，包含普通记录和关键记录。
    std::uint64_t Dropped() const { return dropped_.load(); }
    /// @brief 关键记录的入队失败次数，是 Dropped() 的子集。
    std::uint64_t CriticalDropped() const { return critical_dropped_.load(); }

 private:
    // 消费线程主循环：等待、交换缓冲、解锁消费，停止时完成最后一批后退出。
    void Run() {
        for (;;) {
            bool done;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                // 正常运行按 20 ms 等待周期批处理；TryPush() 不通知，Stop() 提前唤醒。
                // 消费回调和线程调度也会影响实际处理间隔。
                wake_.wait_for(lock, std::chrono::milliseconds(20), [this] { return stopping_; });
                // 交换两个 vector 的存储，无需逐条搬运；生产者随后可填充新的 pending_。
                pending_.swap(batch_);
                done = stopping_;
            }
            // 回调在锁外执行，格式化或 I/O 不占用队列锁。
            consumer_(batch_);
            // 清除本批记录并保留容量，供下次交换后继续接收。
            batch_.clear();
            if (done) return;
        }
    }

    const std::size_t capacity_, reserve_;  ///< 待接收缓冲区总容量，以及关键记录的预留容量。
    Consumer consumer_;  ///< 由唯一消费线程执行的批处理回调。
    std::vector<T> pending_, batch_;  ///< 双缓冲：生产者追加 pending_，消费者处理交换得到的 batch_。
    std::mutex mutex_;  ///< 保护入队、双缓冲交换、序号分配和停止标志。
    std::condition_variable wake_;  ///< 控制周期等待，并在 Stop() 时提前唤醒消费者。
    std::thread worker_;  ///< 后台消费线程，由 Start() 创建、Stop() 等待退出。
    // 依次统计接收尝试、全部入队失败、关键记录入队失败；支持多生产线程并发更新。
    std::atomic<std::uint64_t> received_{0}, dropped_{0}, critical_dropped_{0};
    std::uint64_t next_sequence_ = 0;  ///< 下一个成功入队记录的序号，在队列锁下递增。
    bool stopping_ = false;  ///< 在队列锁下访问；置位后拒绝新记录并结束消费循环。
};
}  // namespace lightning
