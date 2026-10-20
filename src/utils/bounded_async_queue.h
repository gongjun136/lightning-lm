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

// The causal-trace double-buffer design, shared with online text logging.
// T has a uint64_t sequence member. Sequence denotes admission order, not
// cross-thread event time. Producers never wait for the consumer or storage.
// Start before producers; Stop only after producers have been joined.
template <typename T>
class BoundedAsyncQueue {
 public:
    using Consumer = std::function<void(const std::vector<T>&)>;
    BoundedAsyncQueue(std::size_t capacity, Consumer consumer,
                      std::size_t critical_reserve = 0)
        : capacity_(capacity), reserve_(critical_reserve), consumer_(std::move(consumer)) {
        pending_.reserve(capacity_);
        batch_.reserve(capacity_);
    }
    ~BoundedAsyncQueue() { Stop(); }
    void Start() { worker_ = std::thread([this] { Run(); }); }
    bool TryPush(T row, bool critical = false) {
        ++received_;
        std::unique_lock<std::mutex> lock(mutex_, std::try_to_lock);
        if (!lock.owns_lock() || stopping_ ||
            pending_.size() >= (critical ? capacity_ : capacity_ - reserve_)) {
            ++dropped_;
            if (critical) ++critical_dropped_;
            return false;
        }
        // Allocate the sequence while holding the same lock as insertion.
        row.sequence = next_sequence_++;
        pending_.push_back(std::move(row));
        return true;
    }
    void Stop() {
        if (!worker_.joinable()) return;
        { std::lock_guard<std::mutex> lock(mutex_); stopping_ = true; }
        wake_.notify_one();
        worker_.join();
    }
    std::uint64_t Received() const { return received_.load(); }
    std::uint64_t Dropped() const { return dropped_.load(); }
    std::uint64_t CriticalDropped() const { return critical_dropped_.load(); }

 private:
    void Run() {
        for (;;) {
            bool done;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                wake_.wait_for(lock, std::chrono::milliseconds(20), [this] { return stopping_; });
                pending_.swap(batch_);
                done = stopping_;
            }
            // No queue lock is held across formatting, writes, or flushes.
            consumer_(batch_);
            batch_.clear();
            if (done) return;
        }
    }
    const std::size_t capacity_, reserve_;
    Consumer consumer_;
    std::vector<T> pending_, batch_;
    std::mutex mutex_;
    std::condition_variable wake_;
    std::thread worker_;
    std::atomic<std::uint64_t> received_{0}, dropped_{0}, critical_dropped_{0};
    std::uint64_t next_sequence_ = 0;
    bool stopping_ = false;
};
}  // namespace lightning
