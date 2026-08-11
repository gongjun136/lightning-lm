//
// Created by xiang on 2022/2/9.
//

#ifndef ASYNC_MESSAGE_PROCESS_H
#define ASYNC_MESSAGE_PROCESS_H

#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>
#include <utility>

#include <glog/logging.h>

namespace lightning::sys {
using UL = std::unique_lock<std::mutex>;

/**
 * 异步消息处理类
 * 内部有线程和队列机制，保证回调是串行的
 * @tparam T
 *
 * NOTE skip设为1的时候实际不会跳帧。。设为2的时候实际走一帧跳一帧
 */
template <typename T>
class AsyncMessageProcess {
   public:
    using ProcFunc = std::function<void(const T&)>;  // 消息回调函数
    AsyncMessageProcess() = default;
    AsyncMessageProcess(ProcFunc proc_func, std::string name = "");

    /// 设置处理函数
    void SetProcFunc(ProcFunc proc_func) { custom_func_ = proc_func; }

    /// 设置队列最大长度
    void SetMaxSize(size_t size) { max_size_ = size; }

    /// 开始处理消息
    void Start();

    /// 添加一条消息
    void AddMessage(const T& msg);

    /// 退出
    void Quit();

    /// 清空跳帧计数器，下一个数据会立即执行
    void CleanSkipCnt();

    /// Number of accepted messages waiting or currently being processed.
    size_t PendingCount() const;

    /// Number of accepted messages evicted because the bounded queue was full.
    size_t DroppedCount() const { return dropped_count_.load(); }

    /// Number of messages whose callback completed.
    size_t ProcessedCount() const { return processed_count_.load(); }

    void SetName(std::string name) { name_ = std::move(name); }
    void SetSkipParam(bool enable_skip, int skip_num) { enable_skip_ = enable_skip, skip_num_ = skip_num; }

    AsyncMessageProcess(const AsyncMessageProcess&) = delete;
    void operator=(const AsyncMessageProcess&) = delete;

   private:
    void ProcLoop();

    std::thread proc_;
    mutable std::mutex mutex_;
    std::condition_variable cv_msg_;
    std::deque<T> msg_buffer_;
    bool update_flag_ = false;
    bool exit_flag_ = false;
    size_t max_size_ = 100;  // 40
    std::string name_;

    /// 跳帧
    bool enable_skip_ = false;
    int skip_num_ = 0;
    int skip_cnt_ = 0;

    ProcFunc custom_func_;
    std::atomic<size_t> in_flight_count_{0};
    std::atomic<size_t> dropped_count_{0};
    std::atomic<size_t> processed_count_{0};
};

template <typename T>
size_t AsyncMessageProcess<T>::PendingCount() const {
    UL lock(mutex_);
    return msg_buffer_.size() + in_flight_count_.load();
}

template <typename T>
void AsyncMessageProcess<T>::CleanSkipCnt() {
    UL lock(mutex_);
    skip_cnt_ = 0;
}

template <typename T>
AsyncMessageProcess<T>::AsyncMessageProcess(AsyncMessageProcess::ProcFunc proc_func, std::string name) {
    custom_func_ = std::move(proc_func);
    name_ = name;
}

template <typename T>
void AsyncMessageProcess<T>::Start() {
    {
        UL lock(mutex_);
        exit_flag_ = false;
        update_flag_ = false;
        in_flight_count_ = 0;
        dropped_count_ = 0;
        processed_count_ = 0;
    }
    proc_ = std::thread([this]() { ProcLoop(); });
}

template <typename T>
void AsyncMessageProcess<T>::ProcLoop() {
    while (true) {
        UL lock(mutex_);
        cv_msg_.wait(lock, [this]() { return update_flag_ || exit_flag_; });
        if (exit_flag_ && msg_buffer_.empty()) break;

        // Take every message currently available. If Quit() arrives while the
        // batch is being processed, the next loop drains messages accumulated
        // in the meantime before the worker exits.
        std::deque<T> buffer;
        buffer.swap(msg_buffer_);
        in_flight_count_ = buffer.size();
        update_flag_ = false;
        lock.unlock();

        // 处理之
        for (const auto& msg : buffer) {
            custom_func_(msg);
            --in_flight_count_;
            ++processed_count_;
        }
    }
}

template <typename T>
void AsyncMessageProcess<T>::AddMessage(const T& msg) {
    UL lock(mutex_);
    if (exit_flag_) return;
    if (enable_skip_) {
        if (skip_cnt_ != 0) {
            skip_cnt_++;
            skip_cnt_ = skip_cnt_ % skip_num_;
            return;
        }

        skip_cnt_++;
        skip_cnt_ = skip_cnt_ % skip_num_;
    }

    msg_buffer_.push_back(msg);
    while (msg_buffer_.size() > max_size_) {
        LOG_EVERY_N(INFO, 100) << name_ << " exceeds largest size: " << max_size_;
        msg_buffer_.pop_front();
        ++dropped_count_;
    }

    update_flag_ = true;
    cv_msg_.notify_one();
}

template <typename T>
void AsyncMessageProcess<T>::Quit() {
    {
        UL lock(mutex_);
        update_flag_ = true;
        exit_flag_ = true;
    }
    cv_msg_.notify_one();

    if (proc_.joinable()) {
        proc_.join();
    }
}

}  // namespace lightning::sys

#endif
