//
// Created by xiang on 2022/2/9.
//

#pragma once

#include <glog/logging.h>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>
#include <utility>

#include "common/eigen_types.h"
#include "common/std_types.h"

namespace lightning {

/**
 * @brief 单线程异步消息处理队列。
 *
 * AsyncMessageProcess内部维护一个消息队列和一个后台线程。
 * AddMessage()只负责入队并唤醒线程，ProcLoop()在后台串行调用用户注册的处理函数。
 * 该类主要用于把耗时模块从前端调用链中拆开，避免在业务代码里重复编写锁和条件变量。
 *
 * @tparam T 队列中传递的消息类型。
 *
 * @note skip_num_设为1时不会实际跳帧；设为2时表现为处理一帧、跳过一帧。
 */
template <typename T>
class AsyncMessageProcess {
   public:
    using ProcFunc = std::function<void(const T&)>;  ///< 消息处理回调函数。

    AsyncMessageProcess() = default;

    /**
     * @brief 构造异步处理队列并设置处理函数。
     * @param proc_func 后台线程中执行的消息处理函数。
     * @param name 队列名称，用于日志输出。
     */
    AsyncMessageProcess(ProcFunc proc_func, std::string name = "");

    /**
     * @brief 析构时停止后台线程。
     */
    ~AsyncMessageProcess() { Quit(); }

    /**
     * @brief 设置或替换消息处理函数。
     * @param proc_func 后台线程中执行的消息处理函数。
     */
    void SetProcFunc(ProcFunc proc_func) { custom_func_ = proc_func; }

    /**
     * @brief 设置队列最大长度。
     * @param size 队列允许保留的最大消息数，超出后丢弃最旧消息。
     */
    void SetMaxSize(size_t size) { max_size_ = size; }

    /**
     * @brief 启动后台处理线程。
     */
    void Start();

    /**
     * @brief 添加一条待处理消息。
     * @param msg 待入队消息。
     */
    void AddMessage(const T& msg);

    /**
     * @brief 请求后台线程退出并等待线程结束。
     */
    void Quit();

    /**
     * @brief 清空跳帧计数器，使下一条消息立即进入队列。
     */
    void CleanSkipCnt();

    /**
     * @brief 设置队列名称，用于日志输出。
     * @param name 队列名称。
     */
    void SetName(std::string name) { name_ = std::move(name); }

    /**
     * @brief 设置跳帧参数。
     * @param enable_skip 是否启用跳帧。
     * @param skip_num 跳帧周期。
     */
    void SetSkipParam(bool enable_skip, int skip_num) { enable_skip_ = enable_skip, skip_num_ = skip_num; }

    AsyncMessageProcess(const AsyncMessageProcess&) = delete;
    void operator=(const AsyncMessageProcess&) = delete;

   private:
    /**
     * @brief 后台线程主循环。
     */
    void ProcLoop();

    std::thread proc_;               ///< 后台处理线程。
    std::mutex mutex_;               ///< 保护消息队列和状态标志的互斥锁。
    std::condition_variable cv_msg_;  ///< 用于唤醒后台线程的条件变量。
    std::deque<T> msg_buffer_;        ///< 待处理消息队列。
    bool update_flag_ = false;        ///< 是否有新消息或退出请求。
    bool exit_flag_ = false;          ///< 后台线程退出标志。
    size_t max_size_ = 40;            ///< 队列最大长度。
    std::string name_;                ///< 队列名称，用于日志输出。

    bool enable_skip_ = false;  ///< 是否启用跳帧。
    int skip_num_ = 0;          ///< 跳帧周期。
    int skip_cnt_ = 0;          ///< 当前跳帧计数。

    ProcFunc custom_func_;  ///< 用户注册的消息处理函数。
};

template <typename T>
void AsyncMessageProcess<T>::CleanSkipCnt() {
    // 重置计数后，下一条AddMessage()收到的消息会被放入队列。
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
    // 每次启动都重置状态标志，后台线程随后进入ProcLoop()等待消息。
    exit_flag_ = false;
    update_flag_ = false;
    proc_ = std::thread([this]() { ProcLoop(); });
}

template <typename T>
void AsyncMessageProcess<T>::ProcLoop() {
    while (!exit_flag_) {
        UL lock(mutex_);
        cv_msg_.wait(lock, [this]() { return update_flag_; });

        // 将当前队列复制到局部变量后清空队列，缩短持锁时间。
        // 这样custom_func_执行较慢时，AddMessage()仍能继续入队新消息。
        auto buffer = msg_buffer_;
        msg_buffer_.clear();
        update_flag_ = false;
        lock.unlock();

        // 同一批消息在后台线程中串行处理，避免回调并发执行。
        for (const auto& msg : buffer) {
            custom_func_(msg);
        }
    }
}

template <typename T>
void AsyncMessageProcess<T>::AddMessage(const T& msg) {
    UL lock(mutex_);
    if (enable_skip_) {
        // skip_cnt_为0的消息会进入队列，其余消息按skip_num_周期丢弃。
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
        // 后台线程处理不过来时丢弃最旧消息，保证队列长度有界。
        LOG(ERROR) << name_ << " exceeds largest size: " << max_size_;
        msg_buffer_.pop_front();
    }

    // 标记有新消息并唤醒后台线程。
    update_flag_ = true;
    cv_msg_.notify_one();
}

template <typename T>
void AsyncMessageProcess<T>::Quit() {
    // 即使线程正阻塞在条件变量上，也通过update_flag_唤醒它检查exit_flag_。
    update_flag_ = true;
    exit_flag_ = true;
    cv_msg_.notify_one();

    if (proc_.joinable()) {
        proc_.join();
    }
}

}  // namespace lightning
