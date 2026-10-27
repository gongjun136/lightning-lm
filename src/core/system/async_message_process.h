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
 * @brief 使用有界队列和一个工作线程串行处理消息。
 * @details 定位系统用它分别组织传感器、地图匹配和高频输出；同一实例的回调串行执行。
 * @par 调用关系
 * 1. <b>接收消息</b>：AddMessage() 根据跳帧策略入队，超出容量时淘汰最旧的待处理消息。
 * 2. <b>执行回调</b>：ProcLoop() 每次取一条消息，释放队列锁后执行处理函数。
 * 3. <b>结束处理</b>：Quit() 停止接收新消息，处理完已接收队列并等待线程退出。
 * @tparam T 队列中保存并传给处理回调的消息类型。
 * @note 处理函数、容量和跳帧策略应在 Start() 前配置；对象销毁前由所有者调用 Quit()。
 * @see Start(), AddMessage(), Quit(), PendingCount()
 */
template <typename T>
class AsyncMessageProcess {
   public:
    using ProcFunc = std::function<void(const T&)>;  // 消息回调函数
    AsyncMessageProcess() = default;
    AsyncMessageProcess(ProcFunc proc_func, std::string name = "");

    /// @brief 设置工作线程使用的处理函数，应在 Start() 前完成。
    void SetProcFunc(ProcFunc proc_func) { custom_func_ = proc_func; }
    /// Run once on the worker before processing messages; configure before Start().
    void SetThreadInit(std::function<void()> init) { thread_init_ = std::move(init); }

    /**
     * @brief 设置后续入队使用的容量上限，不包括正在执行的那一条消息。
     * @pre size 大于 0，且在 Start() 前设置。
     * @see AddMessage(), PendingCount()
     */
    void SetMaxSize(size_t size) { max_size_ = size; }

    /**
     * @brief 重置运行计数并启动一个消费线程。
     * @pre 已设置有效处理函数；再次启动前先调用 Quit() 等待上一次线程退出。
     * @see SetProcFunc(), ProcLoop(), Quit()
     */
    void Start();

    /**
     * @brief 按跳帧策略接收消息，必要时淘汰旧数据，并唤醒消费线程。
     * @details 入队不会等待回调完成；开始退出后不再接收新消息。
     * @par 主要流程
     * 1. <b>检查退出与跳帧条件</b>：在队列锁下拒收退出后的输入，按计数决定是否跳过。
     * 2. <b>入队并限制等待容量</b>：加入新消息，超限时淘汰队头的旧消息。
     * 3. <b>唤醒消费线程</b>：设置更新标志并通知 worker，处理函数由 ProcLoop() 执行。
     * @note 阶段名称对应函数体；容量不包含已出队、正在执行的消息。
     * @htmlonly[block]
     * <details><summary>关键代码：队列超出容量时保留较新的输入</summary>
     * @endhtmlonly
     * @snippet core/system/async_message_process.h async-message-capacity
     * @htmlonly[block]
     * </details>
     * @endhtmlonly
     * @see SetMaxSize(), SetSkipParam(), DroppedCount()
     */
    void AddMessage(const T& msg);

    /**
     * @brief 拒绝新输入，等待已接收消息处理完毕，再回收工作线程。
     * @note 此调用会阻塞；由工作线程之外的所有者在停止生产输入后调用。
     * @see Start(), PendingCount()
     */
    void Quit();

    /// @brief 清空跳帧计数，使下一条输入不因跳帧计数而被跳过。
    void CleanSkipCnt();

    /// @brief 返回待处理队列长度与正在执行的消息数之和。
    size_t PendingCount() const;

    /**
     * @brief 返回容量淘汰和 RecordDroppedMessage() 显式记录的累计丢弃数。
     * @note 不包含跳帧策略跳过的消息或退出后拒收的消息；Start() 会重置计数。
     */
    size_t DroppedCount() const { return dropped_count_.load(); }

    /// @brief 由队列所有者记录入队前的主动丢弃。
    void RecordDroppedMessage() { ++dropped_count_; }

    /// @brief 返回处理回调已经完成的消息数，Start() 会重置计数。
    size_t ProcessedCount() const { return processed_count_.load(); }

    void SetName(std::string name) { name_ = std::move(name); }
    /**
     * @brief 设置周期跳帧；开启时每 skip_num 条输入接收一条。
     * @pre 开启时 skip_num 大于 0，且在 Start() 前配置；1 表示不跳帧。
     */
    void SetSkipParam(bool enable_skip, int skip_num) { enable_skip_ = enable_skip, skip_num_ = skip_num; }

    AsyncMessageProcess(const AsyncMessageProcess&) = delete;
    void operator=(const AsyncMessageProcess&) = delete;

   private:
    std::function<void()> thread_init_;  ///< Optional worker affinity/name initialization.
    /**
     * @brief 等待消息，逐条执行回调，直到退出请求已提出且队列为空。
     * @details 回调执行期间不持有队列锁；容量限制仅作用于等待中的消息。
     * @par 主要流程
     * 1. <b>初始化工作线程</b>：执行可选的线程初始化回调。
     * 2. <b>等待消息或退出</b>：等待更新/退出信号，仅在已请求退出且队列为空时结束。
     * 3. <b>取出一条消息并解锁</b>：从队头取消息，更新在途计数后释放队列锁。
     * 4. <b>执行回调并更新统计</b>：同步执行处理函数，完成后清除在途计数并累计完成数。
     * @note 阶段名称对应函数体；阶段 2–4 在循环中重复，退出请求不会取消正在执行的回调。
     * @see AddMessage(), Quit(), ProcessedCount()
     */
    void ProcLoop();

    std::thread proc_;  ///< 执行消息处理回调的单个消费线程。
    mutable std::mutex mutex_;  ///< 保护队列、退出标志和跳帧状态。
    std::condition_variable cv_msg_;  ///< 在新消息或退出请求到达时唤醒消费线程。
    std::deque<T> msg_buffer_;  ///< 按入队顺序保存尚未执行的消息。
    bool update_flag_ = false;  ///< 标记队列更新，供消费线程唤醒条件使用。
    bool exit_flag_ = false;  ///< 退出请求标志，置位后停止接收并排空队列。
    size_t max_size_ = 100;  ///< 待处理队列容量上限，超限时淘汰最旧消息。
    std::string name_;  ///< 用于线程标识及日志的队列名称。

    /// 跳帧
    bool enable_skip_ = false;
    int skip_num_ = 0;  ///< 跳帧周期参数，配合 enable_skip_ 使用。
    int skip_cnt_ = 0;  ///< 当前跳帧计数，用于决定是否接收入站消息。

    ProcFunc custom_func_;  ///< 由所有者设置的消息处理回调。
    std::atomic<size_t> in_flight_count_{0};  ///< 已出队且回调尚未完成的消息数。
    std::atomic<size_t> dropped_count_{0};  ///< 容量淘汰与所有者显式记录的丢弃总数。
    std::atomic<size_t> processed_count_{0};  ///< 处理回调已经完成的消息数。
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
    // 阶段 1：初始化工作线程。
    if (thread_init_) thread_init_();
    while (true) {
        // 阶段 2：等待消息或退出。
        // [async-message-exit-condition]
        UL lock(mutex_);
        cv_msg_.wait(lock, [this]() { return update_flag_ || exit_flag_; });
        if (exit_flag_ && msg_buffer_.empty()) break;
        // [async-message-exit-condition]

        // 阶段 3：取出一条消息并解锁。
        // Keep only the message currently executing outside the bounded
        // buffer. Swapping the complete queue into an unbounded in-flight
        // batch defeats max_size_ while a slow callback is running.
        // [async-message-execute-one]
        T msg = std::move(msg_buffer_.front());
        msg_buffer_.pop_front();
        in_flight_count_ = 1;
        update_flag_ = !msg_buffer_.empty();
        lock.unlock();

        // 阶段 4：执行回调并更新统计。
        custom_func_(msg);
        in_flight_count_ = 0;
        ++processed_count_;
        // [async-message-execute-one]
    }
}

template <typename T>
void AsyncMessageProcess<T>::AddMessage(const T& msg) {
    // 阶段 1：检查退出与跳帧条件。
    // [async-message-refuse-after-exit]
    UL lock(mutex_);
    if (exit_flag_) return;
    // [async-message-refuse-after-exit]
    if (enable_skip_) {
        if (skip_cnt_ != 0) {
            skip_cnt_++;
            skip_cnt_ = skip_cnt_ % skip_num_;
            return;
        }

        skip_cnt_++;
        skip_cnt_ = skip_cnt_ % skip_num_;
    }

    // 阶段 2：入队并限制等待容量。
    // [async-message-capacity]
    msg_buffer_.push_back(msg);
    while (msg_buffer_.size() > max_size_) {
        LOG_EVERY_N(INFO, 100) << name_ << " exceeds largest size: " << max_size_;
        msg_buffer_.pop_front();
        ++dropped_count_;
    }
    // [async-message-capacity]

    // 阶段 3：唤醒消费线程。
    update_flag_ = true;
    cv_msg_.notify_one();
}

template <typename T>
void AsyncMessageProcess<T>::Quit() {
    // [async-message-quit-join]
    {
        UL lock(mutex_);
        update_flag_ = true;
        exit_flag_ = true;
    }
    cv_msg_.notify_one();

    if (proc_.joinable()) {
        proc_.join();
    }
    // [async-message-quit-join]
}

}  // namespace lightning::sys

#endif
