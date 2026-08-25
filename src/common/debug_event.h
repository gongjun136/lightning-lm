#pragma once

#include <chrono>
#include <functional>
#include <map>
#include <mutex>
#include <string>
#include <string_view>
#include <utility>

namespace lightning::debug_event {

using Sink = std::function<void(const std::string&)>;

namespace detail {

using Clock = std::chrono::steady_clock;

struct TransitionState {
    bool initialized = false;
    bool published = false;
    bool pending = false;
    Clock::time_point pending_since{};
    std::string active_message;
    std::string recovery_message;
};

struct Dispatcher {
    std::mutex mutex;
    Sink sink;
    std::map<std::string, Clock::time_point> last_emit;
    std::map<std::string, TransitionState> transitions;
};

inline Dispatcher& GetDispatcher() {
    static Dispatcher dispatcher;
    return dispatcher;
}

inline void Dispatch(Sink sink, std::string message) {
    if (sink && !message.empty()) sink(message);
}

}  // namespace detail

/// Install the process-wide sink used by the online localization node.
/// Setting a new sink also starts a fresh debounce epoch.
inline void SetSink(Sink sink) {
    auto& dispatcher = detail::GetDispatcher();
    std::lock_guard<std::mutex> lock(dispatcher.mutex);
    dispatcher.sink = std::move(sink);
    dispatcher.last_emit.clear();
    dispatcher.transitions.clear();
}

inline void ClearSink() {
    auto& dispatcher = detail::GetDispatcher();
    std::lock_guard<std::mutex> lock(dispatcher.mutex);
    dispatcher.sink = nullptr;
    dispatcher.last_emit.clear();
    dispatcher.transitions.clear();
}

/// Emit a one-shot event immediately.
inline void Emit(const std::string& message) {
    Sink sink;
    {
        auto& dispatcher = detail::GetDispatcher();
        std::lock_guard<std::mutex> lock(dispatcher.mutex);
        sink = dispatcher.sink;
    }
    detail::Dispatch(std::move(sink), message);
}

/// Emit at most once per interval for a noisy event key.
inline void EmitThrottled(const std::string& key, const std::string& message,
                          std::chrono::milliseconds interval) {
    Sink sink;
    bool should_emit = false;
    {
        auto& dispatcher = detail::GetDispatcher();
        const auto now = detail::Clock::now();
        std::lock_guard<std::mutex> lock(dispatcher.mutex);
        const auto found = dispatcher.last_emit.find(key);
        if (found == dispatcher.last_emit.end() || now - found->second >= interval) {
            dispatcher.last_emit[key] = now;
            sink = dispatcher.sink;
            should_emit = true;
        }
    }
    if (should_emit) detail::Dispatch(std::move(sink), message);
}

/// Publish only after an active/recovered state remains stable for debounce.
/// The first normal observation establishes the baseline without publishing.
inline void ReportState(const std::string& key, bool active,
                        std::string_view active_message,
                        std::string_view recovery_message,
                        std::chrono::milliseconds debounce) {
    Sink sink;
    std::string message;
    {
        auto& dispatcher = detail::GetDispatcher();
        const auto now = detail::Clock::now();
        std::lock_guard<std::mutex> lock(dispatcher.mutex);
        auto& state = dispatcher.transitions[key];

        if (!state.initialized) {
            state.active_message.assign(active_message.begin(), active_message.end());
            state.recovery_message.assign(recovery_message.begin(), recovery_message.end());
            state.initialized = true;
            state.pending = active;
            state.pending_since = now;
            if (!active) return;
        } else if (state.pending != active) {
            state.active_message.assign(active_message.begin(), active_message.end());
            state.recovery_message.assign(recovery_message.begin(), recovery_message.end());
            state.pending = active;
            state.pending_since = now;
        }

        if (state.pending == state.published || now - state.pending_since < debounce) return;
        state.published = state.pending;
        message = state.published ? state.active_message : state.recovery_message;
        sink = dispatcher.sink;
    }
    detail::Dispatch(std::move(sink), std::move(message));
}

}  // namespace lightning::debug_event
