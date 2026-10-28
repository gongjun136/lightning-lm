#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <rclcpp/serialization.hpp>
#include <rclcpp/serialized_message.hpp>
#include <rosidl_runtime_cpp/traits.hpp>
#include <yaml-cpp/yaml.h>

namespace lightning::diagnostics {

// No subscriptions or processes: retain existing immutable callback messages.
// Serialization, filesystem access, retention and message destruction run off
// the algorithm thread. Admission uses try_lock and reports every dropped item.
class IncidentRecorder {
 public:
    enum class Kind { kLidar, kImu, kWheelSpeed, kOutput, kDiagnostic };
    struct Options {
        bool enabled = false;
        double pre_seconds = 30, post_seconds = 10, max_event_seconds = 60;
        double startup_timeout_seconds = 8, output_timeout_seconds = .5;
        double input_timeout_seconds = .5, cooldown_seconds = 10;
        double write_mb_per_second = 20;
        std::uint64_t memory_mb = 2048, max_message_mb = 64, event_mb = 2048;
        std::uint64_t disk_mb = 8192, min_free_mb = 2048, segment_mb = 32;
        unsigned max_events = 10;
        static Options FromYaml(const YAML::Node& node);
        void Validate() const;
    };
    struct Safety {
        int raw = 0, selected = 0, published = 0, strict = 0;
        int consistency = 0, uncertainty = 0;
        bool covariance_valid = false, unobservable = false;
        double major_sigma = -1, lateral_sigma = -1, yaw_sigma = -1;
    };
    struct Pipeline {
        std::uint64_t sensor_pending = 0, localization_pending = 0, high_frequency_pending = 0;
        std::uint64_t sensor_dropped = 0, localization_dropped = 0;
        double sensor_lag = -1, match_age = -1, lio_ms = -1;
        bool match_stale = false, output_gate = false;
        unsigned missing_lidars = 0;
    };
    struct Stats {
        std::uint64_t held_bytes = 0, peak_bytes = 0, received = 0, admitted = 0;
        std::uint64_t dropped_busy = 0, dropped_budget = 0, dropped_oversize = 0;
        std::uint64_t events = 0, merged = 0, suppressed = 0, writer_errors = 0, written = 0;
    };
    using Serializer = std::function<std::shared_ptr<rclcpp::SerializedMessage>(const void*)>;
    IncidentRecorder(Options options, std::string directory, std::string config_path);
    ~IncidentRecorder();
    IncidentRecorder(const IncidentRecorder&) = delete;
    IncidentRecorder& operator=(const IncidentRecorder&) = delete;

    int Register(const std::string& topic, const std::string& type, Kind kind, Serializer serializer);
    template <class Message>
    int Register(const std::string& topic, Kind kind) {
        return Register(topic, rosidl_generator_traits::name<Message>(), kind, [](const void* value) {
            auto serialized = std::make_shared<rclcpp::SerializedMessage>();
            rclcpp::Serialization<Message> serializer;
            serializer.serialize_message(static_cast<const Message*>(value), serialized.get());
            return serialized;
        });
    }
    bool Start();
    bool Capture(int stream, std::shared_ptr<const void> message, std::uint64_t retained_bytes,
                 std::int64_t sensor_stamp_ns = 0);
    void ObserveSafety(const Safety& snapshot);
    void ObservePipeline(const Pipeline& snapshot);
    void PosePrepared();
    void PublishBegin();
    void PublishEnd();
    void BeginShutdown();
    void Stop();
    Stats GetStats() const;

 private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};
}  // namespace lightning::diagnostics
