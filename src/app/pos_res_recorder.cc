#include <rclcpp/rclcpp.hpp>

#include <geosun_msgs/msg/pos_res.hpp>
#include <tf2/LinearMath/Quaternion.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>

namespace {

namespace fs = std::filesystem;
using SteadyClock = std::chrono::steady_clock;

double ToSeconds(const builtin_interfaces::msg::Time& stamp) {
    return static_cast<double>(stamp.sec) + static_cast<double>(stamp.nanosec) * 1e-9;
}

double ElapsedSeconds(const SteadyClock::time_point& begin, const SteadyClock::time_point& end) {
    return std::chrono::duration<double>(end - begin).count();
}

void OpenOutput(std::ofstream& stream, const fs::path& path, const std::string& description) {
    if (path.has_parent_path()) fs::create_directories(path.parent_path());
    stream.open(path, std::ios::out | std::ios::trunc);
    if (!stream) throw std::runtime_error("failed to open " + description + ": " + path.string());
}

class PosResRecorder final : public rclcpp::Node {
   public:
    PosResRecorder() : Node("pos_res_recorder") {
        const std::string topic = declare_parameter<std::string>("topic", "/PosRes");
        const fs::path output_dir = declare_parameter<std::string>("output_dir", ".");
        const fs::path tum_path = output_dir / declare_parameter<std::string>(
                                                   "tum_file", "posres_trajectory.tum");
        const fs::path timing_path = output_dir / declare_parameter<std::string>(
                                                      "timing_file", "posres_timing.csv");
        const fs::path rate_path = output_dir / declare_parameter<std::string>(
                                                    "rate_file", "posres_rate.csv");

        OpenOutput(tum_stream_, tum_path, "TUM trajectory");
        OpenOutput(timing_stream_, timing_path, "timing CSV");
        OpenOutput(rate_stream_, rate_path, "rate CSV");

        tum_stream_ << "# timestamp tx ty tz qx qy qz qw\n";
        timing_stream_
            << "sequence,msg_stamp_s,receive_stamp_s,receive_elapsed_s,latency_ms,"
               "interarrival_receive_ms,interarrival_message_ms\n";
        rate_stream_ << "window_start_elapsed_s,window_end_elapsed_s,window_duration_s,"
                        "message_count,rate_hz,mean_interarrival_ms\n";
        tum_stream_ << std::fixed;
        timing_stream_ << std::fixed;
        rate_stream_ << std::fixed;

        start_steady_ = SteadyClock::now();
        window_start_steady_ = start_steady_;
        const auto qos = rclcpp::QoS(rclcpp::KeepLast(1000)).reliable().durability_volatile();
        subscription_ = create_subscription<geosun_msgs::msg::PosRes>(
            topic, qos, [this](const geosun_msgs::msg::PosRes::SharedPtr message) { Record(*message); });
        rate_timer_ = create_wall_timer(std::chrono::seconds(1), [this]() { CloseRateWindow(false); });

        RCLCPP_INFO(get_logger(),
                    "recording %s to TUM=%s, timing=%s, rate=%s",
                    topic.c_str(), tum_path.c_str(), timing_path.c_str(), rate_path.c_str());
    }

    ~PosResRecorder() override {
        CloseRateWindow(true);
        RCLCPP_INFO(get_logger(), "received %llu messages and wrote %llu TUM poses",
                    static_cast<unsigned long long>(sequence_),
                    static_cast<unsigned long long>(tum_pose_count_));
    }

   private:
    void Record(const geosun_msgs::msg::PosRes& message) {
        const auto receive_steady = SteadyClock::now();
        const double receive_stamp = now().seconds();
        const double receive_elapsed = ElapsedSeconds(start_steady_, receive_steady);
        const double message_stamp = ToSeconds(message.header.stamp);
        const double latency_ms = message_stamp > 0.0
                                      ? (receive_stamp - message_stamp) * 1000.0
                                      : std::numeric_limits<double>::quiet_NaN();
        const double receive_interval_ms =
            has_previous_ ? ElapsedSeconds(previous_receive_steady_, receive_steady) * 1000.0
                          : std::numeric_limits<double>::quiet_NaN();
        const double message_interval_ms =
            has_previous_ ? (message_stamp - previous_message_stamp_) * 1000.0
                          : std::numeric_limits<double>::quiet_NaN();

        timing_stream_ << sequence_ << ',' << std::setprecision(9) << message_stamp << ','
                       << receive_stamp << ',' << std::setprecision(6) << receive_elapsed << ','
                       << latency_ms << ',' << receive_interval_ms << ',' << message_interval_ms
                       << '\n';

        if (has_previous_ && std::isfinite(receive_interval_ms)) {
            window_interval_sum_ms_ += receive_interval_ms;
            ++window_interval_count_;
        }
        previous_receive_steady_ = receive_steady;
        previous_message_stamp_ = message_stamp;
        has_previous_ = true;
        ++sequence_;
        ++window_message_count_;

        if (!WriteTumPose(message, message_stamp)) return;
        // Flush per message so Ctrl-C leaves useful output on disk.
        tum_stream_.flush();
        timing_stream_.flush();
    }

    bool WriteTumPose(const geosun_msgs::msg::PosRes& message, double stamp) {
        const bool finite_pose = std::all_of(message.f8enh.begin(), message.f8enh.end(),
                                             [](double value) { return std::isfinite(value); }) &&
                                 std::all_of(message.f8pry.begin(), message.f8pry.end(),
                                             [](double value) { return std::isfinite(value); });
        if (!(stamp > 0.0) || !finite_pose) {
            RCLCPP_WARN_THROTTLE(get_logger(), *get_clock(), 5000,
                                 "skipping TUM pose with invalid timestamp or non-finite pose");
            return false;
        }
        if (has_last_tum_stamp_ && stamp <= last_tum_stamp_) {
            RCLCPP_WARN_THROTTLE(get_logger(), *get_clock(), 5000,
                                 "skipping non-monotonic TUM timestamp %.9f (last %.9f)", stamp,
                                 last_tum_stamp_);
            return false;
        }

        constexpr double kRadiansPerDegree = 3.14159265358979323846 / 180.0;
        tf2::Quaternion quaternion;
        quaternion.setRPY(message.f8pry[0] * kRadiansPerDegree,
                          message.f8pry[1] * kRadiansPerDegree,
                          message.f8pry[2] * kRadiansPerDegree);
        quaternion.normalize();
        tum_stream_ << std::setprecision(9) << stamp << ' ' << std::setprecision(12)
                    << message.f8enh[0] << ' ' << message.f8enh[1] << ' ' << message.f8enh[2] << ' '
                    << quaternion.x() << ' ' << quaternion.y() << ' ' << quaternion.z() << ' '
                    << quaternion.w() << '\n';
        last_tum_stamp_ = stamp;
        has_last_tum_stamp_ = true;
        ++tum_pose_count_;
        return true;
    }

    void CloseRateWindow(bool final_window) {
        const auto now_steady = SteadyClock::now();
        const double duration = ElapsedSeconds(window_start_steady_, now_steady);
        if (duration <= 0.0 || (final_window && window_message_count_ == 0)) return;
        const double start_elapsed = ElapsedSeconds(start_steady_, window_start_steady_);
        const double end_elapsed = ElapsedSeconds(start_steady_, now_steady);
        const double rate = static_cast<double>(window_message_count_) / duration;
        const double mean_interval = window_interval_count_ > 0
                                         ? window_interval_sum_ms_ /
                                               static_cast<double>(window_interval_count_)
                                         : std::numeric_limits<double>::quiet_NaN();

        rate_stream_ << std::setprecision(6) << start_elapsed << ',' << end_elapsed << ','
                     << duration << ',' << window_message_count_ << ',' << rate << ','
                     << mean_interval << '\n';
        rate_stream_.flush();
        RCLCPP_INFO(get_logger(), "window %.3f-%.3f s: %llu messages, %.2f Hz",
                    start_elapsed, end_elapsed,
                    static_cast<unsigned long long>(window_message_count_), rate);

        window_start_steady_ = now_steady;
        window_message_count_ = 0;
        window_interval_sum_ms_ = 0.0;
        window_interval_count_ = 0;
    }

    rclcpp::Subscription<geosun_msgs::msg::PosRes>::SharedPtr subscription_;
    rclcpp::TimerBase::SharedPtr rate_timer_;
    std::ofstream tum_stream_;
    std::ofstream timing_stream_;
    std::ofstream rate_stream_;
    SteadyClock::time_point start_steady_;
    SteadyClock::time_point window_start_steady_;
    SteadyClock::time_point previous_receive_steady_;
    bool has_previous_ = false;
    bool has_last_tum_stamp_ = false;
    double previous_message_stamp_ = 0.0;
    double last_tum_stamp_ = 0.0;
    double window_interval_sum_ms_ = 0.0;
    std::uint64_t sequence_ = 0;
    std::uint64_t tum_pose_count_ = 0;
    std::uint64_t window_message_count_ = 0;
    std::uint64_t window_interval_count_ = 0;
};

}  // namespace

int main(int argc, char** argv) {
    rclcpp::init(argc, argv);
    try {
        auto node = std::make_shared<PosResRecorder>();
        rclcpp::spin(node);
        node.reset();
    } catch (const std::exception& error) {
        RCLCPP_FATAL(rclcpp::get_logger("pos_res_recorder"), "%s", error.what());
        rclcpp::shutdown();
        return 1;
    }
    rclcpp::shutdown();
    return 0;
}
