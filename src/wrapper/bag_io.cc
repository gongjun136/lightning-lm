//
// Created by xiang on 23-12-14.
//

#include "bag_io.h"

#include <glog/logging.h>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <thread>
#include <rosbag2_cpp/reader.hpp>
#include <rosbag2_cpp/readers/sequential_reader.hpp>

namespace lightning {
namespace {

std::string StorageIdForBag(const std::string& bag_path) {
    const std::filesystem::path path(bag_path);
    if (path.extension() == ".mcap") return "mcap";
    if (std::filesystem::is_directory(path)) {
        for (const auto& entry : std::filesystem::directory_iterator(path)) {
            if (entry.path().extension() == ".mcap") return "mcap";
        }
    }
    return "sqlite3";
}

}  // namespace

void RosbagIO::Go(int sleep_usec) {
    rosbag2_cpp::Reader reader(std::make_unique<rosbag2_cpp::readers::SequentialReader>());
    rosbag2_cpp::ConverterOptions cv_options{"cdr", "cdr"};
    reader.open({bag_file_, StorageIdForBag(bag_file_)}, cv_options);

    while (reader.has_next()) {
        auto msg = reader.read_next();
        auto iter = process_func_.find(msg->topic_name);
        if (iter != process_func_.end()) {
            iter->second(msg);
        }

        if (sleep_usec > 0) {
            usleep(sleep_usec);
        }

        if (lightning::debug::flg_exit) {
            LOG(INFO) << "bag " << bag_file_ << " exit.";
            return;
        }
    }

    LOG(INFO) << "bag " << bag_file_ << " finished.";
}

void RosbagIO::GoRealtime(double playback_rate) {
    if (!std::isfinite(playback_rate) || playback_rate <= 0.0) {
        throw std::invalid_argument("playback rate must be finite and positive");
    }

    rosbag2_cpp::Reader reader(std::make_unique<rosbag2_cpp::readers::SequentialReader>());
    rosbag2_cpp::ConverterOptions cv_options{"cdr", "cdr"};
    reader.open({bag_file_, StorageIdForBag(bag_file_)}, cv_options);

    std::optional<rcutils_time_point_value_t> first_bag_time;
    std::chrono::steady_clock::time_point first_wall_time;
    while (reader.has_next()) {
        auto msg = reader.read_next();
        if (!first_bag_time) {
            first_bag_time = msg->time_stamp;
            first_wall_time = std::chrono::steady_clock::now();
        }
        const auto elapsed_ns = static_cast<int64_t>((msg->time_stamp - *first_bag_time) / playback_rate);
        std::this_thread::sleep_until(first_wall_time + std::chrono::nanoseconds(elapsed_ns));

        const auto iter = process_func_.find(msg->topic_name);
        if (iter != process_func_.end()) iter->second(msg);
        if (lightning::debug::flg_exit) {
            LOG(INFO) << "realtime bag " << bag_file_ << " exit.";
            return;
        }
    }
    LOG(INFO) << "realtime bag " << bag_file_ << " finished.";
}

}  // namespace lightning
