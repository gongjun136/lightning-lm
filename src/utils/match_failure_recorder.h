#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <yaml-cpp/yaml.h>

namespace lightning::diagnostics {

// Captures the actual tracking acceptance input, before rejected poses fall back
// to odometry. The callback only performs bounded copying into preallocated slots.
// The immutable map must outlive Stop(); all filesystem work runs on the writer.
class MatchFailureRecorder {
 public:
    struct Point { float x = 0, y = 0, z = 0; std::uint32_t lidar_id = 0; };
    struct CloudView {
        const void* data = nullptr;
        std::size_t size = 0;
        Point (*point)(const void*, std::size_t) = nullptr;
    };
    using Pose = std::array<double, 16>;  // row-major T_map_lidar
    enum class Reason { kOverlap, kBounds, kGravity, kOdometry, kRegistration };
    struct Frame {
        std::int64_t sensor_stamp_ns = 0;
        Pose guess{}, ndt{}, check{};
        Reason reason = Reason::kRegistration;
        bool evaluated = false;
        std::uint64_t validation_stride = 1, validation_points = 0;
        double overlap = 0, inside_xy = 0, inside_xyz = 0, gravity_cos = 0;
        double min_overlap = .75, min_inside_xy = .9, min_gravity_cos = .95;
        double nearest_neighbor_distance_m = .5, bounds_margin_m = .5;
        double confidence = 0, raw_odom_delta = 0;
        int iterations = 0;
        bool raw_odom_valid = false;
    };
    struct Options {
        bool enabled = false;
        double cooldown_seconds = 30, recovery_seconds = 2, write_mb_per_second = 10;
        unsigned protect_days = 7;
        std::uint64_t disk_mb = 1024, min_free_mb = 2048, event_mb = 64;
        unsigned max_scan_points = 50000;
        static Options FromYaml(const YAML::Node& node);
        void Validate() const;
    };
    struct Context {
        std::string directory, config_path, config_snapshot, run_dir, map_path;
        std::string build_commit, build_dirty, build_utc;
        bool synthetic_fixture = false;
        std::array<double, 3> map_min{}, map_max{};
    };
    struct Stats {
        std::uint64_t received = 0, admitted = 0, saved = 0, suppressed = 0;
        std::uint64_t dropped_busy = 0, dropped_queue = 0, dropped_oversize = 0;
        std::uint64_t rejected_disk = 0, rejected_map = 0, writer_errors = 0;
    };
    MatchFailureRecorder(Options options, Context context, CloudView immutable_map);
    ~MatchFailureRecorder();
    MatchFailureRecorder(const MatchFailureRecorder&) = delete;
    MatchFailureRecorder& operator=(const MatchFailureRecorder&) = delete;
    bool Start();
    void ObserveTracking(bool accepted);
    bool Capture(const CloudView& input, const Frame& frame);
    void Stop();
    Stats GetStats() const;

 private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};
}  // namespace lightning::diagnostics
