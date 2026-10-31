#include "utils/incident_recorder.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>
#include <unistd.h>
#include <rosbag2_cpp/reader.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <std_msgs/msg/string.hpp>

namespace {
namespace fs = std::filesystem;
using Recorder = lightning::diagnostics::IncidentRecorder;
using Clock = std::chrono::steady_clock;
void Check(bool value, const std::string& message) { if (!value) throw std::runtime_error(message); }
void Sleep(double seconds) { std::this_thread::sleep_for(std::chrono::duration<double>(seconds)); }
struct Fixture {
    fs::path root;
    Fixture() {
        char pattern[] = "/tmp/lightning_incident_test.XXXXXX";
        auto path = mkdtemp(pattern); Check(path != nullptr, "mkdtemp failed"); root = path;
        std::ofstream(root / "config.yaml") << "diagnostics: {safety_fault_reporting_enabled: false}\n";
    }
    ~Fixture() { std::error_code ec; fs::remove_all(root, ec); }
};
Recorder::Options Options() {
    Recorder::Options o; o.enabled = true; o.pre_seconds = .5; o.post_seconds = .3;
    o.max_event_seconds = 1; o.cooldown_seconds = .01;
    o.startup_timeout_seconds = 30; o.output_timeout_seconds = 30; o.input_timeout_seconds = 30;
    o.memory_mb = 16; o.max_message_mb = 1; o.event_mb = 8; o.disk_mb = 32;
    o.min_free_mb = 0; o.segment_mb = 1; o.write_mb_per_second = 1000;
    return o;
}
std::vector<fs::path> EventPaths(const fs::path& archive) {
    std::vector<fs::path> paths;
    if (!fs::exists(archive)) return paths;
    for (const auto& p : fs::directory_iterator(archive))
        if (p.is_directory() && fs::exists(p.path() / "summary.json")) paths.push_back(p.path());
    std::sort(paths.begin(), paths.end()); return paths;
}
fs::path WaitComplete(const fs::path& archive, double timeout = 8, const fs::path& previous = {}) {
    const auto deadline = Clock::now() + std::chrono::duration<double>(timeout);
    while (Clock::now() < deadline) {
        auto paths = EventPaths(archive);
        if (!paths.empty() && paths.back() != previous && YAML::LoadFile((paths.back() / "summary.json").string())["status"].as<std::string>() != "capturing")
            return paths.back();
        Sleep(.02);
    }
    throw std::runtime_error("incident did not finish: " + archive.string());
}
std::string Read(const fs::path& path) {
    std::ifstream file(path); return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}
void MaskedFaultAndReadback() {
    Fixture f; const auto archive = f.root / "incidents";
    Recorder r(Options(), archive.string(), (f.root / "config.yaml").string());
    int stream = r.Register<std_msgs::msg::String>("/test/raw", Recorder::Kind::kDiagnostic);
    Check(r.Start(), "start failed");
    for (int i = 0; i < 8; ++i) {
        auto msg = std::make_shared<std_msgs::msg::String>(); msg->data = "before-" + std::to_string(i);
        r.Capture(stream, msg, sizeof(*msg) + msg->data.capacity()); Sleep(.025);
    }
    Recorder::Safety fault; fault.raw = fault.selected = 4; fault.published = 0; fault.uncertainty = 2;
    fault.yaw_condition = 2; fault.confirmed_failure = true; fault.error_yaw_rad = .0428282;
    fault.policy_version = "localization-dtc-v5-confirmed-loss";
    fault.startup_age_sec = 12.5; fault.valid_match_age_sec = .51;
    fault.match_stale = true; fault.latest_match_valid = true;
    r.ObserveSafety(fault); Sleep(.01); r.ObserveSafety({});  // shorter than a monitor period
    for (int i = 0; i < 6; ++i) {
        auto msg = std::make_shared<std_msgs::msg::String>(); msg->data = "after-" + std::to_string(i);
        r.Capture(stream, msg, sizeof(*msg) + msg->data.capacity()); Sleep(.025);
    }
    const auto path = WaitComplete(archive);
    auto summary = YAML::LoadFile((path / "summary.json").string());
    Check(summary["status"].as<std::string>() == "complete", "capture status");
    Check(summary["reasons"][0].as<std::string>() == "raw_4", "masked raw fault did not trigger");
    Check(summary["observed_first_stage"].as<std::string>() == "uncertainty_check", "short raw transition was lost");
    Check(summary["streams"][0]["first_relative_sec"].as<double>() < 0 &&
          summary["streams"][0]["last_relative_sec"].as<double>() > 0, "per-topic capture coverage missing");
    Check(Read(path / "timeline.csv").find(",4,4,0,") != std::string::npos, "masked transition missing in timeline");
    const auto timeline = Read(path / "timeline.csv");
    Check(timeline.find("consistency,uncertainty,translation_condition,yaw_condition,confirmed_failure") != std::string::npos &&
          timeline.find("0.0428282") != std::string::npos, "standalone event lost detailed Safety evidence");
    auto columns = [](const std::string& line) {
        std::vector<std::string> values; std::istringstream input(line); std::string value;
        while (std::getline(input, value, ',')) values.push_back(value);
        return values;
    };
    std::istringstream csv(timeline); std::string line; std::getline(csv, line);
    const auto header = columns(line);
    const auto yaw_column = std::find(header.begin(), header.end(), "yaw_condition") - header.begin();
    const auto error_column = std::find(header.begin(), header.end(), "error_yaw_rad") - header.begin();
    const auto policy_column = std::find(header.begin(), header.end(), "policy_version") - header.begin();
    const auto lost_column = std::find(header.begin(), header.end(), "confirmed_tracking_lost") - header.begin();
    const auto match_column = std::find(header.begin(), header.end(), "valid_match_age_sec") - header.begin();
    const auto stale_column = std::find(header.begin(), header.end(), "match_stale") - header.begin();
    const auto suppression_column = std::find(header.begin(), header.end(), "partial_radar_suppressed") - header.begin();
    Check(policy_column < header.size() && lost_column < header.size() && match_column < header.size() &&
          stale_column < header.size() && suppression_column < header.size(), "policy evidence columns missing");
    bool detailed_fault_seen = false;
    while (std::getline(csv, line)) {
        const auto values = columns(line);
        Check(values.size() == header.size(), "timeline header/row column mismatch");
        if (values[2] == "4") {
            Check(values[yaw_column] == "2" && std::abs(std::stod(values[error_column]) - .0428282) < 1e-9,
                  "fault details assigned to wrong CSV columns");
            Check(values[policy_column] == "localization-dtc-v5-confirmed-loss" && values[lost_column] == "0" &&
                  values[stale_column] == "1" && std::abs(std::stod(values[match_column]) - .51) < 1e-9 &&
                  values[suppression_column] == "0", "stale-vs-confirmed-loss policy evidence assigned incorrectly");
            detailed_fault_seen = true;
        }
    }
    Check(detailed_fault_seen, "standalone event omitted masked fault details");
    Check(fs::file_size(path / "summary.html") > 1000, "HTML report missing");
    unsigned before = 0, after = 0;
    for (auto segment : summary["closed_segments"]) {
        rosbag2_cpp::Reader reader;
        reader.open((path / segment.as<std::string>()).string());
        while (reader.has_next()) {
            auto stored = reader.read_next();
            rclcpp::SerializedMessage bytes(*stored->serialized_data);
            std_msgs::msg::String value;
            rclcpp::Serialization<std_msgs::msg::String>().deserialize_message(&bytes, &value);
            if (value.data.find("before-") == 0) ++before;
            if (value.data.find("after-") == 0) ++after;
        }
    }
    Check(before == 8 && after == 6, "pre/post message readback incomplete");
    r.Stop(); Check(r.GetStats().held_bytes == 0, "retained messages leaked on Stop");
    std::cout << "masked short fault: " << before << " pre + " << after << " post MCAP messages\n";
}
void BoundedAdmissionAndSlowWriter() {
    Fixture f; auto o = Options(); o.memory_mb = 1; o.write_mb_per_second = .2;
    Recorder r(o, (f.root / "incidents").string(), (f.root / "config.yaml").string());
    int stream = r.Register<sensor_msgs::msg::PointCloud2>("/test/cloud", Recorder::Kind::kLidar);
    Check(r.Start(), "start failed");
    auto cloud = std::make_shared<sensor_msgs::msg::PointCloud2>(); cloud->data.resize(300000);
    r.Capture(stream, cloud, sizeof(*cloud) + cloud->data.capacity(), 1);
    Recorder::Safety fault; fault.raw = -9; r.ObserveSafety(fault); Sleep(.15);
    std::vector<double> latency;
    for (int i = 0; i < 2000; ++i) {
        const auto start = Clock::now();
        r.Capture(stream, cloud, sizeof(*cloud) + cloud->data.capacity(), i + 2);
        latency.push_back(std::chrono::duration<double, std::milli>(Clock::now() - start).count());
    }
    auto stats = r.GetStats();
    Check(stats.peak_bytes <= o.memory_mb * 1024 * 1024, "retained byte budget exceeded");
    Check(stats.dropped_budget + stats.dropped_busy > 0, "pressure drops not reported");
    std::sort(latency.begin(), latency.end());
    Check(latency[1980] < 5, "Capture blocked behind the slow writer");
    std::cout << "slow writer: Capture p99_ms=" << latency[1980] << " peak_bytes=" << stats.peak_bytes
              << " dropped=" << stats.dropped_budget + stats.dropped_busy << '\n';
    r.Stop(); Check(r.GetStats().held_bytes == 0, "slow-writer retention leak");
}
void IndependentTimeouts() {
    Fixture f; auto o = Options(); o.startup_timeout_seconds = .15; o.post_seconds = .15;
    {
        Recorder r(o, (f.root / "startup").string(), (f.root / "config.yaml").string());
        Check(r.Start(), "start failed");
        auto summary = YAML::LoadFile((WaitComplete(f.root / "startup") / "summary.json").string());
        Check(summary["reasons"][0].as<std::string>() == "startup_without_output", "startup watchdog missing");
        r.Stop();
    }
    o.startup_timeout_seconds = 30; o.output_timeout_seconds = .15;
    Recorder r(o, (f.root / "publish").string(), (f.root / "config.yaml").string());
    Check(r.Start(), "start failed"); r.PublishEnd(); r.PosePrepared(); r.PublishBegin();
    auto summary = YAML::LoadFile((WaitComplete(f.root / "publish") / "summary.json").string());
    Check(summary["observed_first_stage"].as<std::string>() == "publish_call_pending", "publish begin/end not distinguished");
    r.Stop();

    o.output_timeout_seconds = 30; o.input_timeout_seconds = .15;
    Recorder frozen(o, (f.root / "frozen").string(), "");
    const auto imu = frozen.Register<std_msgs::msg::String>("/test/imu", Recorder::Kind::kImu);
    Check(frozen.Start(), "frozen stamp fixture start failed"); frozen.PublishEnd();
    for (int i = 0; i < 20; ++i) {
        frozen.Capture(imu, std::make_shared<std_msgs::msg::String>(), 100, 1000);
        Sleep(.03);
    }
    summary = YAML::LoadFile((WaitComplete(f.root / "frozen") / "summary.json").string());
    Check(summary["reasons"][0].as<std::string>() == "input_stamp_not_advancing", "frozen sensor clock not detected");
    Check(summary["observed_first_stage"].as<std::string>() == "input_stamp_not_advancing", "arrival silence confused with frozen stamp");
    frozen.Stop();
}

void ExportBudgetAndMergedFaults() {
    Fixture f; auto o = Options(); o.event_mb = 1; o.disk_mb = 4;
    Recorder r(o, (f.root / "incidents").string(), "");
    const auto id = r.Register<sensor_msgs::msg::PointCloud2>("/test/cloud", Recorder::Kind::kLidar);
    Check(r.Start(), "budget fixture start failed");
    for (int i = 0; i < 10; ++i) {
        auto cloud = std::make_shared<sensor_msgs::msg::PointCloud2>(); cloud->data.resize(300000);
        Check(r.Capture(id, cloud, sizeof(*cloud) + cloud->data.capacity(), i + 1), "budget fixture admission failed");
        Sleep(.005);
    }
    Recorder::Safety fault; fault.raw = 2; r.ObserveSafety(fault); Sleep(.15);
    fault.raw = 4; r.ObserveSafety(fault); Sleep(.15);
    fault.raw = -9; r.ObserveSafety(fault);
    const auto summary = YAML::LoadFile((WaitComplete(f.root / "incidents") / "summary.json").string());
    Check(summary["status"].as<std::string>() == "partial", "capacity loss must be reported");
    Check(summary["serialized_bytes"].as<std::uint64_t>() <= 1024 * 1024, "serialized event budget exceeded");
    Check(summary["omitted_messages"].as<unsigned>() > 0, "export omissions hidden");
    Check(summary["reasons"].size() == 3 && r.GetStats().events == 1, "faults did not merge");
    r.Stop();
}
void RetentionAndFailureIsolation() {
    Fixture f; auto o = Options(); o.max_events = 1;
    auto archive = f.root / "archive";
    fs::create_directories(archive / "event_user_files");
    std::ofstream(archive / "event_user_files/keep.txt") << "unowned";
    fs::path previous;
    for (int i = 0; i < 2; ++i) {
        Recorder r(o, archive.string(), (f.root / "config.yaml").string());
        Check(r.Start(), "archive lock not released");
        if (!i) {
            Recorder competing(o, archive.string(), "");
            Check(!competing.Start(), "two writers acquired the same archive");
        }
        Recorder::Safety fault; fault.raw = 2; r.ObserveSafety(fault);
        previous = WaitComplete(archive, 8, previous); r.Stop();
    }
    Check(EventPaths(archive).size() == 1, "owned retention count exceeded");
    Check(fs::exists(archive / "event_user_files/keep.txt"), "unowned files deleted");
    std::ofstream(f.root / "not_directory") << "file";
    Recorder invalid(o, (f.root / "not_directory").string(), "");
    Check(!invalid.Start() && invalid.GetStats().writer_errors == 1, "filesystem failure not isolated");
    Recorder disabled(Recorder::Options{}, (f.root / "disabled").string(), "");
    Check(!disabled.Start() && !fs::exists(f.root / "disabled"), "disabled collector touched disk");

    Recorder broken(o, (f.root / "broken_writer").string(), "");
    const auto id = broken.Register("/test/error", "std_msgs/msg/String", Recorder::Kind::kDiagnostic,
        [](const void*) -> std::shared_ptr<rclcpp::SerializedMessage> { throw std::runtime_error("injected export error"); });
    Check(broken.Start(), "failure fixture start failed");
    broken.Capture(id, std::make_shared<std_msgs::msg::String>(), 100);
    Recorder::Safety fault; fault.raw = 4; broken.ObserveSafety(fault);
    const auto failed = YAML::LoadFile((WaitComplete(f.root / "broken_writer") / "summary.json").string());
    Check(failed["status"].as<std::string>() == "error" && broken.GetStats().writer_errors == 1, "writer error not isolated/reported");
    broken.Stop();
}

void SuppressedRadarEvidence() {
    Fixture f;
    const auto archive = f.root / "incidents";
    Recorder r(Options(), archive.string(), (f.root / "config.yaml").string());
    Check(r.Start(), "suppressed-radar fixture start failed");
    Recorder::Safety fault;
    fault.raw = fault.strict = 2;
    fault.policy_version = "localization-dtc-v5-confirmed-loss";
    fault.output_gate_open = fault.latest_match_valid = fault.output_continuous = true;
    fault.localization_healthy_for_degraded = fault.partial_radar_suppressed = true;
    fault.valid_match_age_sec = .02;
    r.ObserveSafety(fault);
    const auto path = WaitComplete(archive);
    const auto summary = YAML::LoadFile((path / "summary.json").string());
    Check(summary["reasons"][0].as<std::string>() == "raw_2", "suppressed radar fault did not trigger capture");
    const auto csv = Read(path / "timeline.csv");
    Check(csv.find(",2,0,0,2,") != std::string::npos &&
          csv.find(",1,1,healthy_localization") != std::string::npos,
          "suppressed radar capture lost raw/selected/strict or suppression reason");
    r.Stop();
}

void CrashFixture(const fs::path& path) {
    Check(!fs::exists(path) || fs::is_empty(path), "crash fixture directory must be empty");
    auto o = Options(); o.pre_seconds = 1; o.post_seconds = o.max_event_seconds = 60;
    o.memory_mb = 16; o.event_mb = 32; o.disk_mb = 64;
    Recorder r(o, fs::absolute(path).string(), "");
    const auto id = r.Register<sensor_msgs::msg::PointCloud2>("/test/cloud", Recorder::Kind::kLidar);
    Check(r.Start(), "crash fixture start failed");
    Recorder::Safety fault; fault.raw = -9; r.ObserveSafety(fault);
    const auto deadline = Clock::now() + std::chrono::seconds(30);
    while (Clock::now() < deadline) {
        auto cloud = std::make_shared<sensor_msgs::msg::PointCloud2>(); cloud->data.resize(100000, 42);
        r.Capture(id, cloud, sizeof(*cloud) + cloud->data.capacity(), 1);
        Sleep(.02);
    }
    throw std::runtime_error("crash fixture was not killed by harness");
}
}  // namespace
int main(int argc, char** argv) {
    try {
        using lightning::diagnostics::detail::IncidentAgeSeconds;
        Check(IncidentAgeSeconds(1000000000, 1000020000) == 0.0, "concurrent publication was mistaken for never published");
        Check(IncidentAgeSeconds(1000000000, 0) == -1.0, "missing publication sentinel changed");
        Check(std::abs(IncidentAgeSeconds(1000000000, 800000000) - .2) < 1e-9, "stale observation age changed");
        if (argc == 3 && std::string(argv[1]) == "--crash-fixture") { CrashFixture(argv[2]); return 1; }
        MaskedFaultAndReadback(); SuppressedRadarEvidence(); BoundedAdmissionAndSlowWriter(); IndependentTimeouts();
        ExportBudgetAndMergedFaults(); RetentionAndFailureIsolation();
        std::cout << "incident_recorder_test passed\n"; return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
