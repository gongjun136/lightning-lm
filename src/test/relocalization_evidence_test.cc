#include "core/localization/lidar_loc/lidar_loc.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <thread>
#include <unistd.h>

namespace fs = std::filesystem;
using Recorder = lightning::diagnostics::MatchFailureRecorder;

namespace lightning::loc {
// Controlled retrieval responses exercise the real candidate rejection path,
// map-consistency calculation and asynchronous evidence writer together.
class EvidenceRetrievalFixture final : public GlobalRelocalizer {
 public:
    bool ready = true, candidate = false;
    bool Init(const std::string&, const std::string&, const SE3&) override { return true; }
    void ResetQuery() override {}
    bool IsReady() const override { return ready; }
    std::size_t DatabaseSize() const override { return 1; }
    const char* Name() const override { return "synthetic_evidence_fixture"; }
    std::optional<RelocalizationResult> AddFrame(const CloudPtr&, const SE3&, double stamp) override {
        RelocalizationResult result;
        result.attempted = true;
        result.timestamp = stamp;
        result.reason = candidate ? "candidate_found" : "no_candidate";
        if (candidate) {
            result.candidate_found = result.accepted = true;
            result.candidate_id = 17;
            RelocalizationCandidate item;
            item.candidate_id = 17;
            item.score = 0.8;
            item.T_world_imu = SE3(SO3(), Vec3d(100, 0, 0));
            result.candidates.push_back(item);
        }
        return result;
    }
};

struct LidarLocRegressionAccess {
    static EvidenceRetrievalFixture* Prepare(LidarLoc& loc, const CloudPtr& map,
                                             const fs::path& folder) {
        loc.current_lo_pose_set_ = true;
        loc.current_lo_pose_ = SE3();
        loc.current_timestamp_ = 1000;
        loc.options_.relocalization_validation_workers_ = 1;
        loc.relocalization_static_map_ = map;
        loc.relocalization_map_min_ = Vec3d(-2, -2, -2);
        loc.relocalization_map_max_ = Vec3d(2, 2, 2);
        auto tree = std::make_unique<pcl::KdTreeFLANN<PointType>>();
        tree->setInputCloud(map);
        loc.relocalization_kdtrees_.push_back(std::move(tree));
        auto backend = std::make_unique<EvidenceRetrievalFixture>();
        auto* pointer = backend.get();
        loc.global_relocalizer_ = std::move(backend);
        Recorder::Options options;
        options.enabled = true;
        options.min_free_mb = 0;
        Recorder::Context context;
        context.directory = (folder / "evidence").string();
        context.run_dir = folder.string();
        context.config_path = (folder / "config.yaml").string();
        context.config_snapshot = "diagnostics:\n  output_timeout_sec: 0.5\n  startup_grace_sec: 8.0\n";
        context.map_path = "synthetic_in_memory_map";
        context.map_min = {{-2, -2, -2}};
        context.map_max = {{2, 2, 2}};
        context.synthetic_fixture = true;
        context.build_commit = "test_fixture";
        Recorder::CloudView view{map.get(), map->size(), [](const void* data, std::size_t index) {
            const auto& p = (*static_cast<const PointCloudType*>(data))[index];
            return Recorder::Point{p.x, p.y, p.z, p.lidar_id};
        }};
        loc.match_failure_recorder_ = std::make_unique<Recorder>(options, context, view);
        if (!loc.match_failure_recorder_->Start()) throw std::runtime_error("fixture recorder start failed");
        return pointer;
    }
    static bool Attempt(LidarLoc& loc, const CloudPtr& scan) {
        loc.current_timestamp_ += 0.1;
        return loc.TryGlobalRelocalization(scan);
    }
    static void AfterPreviousLocalization(LidarLoc& loc) { loc.evidence_ever_localized_ = true; }
    static Recorder::Stats Stats(LidarLoc& loc) { return loc.match_failure_recorder_->GetStats(); }
    static void Stop(LidarLoc& loc) { loc.match_failure_recorder_->Stop(); }
};
}  // namespace lightning::loc

namespace {
void Require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
void WaitSaved(lightning::loc::LidarLoc& loc, std::uint64_t expected) {
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
    while (std::chrono::steady_clock::now() < deadline) {
        if (lightning::loc::LidarLocRegressionAccess::Stats(loc).saved >= expected) return;
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    throw std::runtime_error("real rejection hook did not save expected evidence");
}
}

int main(int argc, char** argv) {
    const bool exporting = argc == 3 && std::string(argv[1]) == "--export-fixture";
    if (argc != 1 && !exporting) return 2;
    const fs::path folder = exporting ? fs::absolute(argv[2]) :
        fs::temp_directory_path() / ("lightning_relocalization_evidence_" + std::to_string(getpid()));
    if (fs::exists(folder)) {
        std::cerr << "refusing to reuse fixture directory: " << folder << '\n';
        return 2;
    }
    try {
        fs::create_directories(folder);
        using namespace lightning;
        using Access = loc::LidarLocRegressionAccess;
        CloudPtr cloud(new PointCloudType);
        for (int x = -5; x <= 5; ++x) for (int y = -5; y <= 5; ++y) {
            PointType p{};
            p.x = x * .2f; p.y = y * .2f; p.z = .2f;
            p.lidar_id = 0; p.data[3] = 1;
            cloud->push_back(p);
        }
        loc::LidarLoc localization;
        auto* backend = Access::Prepare(localization, cloud, folder);
        Require(!Access::Attempt(localization, cloud), "no candidate must remain unlocalized");
        WaitSaved(localization, 1);
        for (int i = 0; i < 100; ++i)
            Require(!Access::Attempt(localization, cloud), "repeat no-candidate unexpectedly accepted");
        Require(Access::Stats(localization).admitted == 1, "repeated no-candidate wrote more point clouds");
        Access::AfterPreviousLocalization(localization);
        backend->candidate = true;
        Require(!Access::Attempt(localization, cloud), "out-of-map candidate must be rejected by real precheck");
        WaitSaved(localization, 2);
        backend->ready = false;
        Require(!Access::Attempt(localization, cloud), "unavailable backend unexpectedly accepted");
        WaitSaved(localization, 3);
        for (int i = 0; i < 30; ++i) {
            backend->ready = i % 2 == 0;
            backend->candidate = false;
            Access::Attempt(localization, cloud);
        }
        Access::Stop(localization);
        const auto stats = Access::Stats(localization);
        Require(stats.admitted == 3 && stats.saved == 3, "continuous failure episode exceeded three packages");
        Require(stats.writer_errors == 0, "fixture writer errors");
        bool saw_startup = false, saw_global = false, saw_unavailable = false;
        for (const auto& entry : fs::directory_iterator(folder / "evidence")) {
            if (!entry.is_directory() || !fs::exists(entry.path() / "summary.json")) continue;
            const auto summary = YAML::LoadFile((entry.path() / "summary.json").string());
            Require(summary["synthetic_fixture"].as<bool>(), "synthetic hook test must be labelled");
            Require(summary["scan"]["saved_points"].as<int>() == 121, "hook lost original scan points");
            const auto stage = summary["stage"].as<std::string>();
            if (stage == "search") {
                saw_startup = summary["phase"].as<std::string>() == "startup";
                Require(summary["poses"]["guess"].IsNull() && summary["poses"]["ndt"].IsNull(),
                        "no candidate must not fabricate a map pose");
            } else if (stage == "map_precheck") {
                saw_global = summary["phase"].as<std::string>() == "global_relocalization";
                Require(summary["poses"]["ndt"].IsNull(), "unattempted refinement has no pose");
                Require(summary["poses"]["guess"][3].as<double>() == 100, "saved candidate pose mismatch");
                Require(summary["criterion"]["evaluated"].as<bool>(), "missing actual precheck result");
                Require(summary["criterion"]["overlap"].as<double>() == 0, "out-of-map overlap must be zero");
                Require(fs::file_size(entry.path() / "map_local.bin") == 0, "out-of-map crop should be empty");
            } else if (stage == "unavailable") saw_unavailable = true;
        }
        Require(saw_startup && saw_global && saw_unavailable, "missing real hook phase coverage");
        std::cout << "relocalization evidence hooks passed: startup search, global precheck, unavailable, episode cap\n";
        if (exporting) std::cout << "SYNTHETIC_HOOK_FIXTURE=" << folder << '\n';
        else fs::remove_all(folder);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAILED: " << e.what() << "; retained evidence=" << folder << '\n';
        return 1;
    }
}
