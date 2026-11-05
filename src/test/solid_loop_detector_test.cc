#include <cmath>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "core/backend/backend_pipeline.h"

namespace {
using namespace lightning;

CloudPtr Scene() {
    auto cloud = std::make_shared<PointCloudType>();
    auto add = [&](double x, double y, double z) {
        PointType p; p.x = x; p.y = y; p.z = z; p.intensity = 1;
        cloud->push_back(p);
    };
    for (int x = -30; x <= 30; ++x)
        for (int y = -30; y <= 30; ++y) add(x * .2, y * .2, -1.0);
    for (int a = -30; a <= 30; ++a) {
        for (int h = -5; h <= 20; ++h) {
            add(-6, a * .2, h * .2); add(a * .2, -6, h * .2);
            add(6, a * .2, h * .2);
            if (a < 5) add(a * .2, 6, h * .2);
        }
    }
    for (int a = 0; a < 72; ++a)
        for (int h = -5; h <= 20; ++h)
            add(2 + .45 * std::cos(a * M_PI / 36),
                -1.5 + .45 * std::sin(a * M_PI / 36), h * .2);
    return cloud;
}

Keyframe::Ptr Frame(unsigned long id, const SE3& truth, const SE3& odometry) {
    auto cloud = Scene();
    for (auto& p : cloud->points) {
        const Vec3d local = truth.inverse() * Vec3d(p.x, p.y, p.z);
        p.x = local.x(); p.y = local.y(); p.z = local.z();
    }
    NavState state;
    state.timestamp_ = id * 40.;
    state.SetPose(odometry);
    return std::make_shared<Keyframe>(id, cloud, state);
}

bool Require(bool condition, const char* message) {
    if (!condition) std::cerr << message << '\n';
    return condition;
}
}  // namespace

int main() {
    using namespace lightning;
    using namespace lightning::backend;
    BtcLoopDetectorOptions options;
    options.retrieval = LoopRetrieval::kSolid;
    options.descriptor_submap_size = 1;
    options.max_points_per_submap = 30000;
    options.downsample_leaf_size = .15;
    options.min_loop_score = .8;
    options.solid_skip_near_descriptors = 1;
    options.solid_min_time_separation_sec = 30;
    options.confirmation_count = 1;
    options.loop_cooldown_descriptors = 0;
    options.allow_degenerate_plane_icp = false;
    options.plane_icp_min_observability = 3;
    options.max_drift_ratio = 1;
    options.min_optimization_translation = .1;
    // 故意让 BTC 三角形不可用，验证 SOLiD 检索不依赖 BTC 描述子。
    options.descriptor.descriptor_min_len_ = 10000;
    const Quatd mount(AngAxisd(M_PI / 6, Vec3d::UnitY()));
    const SE3 first(mount, Vec3d::Zero());
    const SE3 far(mount, Vec3d(10, 0, 0));
    const SE3 revisited(Quatd(AngAxisd(M_PI / 3, Vec3d::UnitZ())) * mount, Vec3d::Zero());
    const SE3 drifted(revisited.unit_quaternion(), Vec3d(.4, 0, 0));
    auto history = Frame(0, first, first);
    auto middle = Frame(1, far, far);
    auto current = Frame(2, revisited, drifted);
    BtcLoopDetector detector(options);
    const auto initial = detector.AddKeyframe(history, SE3());
    const auto recent = detector.AddKeyframe(middle, SE3());
    const auto loop = detector.AddKeyframe(current, SE3());
    if (!Require(initial && !initial->candidate_found && recent && !recent->candidate_found,
                 "SOLiD must exclude recent history") ||
        !Require(loop && loop->candidate_source == "solid" && loop->accepted && loop->optimization_warranted,
                 "SOLiD revisit failed geometry verification") ||
        !Require(loop->history_descriptor_id == 0 && loop->descriptor_count == 1 &&
                 loop->plane_icp_matches >= 20, "unexpected SOLiD candidate")) return 1;
    const SE3 error = (first.inverse() * revisited).inverse() * loop->T_history_lidar_current_lidar;
    if (!Require(error.translation().norm() < .05 && error.so3().log().norm() < .03,
                 "pitched-LiDAR candidate-from-query transform direction is wrong")) return 2;
    auto negative_options = options;
    negative_options.min_loop_score = -1;
    BtcLoopDetector negative(negative_options);
    negative.AddKeyframe(history, SE3()); negative.AddKeyframe(middle, SE3());
    const auto wrong_scene = negative.AddKeyframe(
        Frame(2, SE3(mount, Vec3d(-30, 0, 0)), drifted), SE3());
    if (!Require(wrong_scene && wrong_scene->candidate_found && !wrong_scene->accepted &&
                 wrong_scene->rejection_reason == "plane_icp_rejected",
                 "a SOLiD score candidate bypassed point-cloud geometry")) return 2;
    const auto directory = std::filesystem::temp_directory_path() /
                           ("solid_loop_test_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directory(directory);
    const auto config = directory / "config.yaml";
    {
        std::ofstream file(config);
        file << "fasterlio:\n  extrinsic_T: [0, 0, 0]\n  extrinsic_R: [1, 0, 0, 0, 1, 0, 0, 0, 1]\n"
             << "backend:\n  mode: ba_solid_hba\n  local_ba: {enabled: false}\n  hba: {enabled: false}\n"
             << "  pose_graph: {loop_translation_noise: 0.2}\n  solid:\n"
             << "    descriptor_submap_size: 1\n    downsample_leaf_size: 0.15\n"
             << "    min_loop_score: 0.8\n    skip_near_descriptors: 1\n    confirmation_count: 1\n"
             << "    loop_cooldown_descriptors: 0\n    allow_degenerate_plane_icp: false\n"
             << "    plane_icp_min_observability: 3\n    max_drift_ratio: 1\n"
             << "    min_optimization_translation: 0.1\n";
    }
    BackendPipeline pipeline;
    if (!Require(ReadBackendMode(config.string()) == BackendMode::kBaSolidHba &&
                 pipeline.Init(config.string(), false), "SOLiD backend dispatch failed")) return 3;
    pipeline.AddKeyframe(history); pipeline.AddKeyframe(middle); pipeline.AddKeyframe(current);
    pipeline.WaitUntilIdle(false);
    const auto runtime = pipeline.GetRuntimeSummary();
    if (!Require(runtime.btc_descriptors == 3 && runtime.loops_applied == 1 && runtime.loops_graph_inliers == 1 &&
                 current->GetOptPose().translation().norm() < .35,
                 "SOLiD loop did not correct the pose graph") ||
        !Require(pipeline.SaveDiagnostics(directory.string()) &&
                 std::filesystem::exists(directory / "solid_loop_candidates.csv") &&
                 !std::filesystem::exists(directory / "btc_loop_candidates.csv"),
                 "SOLiD diagnostics are mislabeled")) return 4;
    pipeline.Shutdown();
    std::ifstream input(config);
    std::string rejected_config((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    const std::string graph_config = "pose_graph: {loop_translation_noise: 0.2}";
    rejected_config.replace(rejected_config.find(graph_config), graph_config.size(),
                            "pose_graph: {loop_translation_noise: 0.2, loop_outlier_chi2: 0.000000000001}");
    const auto rejected_path = directory / "reject_all.yaml";
    { std::ofstream file(rejected_path); file << rejected_config; }
    BackendPipeline rejected;
    auto rejected_history = Frame(0, first, first);
    auto rejected_middle = Frame(1, far, far);
    auto rejected_current = Frame(2, revisited, drifted);
    if (!Require(rejected.Init(rejected_path.string(), false), "outlier regression setup failed")) return 5;
    rejected.AddKeyframe(rejected_history); rejected.AddKeyframe(rejected_middle); rejected.AddKeyframe(rejected_current);
    const auto rejected_runtime = rejected.GetRuntimeSummary();
    if (!Require(rejected_runtime.loops_applied == 1 && rejected_runtime.loops_graph_inliers == 0 &&
                 (rejected_current->GetOptPose().inverse() * drifted).log().norm() < 1e-12 &&
                 (rejected_middle->GetOptPose().inverse() * far).log().norm() < 1e-12,
                 "zero graph inliers still changed keyframe poses")) return 5;
    rejected.Shutdown();
    // 本测试创建的临时目录仅含自身输出。
    std::filesystem::remove_all(directory);
    return 0;
}
