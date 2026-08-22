#include "core/lio/multi_lidar_fusion.h"
#include "core/lightning_math.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

using lightning::CloudPtr;
using lightning::AdaptiveLidarLoadController;
using lightning::AdaptiveLidarSelection;
using lightning::FusedLidarFrame;
using lightning::Mat3d;
using lightning::MultiLidarConfig;
using lightning::MultiLidarFrameAssembler;
using lightning::MultiLidarSensorConfig;
using lightning::PointCloudType;
using lightning::PointType;
using lightning::SelfPointFilterConfig;
using lightning::Vec3d;

void Require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << std::endl;
        std::exit(1);
    }
}

CloudPtr MakeCloud(double x, double time_ms = 99.0) {
    CloudPtr cloud(new PointCloudType);
    PointType point;
    point.x = x;
    point.y = 0.0;
    point.z = 0.0;
    point.intensity = 1.0;
    point.time = time_ms;
    cloud->push_back(point);
    return cloud;
}

MultiLidarConfig MakeConfig() {
    MultiLidarConfig config;
    config.enabled = true;
    config.primary_lidar_id = 0;
    config.frame_period = 0.1;
    config.match_tolerance = 0.002;
    config.reorder_window = 0.5;
    config.min_lidars = 1;
    for (int id = 0; id < 4; ++id) {
        MultiLidarSensorConfig sensor;
        sensor.id = id;
        sensor.lidar_topic = "/lidar" + std::to_string(id);
        sensor.imu_topic = "/imu" + std::to_string(id);
        sensor.t_lidar_to_primary = Vec3d(id, 0.0, 0.0);
        config.lidars.push_back(sensor);
    }
    return config;
}

void TestCompleteFrameAndTransform() {
    MultiLidarFrameAssembler assembler(MakeConfig());
    Require(assembler.AddCloud(0, 100.0000, MakeCloud(1.0)), "add lidar0");
    Require(assembler.AddCloud(1, 100.0004, MakeCloud(1.0)), "add lidar1");
    Require(assembler.AddCloud(2, 100.0002, MakeCloud(1.0)), "add lidar2");
    Require(assembler.AddCloud(3, 100.0001, MakeCloud(1.0)), "add lidar3");
    FusedLidarFrame frame;
    Require(assembler.PopReady(frame), "complete frame is ready");
    Require(frame.cloud->size() == 4, "complete frame has four points");
    Require(!frame.stats.partial, "complete frame is not partial");
    for (const auto& point : frame.cloud->points) {
        Require(std::abs(point.x - (1.0 + point.lidar_id)) < 1e-6, "point transformed to primary");
    }
    Require(std::abs(frame.stats.end_time - 100.0994) < 1e-6, "end timestamp includes header offset");
}

void TestThreeCycleLateFrameIsNotPrematurelyReleased() {
    MultiLidarFrameAssembler assembler(MakeConfig());
    for (int id = 0; id < 3; ++id) assembler.AddCloud(id, 200.0, MakeCloud(id));
    for (int cycle = 1; cycle <= 3; ++cycle) {
        for (int id = 0; id < 4; ++id) assembler.AddCloud(id, 200.0 + 0.1 * cycle, MakeCloud(id));
    }
    FusedLidarFrame frame;
    Require(!assembler.PopReady(frame), "old partial frame waits inside reorder window");
    Require(assembler.AddCloud(3, 200.0003, MakeCloud(3)), "late lidar is accepted");
    Require(assembler.PopReady(frame), "completed old frame is released first");
    Require(!frame.stats.partial && std::abs(frame.stats.begin_time - 200.0) < 1e-6,
            "late frame remains complete and ordered");
}

void TestTimeoutDegradeAndRecovery() {
    MultiLidarFrameAssembler assembler(MakeConfig());
    for (int id = 0; id < 3; ++id) assembler.AddCloud(id, 300.0, MakeCloud(id));
    for (int cycle = 1; cycle <= 6; ++cycle) {
        for (int id = 0; id < 4; ++id) assembler.AddCloud(id, 300.0 + 0.1 * cycle, MakeCloud(id));
    }
    FusedLidarFrame frame;
    Require(assembler.PopReady(frame), "expired partial frame is released");
    Require(frame.stats.partial && frame.stats.missing_lidar_ids.size() == 1 &&
                frame.stats.missing_lidar_ids.front() == 3,
            "missing lidar is recorded");
    Require(assembler.PopReady(frame), "next complete frame follows partial");
    Require(!frame.stats.partial, "recovered full frame is complete");
    Require(!assembler.AddCloud(3, 300.0002, MakeCloud(3)), "late point for emitted bucket is dropped");
    Require(assembler.LateDropCount() == 1, "late drop is counted");
}

void TestFlushAndDuplicate() {
    MultiLidarFrameAssembler assembler(MakeConfig());
    Require(assembler.AddCloud(0, 400.0, MakeCloud(0)), "first cloud accepted");
    Require(!assembler.AddCloud(0, 400.0001, MakeCloud(0)), "duplicate sensor frame rejected");
    Require(assembler.DuplicateDropCount() == 1, "duplicate is counted");
    assembler.Flush();
    FusedLidarFrame frame;
    Require(assembler.PopReady(frame), "EOF flush emits partial frame");
    Require(frame.stats.partial && frame.cloud->size() == 1, "flushed frame content");
}

void TestPairwiseToleranceAndCommonPhaseDrift() {
    MultiLidarFrameAssembler assembler(MakeConfig());
    Require(assembler.AddCloud(0, 500.0030, MakeCloud(0)), "common phase drift anchor accepted");
    Require(assembler.AddCloud(1, 500.0049, MakeCloud(1)), "spread below tolerance accepted");
    Require(!assembler.AddCloud(2, 500.0051, MakeCloud(2)), "spread above tolerance rejected");
    Require(assembler.ToleranceDropCount() == 1, "tolerance rejection is counted");
    Require(assembler.AddCloud(2, 500.0040, MakeCloud(2)), "replacement inside tolerance accepted");
    Require(assembler.AddCloud(3, 500.0035, MakeCloud(3)), "complete drifted group accepted");
    FusedLidarFrame frame;
    Require(assembler.PopReady(frame), "drifted complete group is emitted");
    Require(!frame.stats.partial, "drifted group remains complete");
    Require(std::abs(lightning::math::ToSec(frame.cloud->header.stamp) - frame.stats.begin_time) < 1e-9,
            "PCL header uses nanoseconds expected by math::ToSec");
}

void TestDownsampleKeepsRealIds() {
    CloudPtr cloud(new PointCloudType);
    auto a = MakeCloud(0.01)->front();
    auto b = MakeCloud(0.02)->front();
    a.lidar_id = 0;
    b.lidar_id = 3;
    cloud->push_back(a);
    cloud->push_back(b);
    const auto filtered = lightning::DownsamplePreservingSource(cloud, 1.0);
    Require(filtered->size() == 1, "same voxel has one representative");
    Require(filtered->front().lidar_id == 0 || filtered->front().lidar_id == 3, "source id is not averaged");
}

void TestBodyAlignedSelfPointFilter() {
    const YAML::Node root = YAML::Load(R"(
self_point_filter:
  enabled: true
  front: 4.2218
  back: 4.4582
  left: 1.543
  right: 1.543
  bottom: 3.577
  top: 0.05
  padding: 0.10
)");
    SelfPointFilterConfig config;
    std::string error;
    Require(lightning::LoadSelfPointFilterConfig(root, config, &error), "load self-point filter");
    Require((config.min_body - Vec3d(-4.5582, -1.643, -3.677)).norm() < 1e-9,
            "padding expands negative box faces");
    Require((config.max_body - Vec3d(4.3218, 1.643, 0.15)).norm() < 1e-9,
            "padding expands positive box faces");

    PointCloudType cloud;
    PointType inside;
    inside.x = 1.0F;
    inside.y = 0.0F;
    inside.z = 0.0F;
    inside.lidar_id = 2;
    inside.time = 12.0;
    cloud.push_back(inside);
    PointType outside = inside;
    outside.x = 5.0F;
    outside.lidar_id = 3;
    cloud.push_back(outside);

    const double half_pi = 0.5 * std::acos(-1.0);
    const Mat3d R_primary_to_body =
        Eigen::AngleAxisd(half_pi, Vec3d::UnitZ()).toRotationMatrix();
    const std::size_t removed = lightning::FilterSelfPoints(cloud, config, R_primary_to_body);
    Require(removed == 1, "point inside rotated body box is removed");
    Require(cloud.size() == 1 && cloud.front().lidar_id == 3,
            "point outside body box and source id are preserved");

    const YAML::Node export_root = YAML::Load(R"(
map_export:
  self_point_filter:
    enabled: true
    min_body: [1.5, -4.8, -0.75]
    max_body: [5.4, 4.8, 1.0]
    padding: 0.10
)");
    SelfPointFilterConfig export_config;
    Require(lightning::LoadSelfPointFilterConfig(
                export_root["map_export"], export_config, &error),
            "load explicit map-export self-point filter");
    Require((export_config.min_body - Vec3d(1.4, -4.9, -0.85)).norm() < 1e-9,
            "explicit export minimum is padded");
    Require((export_config.max_body - Vec3d(5.5, 4.9, 1.1)).norm() < 1e-9,
            "explicit export maximum is padded");
}

lightning::MultiLidarFrameStats MakeFrameStats(std::vector<int> ids) {
    lightning::MultiLidarFrameStats stats;
    stats.present_lidar_ids = std::move(ids);
    for (const int id : stats.present_lidar_ids) {
        stats.points_by_lidar[id] = 1000 - static_cast<std::size_t>(id) * 100;
    }
    return stats;
}

void TestAdaptiveLoadOrderAndQualityGate() {
    MultiLidarConfig config = MakeConfig();
    auto& load = config.adaptive_load;
    load.enabled = true;
    load.tracking_min_lidars = 3;
    load.relocalization_min_lidars = 3;
    load.cloud_publish_min_lidars = 3;
    load.degrade_consecutive_frames = 2;
    load.recover_consecutive_frames = 3;
    load.point_strides = {1, 2, 3};

    AdaptiveLidarLoadController controller;
    controller.Reset(config);
    controller.SetLocalizationGood(true);
    const auto full = MakeFrameStats({0, 1, 2, 3});
    auto selection = controller.Select(full);
    Require(selection.lidar_ids.size() == 4 && selection.point_stride == 1,
            "healthy baseline uses every lidar and point");

    controller.Observe(0.17, 0.21, true);
    controller.Observe(0.17, 0.21, true);
    selection = controller.Select(full);
    Require(selection.lidar_ids.size() == 4 && selection.point_stride == 2,
            "first overload step increases point stride before removing lidar");
    controller.Observe(0.17, 0.21, true);
    controller.Observe(0.17, 0.21, true);
    selection = controller.Select(full);
    Require(selection.lidar_ids.size() == 4 && selection.point_stride == 3,
            "second overload step keeps every lidar with stronger sampling");
    controller.Observe(0.17, 0.21, true);
    controller.Observe(0.17, 0.21, true);
    selection = controller.Select(full);
    Require(selection.lidar_ids.size() == 3 && selection.lidar_ids.front() == 0,
            "only after point reduction does overload remove a secondary lidar");

    controller.Observe(0.05, 0.05, false);
    controller.Observe(0.05, 0.05, false);
    Require(controller.DegradationStep() == 2,
            "bad tracking restores data instead of degrading further");

    for (int i = 0; i < 8; ++i) {
        controller.Observe(0.17, 0.21, true);
    }
    Require(controller.TargetLidarCount() == 3,
            "sustained overload never reduces localization below three lidars");
    controller.SetLocalizationGood(false);
    Require(controller.TargetLidarCount() == 3,
            "initialization and relocalization clamp the minimum lidar count");
}

void TestAdaptiveCloudPublicationAndPointSelection() {
    MultiLidarConfig config = MakeConfig();
    config.adaptive_load.enabled = true;
    config.adaptive_load.relocalization_min_lidars = 3;
    config.adaptive_load.cloud_publish_min_lidars = 3;
    AdaptiveLidarLoadController controller;
    controller.Reset(config);
    Require(controller.CanPublishCloud(MakeFrameStats({0, 1, 2})),
            "three-lidar cloud containing primary is publishable");
    Require(!controller.CanPublishCloud(MakeFrameStats({1, 2, 3})),
            "cloud without primary is suppressed");
    Require(!controller.CanPublishCloud(MakeFrameStats({0, 1})),
            "cloud below publication minimum is suppressed");
    Require(controller.IsHardStale(0.31) && !controller.IsHardStale(0.29),
            "hard stale deadline is explicit");

    CloudPtr cloud(new PointCloudType);
    for (int id = 0; id < 3; ++id) {
        for (int i = 0; i < 4; ++i) {
            auto point = MakeCloud(static_cast<double>(i))->front();
            point.lidar_id = static_cast<std::uint8_t>(id);
            cloud->push_back(point);
        }
    }
    AdaptiveLidarSelection selection;
    selection.lidar_ids = {0, 2};
    selection.point_stride = 2;
    const auto selected = lightning::SelectLidarPoints(cloud, selection);
    Require(selected->size() == 4, "stride is applied independently to each selected lidar");
    for (const auto& point : selected->points) {
        Require(point.lidar_id == 0 || point.lidar_id == 2, "unselected lidar points are excluded");
    }
}

void TestFormalSanyAdaptiveConfigs() {
    const auto check_config = [](const std::string& relative_path,
                                 int expected_lidar_count,
                                 int expected_relocalization_min) {
        const YAML::Node root = YAML::LoadFile(std::string(ROOT_DIR) + relative_path);
        MultiLidarConfig config;
        std::string error;
        Require(lightning::LoadMultiLidarConfig(root, config, &error),
                "formal SANY multi-lidar config parses");
        Require(config.enabled && config.adaptive_load.enabled,
                "formal SANY adaptive load is enabled");
        Require(static_cast<int>(config.lidars.size()) == expected_lidar_count,
                "formal SANY lidar count");
        Require(config.primary_lidar_id == 0 && config.reorder_window == 0.1 &&
                    config.min_lidars == 3,
                "formal SANY primary, three-lidar minimum and bounded reorder window");
        Require(config.adaptive_load.tracking_min_lidars == 3 &&
                    config.adaptive_load.relocalization_min_lidars ==
                        expected_relocalization_min &&
                    config.adaptive_load.cloud_publish_min_lidars == 3 &&
                    config.adaptive_load.cloud_publish_require_primary,
                "formal SANY localization and publication minima");
        Require(config.adaptive_load.target_latency_sec == 0.2 &&
                    config.adaptive_load.hard_latency_sec == 0.3 &&
                    config.adaptive_load.point_strides == std::vector<int>({1, 2, 3}),
                "formal SANY latency and point-stride policy");
        const YAML::Node system = root["system"];
        Require(root["fasterlio"]["skip_lidar_num"].as<int>() == 0 && system &&
                    system["enable_lidar_loc_skip"] &&
                    !system["enable_lidar_loc_skip"].as<bool>() &&
                    system["pub_tf"] && !system["pub_tf"].as<bool>() &&
                    system["enable_wheel_speed_dr_observation"].as<bool>() &&
                    system["localization_output_max_lidar_age_sec"].as<double>() == 0.5 &&
                    system["wheel_speed_dr_max_age_sec"].as<double>() == 0.25 &&
                    system["wheel_speed_dr_max_velocity_step_mps"].as<double>() == 0.35,
                "formal SANY fixed-skip, output freshness and wheel-speed DR gates are explicit");
    };

    check_config(
        "config/reproduction/multi_lidar/sany_3livox/"
        "sany_3lidar_localization_solid.yaml",
        3, 3);
    check_config(
        "config/reproduction/multi_lidar/sany_4livox/"
        "sany_4lidar_localization_solid.yaml",
        4, 3);
}

}  // namespace

int main() {
    TestCompleteFrameAndTransform();
    TestThreeCycleLateFrameIsNotPrematurelyReleased();
    TestTimeoutDegradeAndRecovery();
    TestFlushAndDuplicate();
    TestPairwiseToleranceAndCommonPhaseDrift();
    TestDownsampleKeepsRealIds();
    TestBodyAlignedSelfPointFilter();
    TestAdaptiveLoadOrderAndQualityGate();
    TestAdaptiveCloudPublicationAndPointSelection();
    TestFormalSanyAdaptiveConfigs();
    std::cout << "multi_lidar_fusion_test passed" << std::endl;
    return 0;
}
