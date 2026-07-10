#include "core/lio/multi_lidar_fusion.h"
#include "core/lightning_math.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

using lightning::CloudPtr;
using lightning::FusedLidarFrame;
using lightning::Mat3d;
using lightning::MultiLidarConfig;
using lightning::MultiLidarFrameAssembler;
using lightning::MultiLidarSensorConfig;
using lightning::PointCloudType;
using lightning::PointType;
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

}  // namespace

int main() {
    TestCompleteFrameAndTransform();
    TestThreeCycleLateFrameIsNotPrematurelyReleased();
    TestTimeoutDegradeAndRecovery();
    TestFlushAndDuplicate();
    TestPairwiseToleranceAndCommonPhaseDrift();
    TestDownsampleKeepsRealIds();
    std::cout << "multi_lidar_fusion_test passed" << std::endl;
    return 0;
}
