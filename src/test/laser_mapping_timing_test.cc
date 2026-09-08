#include "core/lio/laser_mapping.h"

#include <cmath>
#include <iostream>
#include <initializer_list>
#include <limits>
#include <stdexcept>

namespace lightning {
class LaserMappingTimingTestPeer {
 public:
    static IMUPtr Imu(double timestamp) {
        auto imu = std::make_shared<IMU>();
        imu->timestamp = timestamp;
        imu->linear_acceleration = Vec3d(1.0, 0.0, 9.81);
        return imu;
    }
    static NavState Rebuild(LaserMapping& mapping, double epoch, double tail,
                            std::initializer_list<double> buffered) {
        NavState state;
        state.timestamp_ = epoch;
        mapping.kf_.ChangeX(state);
        mapping.kf_.SetPropagateVelocity(true);
        mapping.measures_.imu_ = {Imu(tail)};
        mapping.imu_buffer_.clear();
        for (double stamp : buffered) mapping.imu_buffer_.push_back(Imu(stamp));
        mapping.RebuildHighFrequencyState();
        return mapping.kf_imu_.GetX();
    }
    static void InitializeImu(LaserMapping& mapping) {
        ImuProcess::InitializationOptions options;
        options.min_duration = 0.1;
        options.min_samples = 20;
        mapping.p_imu_->SetInitializationOptions(options);
        MeasureGroup init;
        for (int i = 0; i <= 20; ++i) {
            auto sample = Imu(0.01 + 0.01 * i);
            sample->linear_acceleration = Vec3d(0.0, 0.0, 9.81);
            init.imu_.push_back(sample);
        }
        CloudPtr cloud;
        mapping.p_imu_->Process(init, mapping.kf_, cloud);
    }
    static std::size_t BufferedImus(const LaserMapping& mapping) { return mapping.imu_buffer_.size(); }
};
}  // namespace lightning

void Require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void TestScanEpoch(double tail) {
    using namespace lightning;
    ImuProcess imu;
    ImuProcess::InitializationOptions options;
    options.min_duration = 0.1;
    options.min_samples = 20;
    imu.SetInitializationOptions(options);
    imu.SetAccCov(Vec3d::Ones());
    imu.SetGyrCov(Vec3d::Ones());
    ESKF filter;
    MeasureGroup init;
    for (int i = 0; i <= 20; ++i) {
        auto sample = LaserMappingTimingTestPeer::Imu(10.0 + 0.01 * i);
        sample->linear_acceleration = Vec3d(0.0, 0.0, 9.81);
        init.imu_.push_back(sample);
    }
    CloudPtr output;
    imu.Process(init, filter, output);
    Require(imu.IsIMUInited(), "test IMU initialization failed");
    auto state = filter.GetX();
    state.timestamp_ = 10.205;
    filter.ChangeX(state);
    MeasureGroup scan;
    scan.lidar_begin_time_ = 10.2;
    scan.lidar_end_time_ = 10.3;
    scan.scan_.reset(new PointCloudType());
    scan.imu_ = {LaserMappingTimingTestPeer::Imu(10.210), LaserMappingTimingTestPeer::Imu(tail)};
    imu.Process(scan, filter, output);
    Require(std::abs(filter.GetX().timestamp_ - 10.3) < 1e-12,
            "scan prediction must integrate from state epoch and stop at scan end");
}

int main() {
    using namespace lightning;
    LaserMapping mapping;
    const auto first = LaserMappingTimingTestPeer::Rebuild(mapping, 1.000, 0.995, {1.100});
    Require(std::abs(first.timestamp_ - 1.100) < 1e-12,
            "rebuild must start at copied filter epoch, not the frame-tail IMU");
    Require(std::abs(first.vel_.x() - 0.100) < 1e-12,
            "already integrated frame-tail interval must not add velocity twice");
    const auto next = LaserMappingTimingTestPeer::Rebuild(mapping, 1.100, 1.099, {1.101});
    Require(next.timestamp_ > first.timestamp_, "monotonic IMU must not produce a reconstructed time rollback");
    const auto duplicates = LaserMappingTimingTestPeer::Rebuild(mapping, 2.000, 1.995,
                                                               {1.999, 2.000, 2.010, 2.010, 2.009, 2.020});
    Require(std::abs(duplicates.timestamp_ - 2.020) < 1e-12,
            "covered, duplicate and out-of-order samples must not change final epoch");
    Require(std::abs(duplicates.vel_.x() - 0.020) < 1e-12,
            "only uncovered positive time intervals may be integrated");
    const auto empty = LaserMappingTimingTestPeer::Rebuild(mapping, 3.000, 2.995, {});
    Require(empty.timestamp_ == 3.000, "empty replay must preserve copied state epoch");
    TestScanEpoch(10.295);
    TestScanEpoch(10.305);

    LaserMappingTimingTestPeer::InitializeImu(mapping);
    LaserMappingTimingTestPeer::Rebuild(mapping, 4.000, 3.995, {4.020});
    mapping.ProcessIMU(LaserMappingTimingTestPeer::Imu(4.030));
    const auto incremental = mapping.GetIMUState();
    Require(std::abs(incremental.timestamp_ - 4.030) < 1e-12 &&
                std::abs(incremental.vel_.x() - 0.030) < 1e-12,
            "incremental IMU prediction must continue from the reconstructed epoch");
    const auto buffer_size = LaserMappingTimingTestPeer::BufferedImus(mapping);
    mapping.ProcessIMU(LaserMappingTimingTestPeer::Imu(4.030));
    mapping.ProcessIMU(LaserMappingTimingTestPeer::Imu(4.025));
    mapping.ProcessIMU(LaserMappingTimingTestPeer::Imu(std::numeric_limits<double>::quiet_NaN()));
    Require(mapping.GetIMUState().timestamp_ == incremental.timestamp_ &&
                LaserMappingTimingTestPeer::BufferedImus(mapping) == buffer_size,
            "duplicate, rollback and invalid IMU inputs must not clear or advance the bridge");
    ESKF filter;
    const auto covariance = filter.GetP();
    const auto noise = ESKF::ProcessNoiseType::Identity();
    Require(!filter.PredictTo(-1.0, noise, Vec3d::Zero(), Vec3d::Zero()) &&
                !filter.PredictTo(0.0, noise, Vec3d::Zero(), Vec3d::Zero()) &&
                !filter.PredictTo(std::numeric_limits<double>::quiet_NaN(), noise,
                                   Vec3d::Zero(), Vec3d::Zero()) &&
                (filter.GetP() - covariance).norm() == 0.0,
            "rejected prediction must preserve covariance as well as state time");
    std::cout << "laser mapping timing tests passed\n";
}
