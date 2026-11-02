#include "utils/input_timing_monitor.h"
#include <cstdlib>
#include <iostream>
#include <limits>

using Monitor = lightning::diagnostics::InputTimingMonitor;
auto At(double sec) {
    return Monitor::Clock::time_point(std::chrono::duration_cast<Monitor::Clock::duration>(std::chrono::duration<double>(sec)));
}
void Require(bool good, const char* message) {
    if (!good) { std::cerr << message << '\n'; std::exit(1); }
}
int main() {
    Monitor imu;
    Monitor::Record r;
    imu.Observe(100, At(1));
    Require(imu.Pop(r) && r.kind == Monitor::Kind::kStart, "initial source sample");
    imu.Poll(At(1.2));
    Require(imu.Pop(r) && r.kind == Monitor::Kind::kSilence && r.source_gap == -1, "silence has no invented source gap");
    imu.Observe(100.005, At(1.21));
    Require(imu.Pop(r) && r.kind == Monitor::Kind::kInterval && r.source_gap < .006 && r.arrival_gap > .2,
            "continuous source delayed at callback");
    imu.Observe(100.3, At(1.215));
    Require(imu.Pop(r) && r.source_gap > .29 && r.arrival_gap < .006, "source gap also detected during fast catchup");
    imu.Observe(100.3, At(1.22));
    Require(imu.Pop(r) && r.source_gap == 0, "duplicate is not fresh source progress");
    imu.Observe(100.2, At(1.225));
    Require(imu.Pop(r) && r.source_gap < 0, "source rollback retained");
    imu.Observe(std::numeric_limits<double>::quiet_NaN(), At(1.23));
    Require(imu.Pop(r) && std::isnan(r.stamp) && r.arrival_gap < .006, "invalid CAN header still records real arrival");
    Require(imu.Summary().invalid == 1 && imu.Summary().duplicates == 1 && imu.Summary().rollbacks == 1,
            "quality counters retained");

    Monitor radar(.2);
    for (int i = 0; i < 100; ++i) radar.Observe(1000 + i*.1001, At(10+i*.1001));
    Require(radar.Pop(r) && r.kind == Monitor::Kind::kStart && !radar.Pop(r), "10Hz radar jitter does not flood log");
    Require(radar.Summary().arrival_gaps == 0 && radar.Summary().source_gaps == 0, "radar threshold is independent");

    Monitor normal;
    for (int i = 0; i < 1000000; ++i) normal.Observe(1000+i*.005, At(100+i*.005));
    Require(normal.Pop(r) && !normal.Pop(r) && normal.Summary().sample == 1000000,
            "one million regular samples need only the start event");
    Monitor bounded;
    for (int i = 0; i < 100; ++i) bounded.Observe(100, At(i*.005));
    int retained = 0;
    while (bounded.Pop(r)) ++retained;
    Require(retained == 16 && bounded.Summary().lost == 84 && bounded.Summary().duplicates == 99,
            "overflow is bounded and explicitly counted");
    std::cout << "input timing monitor OK\n";
}
