#include "core/localization/radar_input_frame_monitor.h"
#include <iostream>
#include <stdexcept>
using Monitor = lightning::loc::RadarInputFrameMonitor;
using Clock = Monitor::Clock;
static Clock::time_point At(int ms) { return Clock::time_point(std::chrono::milliseconds(ms)); }
static void Check(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
static lightning::MultiLidarConfig Config() {
    lightning::MultiLidarConfig c;
    c.enabled = true; c.frame_period = .1; c.match_tolerance = .02; c.reorder_window = .1; c.min_lidars = 1;
    for (int id = 0; id < 3; ++id) { lightning::MultiLidarSensorConfig sensor; sensor.id = id; c.lidars.push_back(sensor); }
    return c;
}
static void Frame(Monitor& monitor, int ms) {
    for (int id = 0; id < 3; ++id) monitor.ObserveInput(id, 10 + ms * .001, true, At(ms));
}
int main() {
    try {
        const auto config = Config();
        Monitor::Watchdog::Config timing; timing.startup_grace_sec = 1;
        lightning::MultiLidarFrameAssembler metadata;
        metadata.Reset(config, lightning::MultiLidarFrameAssembler::Mode::kInputStamps);
        for (int id = 0; id < 3; ++id) metadata.AddStamp(id, 10);
        lightning::FusedLidarFrame grouped;
        Check(metadata.PopReady(grouped) && !grouped.cloud && grouped.stats.present_lidar_ids.size() == 3,
              "input grouping must not merge or allocate point data");

        Monitor healthy(config, timing); healthy.Start(At(0));
        // No LIO consumer, preprocessing, matching or algorithm progress is run.
        // Arbitrarily delayed/missing algorithm results cannot become radar 2.
        for (int ms = 100; ms <= 2000; ms += 100) {
            Frame(healthy, ms);
            Check(!healthy.Evaluate(At(ms)).fault, "healthy radar input must not depend on algorithm completion");
        }
        Check(healthy.Evaluate(At(2299)).warning && !healthy.Evaluate(At(2299)).fault,
              "200 ms input-group warning is not yet radar 2");
        Check(healthy.Evaluate(At(2300)).fault, "300 ms with no eligible input group is radar 2");
        for (int i = 0; i < 3; ++i) {
            Frame(healthy, 2400 + i * 100);
            Check(healthy.Evaluate(At(2400 + i * 100)).fault == (i < 2), "three live input groups recover radar 2");
        }
        for (int ms = 2700; ms <= 3300; ms += 100)
            for (int id = 0; id < 3; ++id) healthy.ObserveInput(id, 12.6, true, At(ms));
        Check(healthy.Evaluate(At(3300)).fault, "fresh callbacks with frozen source stamps cannot mask radar failure");

        Monitor partial(config, timing); partial.Start(At(0)); Frame(partial, 100);
        for (int ms = 200; ms <= 1500; ms += 100) {
            partial.ObserveInput(0, 10 + ms * .001, true, At(ms));
            Check(!partial.Evaluate(At(ms)).fault, "usable primary-only input is allowed by min_lidars=1");
        }
        Monitor missing_primary(config, timing); missing_primary.Start(At(0)); Frame(missing_primary, 100);
        for (int ms = 200; ms <= 1000; ms += 100)
            for (int id : {1, 2}) missing_primary.ObserveInput(id, 10 + ms * .001, true, At(ms));
        Check(missing_primary.Evaluate(At(1000)).fault, "auxiliary traffic cannot replace a required primary");

        auto strict = config; strict.min_lidars = 3;
        Monitor misaligned(strict, timing); misaligned.Start(At(0)); Frame(misaligned, 100);
        for (int ms = 200; ms <= 1000; ms += 100) {
            misaligned.ObserveInput(0, 10 + ms * .001, true, At(ms));
            misaligned.ObserveInput(1, 10 + ms * .001 + .03, true, At(ms));
            misaligned.ObserveInput(2, 10 + ms * .001 + .04, true, At(ms));
        }
        Check(misaligned.Evaluate(At(1000)).fault, "required inputs outside the shared 20 ms tolerance cannot form a group");
        Monitor empty(config, timing); empty.Start(At(0)); Frame(empty, 100);
        for (int ms = 200; ms <= 1000; ms += 100)
            for (int id = 0; id < 3; ++id) empty.ObserveInput(id, 10 + ms * .001, false, At(ms));
        Check(empty.Evaluate(At(1000)).fault, "empty raw input cannot refresh grouping progress");

        auto delayed = config; delayed.reorder_window = 2;
        Monitor old_input(delayed, timing); old_input.Start(At(0));
        for (int ms : {100, 200, 300}) old_input.ObserveInput(0, 10 + ms * .001, true, At(ms));
        Check(old_input.Evaluate(At(1000)).fault, "no first eligible group after finite startup grace");
        old_input.ObserveInput(1, 20, true, At(1100));
        Check(old_input.Evaluate(At(1100)).fault && old_input.Evaluate(At(1100)).accepted == 0,
              "late auxiliary input must not revive old primary metadata");
        auto single = config; single.enabled = false;
        Monitor one(single, timing); one.Start(At(0));
        for (int ms = 100; ms <= 1000; ms += 100) one.ObserveInput(0, 10 + ms * .001, true, At(ms));
        Check(!one.Evaluate(At(1000)).fault, "single-radar raw input has no algorithm-progress dependency");
        std::cout << "radar_input_frame_monitor_test passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
