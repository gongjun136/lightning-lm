#include "core/localization/radar_input_frame_monitor.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
using Monitor = lightning::loc::RadarInputFrameMonitor;
using Clock = Monitor::Clock;
static Clock::time_point At(int ms) { return Clock::time_point(std::chrono::milliseconds(ms)); }
static void Check(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
static lightning::MultiLidarConfig Config() {
    lightning::MultiLidarConfig c;
    c.enabled = true; c.frame_period = .1; c.match_tolerance = .02; c.reorder_window = .1; c.min_lidars = 1;
    for (int id = 0; id < 3; ++id) { lightning::MultiLidarSensorConfig sensor; sensor.id = id;
        sensor.lidar_topic = "/lidar/" + std::to_string(id); c.lidars.push_back(sensor); }
    return c;
}
static void Frame(Monitor& monitor, int ms) {
    for (int id = 0; id < 3; ++id) monitor.ObserveInput(id, 10 + ms * .001, true, At(ms));
}
static void Dump(const char* name, Monitor& m, int ms, int sequence) {
    const auto state = m.Evaluate(At(ms));
    std::ostringstream text;
    m.WriteEvidence(text, m.SnapshotEvidence(state.accepted, At(ms), At(ms)));
    std::cout << "FIXTURE scenario=" << name << " SAFETY_DIAGNOSTIC sequence=" << sequence
        << " policy_version=localization-dtc-v9-match-grace raw_runstate=" << (state.fault ? 2 : 0)
        << " selected_runstate=" << (state.fault ? 2 : 0) << " published_runstate=" << (state.fault ? 2 : 0)
        << " strict_runstate=0 lidar_frame_started=1 lidar_frame_seen=" << state.seen
        << " lidar_frame_fault=" << state.fault << " lidar_frame_warning=" << state.warning
        << " lidar_frame_age_sec=" << state.age_sec << " lidar_frame_recovery_count=" << state.recovery_count
        << " radar_frame_timeout_sec=0.3 radar_frame_recovery_frames=1" << text.str() << '\n';
}
int main(int argc, char**) {
    try {
        const auto config = Config();
        Monitor::Watchdog::Config timing; timing.startup_grace_sec = 1; timing.recovery_frames = 3;
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

        Monitor::Watchdog::Config fast; fast.startup_grace_sec = 1;
        Monitor evidence(config, fast); evidence.Start(At(0));
        for (int id = 0; id < 3; ++id) evidence.ObserveInput(id, 10.1, true, At(100), 10.15);
        const auto first = evidence.Evaluate(At(110));
        auto snapshot = evidence.SnapshotEvidence(first.accepted, At(110), At(110));
        Check(snapshot.available && snapshot.consistent && snapshot.inputs.size() == 3, "all configured routes recorded");
        Check(std::abs(snapshot.inputs[0].lag_sec - .05) < 1e-10 &&
              std::abs(snapshot.inputs[0].callback_age_sec - .01) < 1e-10, "ROS lag and monotonic age use different clocks");
        for (int ms = 200; ms <= 600; ms += 100)
            for (int id : {1, 2}) evidence.ObserveInput(id, 10 + ms * .001, true, At(ms), 10 + ms * .001 + .05);
        auto stopped = evidence.Evaluate(At(600));
        snapshot = evidence.SnapshotEvidence(stopped.accepted, At(600), At(600));
        Check(stopped.fault && snapshot.inputs[0].callback_age_sec == .5 && snapshot.inputs[1].callback_age_sec == 0,
              "primary silence is distinguished from healthy auxiliary callbacks");
        Check(snapshot.pending && snapshot.oldest.missing_lidar_ids == std::vector<int>{0}, "pending bucket exposes missing primary");
        if (argc > 1) Dump("primary_silence", evidence, 600, 1);
        evidence.ObserveInput(0, 10.6, true, At(610), 10.65);
        // 早桶尚可能等待源时间推进；主雷达来包不等于已经输出合格组帧。
        // 下一次源时间推进释放一个完整的 10.6 桶，10.7 桶仍仅含主雷达。
        evidence.ObserveInput(0, 10.7, true, At(700), 10.75);
        const auto recovered = evidence.Evaluate(At(700));
        Check(!recovered.fault && recovered.accepted == stopped.accepted + 1,
              "one eligible multi-lidar frame clears the default fault");
        auto changed = evidence.SnapshotEvidence(stopped.accepted, At(600), At(700));
        Check(!changed.consistent, "concurrent recovery cannot be described as the triggering snapshot");
        if (argc > 1) Dump("recovered", evidence, 700, 2);

        Monitor delayed_aux(config, fast); delayed_aux.Start(At(0)); Frame(delayed_aux, 100);
        for (int ms = 200; ms <= 600; ms += 100)
            for (int id : {0, 1}) delayed_aux.ObserveInput(id, 10 + ms * .001, true, At(ms), 10 + ms * .001 + .05);
        delayed_aux.ObserveInput(2, 10.2, true, At(610), 10.65);
        snapshot = delayed_aux.SnapshotEvidence(delayed_aux.Evaluate(At(610)).accepted, At(610), At(610));
        Check(!delayed_aux.Evaluate(At(610)).fault && snapshot.inputs[2].assembly == std::string("late_bucket"),
              "late auxiliary input need not cause fault with usable primary-only grouping");
        Check(std::abs(snapshot.inputs[2].lag_sec - .45) < 1e-9, "constant-clock delayed input retains receive lag");
        if (argc > 1) Dump("late_auxiliary_without_fault", delayed_aux, 610, 3);

        Monitor skew(strict, fast); skew.Start(At(0)); Frame(skew, 100);
        for (int ms = 200; ms <= 600; ms += 100) {
            skew.ObserveInput(0, 10 + ms * .001, true, At(ms), 10 + ms * .001 + .05);
            skew.ObserveInput(1, 10 + ms * .001 + .03, true, At(ms), 10 + ms * .001 + .05);
            skew.ObserveInput(2, 10 + ms * .001 + .04, true, At(ms), 10 + ms * .001 + .05);
        }
        snapshot = skew.SnapshotEvidence(skew.Evaluate(At(600)).accepted, At(600), At(600));
        Check(skew.Evaluate(At(600)).fault && snapshot.inputs[1].assembly == std::string("tolerance") &&
              snapshot.inputs[1].callback_age_sec == 0, "fresh callbacks can fail timestamp grouping rather than receive delivery");
        if (argc > 1) Dump("timestamp_tolerance", skew, 600, 4);
        if (argc > 1) Dump("all_silent", skew, 1000, 5);
        skew.ObserveInput(1, 11, true, At(1010), 10.5);
        snapshot = skew.SnapshotEvidence(skew.Evaluate(At(1010)).accepted, At(1010), At(1010));
        Check(snapshot.inputs[1].lag_sec == -.5, "negative clock difference must not be clamped to healthy zero");
        if (argc > 1) Dump("clock_mismatch", skew, 1010, 6);
        std::cout << "radar_input_frame_monitor_test passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
