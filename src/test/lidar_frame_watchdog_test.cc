#include "utils/lidar_frame_watchdog.h"
#include <iostream>
#include <stdexcept>
using Watchdog = lightning::functional_safety::LidarFrameWatchdog;
static Watchdog::Clock::time_point At(int ms) { return Watchdog::Clock::time_point(std::chrono::milliseconds(ms)); }
static void Require(bool ok, const char* text) { if (!ok) throw std::runtime_error(text); }
int main() {
    try {
        Watchdog w(Watchdog::Config{});
        Require(!w.Evaluate(At(9000)).fault, "unarmed watchdog must not fault");
        w.Start(At(1000));
        Require(!w.Evaluate(At(8999)).fault && w.Evaluate(At(9000)).fault, "finite no-input startup grace");
        w.Start(At(9000));
        Require(w.Evaluate(At(9000)).fault, "repeated start cannot hide fault");
        for (int i = 0; i < 3; ++i) {
            Require(w.Observe(10 + i * .1, At(9100 + i * 100), true, At(9100 + i * 100)), "fresh recovery frame");
            Require(w.Evaluate(At(9100 + i * 100)).fault == (i < 2), "three fresh frames required");
        }
        Require(!w.Evaluate(At(9499)).warning && w.Evaluate(At(9500)).warning, "200 ms warning boundary");
        Require(!w.Evaluate(At(9599)).fault && w.Evaluate(At(9600)).fault, "300 ms independent silence detection");
        Require(!w.Observe(10.3, At(9400), true, At(9800)), "old queued frame cannot revive watchdog");
        Require(!w.Observe(10.2, At(9800), true, At(9800)), "repeated source stamp cannot revive watchdog");
        Require(!w.Observe(10.3, At(9801), false, At(9801)), "ineligible partial frame cannot revive watchdog");
        Require(!w.Observe(10.3, At(9999), true, At(9801)), "future arrival is rejected");
        Require(w.Observe(10.3, At(9900), true, At(9900)), "first valid recovery");
        Require(!w.Observe(10.4, At(9950), false, At(9950)), "invalid bucket interrupts recovery");
        Require(w.Observe(10.4, At(10000), true, At(10000)), "recovery restarts");
        Require(w.Evaluate(At(10000)).recovery_count == 1, "invalid frame resets streak");
        w.Observe(10.5, At(10250), true, At(10250));
        Require(w.Evaluate(At(10250)).recovery_count == 1, "irregular recovery interval resets streak");
        w.Observe(10.6, At(10350), true, At(10350));
        w.Observe(10.7, At(10450), true, At(10450));
        Require(!w.Evaluate(At(10450)).fault, "continuous fresh recovery succeeds");
        Watchdog batch(Watchdog::Config{});
        batch.Start(At(1000)); batch.Observe(1, At(1100), true, At(1100));
        Require(batch.Evaluate(At(1400)).fault, "batch fixture faults independently");
        batch.Observe(1.1, At(1410), true, At(1600));
        batch.Observe(1.2, At(1510), true, At(1600));
        batch.Observe(1.3, At(1590), true, At(1600));
        Require(batch.Evaluate(At(1600)).fault && batch.Evaluate(At(1600)).recovery_count == 1,
                "three cached arrivals processed together cannot clear a fault");
        batch.Observe(1.4, At(1700), true, At(1700)); batch.Observe(1.5, At(1800), true, At(1800));
        Require(!batch.Evaluate(At(1800)).fault, "live arrivals finish recovery after a cached batch");
        bool invalid = false;
        try { auto c = Watchdog::Config{}; c.timeout_sec = c.warning_sec; Watchdog bad(c); }
        catch (const std::invalid_argument&) { invalid = true; }
        Require(invalid, "invalid thresholds fail configuration");
        std::cout << "lidar_frame_watchdog_test PASS\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
