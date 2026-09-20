#include "utils/causal_trace.h"
#include <string>
#include <thread>
#include <vector>

int main(int argc, char** argv) {
    const bool enabled = argc > 1 && std::string(argv[1]) == "on";
    if (lightning::profiling::CausalTraceEnabled() != enabled) return 1;
    if (!enabled) return 0;
    std::vector<std::thread> threads;
    for (int i=0;i<4;++i) threads.emplace_back([] {
        for (int j=0;j<1000;++j) {
            lightning::profiling::CausalEvent event;
            event.kind="can_update";event.source="test";event.before_v=1;event.after_v=2;
            lightning::profiling::RecordCausalEvent(event);
        }
    });
    for (auto& thread:threads) thread.join();
    return 0;
}
