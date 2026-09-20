#include "utils/thread_scheduling.h"

#include <algorithm>
#include <cerrno>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <pthread.h>
#include <sched.h>
#include <sys/resource.h>
#include <sys/syscall.h>
#include <unistd.h>

namespace {

void Require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

int main() {
    using lightning::threading::ConfigureCurrentThread;
    using lightning::threading::FormatCpuList;
    using lightning::threading::ParseCpuList;

    std::vector<int> cpus;
    std::string error;
    Require(ParseCpuList("1-3,5,3", cpus, &error), "valid CPU list parses");
    Require(cpus == std::vector<int>({1, 2, 3, 5}),
            "CPU list is expanded, sorted and deduplicated");
    Require(FormatCpuList(cpus) == "1-3,5", "CPU list is compactly formatted");
    Require(!ParseCpuList("3-1", cpus, &error), "descending range is rejected");
    Require(!ParseCpuList("1,,2", cpus, &error), "empty list item is rejected");

    bool thread_success = false;
    std::thread worker([&]() {
        cpu_set_t inherited;
        CPU_ZERO(&inherited);
        Require(sched_getaffinity(0, sizeof(inherited), &inherited) == 0,
                "thread affinity can be read");
        int selected_cpu = -1;
        for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
            if (CPU_ISSET(cpu, &inherited)) {
                selected_cpu = cpu;
                break;
            }
        }
        Require(selected_cpu >= 0, "at least one CPU is available");

        errno = 0;
        const auto thread_id = static_cast<id_t>(syscall(SYS_gettid));
        const int inherited_nice = getpriority(PRIO_PROCESS, thread_id);
        Require(!(inherited_nice == -1 && errno != 0), "thread nice can be read");
        const int requested_nice = std::min(19, std::max(1, inherited_nice + 1));
        Require(ConfigureCurrentThread("solid_sched_t", requested_nice,
                                       {selected_cpu}, &error),
                "thread scheduling policy applies");

        char name[16] = {};
        Require(pthread_getname_np(pthread_self(), name, sizeof(name)) == 0,
                "thread name can be read");
        Require(std::string(name) == "solid_sched_t", "thread name is applied");

        cpu_set_t applied;
        CPU_ZERO(&applied);
        Require(sched_getaffinity(0, sizeof(applied), &applied) == 0,
                "applied affinity can be read");
        Require(CPU_COUNT(&applied) == 1 && CPU_ISSET(selected_cpu, &applied),
                "thread is restricted to the requested CPU");

        const int outside_cpu = selected_cpu == 0 ? 1 : 0;
        Require(!ConfigureCurrentThread("solid_sched_t", 0, {outside_cpu}, &error),
                "thread cannot escape its inherited process CPU mask");
        CPU_ZERO(&applied);
        Require(sched_getaffinity(0, sizeof(applied), &applied) == 0,
                "affinity remains readable after rejected widening");
        Require(CPU_COUNT(&applied) == 1 && CPU_ISSET(selected_cpu, &applied),
                "rejected widening leaves the inherited mask unchanged");

        errno = 0;
        const int applied_nice = getpriority(PRIO_PROCESS, thread_id);
        Require(!(applied_nice == -1 && errno != 0), "applied nice can be read");
        Require(applied_nice >= requested_nice,
                "thread priority is no higher than the configured nice floor");
        thread_success = true;
    });
    worker.join();
    Require(thread_success, "thread scheduling test completes");
    return 0;
}
