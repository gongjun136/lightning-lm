#include "utils/thread_scheduling.h"

#include <algorithm>
#include <cerrno>
#include <charconv>
#include <cstring>
#include <string>

#include <pthread.h>
#include <sched.h>
#include <sys/resource.h>
#include <sys/syscall.h>
#include <unistd.h>

namespace lightning::threading {
namespace {

constexpr int kMaxNice = 19;

std::string_view Trim(std::string_view text) {
    while (!text.empty() && (text.front() == ' ' || text.front() == '\t')) {
        text.remove_prefix(1);
    }
    while (!text.empty() && (text.back() == ' ' || text.back() == '\t')) {
        text.remove_suffix(1);
    }
    return text;
}

bool ParseCpu(std::string_view text, int& cpu, std::string* error) {
    text = Trim(text);
    if (text.empty()) {
        if (error) *error = "CPU index is empty";
        return false;
    }
    int parsed = 0;
    const auto result =
        std::from_chars(text.data(), text.data() + text.size(), parsed);
    if (result.ec != std::errc() || result.ptr != text.data() + text.size() ||
        parsed < 0 || parsed >= CPU_SETSIZE) {
        if (error) {
            *error = "CPU index must be an integer in [0, " +
                     std::to_string(CPU_SETSIZE - 1) + "]";
        }
        return false;
    }
    cpu = parsed;
    return true;
}

void SetError(std::string* error, const std::string& message) {
    if (error && error->empty()) *error = message;
}

}  // namespace

bool ParseCpuList(std::string_view text, std::vector<int>& cpus,
                  std::string* error) {
    std::vector<int> parsed;
    text = Trim(text);
    if (text.empty()) {
        cpus.clear();
        return true;
    }

    std::size_t offset = 0;
    while (offset <= text.size()) {
        const std::size_t comma = text.find(',', offset);
        std::string_view item = text.substr(
            offset, comma == std::string_view::npos ? text.size() - offset
                                                    : comma - offset);
        item = Trim(item);
        if (item.empty()) {
            if (error) *error = "CPU list contains an empty item";
            return false;
        }

        const std::size_t dash = item.find('-');
        int first = 0;
        int last = 0;
        if (dash == std::string_view::npos) {
            if (!ParseCpu(item, first, error)) return false;
            last = first;
        } else {
            if (item.find('-', dash + 1) != std::string_view::npos ||
                !ParseCpu(item.substr(0, dash), first, error) ||
                !ParseCpu(item.substr(dash + 1), last, error)) {
                return false;
            }
            if (first > last) {
                if (error) *error = "CPU range start must not exceed its end";
                return false;
            }
        }
        for (int cpu = first; cpu <= last; ++cpu) parsed.push_back(cpu);

        if (comma == std::string_view::npos) break;
        offset = comma + 1;
    }

    std::sort(parsed.begin(), parsed.end());
    parsed.erase(std::unique(parsed.begin(), parsed.end()), parsed.end());
    cpus = std::move(parsed);
    return true;
}

std::string FormatCpuList(const std::vector<int>& cpus) {
    if (cpus.empty()) return "unrestricted";
    std::string output;
    for (std::size_t index = 0; index < cpus.size();) {
        std::size_t end = index;
        while (end + 1 < cpus.size() && cpus[end + 1] == cpus[end] + 1) ++end;
        if (!output.empty()) output += ',';
        output += std::to_string(cpus[index]);
        if (end != index) output += '-' + std::to_string(cpus[end]);
        index = end + 1;
    }
    return output;
}

bool ConfigureCurrentThread(const std::string& name, int nice_floor,
                            const std::vector<int>& cpus, std::string* error) {
    if (error) error->clear();
    bool success = true;

    std::string thread_name = name.substr(0, 15);
    const int name_status = pthread_setname_np(pthread_self(), thread_name.c_str());
    if (name_status != 0) {
        SetError(error, "pthread_setname_np failed: " +
                            std::string(std::strerror(name_status)));
        success = false;
    }

    if (!cpus.empty()) {
        cpu_set_t inherited;
        CPU_ZERO(&inherited);
        if (sched_getaffinity(0, sizeof(inherited), &inherited) != 0) {
            SetError(error, "sched_getaffinity failed: " +
                                std::string(std::strerror(errno)));
            success = false;
        } else {
            cpu_set_t requested;
            CPU_ZERO(&requested);
            bool outside_inherited_mask = false;
            for (const int cpu : cpus) {
                if (!CPU_ISSET(cpu, &inherited)) outside_inherited_mask = true;
                CPU_SET(cpu, &requested);
            }
            if (outside_inherited_mask) {
                SetError(error,
                         "requested CPU affinity is outside the inherited process mask");
                success = false;
            } else if (sched_setaffinity(0, sizeof(requested), &requested) != 0) {
                SetError(error, "sched_setaffinity failed: " +
                                    std::string(std::strerror(errno)));
                success = false;
            }
        }
    }

    if (nice_floor > 0) {
        const auto thread_id = static_cast<id_t>(syscall(SYS_gettid));
        errno = 0;
        const int current_nice = getpriority(PRIO_PROCESS, thread_id);
        if (current_nice == -1 && errno != 0) {
            SetError(error, "getpriority failed: " +
                                std::string(std::strerror(errno)));
            success = false;
        } else {
            const int target_nice = std::min(kMaxNice, std::max(current_nice, nice_floor));
            if (target_nice != current_nice &&
                setpriority(PRIO_PROCESS, thread_id, target_nice) != 0) {
                SetError(error, "setpriority failed: " +
                                    std::string(std::strerror(errno)));
                success = false;
            }
        }
    }
    return success;
}

}  // namespace lightning::threading
