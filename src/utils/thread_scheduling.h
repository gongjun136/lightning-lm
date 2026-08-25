#pragma once

#include <string>
#include <string_view>
#include <vector>

namespace lightning::threading {

// Parse Linux taskset-style CPU lists such as "1-5,7-11". An empty value
// means unrestricted and produces an empty vector.
bool ParseCpuList(std::string_view text, std::vector<int>& cpus,
                  std::string* error = nullptr);

std::string FormatCpuList(const std::vector<int>& cpus);

// Configure only the calling Linux thread. A nice floor of zero leaves its
// scheduling priority unchanged. A non-empty CPU list must be a subset of the
// process mask inherited by the thread; requesting a CPU outside that mask is
// treated as an error instead of silently weakening the requested isolation.
bool ConfigureCurrentThread(const std::string& name, int nice_floor,
                            const std::vector<int>& cpus,
                            std::string* error = nullptr);

}  // namespace lightning::threading
