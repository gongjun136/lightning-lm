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

// Initialize once, before creating threads. Optional environment lists:
// LIGHTNING_LM_BACKGROUND_CPUS, LIGHTNING_LM_SENSOR_CPUS,
// LIGHTNING_LM_OUTPUT_CPUS, LIGHTNING_LM_EXECUTOR_CPUS, LIGHTNING_LM_MATCH_CPUS,
// LIGHTNING_LM_LIO_WORKER_CPUS.
// All lists must be within the startup mask; background must be disjoint
// from sensor/output/executor. Match defaults to background for old launchers;
// an explicit match mask may share output/executor, but not the sensor core.
// LIO workers default to background for old launchers. An explicit worker mask
// may share output/executor, but must exclude sensor, background and matching.
// Threads initially inherit the background mask.
bool InitializeOnlineThreadPolicy(std::string* error = nullptr);
enum class OnlineThreadRole { kBackground, kSensor, kOutput, kExecutor, kMatch, kLioWorker };
// Explicitly restore only a validated role mask from the startup policy.
// Throws on an application/readback error, never silently runs unpinned.
void ApplyOnlineThreadRole(OnlineThreadRole role, const std::string& name);
// Call inside LIO OpenMP loops; thread zero retains the sensor role, additional
// threads use the worker role instead of inheriting the sensor core.
void ApplyOnlineParallelWorkerPolicy();

}  // namespace lightning::threading
