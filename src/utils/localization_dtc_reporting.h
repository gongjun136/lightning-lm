#pragma once
#include "core/localization/localization_dtc_monitor.h"
#include "utils/functional_safety_heartbeat.h"

namespace lightning::functional_safety {
inline void UpdateLocalizationFaults(HeartbeatPublisher& heartbeat, const loc::dtc::Snapshot& snapshot) {
    std::vector<uint32_t> failed, passed;
    for (const auto& check : snapshot.observations) {
        if (check.condition == loc::dtc::kFail) failed.push_back(check.code);
        if (check.condition == loc::dtc::kPass) passed.push_back(check.code);
    }
    heartbeat.UpdateFaults(failed, passed);
}
}
