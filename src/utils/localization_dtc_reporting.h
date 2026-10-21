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
    // The monitor tracks translation/yaw failures separately. This flag is
    // allowed only after required yaw evidence passes and no failure remains;
    // a CAN dropout alone can never clear a confirmed translation failure.
    if (snapshot.optional_can_unavailable) passed.push_back(loc::dtc::kConsistency);
    heartbeat.UpdateFaults(failed, passed);
}
}
