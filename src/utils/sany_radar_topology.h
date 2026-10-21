#pragma once

#include <array>
#include <set>
#include <stdexcept>
#include <string>
#include <yaml-cpp/yaml.h>

namespace lightning::functional_safety {

/// Use the immutable run snapshot resolved from vehicle_config.yaml, not topic discovery.
inline std::array<bool, 4> RequiredSanyRadarRoles(const YAML::Node& root) {
    if (root["sany_vehicle_config"])
        throw std::invalid_argument("resolve SANY vehicle configuration before starting localization");
    const auto vehicle = root["vehicle_configuration"];
    // Preserve the legacy contract for non-SANY configurations without a vehicle snapshot.
    if (!vehicle) return {{true, true, true, true}};
    const auto layout = vehicle["LIDAR"].as<std::string>();
    const int count = layout == "MASTER AND SLAVE" ? 3 : layout == "MASTER" ? 4 : 0;
    if (count == 0 || vehicle["lidar_count"].as<int>() != count)
        throw std::invalid_argument("inconsistent SANY LIDAR layout/count");
    const auto multi = root["multi_lidar"];
    if (!multi || !multi["enabled"].as<bool>() || multi["primary_lidar_id"].as<int>() != 0)
        throw std::invalid_argument("SANY radar diagnostics require enabled multi_lidar with primary 0");
    const auto topics = multi["topics"];
    if (!topics.IsMap() || topics.size() != static_cast<size_t>(2 * count))
        throw std::invalid_argument("SANY radar topics do not match vehicle layout");
    std::set<std::string> lidar_topics;
    for (int id = 0; id < count; ++id) {
        const auto lidar = topics["lidar_" + std::to_string(id)].as<std::string>();
        const auto imu = topics["imu_" + std::to_string(id)].as<std::string>();
        if (lidar.empty() || imu.empty() || !lidar_topics.insert(lidar).second)
            throw std::invalid_argument("SANY radar topics must be nonempty and distinct");
    }
    // Localization IDs: front=0,left=1,right=2,rear=3. DTC roles: front,rear,left,right.
    return {{true, count == 4, true, true}};
}

}  // namespace lightning::functional_safety
