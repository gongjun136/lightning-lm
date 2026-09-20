#pragma once

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "common/point_def.h"
#include "core/maps/map_frame.h"

namespace lightning::navigation_map {

struct RaycastFrame {
    CloudPtr cloud;
    SE3 T_map_primary;
};

using SensorOrigins = std::map<std::uint8_t, Vec3d>;

struct ExportResult {
    int width = 0;
    int height = 0;
    double origin_x = 0.0;
    double origin_y = 0.0;
    std::size_t occupied_cells = 0;
    std::size_t free_cells = 0;
    std::size_t unknown_cells = 0;
    std::size_t ray_count = 0;
};

bool ExportPgmAndYaml(const CloudPtr& map,
                      const std::vector<RaycastFrame>& raycast_frames,
                      const SensorOrigins& sensor_origins_primary,
                      const std::string& map_directory,
                      const map_frame::ExportOptions& options,
                      ExportResult& result, std::string& error);

}  // namespace lightning::navigation_map
