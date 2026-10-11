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
    CloudPtr cloud;  ///< 该帧用于射线投射的主雷达坐标系点云。
    SE3 T_map_primary;  ///< 主雷达坐标系到导出地图坐标系的该帧位姿。
};

using SensorOrigins = std::map<std::uint8_t, Vec3d>;

struct ExportResult {
    int width = 0;  ///< 导出栅格宽度，单位像素。
    int height = 0;  ///< 导出栅格高度，单位像素。
    double origin_x = 0.0;  ///< 导航栅格原点的地图 x 坐标（m）。
    double origin_y = 0.0;  ///< 导航栅格原点的地图 y 坐标（m）。
    std::size_t occupied_cells = 0;  ///< 导出栅格中障碍单元的数量。
    std::size_t free_cells = 0;  ///< 导出栅格中射线确认自由单元的数量。
    std::size_t unknown_cells = 0;  ///< 导出栅格中未观测单元的数量。
    std::size_t ray_count = 0;  ///< 用于清理自由空间的射线数量。
};

bool ExportPgmAndYaml(const CloudPtr& map,
                      const std::vector<RaycastFrame>& raycast_frames,
                      const SensorOrigins& sensor_origins_primary,
                      const std::string& map_directory,
                      const map_frame::ExportOptions& options,
                      ExportResult& result, std::string& error);

}  // namespace lightning::navigation_map
