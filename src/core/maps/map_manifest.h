#pragma once

#include <string>
#include <yaml-cpp/yaml.h>

#include "common/eigen_types.h"

namespace lightning::map_manifest {

inline constexpr char kFilename[] = "map_manifest.yaml";

struct ExportOptions {
    std::string site_id;
    std::string map_id;
    std::string alignment_status = "unknown";
};

/// New maps cannot inherit a previous map's alignment from the localization YAML.
/// Only unknown/not_required may be declared before a new alignment is computed.
bool ReadExportOptions(const YAML::Node& root, ExportOptions& options, std::string& error);

/// Write portable geometry checksums, ground-frame provenance and output policy.
/// Call after map_frame.yaml and all map geometry have been finalized.
bool SaveExportManifest(const std::string& directory, const ExportOptions& options,
                        std::string& error);

/**
 * @brief 将调用方接受的新图相对旧图的对齐关系绑定到确切的地图几何文件。
 * T_old_new 将新图坐标转换到旧图坐标；already_transformed 记录该转换是否已作用于文件。
 *
 * @par 主要流程
 * 1. <b>取得源图和参考图身份</b>：验证已有清单，或为没有清单的地图生成身份资料。
 * 2. <b>记录对齐变换及溯源</b>：校验位姿，记录参考 global.pcd 哈希、变换是否已作用于文件以及配准报告。
 * 3. <b>写入输出策略与清单</b>：按文件是否已变换确定坐标身份和输出补偿，保存清单。
 * @note already_transformed 只记录调用方已完成的文件转换；本函数不修改点云，也不自行判定配准质量是否合格。
 */
bool SaveAlignmentManifest(const std::string& directory, const std::string& old_directory,
                           const SE3& T_old_new, const std::string& report, double quality,
                           bool already_transformed, std::string& error);

/**
 * @brief 验证分块地图清单和几何文件，返回权威输出配置。
 * 返回 output.map_frame/output.fixed_map_transform；清单缺失或对齐状态 unknown 时拒绝。
 *
 * @par 主要流程
 * 1. <b>验证地图身份与可用状态</b>：校验几何清单，要求 tiled 地图且对齐状态已明确。
 * 2. <b>检查变换方向和应用状态</b>：验证源/目标身份、单位、位姿和文件坐标系，确定是否仍需输出补偿。
 * 3. <b>交叉核对输出声明</b>：核对 output 的开关、方向及各位姿表示与对齐结果一致。
 * 4. <b>返回权威输出配置</b>：生成输出块并写入清单路径及地图身份；异常转为 error 和 false。
 * @note 此处验证声明的一致性和文件身份，不重新计算配准；unknown 状态拒绝启动，已作用于文件的变换不会再次用于输出。
 */
bool LoadOutputConfig(const std::string& directory, YAML::Node& output, std::string& error);

/// Streaming SHA-256; throws on unreadable files or digest failure.
std::string FileSha256(const std::string& filename);

}  // namespace lightning::map_manifest
