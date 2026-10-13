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

/// Bind an accepted old-from-new transform to the exact source/reference geometry.
/// already_transformed records a physical conversion and disables its reapplication.
bool SaveAlignmentManifest(const std::string& directory, const std::string& old_directory,
                           const SE3& T_old_new, const std::string& report, double quality,
                           bool already_transformed, std::string& error);

/// Require a complete tiled-map manifest, verify geometry, and return the authoritative
/// output.map_frame/output.fixed_map_transform block. Missing/unknown states fail closed.
bool LoadOutputConfig(const std::string& directory, YAML::Node& output, std::string& error);

/// Streaming SHA-256; throws on unreadable files or digest failure.
std::string FileSha256(const std::string& filename);

}  // namespace lightning::map_manifest
