#include "core/maps/map_manifest.h"

#include <openssl/evp.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace lightning::map_manifest {
namespace {
namespace fs = std::filesystem;

// [map-manifest-geometry-scope]
bool GeometryFile(const fs::path& filename) {
    const std::string name = filename.string();
    if (name == "global.pcd" || name == "index.txt" || name == "map_frame.yaml") return true;
    const std::string stem = filename.stem().string();
    return filename.extension() == ".pcd" && !stem.empty() &&
           std::all_of(stem.begin(), stem.end(), [](unsigned char c) { return std::isdigit(c); });
}
// [map-manifest-geometry-scope]

std::string DefaultMapId(const fs::path& directory) {
    const fs::path path = fs::weakly_canonical(directory);
    if (path.filename() == "new_map" && path.parent_path().filename() == "data") {
        return path.parent_path().parent_path().filename().string();
    }
    return path.filename().string();
}

YAML::Node PoseNode(const SE3& transform) {
    YAML::Node node;
    const auto t = transform.translation();
    const auto q = transform.unit_quaternion();
    node["translation_xyz"] = std::vector<double>{t.x(), t.y(), t.z()};
    node["quaternion_xyzw"] = std::vector<double>{q.x(), q.y(), q.z(), q.w()};
    const auto matrix = transform.matrix();
    for (int row = 0; row < 4; ++row) {
        node["matrix_4x4"].push_back(std::vector<double>{
            matrix(row, 0), matrix(row, 1), matrix(row, 2), matrix(row, 3)});
    }
    return node;
}

SE3 ReadPose(const YAML::Node& node) {
    const auto t = node["translation_xyz"].as<std::vector<double>>();
    const auto q = node["quaternion_xyzw"].as<std::vector<double>>();
    if (t.size() != 3 || q.size() != 4) throw std::runtime_error("invalid manifest transform dimensions");
    const Vec3d translation(t[0], t[1], t[2]);
    Quatd quaternion(q[3], q[0], q[1], q[2]);
    if (!translation.allFinite() || !quaternion.coeffs().allFinite() ||
        std::abs(quaternion.norm() - 1.0) > 1e-6) {
        throw std::runtime_error("manifest transform must be a finite unit-quaternion rigid transform");
    }
    quaternion.normalize();
    const SE3 transform(quaternion, translation);
    const auto rows = node["matrix_4x4"].as<std::vector<std::vector<double>>>();
    if (rows.size() != 4) throw std::runtime_error("manifest requires a 4x4 transform matrix");
    for (int row = 0; row < 4; ++row) {
        if (rows[row].size() != 4) throw std::runtime_error("manifest requires a 4x4 transform matrix");
        for (int col = 0; col < 4; ++col) {
            if (!std::isfinite(rows[row][col]) ||
                std::abs(rows[row][col] - transform.matrix()(row, col)) > 1e-9) {
                throw std::runtime_error("manifest matrix and translation/quaternion disagree");
            }
        }
    }
    return transform;
}

YAML::Node OutputNode(const SE3& transform, bool enabled) {
    // [map-manifest-effective-output]
    YAML::Node output;
    output["map_frame"] = "map";
    YAML::Node fixed = PoseNode(enabled ? transform : SE3());
    fixed.remove("matrix_4x4");
    fixed["enabled"] = enabled;
    fixed["convention"] = "target_from_localization";
    fixed["source_frame"] = "localization_map";
    fixed["target_frame"] = "map";
    output["fixed_map_transform"] = fixed;
    // [map-manifest-effective-output]
    return output;
}

void Write(const fs::path& directory, const YAML::Node& root) {
    YAML::Emitter emitter;
    emitter.SetDoublePrecision(17);
    emitter << root;
    if (!emitter.good()) throw std::runtime_error("failed to serialize map manifest");
    const fs::path temporary = directory / (std::string(kFilename) + ".tmp");
    std::ofstream stream(temporary);
    stream << emitter.c_str() << '\n';
    stream.close();
    if (!stream) throw std::runtime_error("failed to write map manifest: " + temporary.string());
    fs::rename(temporary, directory / kFilename);
}

YAML::Node MakeRoot(const fs::path& directory, const ExportOptions& options) {
    if (!fs::is_directory(directory)) throw std::runtime_error("map directory does not exist");
    if (!fs::exists(directory / "index.txt") && !fs::exists(directory / "global.pcd")) {
        throw std::runtime_error("map requires index.txt or global.pcd");
    }
    YAML::Node root;
    root["schema_version"] = 1;
    root["map_id"] = options.map_id.empty() ? DefaultMapId(directory) : options.map_id;
    root["site_id"] = options.site_id;
    root["map_format"] = fs::exists(directory / "index.txt") ? "tiled" : "global_pcd_only";
    root["coordinate_frame_id"] = root["map_id"].as<std::string>();
    std::vector<fs::path> files;
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.is_regular_file() && GeometryFile(entry.path().filename())) files.push_back(entry.path());
    }
    std::sort(files.begin(), files.end());
    for (const auto& file : files) root["files_sha256"][file.filename().string()] = FileSha256(file.string());
    const fs::path ground = directory / "map_frame.yaml";
    root["ground_normalization"]["applied"] = fs::exists(ground);
    if (fs::exists(ground)) root["ground_normalization"]["transform"] = YAML::LoadFile(ground.string())["map_frame"];
    YAML::Node alignment;
    alignment["status"] = options.alignment_status;
    alignment["source_map_id"] = root["map_id"].as<std::string>();
    alignment["target_map_id"] = root["map_id"].as<std::string>();
    alignment["convention"] = "target_from_localization";
    alignment["translation_unit"] = "m";
    alignment["applied_to_map_files"] = false;
    if (options.alignment_status == "unknown") alignment["needs_transform_to_old"] = YAML::Node(YAML::NodeType::Null);
    else alignment["needs_transform_to_old"] = false;
    alignment["transform"] = PoseNode(SE3());
    root["legacy_alignment"] = alignment;
    root["output"] = OutputNode(SE3(), false);
    return root;
}

YAML::Node ReadAndVerify(const fs::path& directory) {
    const fs::path filename = directory / kFilename;
    if (!fs::exists(filename)) throw std::runtime_error("required map_manifest.yaml is missing: " + filename.string());
    YAML::Node root = YAML::LoadFile(filename.string());
    if (root["schema_version"].as<int>() != 1 || root["map_id"].as<std::string>().empty()) {
        throw std::runtime_error("invalid map manifest schema or map_id");
    }
    const YAML::Node files = root["files_sha256"];
    if (!files.IsMap() || files.size() == 0 ||
        (root["map_format"].as<std::string>() == "tiled" && !files["index.txt"])) {
        throw std::runtime_error("manifest is missing geometry checksums");
    }
    // [map-manifest-file-verification]
    for (const auto& item : files) {
        const std::string name = item.first.as<std::string>();
        if (fs::path(name).filename().string() != name || !GeometryFile(name)) {
            throw std::runtime_error("invalid map geometry filename in manifest");
        }
        if (FileSha256((directory / name).string()) != item.second.as<std::string>()) {
            throw std::runtime_error("map geometry checksum mismatch: " + name);
        }
    }
    // A newly added or renamed tile must not silently escape the identity check.
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.is_regular_file() && GeometryFile(entry.path().filename()) &&
            !files[entry.path().filename().string()]) {
            throw std::runtime_error("map geometry is absent from manifest: " + entry.path().filename().string());
        }
    }
    // [map-manifest-file-verification]
    return root;
}
}  // namespace

std::string FileSha256(const std::string& filename) {
    std::ifstream stream(filename, std::ios::binary);
    if (!stream) throw std::runtime_error("cannot read map geometry: " + filename);
    std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> context(EVP_MD_CTX_new(), EVP_MD_CTX_free);
    if (!context || EVP_DigestInit_ex(context.get(), EVP_sha256(), nullptr) != 1) {
        throw std::runtime_error("cannot initialize SHA-256");
    }
    std::array<char, 1024 * 1024> buffer;
    while (stream) {
        stream.read(buffer.data(), buffer.size());
        if (EVP_DigestUpdate(context.get(), buffer.data(), stream.gcount()) != 1) throw std::runtime_error("SHA-256 update failed");
    }
    if (!stream.eof()) throw std::runtime_error("cannot finish reading map geometry: " + filename);
    std::array<unsigned char, EVP_MAX_MD_SIZE> digest;
    unsigned int length = 0;
    if (EVP_DigestFinal_ex(context.get(), digest.data(), &length) != 1) throw std::runtime_error("SHA-256 finalization failed");
    std::ostringstream text;
    text << std::hex << std::setfill('0');
    for (unsigned int index = 0; index < length; ++index) text << std::setw(2) << static_cast<unsigned int>(digest[index]);
    return text.str();
}

bool ReadExportOptions(const YAML::Node& root, ExportOptions& options, std::string& error) {
    try {
        const auto node = root["map_export"];
        if (node) {
            if (node["site_id"]) options.site_id = node["site_id"].as<std::string>();
            if (node["map_id"]) options.map_id = node["map_id"].as<std::string>();
            if (node["alignment_status"]) options.alignment_status = node["alignment_status"].as<std::string>();
        }
        if (options.alignment_status != "unknown" && options.alignment_status != "not_required") {
            throw std::runtime_error("map_export.alignment_status must be unknown or not_required; compute alignment for each new map");
        }
        return true;
    } catch (const std::exception& e) { error = e.what(); return false; }
}

bool SaveExportManifest(const std::string& directory, const ExportOptions& options, std::string& error) {
    try {
        if (options.alignment_status != "unknown" && options.alignment_status != "not_required") {
            throw std::runtime_error("new map alignment status must be unknown or not_required");
        }
        Write(directory, MakeRoot(directory, options));
        return true;
    } catch (const std::exception& e) { error = e.what(); return false; }
}

bool SaveAlignmentManifest(const std::string& directory, const std::string& old_directory,
                           const SE3& T_old_new, const std::string& report, double quality,
                           bool already_transformed, std::string& error) {
    try {
        YAML::Node root = fs::exists(fs::path(directory) / kFilename)
                              ? ReadAndVerify(directory) : MakeRoot(directory, ExportOptions{});
        const std::string old_id = fs::exists(fs::path(old_directory) / kFilename)
                                      ? ReadAndVerify(old_directory)["map_id"].as<std::string>()
                                      : DefaultMapId(old_directory);
        const auto pose = PoseNode(T_old_new);
        ReadPose(pose);
        YAML::Node alignment = root["legacy_alignment"];
        alignment["status"] = "aligned";
        alignment["target_map_id"] = old_id;
        alignment["target_global_pcd_sha256"] = FileSha256((fs::path(old_directory) / "global.pcd").string());
        // This records the caller's completed conversion; it does not rewrite point files.
        alignment["applied_to_map_files"] = already_transformed;
        alignment["needs_transform_to_old"] = !already_transformed;
        alignment["transform"] = pose;
        alignment["provenance"]["generated_by"] = "align_maps_offline";
        alignment["provenance"]["report"] = report;
        alignment["provenance"]["quality"] = quality;
        root["coordinate_frame_id"] = already_transformed ? old_id : root["map_id"].as<std::string>();
        root["output"] = OutputNode(T_old_new, !already_transformed);
        Write(directory, root);
        return true;
    } catch (const std::exception& e) { error = e.what(); return false; }
}

bool LoadOutputConfig(const std::string& directory, YAML::Node& output, std::string& error) {
    try {
        // [map-manifest-startup-gates]
        const YAML::Node root = ReadAndVerify(directory);
        if (root["map_format"].as<std::string>() != "tiled") throw std::runtime_error("this artifact is a global PCD preview, not a localization map package");
        const auto alignment = root["legacy_alignment"];
        const std::string status = alignment["status"].as<std::string>();
        if (status == "unknown") throw std::runtime_error("map alignment is unknown; align to the old map or explicitly declare not_required before localization");
        if (status != "aligned" && status != "not_required") throw std::runtime_error("invalid map alignment status");
        // [map-manifest-startup-gates]
        if (alignment["source_map_id"].as<std::string>() != root["map_id"].as<std::string>() ||
            alignment["target_map_id"].as<std::string>().empty() ||
            alignment["convention"].as<std::string>() != "target_from_localization" ||
            alignment["translation_unit"].as<std::string>() != "m") {
            throw std::runtime_error("invalid map alignment identity, direction or units");
        }
        // [map-manifest-alignment-once]
        const bool applied = alignment["applied_to_map_files"].as<bool>();
        const bool enabled = status == "aligned" && !applied;
        if (alignment["needs_transform_to_old"].as<bool>() != enabled) throw std::runtime_error("map transform requirement contradicts its application state");
        // [map-manifest-alignment-once]
        const SE3 transform = ReadPose(alignment["transform"]);
        if (status == "not_required" &&
            (applied || transform.translation().norm() > 1e-9 ||
             transform.unit_quaternion().angularDistance(Quatd::Identity()) > 1e-9 ||
             alignment["target_map_id"].as<std::string>() != root["map_id"].as<std::string>())) {
            throw std::runtime_error("not_required must identify the native map frame with an identity transform");
        }
        const std::string expected_frame = applied ? alignment["target_map_id"].as<std::string>() : root["map_id"].as<std::string>();
        if (root["coordinate_frame_id"].as<std::string>() != expected_frame) throw std::runtime_error("map coordinate frame contradicts its transform state");
        if (status == "aligned" && alignment["target_global_pcd_sha256"].as<std::string>().size() != 64) throw std::runtime_error("aligned map is missing its reference-map SHA-256");
        const YAML::Node declared = root["output"];
        const YAML::Node fixed = declared["fixed_map_transform"];
        if (declared["map_frame"].as<std::string>() != "map" || fixed["enabled"].as<bool>() != enabled ||
            fixed["convention"].as<std::string>() != "target_from_localization" ||
            fixed["source_frame"].as<std::string>() != "localization_map" || fixed["target_frame"].as<std::string>() != "map") {
            throw std::runtime_error("manifest output block contradicts its alignment");
        }
        YAML::Node declared_pose = YAML::Clone(fixed);
        declared_pose["matrix_4x4"] = PoseNode(enabled ? transform : SE3())["matrix_4x4"];
        const SE3 effective = ReadPose(declared_pose);
        const SE3 expected = enabled ? transform : SE3();
        if ((effective.translation() - expected.translation()).norm() > 1e-9 ||
            effective.unit_quaternion().angularDistance(expected.unit_quaternion()) > 1e-9) {
            throw std::runtime_error("manifest output pose contradicts its alignment");
        }
        output = OutputNode(transform, enabled);
        output["fixed_map_transform"]["provenance"]["map_manifest"] = (fs::path(directory) / kFilename).string();
        output["fixed_map_transform"]["provenance"]["source_map_id"] = root["map_id"].as<std::string>();
        output["fixed_map_transform"]["provenance"]["target_map_id"] = alignment["target_map_id"].as<std::string>();
        // [map-manifest-startup-result]
        return true;
    } catch (const std::exception& e) { error = e.what(); return false; }
    // [map-manifest-startup-result]
}

}  // namespace lightning::map_manifest
