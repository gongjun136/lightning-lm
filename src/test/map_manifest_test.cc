#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "core/maps/map_manifest.h"

namespace {
void Require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
void Write(const std::filesystem::path& path, const std::string& value) {
    std::ofstream(path) << value;
}
}

int main() {
    namespace fs = std::filesystem;
    namespace manifest = lightning::map_manifest;
    const fs::path root = fs::temp_directory_path() /
        ("lightning_map_manifest_test_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        const fs::path old_map = root / "old", new_map = root / "new";
        fs::create_directories(old_map);
        fs::create_directories(new_map);
        for (const auto& directory : {old_map, new_map}) {
            Write(directory / "index.txt", "tile index\n");
            Write(directory / "0.pcd", "tile geometry\n");
            Write(directory / "global.pcd", directory.filename().string());
        }
        std::string error;
        YAML::Node output;
        Require(!manifest::LoadOutputConfig(new_map.string(), output, error) &&
                    error.find("missing") != std::string::npos, "missing manifests must stop localization");
        Require(manifest::SaveExportManifest(new_map.string(), {}, error), error);
        Require(!manifest::LoadOutputConfig(new_map.string(), output, error) &&
                    error.find("unknown") != std::string::npos, "unassessed maps must stop localization");
        const auto unknown = YAML::LoadFile((new_map / manifest::kFilename).string());
        Require(unknown["legacy_alignment"]["needs_transform_to_old"].IsNull(), "unknown must not imply no transform");

        manifest::ExportOptions native;
        native.alignment_status = "not_required";
        native.site_id = "test_site";
        Require(manifest::SaveExportManifest(old_map.string(), native, error), error);
        Require(manifest::LoadOutputConfig(old_map.string(), output, error), error);
        Require(!output["fixed_map_transform"]["enabled"].as<bool>(), "native reference must disable previous site transforms");
        const lightning::SE3 transform(
            lightning::Quatd(Eigen::AngleAxisd(0.4, lightning::Vec3d::UnitZ())),
            lightning::Vec3d(-2.4, 18.3, -0.13));
        Require(manifest::SaveAlignmentManifest(new_map.string(), old_map.string(), transform,
                    "accepted_report.yaml", 0.99, false, error), error);
        Require(manifest::LoadOutputConfig(new_map.string(), output, error), error);
        Require(output["fixed_map_transform"]["enabled"].as<bool>(), "aligned native geometry requires output conversion");
        auto aligned = YAML::LoadFile((new_map / manifest::kFilename).string());
        Require(aligned["legacy_alignment"]["target_map_id"].as<std::string>() == "old", "alignment must name its exact reference");
        Require(aligned["legacy_alignment"]["target_global_pcd_sha256"].as<std::string>() ==
                    manifest::FileSha256((old_map / "global.pcd").string()), "reference hash is required");
        const auto translation = output["fixed_map_transform"]["translation_xyz"].as<std::vector<double>>();
        Require(std::abs(translation[1] - 18.3) < 1e-12, "output translation lost precision/direction");
        const fs::path renamed = root / "renamed";
        fs::rename(new_map, renamed);
        Require(manifest::LoadOutputConfig(renamed.string(), output, error), "portable manifests must survive relocation: " + error);

        Write(renamed / "0.pcd", "changed tile geometry\n");
        Require(!manifest::LoadOutputConfig(renamed.string(), output, error) &&
                    error.find("checksum mismatch") != std::string::npos, "changed geometry must reject the old matrix");
        Write(renamed / "0.pcd", "tile geometry\n");
        Write(renamed / "1.pcd", "unrecorded tile\n");
        Require(!manifest::LoadOutputConfig(renamed.string(), output, error), "unrecorded tiles must reject stale manifests");
        fs::remove(renamed / "1.pcd");
        aligned["legacy_alignment"]["transform"]["matrix_4x4"][0][3] = 999.0;
        { std::ofstream stream(renamed / manifest::kFilename); stream << aligned; }
        Require(!manifest::LoadOutputConfig(renamed.string(), output, error) &&
                    error.find("disagree") != std::string::npos, "matrix/quaternion disagreement must fail");
        // Regenerate geometry identity, then declare an already converted artifact.
        Require(manifest::SaveExportManifest(renamed.string(), {}, error), error);
        Require(manifest::SaveAlignmentManifest(renamed.string(), old_map.string(), transform,
                    "accepted_report.yaml", 0.99, true, error), error);
        Require(manifest::LoadOutputConfig(renamed.string(), output, error), error);
        Require(!output["fixed_map_transform"]["enabled"].as<bool>(), "physically transformed geometry must not transform twice");
        const auto converted = YAML::LoadFile((renamed / manifest::kFilename).string());
        Require(converted["coordinate_frame_id"].as<std::string>() == "old" &&
                    converted["legacy_alignment"]["applied_to_map_files"].as<bool>(), "converted geometry must identify its actual frame");

        YAML::Node config = YAML::Load("map_export: {alignment_status: aligned}");
        manifest::ExportOptions options;
        Require(!manifest::ReadExportOptions(config, options, error), "new exports must not inherit an old alignment");
        const fs::path preview = root / "preview";
        fs::create_directories(preview);
        Write(preview / "global.pcd", "preview\n");
        Require(manifest::SaveExportManifest(preview.string(), native, error), error);
        Require(!manifest::LoadOutputConfig(preview.string(), output, error) &&
                    error.find("preview") != std::string::npos, "PCD-only exports cannot masquerade as localization packages");
        Write(root / "abc", "abc");
        Require(manifest::FileSha256((root / "abc").string()) ==
                    "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad", "SHA-256 known vector failed");
        fs::remove_all(root);
        std::cout << "map_manifest_test passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << "; retained test artifacts at " << root << '\n';
        return 1;
    }
}
