#include <gflags/gflags.h>
#include <glog/logging.h>
#include <filesystem>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "core/maps/map_manifest.h"

DEFINE_string(map, "", "existing map package directory");
DEFINE_string(site_id, "", "physical site identity");
DEFINE_string(map_id, "", "map version; default derives from the directory");
DEFINE_string(status, "unknown", "new-map status: unknown or explicitly not_required");
DEFINE_string(old_map, "", "old/reference map required when importing an accepted alignment report");
DEFINE_string(alignment_report, "", "accepted align_maps_offline report to import");
DEFINE_bool(files_already_transformed, false, "map geometry already contains the report's transform");
DEFINE_bool(check_only, false, "validate localization readiness without modifying the map");
DEFINE_string(apply_to_config, "", "validate and replace only the map-output block in this run snapshot");

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    try {
        namespace fs = std::filesystem;
        namespace manifest = lightning::map_manifest;
        if (FLAGS_map.empty()) throw std::invalid_argument("--map is required");
        std::string error;
        if (FLAGS_check_only || !FLAGS_apply_to_config.empty()) {
            YAML::Node output;
            if (!manifest::LoadOutputConfig(FLAGS_map, output, error)) throw std::runtime_error(error);
            if (!FLAGS_apply_to_config.empty()) {
                YAML::Node config = YAML::LoadFile(FLAGS_apply_to_config);
                config["output"]["map_frame"] = output["map_frame"];
                config["output"]["fixed_map_transform"] = output["fixed_map_transform"];
                config["map_package"]["manifest"] = (fs::path(FLAGS_map) / manifest::kFilename).string();
                config["map_package"]["manifest_sha256"] = manifest::FileSha256(
                    (fs::path(FLAGS_map) / manifest::kFilename).string());
                YAML::Emitter emitter;
                emitter.SetDoublePrecision(17);
                emitter << config;
                if (!emitter.good()) throw std::runtime_error("cannot serialize run configuration");
                const fs::path temporary = FLAGS_apply_to_config + ".map_manifest.tmp";
                std::ofstream stream(temporary);
                stream << emitter.c_str() << '\n';
                stream.close();
                if (!stream) throw std::runtime_error("cannot write run configuration");
                fs::rename(temporary, FLAGS_apply_to_config);
            }
            std::cout << "MAP_MANIFEST_OK map=" << FLAGS_map
                      << " transform_enabled=" << output["fixed_map_transform"]["enabled"].as<bool>() << '\n';
            return 0;
        }
        manifest::ExportOptions options;
        options.site_id = FLAGS_site_id;
        options.map_id = FLAGS_map_id;
        options.alignment_status = FLAGS_status;
        YAML::Node report;
        lightning::SE3 transform;
        double quality = 0.0;
        if (!FLAGS_alignment_report.empty()) {
            if (FLAGS_old_map.empty()) throw std::invalid_argument("--alignment_report requires --old_map");
            report = YAML::LoadFile(FLAGS_alignment_report);
            if (!report["success"].as<bool>()) throw std::invalid_argument("alignment report did not pass");
            const auto t = report["transform"]["translation_xyz"].as<std::vector<double>>();
            const auto q = report["transform"]["quaternion_xyzw"].as<std::vector<double>>();
            if (t.size() != 3 || q.size() != 4) throw std::invalid_argument("invalid report transform dimensions");
            lightning::Quatd rotation(q[3], q[0], q[1], q[2]);
            const lightning::Vec3d translation(t[0], t[1], t[2]);
            if (!rotation.coeffs().allFinite() || std::abs(rotation.norm() - 1.0) > 1e-6 || !translation.allFinite()) {
                throw std::invalid_argument("invalid report transform");
            }
            rotation.normalize();
            transform = lightning::SE3(rotation, translation);
            quality = report["quality"].as<double>();
            // Importing an old report must not bind its matrix to different geometry.
            const fs::path reported_old = report["old_map"].as<std::string>();
            if (manifest::FileSha256((reported_old / "global.pcd").string()) !=
                manifest::FileSha256((fs::path(FLAGS_old_map) / "global.pcd").string())) {
                throw std::invalid_argument("reference map differs from the accepted report");
            }
            if (!FLAGS_files_already_transformed) {
                const fs::path reported_new = report["new_map"].as<std::string>();
                if (manifest::FileSha256((reported_new / "global.pcd").string()) !=
                    manifest::FileSha256((fs::path(FLAGS_map) / "global.pcd").string())) {
                    throw std::invalid_argument("source map differs from the accepted report");
                }
            }
        } else if (FLAGS_files_already_transformed || !FLAGS_old_map.empty()) {
            throw std::invalid_argument("physical conversion/reference options require an alignment report");
        }
        if (!manifest::SaveExportManifest(FLAGS_map, options, error)) throw std::runtime_error(error);
        if (!FLAGS_alignment_report.empty() && !manifest::SaveAlignmentManifest(
                FLAGS_map, FLAGS_old_map, transform, FLAGS_alignment_report, quality,
                FLAGS_files_already_transformed, error)) throw std::runtime_error(error);
        std::cout << (fs::path(FLAGS_map) / manifest::kFilename) << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Map manifest rejected: " << e.what() << '\n';
        return 2;
    }
}
