#include <cmath>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>

#include "core/backend/btc_loop_detector.h"
#include "core/localization/btc_relocalizer.h"
#include "core/maps/map_frame.h"

namespace {

lightning::CloudPtr MakeStructuredCloud() {
    auto cloud = std::make_shared<lightning::PointCloudType>();
    auto add_point = [&](double x, double y, double z, float intensity) {
        lightning::PointType point;
        point.x = static_cast<float>(x);
        point.y = static_cast<float>(y);
        point.z = static_cast<float>(z);
        point.intensity = intensity;
        cloud->push_back(point);
    };

    for (int x = -30; x <= 30; ++x) {
        for (int y = -30; y <= 30; ++y) add_point(0.2 * x, 0.2 * y, 0.0, 1.0F);
    }
    for (int along = -30; along <= 30; ++along) {
        for (int height = 0; height <= 24; ++height) {
            const double a = 0.2 * along;
            const double z = 0.2 * height;
            add_point(-6.0, a, z, 2.0F);
            add_point(6.0, a, z, 3.0F);
            add_point(a, -6.0, z, 4.0F);
            if (along < 5 || along > 20) add_point(a, 6.0, z, 5.0F);
        }
    }
    for (int angle = 0; angle < 72; ++angle) {
        const double theta = angle * M_PI / 36.0;
        for (int height = 0; height <= 24; ++height) {
            add_point(2.0 + 0.45 * std::cos(theta), -1.5 + 0.45 * std::sin(theta),
                      0.2 * height, 8.0F);
        }
    }
    return cloud;
}

lightning::Keyframe::Ptr MakeKeyframe(unsigned long id, double timestamp, double x,
                                      const lightning::CloudPtr& cloud) {
    lightning::NavState state;
    state.timestamp_ = timestamp;
    state.SetPose(lightning::SE3(lightning::Quatd::Identity(), lightning::Vec3d(x, 0.0, 0.0)));
    return std::make_shared<lightning::Keyframe>(id, cloud, state);
}

}  // namespace

int main() {
    lightning::backend::BtcLoopDetectorOptions options;
    options.descriptor_submap_size = 2;
    options.max_points_per_submap = 30000;
    options.min_points_per_submap = 100;
    options.downsample_leaf_size = 0.15;
    options.min_loop_score = 0.10;
    options.max_drift_ratio = 0.5;
    options.refine_with_plane_icp = false;
    options.confirmation_count = 1;
    options.descriptor.skip_near_num_ = 1;
    options.descriptor.icp_threshold_ = 0.05;
    options.descriptor.similarity_threshold_ = 0.5;
    options.descriptor.summary_min_thre_ = 2.0;
    options.descriptor.non_max_suppression_radius_ = 0.5;
    options.descriptor.line_filter_enable_ = 0;
    options.descriptor.useful_corner_num_ = 200;

    lightning::backend::BtcLoopDetector detector(options);
    const auto cloud = MakeStructuredCloud();
    const double positions[] = {-0.1, 0.0, 1.9, 2.0, -0.1, 0.0};
    std::optional<lightning::backend::BtcLoopResult> result;
    for (unsigned long index = 0; index < 6; ++index) {
        const auto current = detector.AddKeyframe(
            MakeKeyframe(index, static_cast<double>(index), positions[index], cloud),
            lightning::SE3());
        if (index % 2 == 0 && current.has_value()) {
            std::cerr << "descriptor submap triggered before two frames" << std::endl;
            return 1;
        }
        if (current) {
            std::cout << "BTC submap " << current->current_descriptor_id
                      << ": descriptors=" << current->descriptor_count
                      << ", candidate=" << current->history_descriptor_id
                      << ", score=" << current->score
                      << ", reason=" << current->rejection_reason << std::endl;
            result = current;
        }
    }

    if (detector.Entries().size() != 3 || detector.PendingKeyframes() != 0) {
        std::cerr << "unexpected BTC database state: entries=" << detector.Entries().size()
                  << ", pending=" << detector.PendingKeyframes() << std::endl;
        return 2;
    }
    if (!result || !result->descriptor_generated || result->descriptor_count == 0) {
        std::cerr << "BTC failed to generate descriptors" << std::endl;
        return 3;
    }
    if (!result->candidate_found || result->history_descriptor_id != 0) {
        std::cerr << "BTC failed global revisit retrieval: candidate=" << result->history_descriptor_id
                  << ", reason=" << result->rejection_reason << std::endl;
        return 4;
    }
    if (!result->accepted) {
        std::cerr << "BTC revisit was rejected: score=" << result->score
                  << ", drift_ratio=" << result->drift_ratio
                  << ", reason=" << result->rejection_reason << std::endl;
        return 5;
    }

    lightning::backend::BtcLoopDetectorOptions fallback_options = options;
    fallback_options.min_loop_score = 2.0;
    fallback_options.enable_odom_revisit_fallback = true;
    fallback_options.odom_revisit_search_radius = 0.5;
    fallback_options.odom_revisit_min_journey = 1.0;
    lightning::backend::BtcLoopDetector fallback_detector(fallback_options);
    std::optional<lightning::backend::BtcLoopResult> fallback_result;
    for (unsigned long index = 0; index < 6; ++index) {
        const auto current = fallback_detector.AddKeyframe(
            MakeKeyframe(index, static_cast<double>(index), positions[index], cloud),
            lightning::SE3());
        if (current) fallback_result = current;
    }
    if (!fallback_result || !fallback_result->accepted ||
        fallback_result->candidate_source != "odom_revisit" ||
        fallback_result->history_descriptor_id != 0) {
        std::cerr << "odometry revisit fallback failed: source="
                  << (fallback_result ? fallback_result->candidate_source : "missing")
                  << ", candidate="
                  << (fallback_result ? fallback_result->history_descriptor_id : -1)
                  << ", reason="
                  << (fallback_result ? fallback_result->rejection_reason : "no_result")
                  << std::endl;
        return 10;
    }

    const auto unique = std::chrono::steady_clock::now().time_since_epoch().count();
    const std::filesystem::path temporary_root =
        std::filesystem::temp_directory_path() / ("lightning_btc_relocalizer_test_" + std::to_string(unique));
    const std::filesystem::path database_path = temporary_root / "btc_relocalization";
    if (!detector.SaveRelocalizationDatabase(database_path.string(), lightning::SE3())) {
        std::cerr << "failed to save BTC relocalization database" << std::endl;
        return 6;
    }
    const std::filesystem::path config_path = temporary_root / "config.yaml";
    {
        std::ofstream config(config_path);
        config << "relocalization:\n"
                  "  enabled: true\n"
                  "  query_submap_size: 2\n"
                  "  min_points_per_submap: 100\n"
                  "  min_btc_score: 0.10\n";
    }

    lightning::loc::BtcRelocalizer relocalizer;
    if (!relocalizer.Init(config_path.string(), temporary_root.string(), lightning::SE3()) ||
        !relocalizer.IsReady() || relocalizer.DatabaseSize() != detector.Entries().size()) {
        std::cerr << "failed to reload BTC relocalization database" << std::endl;
        return 7;
    }
    if (relocalizer.AddFrame(cloud, lightning::SE3(lightning::Quatd::Identity(),
                                                   lightning::Vec3d(-0.1, 0.0, 0.0)),
                             10.0)) {
        std::cerr << "BTC relocalization query triggered before two frames" << std::endl;
        return 8;
    }
    const auto relocalization = relocalizer.AddFrame(
        cloud, lightning::SE3(lightning::Quatd::Identity(), lightning::Vec3d(0.0, 0.0, 0.0)), 11.0);
    if (!relocalization || !relocalization->candidate_found || !relocalization->accepted ||
        relocalization->candidate_id != 0) {
        std::cerr << "BTC database reload query failed: candidate="
                  << (relocalization ? relocalization->candidate_id : -1)
                  << ", score=" << (relocalization ? relocalization->score : 0.0)
                  << ", reason=" << (relocalization ? relocalization->reason : "no_result") << std::endl;
        return 9;
    }

    lightning::map_frame::Metadata map_metadata;
    map_metadata.normalized = true;
    map_metadata.transform_id = "start_ground_z:1.250000000000";
    map_metadata.T_export_slam = lightning::SE3(
        lightning::Quatd::Identity(), lightning::Vec3d(0.0, 0.0, 1.25));
    map_metadata.ground_z_slam = -1.25;
    std::string map_frame_error;
    if (!lightning::map_frame::SaveMetadata(
            temporary_root.string(), map_metadata, map_frame_error) ||
        !detector.SaveRelocalizationDatabase(
            database_path.string(), lightning::SE3(), &map_metadata)) {
        std::cerr << "failed to save normalized map package: " << map_frame_error
                  << std::endl;
        return 11;
    }
    {
        std::ofstream config(config_path);
        config << "map_export:\n"
                  "  normalize_start_ground_z: true\n"
                  "relocalization:\n"
                  "  enabled: true\n"
                  "  query_submap_size: 2\n"
                  "  min_points_per_submap: 100\n"
                  "  min_btc_score: 0.10\n";
    }
    lightning::loc::BtcRelocalizer normalized_relocalizer;
    if (!normalized_relocalizer.Init(
            config_path.string(), temporary_root.string(), lightning::SE3())) {
        std::cerr << "failed to load a consistent normalized BTC database"
                  << std::endl;
        return 12;
    }

    lightning::map_frame::Metadata mismatched_metadata = map_metadata;
    mismatched_metadata.transform_id = "start_ground_z:2.000000000000";
    mismatched_metadata.T_export_slam = lightning::SE3(
        lightning::Quatd::Identity(), lightning::Vec3d(0.0, 0.0, 2.0));
    mismatched_metadata.ground_z_slam = -2.0;
    if (!lightning::map_frame::SaveMetadata(
            temporary_root.string(), mismatched_metadata, map_frame_error)) {
        std::cerr << "failed to write mismatch test metadata" << std::endl;
        return 13;
    }
    lightning::loc::BtcRelocalizer mismatched_relocalizer;
    if (mismatched_relocalizer.Init(
            config_path.string(), temporary_root.string(), lightning::SE3())) {
        std::cerr << "mismatched global-map and BTC transforms were accepted"
                  << std::endl;
        return 14;
    }
    std::error_code cleanup_error;
    std::filesystem::remove_all(temporary_root, cleanup_error);

    std::cout << "btc_loop_detector_test passed: descriptors=" << result->descriptor_count
              << ", score=" << result->score << ", drift_ratio=" << result->drift_ratio
              << ", reload_score=" << relocalization->score << std::endl;
    return 0;
}
