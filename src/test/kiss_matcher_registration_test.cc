#include <cmath>
#include <iostream>
#include <limits>

#include "core/localization/kiss_matcher_registration.h"

namespace {

using lightning::AngAxisd;
using lightning::CloudPtr;
using lightning::PointCloudType;
using lightning::PointType;
using lightning::Quatd;
using lightning::SE3;
using lightning::Vec3d;
using lightning::loc::KissMatcherRegistration;

using Correspondence = KissMatcherRegistration::Correspondence;
using Correspondences =
    std::vector<Correspondence, Eigen::aligned_allocator<Correspondence>>;

KissMatcherRegistration::Options TestOptions() {
    KissMatcherRegistration::Options options;
    options.minimum_core_degree = 4;
    options.minimum_inliers = 15;
    options.compatibility_tolerance = 0.15;
    options.minimum_pair_distance = 0.25;
    options.inlier_threshold = 0.20;
    options.maximum_rmse = 0.10;
    options.spatial_cell_size = 2.0;
    options.minimum_spatial_cells = 4;
    return options;
}

SE3 TrueTransform() {
    return SE3(Quatd(AngAxisd(70.0 * M_PI / 180.0, Vec3d::UnitZ())),
               Vec3d(11.0, -7.0, 1.5));
}

void AddPoint(PointCloudType& cloud, const Vec3d& position) {
    PointType point;
    point.x = static_cast<float>(position.x());
    point.y = static_cast<float>(position.y());
    point.z = static_cast<float>(position.z());
    point.intensity = 0.0F;
    cloud.push_back(point);
}

CloudPtr MakeAsymmetricIndustrialScene() {
    CloudPtr cloud(new PointCloudType);
    for (double x = -6.0; x <= 6.0; x += 0.25) {
        for (double y = -5.0; y <= 5.0; y += 0.25) {
            const double z = 0.12 * std::sin(0.37 * x) +
                             0.08 * std::cos(0.51 * y) +
                             0.025 * std::sin(0.23 * x * y);
            AddPoint(*cloud, Vec3d(x, y, z));
        }
    }
    for (double y = -5.0; y <= 3.5; y += 0.25) {
        for (double z = 0.25; z <= 4.0; z += 0.25) {
            AddPoint(*cloud,
                     Vec3d(-5.25 + 0.03 * std::sin(0.7 * y + z), y, z));
        }
    }
    for (double x = -3.5; x <= 5.5; x += 0.25) {
        for (double z = 0.25; z <= 3.25; z += 0.25) {
            AddPoint(*cloud,
                     Vec3d(x, 4.25 + 0.02 * std::cos(x + 0.4 * z), z));
        }
    }
    const std::vector<Vec3d> columns = {
        Vec3d(-2.5, -1.5, 0.0), Vec3d(1.75, 2.0, 0.0),
        Vec3d(4.0, -3.0, 0.0)};
    for (std::size_t column = 0; column < columns.size(); ++column) {
        const double radius = 0.25 + 0.07 * static_cast<double>(column);
        for (double z = 0.25; z <= 3.5; z += 0.2) {
            for (int angle = 0; angle < 20; ++angle) {
                const double theta = 2.0 * M_PI * angle / 20.0;
                AddPoint(*cloud,
                         columns[column] +
                             Vec3d(radius * std::cos(theta),
                                   radius * std::sin(theta), z));
            }
        }
    }
    return cloud;
}

CloudPtr TransformToSource(const CloudPtr& target,
                           const SE3& target_from_source) {
    CloudPtr source(new PointCloudType);
    source->reserve(target->size());
    const SE3 source_from_target = target_from_source.inverse();
    for (const auto& point : target->points) {
        AddPoint(*source, source_from_target *
                              Vec3d(point.x, point.y, point.z));
    }
    return source;
}

Correspondences MakeCorrespondences(int inlier_count, int outlier_count,
                                    bool compact_inliers = false) {
    Correspondences correspondences;
    const SE3 transform = TrueTransform();
    for (int i = 0; i < inlier_count; ++i) {
        const double scale = compact_inliers ? 0.03 : 1.0;
        const Vec3d source(scale * static_cast<double>(i % 6),
                           scale * static_cast<double>(i / 6),
                           scale * static_cast<double>((i * 7) % 5) * 0.4);
        Vec3d target = transform * source;
        target += Vec3d(0.005 * static_cast<double>((i % 3) - 1),
                        0.004 * static_cast<double>(((i + 1) % 3) - 1),
                        0.003 * static_cast<double>(((i + 2) % 3) - 1));
        correspondences.push_back(Correspondence{source, target, 0.1});
    }
    for (int i = 0; i < outlier_count; ++i) {
        const Vec3d source(0.37 * i + 1.0, 0.19 * ((i * 11) % 17),
                           0.13 * ((i * 5) % 13));
        const Vec3d target(-0.41 * i + 25.0,
                           0.23 * ((i * 7) % 19) - 12.0,
                           0.17 * ((i * 3) % 11) + 5.0);
        correspondences.push_back(Correspondence{source, target, 0.2});
    }
    return correspondences;
}

bool TestLargeInitialDiscrepancyWithOutliers() {
    const auto options = TestOptions();
    const auto result = KissMatcherRegistration::SolveCorrespondences(
        MakeCorrespondences(30, 70), options);
    const SE3 error = TrueTransform().inverse() * result.T_target_source;
    if (!result.success || !result.converged || result.inliers < 25 ||
        result.core_correspondences < result.inliers ||
        result.spatial_cells <
            static_cast<std::size_t>(options.minimum_spatial_cells) ||
        !std::isfinite(result.rmse) || result.rmse > 0.02 ||
        error.translation().norm() > 0.02 ||
        error.so3().log().norm() > 0.2 * M_PI / 180.0) {
        std::cerr << "global robust solve failed: success=" << result.success
                  << " reason=" << result.reason
                  << " rough=" << result.rough_correspondences
                  << " core=" << result.core_correspondences
                  << " inliers=" << result.inliers
                  << " cells=" << result.spatial_cells
                  << " rmse=" << result.rmse
                  << " translation_error=" << error.translation().norm()
                  << " rotation_error_deg="
                  << error.so3().log().norm() * 180.0 / M_PI << '\n';
        return false;
    }
    return true;
}

bool TestSmallCoherentObjectDoesNotWin() {
    auto correspondences = MakeCorrespondences(28, 0);
    const SE3 wrong(
        Quatd(AngAxisd(-35.0 * M_PI / 180.0, Vec3d::UnitZ())),
        Vec3d(-16.0, 8.0, -0.5));
    for (int i = 0; i < 10; ++i) {
        const Vec3d source(0.3 * (i % 5), 0.3 * (i / 5), 0.1 * (i % 3));
        correspondences.push_back(
            Correspondence{source, wrong * source, 0.05});
    }
    const auto result = KissMatcherRegistration::SolveCorrespondences(
        correspondences, TestOptions());
    const SE3 error = TrueTransform().inverse() * result.T_target_source;
    if (!result.success || error.translation().norm() > 0.03 ||
        error.so3().log().norm() > 0.3 * M_PI / 180.0) {
        std::cerr << "small coherent outlier cluster won: success="
                  << result.success << " reason=" << result.reason
                  << " translation_error=" << error.translation().norm()
                  << " rotation_error_deg="
                  << error.so3().log().norm() * 180.0 / M_PI << '\n';
        return false;
    }
    return true;
}

bool TestSpatialCoverageGate() {
    auto options = TestOptions();
    options.minimum_pair_distance = 0.01;
    const auto result = KissMatcherRegistration::SolveCorrespondences(
        MakeCorrespondences(30, 0, true), options);
    if (result.success || result.reason != "insufficient_spatial_coverage" ||
        result.spatial_cells >=
            static_cast<std::size_t>(options.minimum_spatial_cells)) {
        std::cerr << "compact object passed coverage gate: success="
                  << result.success << " reason=" << result.reason
                  << " cells=" << result.spatial_cells << '\n';
        return false;
    }
    return true;
}

bool TestOptionValidation() {
    auto options = TestOptions();
    std::string error;
    if (!KissMatcherRegistration::ValidateOptions(options, &error)) return false;
    options.gnc_factor = 1.0;
    if (KissMatcherRegistration::ValidateOptions(options, &error)) {
        std::cerr << "invalid GNC factor accepted\n";
        return false;
    }
    options = TestOptions();
    options.maximum_rmse = std::numeric_limits<double>::quiet_NaN();
    if (KissMatcherRegistration::ValidateOptions(options, &error)) {
        std::cerr << "non-finite RMSE accepted\n";
        return false;
    }
    return true;
}

bool TestFeaturePipelineWithoutInitialGuess() {
    auto options = TestOptions();
    options.voxel_size = 0.25;
    options.normal_radius = 0.75;
    options.feature_radius = 1.25;
    options.descriptor_ratio = 0.98;
    options.max_feature_points = 20000;
    options.max_correspondences = 1000;
    options.compatibility_tolerance = 0.30;
    options.minimum_pair_distance = 0.50;
    options.minimum_core_degree = 3;
    options.minimum_inliers = 12;
    options.inlier_threshold = 0.35;
    options.maximum_rmse = 0.15;
    options.spatial_cell_size = 2.0;
    options.minimum_spatial_cells = 3;
    options.feature_threads = 2;

    const CloudPtr target = MakeAsymmetricIndustrialScene();
    const SE3 truth(
        Quatd(AngAxisd(90.0 * M_PI / 180.0, Vec3d::UnitZ())),
        Vec3d(8.0, -5.0, 1.5));
    const CloudPtr source = TransformToSource(target, truth);
    KissMatcherRegistration registration(options);
    std::string error;
    if (!registration.SetTarget(target, &error)) {
        std::cerr << "feature target construction failed: " << error << '\n';
        return false;
    }
    KissMatcherRegistration::PreparedCloud prepared;
    if (!registration.PrepareSource(source, prepared, &error)) {
        std::cerr << "feature source preparation failed: " << error << '\n';
        return false;
    }
    const auto result = registration.Align(prepared);
    const auto direct_result = registration.Align(source);
    const SE3 pose_error = truth.inverse() * result.T_target_source;
    const SE3 path_delta = result.T_target_source.inverse() *
                           direct_result.T_target_source;
    if (!result.success || !direct_result.success ||
        result.rough_correspondences < 20 ||
        result.inliers < 12 || pose_error.translation().norm() > 0.20 ||
        pose_error.so3().log().norm() > 1.0 * M_PI / 180.0 ||
        path_delta.translation().norm() > 0.02 ||
        path_delta.so3().log().norm() > 0.1 * M_PI / 180.0) {
        std::cerr << "initial-guess-free feature alignment failed: success="
                  << result.success << " reason=" << result.reason
                  << " direct_success=" << direct_result.success
                  << " descriptors=" << result.source_descriptors
                  << " rough=" << result.rough_correspondences
                  << " core=" << result.core_correspondences
                  << " inliers=" << result.inliers
                  << " rmse=" << result.rmse
                  << " translation_error=" << pose_error.translation().norm()
                  << " rotation_error_deg="
                  << pose_error.so3().log().norm() * 180.0 / M_PI
                  << " path_translation_delta="
                  << path_delta.translation().norm()
                  << " path_rotation_delta_deg="
                  << path_delta.so3().log().norm() * 180.0 / M_PI << '\n';
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!TestLargeInitialDiscrepancyWithOutliers() ||
        !TestSmallCoherentObjectDoesNotWin() || !TestSpatialCoverageGate() ||
        !TestOptionValidation() || !TestFeaturePipelineWithoutInitialGuess()) {
        return 1;
    }
    std::cout << "KISS-Matcher registration tests passed\n";
    return 0;
}
