#include <cmath>
#include <iostream>
#include <vector>

#include "core/backend/voxel_bundle_adjustment.h"

namespace {

lightning::CloudPtr MakePlaneCloud(double local_x_offset) {
    auto cloud = std::make_shared<lightning::PointCloudType>();
    for (int x = -12; x <= 12; ++x) {
        for (int y = -12; y <= 12; ++y) {
            lightning::PointType point;
            point.x = static_cast<float>(0.15 * x + local_x_offset);
            point.y = static_cast<float>(0.15 * y);
            point.z = 0.0F;
            point.intensity = 1.0F;
            cloud->push_back(point);
        }
    }
    return cloud;
}

}  // namespace

int main() {
    using lightning::Quatd;
    using lightning::SE3;
    using lightning::Vec3d;
    using lightning::backend::BundleAdjustmentOptions;
    using lightning::backend::BundleFrame;
    using lightning::backend::VoxelBundleAdjuster;

    BundleAdjustmentOptions options;
    options.max_iterations = 8;
    options.voxel_size = 1.0;
    options.max_voxel_depth = 2;
    options.max_points_per_frame = 2000;
    options.min_points_per_voxel = 12;
    options.min_frames_per_voxel = 2;
    options.planarity_ratio = 0.2;
    options.correction_prior_rotation = 1e-5;
    options.correction_prior_translation = 1e-5;
    options.correction_smoothness_rotation = 1e-5;
    options.correction_smoothness_translation = 1e-5;
    options.max_translation_step = 0.25;
    options.max_total_translation = 1.0;

    std::vector<BundleFrame> frames(2);
    frames[0].id = 0;
    frames[0].cloud = MakePlaneCloud(0.0);
    frames[0].pose = SE3(Quatd::Identity(), Vec3d::Zero());
    frames[1].id = 1;
    frames[1].cloud = MakePlaneCloud(-0.6);
    frames[1].pose = SE3(Quatd::Identity(), Vec3d(0.6, 0.0, 0.20));

    VoxelBundleAdjuster optimizer(options);
    const auto summary = optimizer.Optimize(frames);
    if (!summary.accepted) {
        std::cerr << "expected accepted optimization, reason=" << summary.reason
                  << ", factors=" << summary.factor_count << std::endl;
        return 1;
    }
    if (!(summary.final_cost < summary.initial_cost)) {
        std::cerr << "cost did not decrease: " << summary.initial_cost << " -> " << summary.final_cost << std::endl;
        return 2;
    }
    if (!(std::abs(frames[1].pose.translation().z()) < 0.10)) {
        std::cerr << "plane offset was not corrected: z=" << frames[1].pose.translation().z() << std::endl;
        return 3;
    }
    if (frames[0].pose.translation().norm() > 1e-12) {
        std::cerr << "anchor pose moved" << std::endl;
        return 4;
    }

    std::cout << "voxel_bundle_adjustment_test passed: factors=" << summary.factor_count
              << ", cost=" << summary.initial_cost << " -> " << summary.final_cost << std::endl;
    return 0;
}
