#include <cmath>
#include <iostream>
#include <stdexcept>

#include "core/localization/solid_descriptor.h"

namespace {

constexpr double kPi = 3.14159265358979323846;

lightning::PointCloudType MakeCloud(double yaw_deg, int copies = 1) {
    lightning::PointCloudType cloud;
    const double yaw = yaw_deg * kPi / 180.0;
    const double cosine = std::cos(yaw);
    const double sine = std::sin(yaw);
    for (int copy = 0; copy < copies; ++copy) {
        for (int index = 0; index < 180; ++index) {
            const double angle = (index * 1.7 + (index % 7) * 13.0) * kPi / 180.0;
            const double radius = 5.0 + (index % 17) * 0.8;
            const double x = radius * std::cos(angle);
            const double y = radius * std::sin(angle);
            lightning::PointType point;
            point.x = static_cast<float>(cosine * x - sine * y);
            point.y = static_cast<float>(sine * x + cosine * y);
            point.z = static_cast<float>(-2.0 + (index % 11) * 0.5);
            cloud.push_back(point);
        }
    }
    return cloud;
}

bool Near(double left, double right, double tolerance) {
    return std::abs(left - right) <= tolerance;
}

}  // namespace

int main() {
    lightning::loc::SolidDescriptorOptions options;
    options.angle_bins = 36;
    options.min_distance = 1.0;
    options.max_distance = 40.0;
    lightning::loc::SolidDescriptorEngine engine(options);

    const auto candidate = engine.Compute(MakeCloud(0.0));
    const auto query = engine.Compute(MakeCloud(30.0));
    const auto denser_candidate = engine.Compute(MakeCloud(0.0, 2));
    if (!candidate || !query || !denser_candidate) {
        std::cerr << "failed to build synthetic SOLiD descriptors\n";
        return 1;
    }
    if (!Near(lightning::loc::SolidDescriptorEngine::Similarity(
                  *candidate, *denser_candidate),
              1.0, 1e-12)) {
        std::cerr << "range descriptor is not density invariant\n";
        return 1;
    }
    const auto yaw = lightning::loc::SolidDescriptorEngine::EstimateCandidateFromQueryYaw(
        *query, *candidate);
    if (!yaw || !Near(*yaw, -30.0 * kPi / 180.0, 1e-9)) {
        std::cerr << "unexpected candidate-from-query yaw: "
                  << (yaw ? *yaw * 180.0 / kPi : 999.0) << " deg\n";
        return 1;
    }

    lightning::PointCloudType invalid;
    lightning::PointType below_minimum_range;
    below_minimum_range.x = 0.5F;
    below_minimum_range.y = 0.0F;
    below_minimum_range.z = 0.0F;
    invalid.push_back(below_minimum_range);
    if (engine.Compute(invalid)) {
        std::cerr << "out-of-FOV cloud unexpectedly produced a descriptor\n";
        return 1;
    }

    lightning::PointCloudType equal_elevation_counts;
    for (int bin = 0; bin < options.elevation_bins; ++bin) {
        const double elevation_deg =
            options.elevation_min_deg +
            (bin + 0.5) * (options.elevation_max_deg - options.elevation_min_deg) /
                options.elevation_bins;
        const double elevation = elevation_deg * kPi / 180.0;
        lightning::PointType point;
        point.x = 10.0F;
        point.y = 0.0F;
        point.z = static_cast<float>(10.0 * std::tan(elevation));
        equal_elevation_counts.push_back(point);
    }
    const auto equal_counts_descriptor = engine.Compute(equal_elevation_counts);
    if (!equal_counts_descriptor || !equal_counts_descriptor->range.allFinite() ||
        !equal_counts_descriptor->angle.allFinite()) {
        std::cerr << "equal elevation counts produced an invalid descriptor\n";
        return 1;
    }

    try {
        auto invalid_options = options;
        invalid_options.elevation_min_deg = 20.0;
        invalid_options.elevation_max_deg = -20.0;
        lightning::loc::SolidDescriptorEngine invalid_engine(invalid_options);
        (void)invalid_engine;
        std::cerr << "invalid options were accepted\n";
        return 1;
    } catch (const std::invalid_argument&) {
    }
    return 0;
}
