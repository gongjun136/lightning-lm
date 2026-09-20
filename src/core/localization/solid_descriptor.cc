#include "core/localization/solid_descriptor.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace lightning::loc {
namespace {

constexpr double kPi = 3.14159265358979323846;

bool FinitePositive(double value) {
    return std::isfinite(value) && value > 0.0;
}

}  // namespace

SolidDescriptorEngine::SolidDescriptorEngine(SolidDescriptorOptions options)
    : options_(std::move(options)) {
    if (options_.range_bins <= 0 || options_.angle_bins <= 0 ||
        options_.elevation_bins <= 0 || options_.min_distance < 0.0 ||
        !FinitePositive(options_.max_distance) ||
        options_.min_distance >= options_.max_distance ||
        !std::isfinite(options_.elevation_min_deg) ||
        !std::isfinite(options_.elevation_max_deg) ||
        options_.elevation_min_deg < -90.0 ||
        options_.elevation_max_deg > 90.0 ||
        options_.elevation_min_deg >= options_.elevation_max_deg) {
        throw std::invalid_argument("invalid SOLiD descriptor options");
    }
}

std::optional<SolidDescriptor> SolidDescriptorEngine::Compute(
    const PointCloudType& cloud) const {
    MatXd range_height = MatXd::Zero(options_.range_bins, options_.elevation_bins);
    MatXd angle_height = MatXd::Zero(options_.angle_bins, options_.elevation_bins);
    const double range_width = options_.max_distance / options_.range_bins;
    const double angle_width = 2.0 * kPi / options_.angle_bins;
    const double elevation_width =
        (options_.elevation_max_deg - options_.elevation_min_deg) /
        options_.elevation_bins;
    std::size_t accepted_points = 0;

    for (const auto& point : cloud.points) {
        if (!std::isfinite(point.x) || !std::isfinite(point.y) ||
            !std::isfinite(point.z)) {
            continue;
        }
        const double x = point.x;
        const double y = point.y;
        const double radial_distance = std::hypot(x, y);
        if (radial_distance < options_.min_distance ||
            radial_distance >= options_.max_distance) {
            continue;
        }
        double angle = std::atan2(y, x);
        if (angle < 0.0) angle += 2.0 * kPi;
        const double elevation_deg =
            std::atan2(static_cast<double>(point.z), radial_distance) * 180.0 / kPi;
        if (elevation_deg < options_.elevation_min_deg ||
            elevation_deg > options_.elevation_max_deg) {
            continue;
        }

        const int range_index = std::clamp(
            static_cast<int>(radial_distance / range_width), 0,
            options_.range_bins - 1);
        const int angle_index = std::clamp(
            static_cast<int>(angle / angle_width), 0, options_.angle_bins - 1);
        const int elevation_index = std::clamp(
            static_cast<int>((elevation_deg - options_.elevation_min_deg) /
                             elevation_width),
            0, options_.elevation_bins - 1);
        range_height(range_index, elevation_index) += 1.0;
        angle_height(angle_index, elevation_index) += 1.0;
        ++accepted_points;
    }

    if (accepted_points == 0) return std::nullopt;
    VecXd elevation_counts = range_height.colwise().sum().transpose();
    const double minimum = elevation_counts.minCoeff();
    const double maximum = elevation_counts.maxCoeff();
    VecXd elevation_weights;
    if (maximum - minimum > std::numeric_limits<double>::epsilon()) {
        elevation_weights = (elevation_counts.array() - minimum) /
                            (maximum - minimum);
    } else {
        // The public implementation divides by zero here.  Equal non-zero
        // counts contain useful geometry, so give each occupied bin unit weight.
        elevation_weights = (elevation_counts.array() > 0.0).cast<double>();
    }

    SolidDescriptor descriptor;
    descriptor.range = range_height * elevation_weights;
    descriptor.angle = angle_height * elevation_weights;
    descriptor.accepted_points = accepted_points;
    if (!descriptor.range.allFinite() || !descriptor.angle.allFinite() ||
        descriptor.range.norm() <= std::numeric_limits<double>::epsilon() ||
        descriptor.angle.norm() <= std::numeric_limits<double>::epsilon()) {
        return std::nullopt;
    }
    return descriptor;
}

double SolidDescriptorEngine::Similarity(const SolidDescriptor& query,
                                         const SolidDescriptor& candidate) {
    if (query.range.size() == 0 || query.range.size() != candidate.range.size() ||
        !query.range.allFinite() || !candidate.range.allFinite()) {
        return -1.0;
    }
    const double denominator = query.range.norm() * candidate.range.norm();
    if (denominator <= std::numeric_limits<double>::epsilon()) return -1.0;
    return std::clamp(query.range.dot(candidate.range) / denominator, -1.0, 1.0);
}

std::optional<double> SolidDescriptorEngine::EstimateCandidateFromQueryYaw(
    const SolidDescriptor& query, const SolidDescriptor& candidate) {
    if (query.angle.size() == 0 || query.angle.size() != candidate.angle.size() ||
        !query.angle.allFinite() || !candidate.angle.allFinite()) {
        return std::nullopt;
    }
    const int count = static_cast<int>(query.angle.size());
    double best_distance = std::numeric_limits<double>::infinity();
    int best_shift = 0;
    for (int shift = 0; shift < count; ++shift) {
        double distance = 0.0;
        for (int index = 0; index < count; ++index) {
            const int shifted_index = (index + shift) % count;
            distance += std::abs(candidate.angle[shifted_index] - query.angle[index]);
        }
        if (distance < best_distance) {
            best_distance = distance;
            best_shift = shift;
        }
    }
    if (best_shift * 2 >= count) best_shift -= count;
    return best_shift * (2.0 * kPi / count);
}

}  // namespace lightning::loc
