#include "core/lio/laser_mapping.h"

#include <iostream>
#include <stdexcept>

namespace {

void Require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

int main() {
    lightning::LaserMapping disabled;
    const auto disabled_cloud = disabled.GetProjCloud();
    Require(disabled_cloud != nullptr, "disabled projection returns the current scan");
    Require(disabled_cloud == disabled.GetProjCloud(),
            "disabled projection does not allocate or replace the current scan");

    lightning::LaserMapping::Options options;
    options.proj_kfs_ = true;
    lightning::LaserMapping enabled(options);
    const auto first_projection = enabled.GetProjCloud();
    lightning::PointType point;
    point.x = 1.0F;
    first_projection->push_back(point);
    const auto second_projection = enabled.GetProjCloud();
    Require(first_projection != second_projection,
            "enabled projection returns an isolated point-cloud copy");
    Require(second_projection->empty(),
            "mutating a projected result cannot mutate the shared current scan");

    std::cout << "laser mapping projection tests passed\n";
    return 0;
}
