#include "core/lio/rear_axle_pose.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

void Require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << std::endl;
        std::exit(1);
    }
}

}  // namespace

int main() {
    using namespace lightning;
    RearAxlePoseTransformer transformer(Mat3d::Identity(), Vec3d::Zero(), Vec3d(2.0, 0.0, 0.0));

    NavState initial;
    initial.pose_is_ok_ = true;
    initial.rot_ = SO3();
    initial.pos_.setZero();
    const SE3 first = transformer.Transform(initial);
    Require(first.translation().norm() < 1e-12, "first rear axle position is zero");
    Require(first.so3().log().norm() < 1e-12, "first rear axle orientation is identity");

    NavState rotated = initial;
    rotated.rot_ = SO3::exp(Vec3d(0.0, 0.0, M_PI_2));
    const SE3 second = transformer.Transform(rotated);
    Require((second.translation() - Vec3d(2.0, -2.0, 0.0)).norm() < 1e-9,
            "lever arm changes rear axle position during rotation");
    Require((second.so3().log() - Vec3d(0.0, 0.0, M_PI_2)).norm() < 1e-9,
            "rear axle orientation follows rigid body");

    std::cout << "rear_axle_pose_test passed" << std::endl;
    return 0;
}
