#include <cstdint>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "core/lio/pointcloud_preprocess.h"

namespace {

using livox_ros_driver2::msg::CustomMsg;
using livox_ros_driver2::msg::CustomPoint;
using lightning::PointCloudPreprocess;
using lightning::PointCloudType;

CustomMsg::SharedPtr MakeTwoPointMessage(const CustomPoint &point) {
    auto message = std::make_shared<CustomMsg>();
    message->points.resize(2);
    message->points[0].x = 0.0F;
    message->points[0].y = 0.0F;
    message->points[0].z = 0.0F;
    message->points[0].tag = 0x10;
    message->points[0].line = 0;
    message->points[1] = point;
    message->point_num = static_cast<std::uint32_t>(message->points.size());
    return message;
}

PointCloudType::Ptr Process(PointCloudPreprocess &preprocess, const CustomMsg::SharedPtr &message) {
    PointCloudType::Ptr output;
    preprocess.Process(message, output);
    if (!output) throw std::runtime_error("preprocessor returned a null cloud");
    return output;
}

}  // namespace

int main() {
#define CHECK(expression)                                                                        \
    do {                                                                                         \
        if (!(expression)) {                                                                     \
            std::cerr << "CHECK failed at line " << __LINE__ << ": " #expression << '\n';     \
            return 1;                                                                            \
        }                                                                                        \
    } while (false)

    PointCloudPreprocess preprocess;
    preprocess.Set(lightning::LidarType::AVIA, 0.5, 1);
    preprocess.NumScans() = 4;
    preprocess.SetHeightROI(100.0F, -100.0F);

    auto empty = std::make_shared<CustomMsg>();
    empty->point_num = 0;
    CHECK(Process(preprocess, empty)->empty());

    auto mismatched = std::make_shared<CustomMsg>();
    mismatched->point_num = 100;
    CHECK(Process(preprocess, mismatched)->empty());

    CustomPoint near_x;
    near_x.x = 0.1F;
    near_x.tag = 0x10;
    near_x.line = 0;
    CHECK(Process(preprocess, MakeTwoPointMessage(near_x))->empty());

    CustomPoint invalid_tag;
    invalid_tag.x = 1.0F;
    invalid_tag.tag = 0x30;
    invalid_tag.line = 0;
    CHECK(Process(preprocess, MakeTwoPointMessage(invalid_tag))->empty());

    CustomPoint invalid_line;
    invalid_line.x = 1.0F;
    invalid_line.tag = 0x10;
    invalid_line.line = 4;
    CHECK(Process(preprocess, MakeTwoPointMessage(invalid_line))->empty());

    CustomPoint valid;
    valid.x = 1.0F;
    valid.y = 0.2F;
    valid.tag = 0x10;
    valid.line = 3;
    valid.reflectivity = 42;
    valid.offset_time = 1000000;
    const auto output = Process(preprocess, MakeTwoPointMessage(valid));
    CHECK(output->size() == 1);
    CHECK(output->front().x == valid.x);
    CHECK(output->front().lidar_id == 0);
    CHECK(output->front().time == 1.0);

    std::cout << "pointcloud_preprocess_test passed\n";
#undef CHECK
    return 0;
}
