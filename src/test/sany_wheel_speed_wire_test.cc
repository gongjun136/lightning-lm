#include "common/sany_wheel_speed_wire.h"

#include <geosun_msgs/msg/spe_thr_can4.hpp>
#include <rclcpp/serialization.hpp>
#include <rclcpp/serialized_message.hpp>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>

void Require(bool condition) {
    if (!condition) throw std::runtime_error("wheel-speed wire regression failed");
}

int main() {
    // CDR1 fixture: Header(1791432657,101875609,"spe_thr_can4"), -1234.5 rpm, 42 Nm.
    const std::vector<uint8_t> legacy = {
        0,1,0,0, 0xd1,0x17,0xc7,0x6a, 0x99,0x7f,0x12,0x06,
        13,0,0,0, 's','p','e','_','t','h','r','_','c','a','n','4',0, 0,0,0,0,0,0,0,
        0,0,0,0,0,0x4a,0x93,0xc0, 0,0,0,0,0,0,0x45,0x40};
    lightning::SanyWheelSpeedWire sample;
    Require(lightning::DecodeSanyWheelSpeed(legacy.data(), legacy.size(), sample));
    Require(!sample.has_comm_header && sample.rpm == -1234.5 && sample.torque == 42);
    Require(std::abs(sample.stamp - (double(0x6ac717d1) + double(0x06127f99) * 1e-9)) < 1e-6);
    auto big_endian = legacy;
    big_endian[1] = 0;
    for (const size_t offset : {4u, 8u, 12u})
        std::reverse(big_endian.begin() + offset, big_endian.begin() + offset + 4);
    for (const size_t offset : {36u, 44u})
        std::reverse(big_endian.begin() + offset, big_endian.begin() + offset + 8);
    Require(lightning::DecodeSanyWheelSpeed(big_endian.data(), big_endian.size(), sample));
    Require(sample.rpm == -1234.5 && sample.torque == 42);
    // The installed ROS generator is an independent check of layout/alignment.
    geosun_msgs::msg::SpeThrCAN4 message;
    message.header.stamp.sec = 1791432913;
    message.header.stamp.nanosec = 101875609;
    message.header.frame_id = "spe_thr_can4";
    message.x = -1234.5;
    message.y = 42;
    rclcpp::Serialization<geosun_msgs::msg::SpeThrCAN4> serializer;
    for (const auto& frame : {std::string(""), std::string("spe_thr_can4"), std::string(31, 'a')}) {
        message.header.frame_id = frame;
        rclcpp::SerializedMessage serialized;
        serializer.serialize_message(&message, &serialized);
        const auto& raw = serialized.get_rcl_serialized_message();
        Require(lightning::DecodeSanyWheelSpeed(raw.buffer, raw.buffer_length, sample));
        Require(sample.rpm == message.x && sample.torque == message.y);
        Require(std::abs(sample.stamp - (1791432913.0 + 0.101875609)) < 1e-6);
        for (size_t size = 0; size < raw.buffer_length; ++size)
            Require(!lightning::DecodeSanyWheelSpeed(raw.buffer, size, sample));
    }
    auto corrupt = legacy;
    corrupt.push_back(0);  // No acceptance of partial deserialization/trailing payload.
    Require(!lightning::DecodeSanyWheelSpeed(corrupt.data(), corrupt.size(), sample));
    corrupt = legacy;
    corrupt[12] = 255;
    Require(!lightning::DecodeSanyWheelSpeed(corrupt.data(), corrupt.size(), sample));
    corrupt = legacy;
    corrupt[8] = corrupt[9] = corrupt[10] = corrupt[11] = 255;
    Require(!lightning::DecodeSanyWheelSpeed(corrupt.data(), corrupt.size(), sample));
    corrupt = legacy;
    std::fill(corrupt.begin() + 36, corrupt.begin() + 44, 0);
    corrupt[42] = 0xf8;
    corrupt[43] = 0x7f;  // NaN motor speed.
    Require(!lightning::DecodeSanyWheelSpeed(corrupt.data(), corrupt.size(), sample));
    std::cout << "CAN legacy/current layouts, frame alignment and malformed payload checks passed\n";
}
