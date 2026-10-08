#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

namespace lightning {

struct SanyWheelSpeedWire {
    double stamp = 0.0;
    double rpm = 0.0;
    double torque = 0.0;
    bool has_comm_header = false;
};

// Two released geosun_msgs/SpeThrCAN4 definitions share a DDS type name:
// Header,x,y and TopicCommHeader(uint16,uint32,uint64),Header,x,y.
// Deserialize the wire layout explicitly: Humble does not reject this mismatch
// and a typed subscription can silently turn frame_id bytes into a timestamp.
inline bool DecodeSanyWheelSpeed(const uint8_t* data, size_t size, SanyWheelSpeedWire& result) {
    if (!data || size < 4 || data[0] != 0 || data[1] > 1 || data[2] != 0 || data[3] != 0) return false;
    const bool little = data[1] == 1;
    const auto parse = [&](bool extended, SanyWheelSpeedWire& value) {
        size_t offset = 4;
        const auto integer = [&](size_t width, uint64_t& out) {
            offset += (width - (offset - 4) % width) % width;
            if (offset > size || width > size - offset) return false;
            out = 0;
            for (size_t i = 0; i < width; ++i) {
                const size_t shift = little ? i : width - i - 1;
                out |= uint64_t(data[offset + i]) << (8 * shift);
            }
            offset += width;
            return true;
        };
        uint64_t ignored, sec, nsec, length, bits;
        if (extended && (!integer(2, ignored) || !integer(4, ignored) || !integer(8, ignored))) return false;
        if (!integer(4, sec) || !integer(4, nsec) || !integer(4, length)) return false;
        if (sec == 0 || sec > INT32_MAX || nsec >= 1000000000 || length == 0 || length > 1024 ||
            length > size - offset || data[offset + length - 1] != 0) return false;
        for (size_t i = 0; i + 1 < length; ++i) {
            if (data[offset + i] == 0) return false;
        }
        offset += length;
        if (!integer(8, bits)) return false;
        std::memcpy(&value.rpm, &bits, sizeof(bits));
        if (!integer(8, bits)) return false;
        std::memcpy(&value.torque, &bits, sizeof(bits));
        if (offset != size || !std::isfinite(value.rpm) || !std::isfinite(value.torque) ||
            std::abs(value.rpm) > 15000.0 || std::abs(value.torque) > 5000.0) return false;
        value.stamp = double(sec) + double(nsec) * 1e-9;
        value.has_comm_header = extended;
        return true;
    };
    SanyWheelSpeedWire legacy, extended;
    const bool legacy_ok = parse(false, legacy);
    const bool extended_ok = parse(true, extended);
    if (legacy_ok == extended_ok) return false;  // Reject malformed or ambiguous payloads.
    result = legacy_ok ? legacy : extended;
    return true;
}

}  // namespace lightning
