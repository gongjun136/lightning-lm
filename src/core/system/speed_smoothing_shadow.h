#pragma once

#include <cmath>
#include <string>

namespace lightning::sany_output {

// Diagnostic-only candidate. Never use this state as an estimator observation.
class SpeedSmoothingShadow {
 public:
    static bool Enabled(const char* mode, const char* requested, bool reduced_io) {
        return mode && requested && std::string(mode) == "diagnostic" &&
               std::string(requested) == "1" && !reduced_io;
    }

    void Reset() { initialized_ = false; }

    double Observe(double stamp, double raw, bool parking) {
        if (!std::isfinite(stamp) || !std::isfinite(raw)) {
            Reset();
            return raw;
        }
        const double dt = stamp - stamp_;
        reset_last_ = !initialized_ || dt <= 0.0 || dt > 0.1;
        if (parking) {
            value_ = 0.0;
        } else if (reset_last_) {
            value_ = raw;
        } else {
            value_ += -std::expm1(-dt / 0.03) * (raw - value_);
        }
        stamp_ = stamp;
        initialized_ = true;
        return value_;
    }

    bool ResetLast() const { return reset_last_; }

 private:
    bool initialized_ = false;
    bool reset_last_ = true;
    double stamp_ = 0.0;
    double value_ = 0.0;
};

}  // namespace lightning::sany_output
