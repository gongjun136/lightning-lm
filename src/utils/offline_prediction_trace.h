#pragma once

// Opt-in, synchronous forensic logging for isolated offline replay only.
// Never enabled by diagnostic/production mode; not suitable for timing benchmarks.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include "common/nav_state.h"

namespace lightning::profiling {
class OfflinePredictionTrace {
 public:
    static OfflinePredictionTrace& Instance() { static OfflinePredictionTrace sink; return sink; }
    bool Enabled() const { return file_ != nullptr; }
    void Record(const char* source, const NavState& before, const NavState& predicted,
                const NavState& final, const Vec3d& raw_acc, const Vec3d& acc,
                const Vec3d& gyro, bool propagate, bool hold) {
        if (!file_) return;
        std::lock_guard<std::mutex> lock(mutex_);
        const double wall = std::chrono::duration<double>(
            std::chrono::system_clock::now().time_since_epoch()).count();
        int rc = std::fprintf(file_, "%s,%.17g,%.17g,%.17g,%d,%d", source, wall,
                              before.timestamp_, predicted.timestamp_, propagate, hold);
        auto number = [&](double value) { if (std::fprintf(file_, ",%.17g", value) < 0) rc = -1; };
        auto vector = [&](const Vec3d& value) { for (int i=0;i<3;++i) number(value[i]); };
        vector(raw_acc); vector(acc); vector(gyro); vector(before.ba_); vector(before.bg_);
        vector(before.grav_); vector(before.vel_); vector(predicted.vel_);
        for (const NavState* state : {&before, &predicted}) {
            const auto q = state->rot_.unit_quaternion();
            number(q.x()); number(q.y()); number(q.z()); number(q.w());
        }
        number(before.rot_.matrix().col(0).dot(before.vel_));
        number(predicted.rot_.matrix().col(0).dot(predicted.vel_));
        number(final.rot_.matrix().col(0).dot(final.vel_));
        if (std::fprintf(file_, "\n") < 0) rc = -1;
        ++rows_; if (rc < 0) ++errors_;
    }
 private:
    OfflinePredictionTrace() {
        const char* path = std::getenv("LIGHTNING_LM_OFFLINE_PREDICTION_TRACE_PATH");
        if (!path || !*path) return;
        file_ = std::fopen(path, "wx");
        if (!file_) { std::perror("offline prediction trace"); return; }
        std::fprintf(file_, "source,wall,prior,stamp,propagate,hold,raw_ax,raw_ay,raw_az,ax,ay,az,gx,gy,gz,bax,bay,baz,bgx,bgy,bgz,grav_x,grav_y,grav_z,v0x,v0y,v0z,v1x,v1y,v1z,q0x,q0y,q0z,q0w,q1x,q1y,q1z,q1w,before_v,predicted_v,final_v\n");
    }
    ~OfflinePredictionTrace() {
        if (!file_) return;
        std::fprintf(file_, "# rows=%llu errors=%llu\n", rows_, errors_);
        if (std::fclose(file_) != 0) std::perror("offline prediction trace close");
    }
    FILE* file_ = nullptr;
    std::mutex mutex_;
    unsigned long long rows_ = 0, errors_ = 0;
};
}  // namespace lightning::profiling
