//
// Created by xiang on 23-2-24.
//

#pragma once

#include <algorithm>
#include <cmath>
#include <deque>
#include <sstream>

#include "common/debug_event.h"
#include "common/eigen_types.h"
#include "common/std_types.h"

namespace lightning::loc {

/**
 * 输出位姿的平滑器
 */
class PoseSmoother {
   public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**
     * 构造器
     * @param factor 越接近于0,则越平滑（但滞后），越接近于1,则平滑效果越弱
     */
    explicit PoseSmoother(double factor = 0.3) { smooth_factor_ = factor; }
    PoseSmoother(const PoseSmoother&) = default;

    bool PushDRPose(const SE3& dr_pose, double timestamp, double speed_mps) {
        UL lock(data_mutex_);

        if (!std::isfinite(timestamp) || timestamp <= 0.0 || !std::isfinite(speed_mps)) {
            LOG(WARNING) << "smoother rejects invalid DR metadata: timestamp="
                         << timestamp << ", speed_mps=" << speed_mps;
            debug_event::EmitThrottled(
                "smoother_invalid_dr_metadata",
                "Pose smoother rejected invalid DR metadata",
                std::chrono::seconds(1));
            motion_effective_ = false;
            return false;
        }

        /// 检查dr pose和末尾那个pose之间的距离
        if (!dr_queue_.empty()) {
            const double dt = timestamp - dr_queue_.back().timestamp;
            if (dt <= 0.0) {
                const double repeated_pose_distance =
                    (dr_pose.translation() - dr_queue_.back().pose.translation()).norm();
                if (std::abs(dt) <= 1e-9 && repeated_pose_distance <= 1e-9) {
                    // PubResult may run again for an LIO update without a new
                    // DR sample. Preserve the former identity-motion behavior.
                    motion_effective_ = true;
                    dr_queue_.push_back({dr_pose, timestamp, speed_mps});
                    while (dr_queue_.size() > max_size_) dr_queue_.pop_front();
                    return true;
                }
                LOG(WARNING) << "smoother rejects non-increasing DR timestamp: "
                             << timestamp << " <= " << dr_queue_.back().timestamp;
                debug_event::EmitThrottled(
                    "smoother_dr_timestamp_rollback",
                    "Pose smoother rejected a non-increasing DR timestamp",
                    std::chrono::seconds(1));
                motion_effective_ = false;
                return false;
            }
            const double distance =
                (dr_pose.translation() - dr_queue_.back().pose.translation()).norm();
            const double reference_speed =
                std::max(std::abs(speed_mps), std::abs(dr_queue_.back().speed_mps));
            const double allowed_distance =
                std::max(smoother_dr_limit_trans_,
                         smoother_dr_speed_scale_ * reference_speed * dt +
                             smoother_dr_distance_margin_);
            if (distance >= allowed_distance) {
                /// 无效
                LOG(WARNING) << "smoother motion is too large: distance=" << distance
                             << ", dt=" << dt << ", speed_mps=" << reference_speed
                             << ", allowed_distance=" << allowed_distance;
                std::ostringstream message;
                message << "Pose smoother rejected a motion jump: distance_m="
                        << distance << ", allowed_distance_m=" << allowed_distance;
                debug_event::EmitThrottled("smoother_motion_jump", message.str(),
                                           std::chrono::seconds(1));
                dr_queue_.clear();
                motion_effective_ = false;
                return false;
            }
        }

        motion_effective_ = true;
        dr_queue_.push_back({dr_pose, timestamp, speed_mps});
        while (dr_queue_.size() > max_size_) {
            dr_queue_.pop_front();
        }
        return true;
    }

    void PushPose(const SE3& pose) {
        UL lock(data_mutex_);
        if (pose_queue_.size() < 3) {
            pose_queue_.push_back(pose);
            output_pose_ = pose;
            return;
        }

        int n = dr_queue_.size();
        SE3 motion, pred;

        if (motion_effective_) {
            if (n >= 2) {
                motion = dr_queue_[n - 2].pose.inverse() * dr_queue_[n - 1].pose;
            } else {
                motion = SE3();
            }

            pred = pose_queue_.back() * motion;  // 预测pose
            output_pose_.translation() =
                pred.translation() * (1.0 - smooth_factor_) + pose.translation() * smooth_factor_;
            Vec3d rot_err = (pred.so3().inverse() * pose.so3()).log();
            output_pose_.so3() = pred.so3() * SO3::exp(rot_err * smooth_factor_);
        } else {
            pred = pose;
            output_pose_ = pose;
        }

        double dn = (output_pose_.translation() - pose.translation()).head<2>().norm();
        if (dn > smoother_trans_limit_) {
            // 平滑器差别太大，直接跳过去并重置
            LOG(WARNING) << "smoother is too large, jump to target: "
                         << (output_pose_.translation() - pose.translation()).transpose()
                         << ", given: " << pose.translation().transpose()
                         << ", output:" << output_pose_.translation().transpose();
            debug_event::EmitThrottled(
                "smoother_output_jump",
                "Pose smoother diverged from its target and was reset",
                std::chrono::seconds(1));
            output_pose_ = pose;
            pose_queue_.clear();
            return;
        } else if (dn > smoother_trans_limit2_) {
            // 快速收敛到正确位姿
            LOG(INFO) << "smoother factor set to 0.2";
            smooth_factor_ = 0.2;
        } else {
            // 回到默认值
            smooth_factor_ = 0.01;
        }

        pose_queue_.emplace_back(output_pose_);

        while (pose_queue_.size() > max_size_) {
            pose_queue_.pop_front();
        }
        motion_effective_ = false;
    }

    SE3 GetPose() const { return output_pose_; }

    void Reset() {
        UL lock(data_mutex_);
        pose_queue_.clear();
        dr_queue_.clear();
        output_pose_ = SE3();
        motion_effective_ = false;
    }

   private:
    std::mutex data_mutex_;

    bool motion_effective_ = true;  // DR的运动是否有效
    size_t max_size_ = 10;
    double smooth_factor_ = 0.1;            // 平滑因子
    double smoother_trans_limit_ = 5.0;     // 平滑器从输入到输出允许的最大平移量
    double smoother_trans_limit2_ = 2.0;    // 平滑器从输入到输出允许的最大平移量
    double smoother_dr_limit_trans_ = 0.3;  // 平滑器允许的DR跳变量
    double smoother_dr_speed_scale_ = 1.5;
    double smoother_dr_distance_margin_ = 0.1;
    SE3 output_pose_;

    struct TimedDrPose {
        SE3 pose;
        double timestamp = 0.0;
        double speed_mps = 0.0;
    };
    std::deque<SE3> pose_queue_;  // 平滑之后的
    std::deque<TimedDrPose> dr_queue_;  // 带时间和速度的DR观测
};

}  // namespace lightning::loc
