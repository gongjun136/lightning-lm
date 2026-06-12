//
// Created by xiang on 25-3-12.
//

#ifndef LIGHTNING_KEYFRAME_H
#define LIGHTNING_KEYFRAME_H

#include "common/eigen_types.h"
#include "common/nav_state.h"
#include "common/point_def.h"
#include "common/std_types.h"

namespace lightning {

/// 关键帧描述
/// NOTE: 在添加后端后，需要加锁
///
/// Keyframe是前端、后端和显示模块之间传递的最小地图单元：
/// - cloud_保存该关键帧对应的去畸变点云。
/// - state_保存创建关键帧时的滤波器状态。
/// - pose_lio_记录前端里程计位姿，pose_opt_记录后端优化后的位姿。
class Keyframe {
   public:
    using Ptr = std::shared_ptr<Keyframe>;

    Keyframe() {}

    /// 用前端生成的点云和状态初始化关键帧。
    /// 初始优化位姿默认等于LIO位姿，后端优化完成后再通过SetOptPose更新。
    Keyframe(unsigned long id, CloudPtr cloud, NavState state)
        : id_(id), cloud_(cloud), state_(state), pose_lio_(state.GetPose()) {
        timestamp_ = state_.timestamp_;
        pose_opt_ = pose_lio_;
    }

    /// 返回关键帧ID；ID由前端创建关键帧时递增分配。
    unsigned long GetID() const { return id_; }

    /// 返回关键帧点云指针；点云本身不在这里做深拷贝。
    CloudPtr GetCloud() const { return cloud_; }

    /// 获取前端LIO位姿。
    SE3 GetLIOPose() {
        UL lock(data_mutex_);
        return pose_lio_;
    }

    /// 更新前端LIO位姿，并同步重置优化位姿。
    /// 适用于前端重新设定位姿的场景；若只想写入后端优化结果，应调用SetOptPose。
    void SetLIOPose(const SE3& pose) {
        UL lock(data_mutex_);
        pose_lio_ = pose;

        // also set opt
        pose_opt_ = pose_lio_;
    }

    /// 获取后端优化后的位姿；未优化前通常与LIO位姿一致。
    SE3 GetOptPose() {
        UL lock(data_mutex_);
        return pose_opt_;
    }

    /// 写入后端优化结果，不改变前端LIO位姿。
    void SetOptPose(const SE3& pose) {
        UL lock(data_mutex_);
        pose_opt_ = pose;
    }

    /// 更新关键帧携带的滤波器状态。
    void SetState(NavState s) {
        UL lock(data_mutex_);
        state_ = s;
    }

    /// 获取关键帧携带的滤波器状态。
    NavState GetState() {
        UL lock(data_mutex_);
        return state_;
    }

   protected:
    /// 关键帧唯一编号，通常按创建顺序递增。
    unsigned long id_ = 0;

    /// 关键帧时间戳，来自创建关键帧时的NavState。
    double timestamp_ = 0;

    /// 关键帧点云：降采样、去畸变后的lidar系点云。
    CloudPtr cloud_ = nullptr;  /// 降采样之后的点云，去畸变后，lidar系

    /// 保护位姿和状态的互斥锁，避免前后端或UI线程并发读写。
    std::mutex data_mutex_;

    /// 前端LIO估计的位姿，从IMU坐标系到世界坐标系。
    SE3 pose_lio_;  // 前端的pose，从IMU坐标系到世界坐标系

    /// 后端优化后的位姿，供回环、位姿图优化和地图拼接使用。
    SE3 pose_opt_;  // 后端优化后的pose

    /// 创建关键帧时的滤波器状态快照。
    NavState state_;  // 卡尔曼滤波器状态
};

}  // namespace lightning

#endif  // LIGHTNING_KEYFRAME_H
