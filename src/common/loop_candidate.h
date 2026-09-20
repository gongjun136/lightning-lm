//
// Created by xiang on 25-3-12.
//

#ifndef LIGHTNING_LOOP_CANDIDATE_H
#define LIGHTNING_LOOP_CANDIDATE_H

#include "common/eigen_types.h"

namespace lightning {

/**
 * @brief 回环检测候选约束。
 *
 * 记录一对可能形成回环的关键帧，以及点云配准后得到的IMU位姿相对约束和匹配分数。
 *
 * @note 关键帧中的点云保存在Lidar坐标系；回环NDT内部会先用Lidar-IMU外参处理点云位姿。
 *       存入Tij_前，配准结果应转换回位姿图顶点使用的IMU/body坐标系。
 */
struct LoopCandidate {
    /// 构造空候选，成员保持默认值。
    LoopCandidate() {}

    /// 使用历史关键帧ID和当前关键帧ID构造候选。
    LoopCandidate(uint64_t id1, uint64_t id2) : idx1_(id1), idx2_(id2) {}

    uint64_t idx1_ = 0;  ///< 历史关键帧ID，作为回环边起点；顶点位姿为T_w_i，表示idx1_帧IMU/body到世界系。
    uint64_t idx2_ = 0;  ///< 当前关键帧ID，作为回环边终点；顶点位姿为T_w_j，表示idx2_帧IMU/body到世界系。

    /// IMU/body坐标系下的相对位姿约束，Tij_ = T_i_j = T_w_i^{-1} * T_w_j。
    /// 从idx2_关键帧到idx1_关键帧的相对位姿约束。
    SE3 Tij_;

    double ndt_score_ = 0.0;  ///< NDT配准评分，用于筛选可信回环候选。
};

}  // namespace lightning

#endif  // LIGHTNING_LOOP_CANDIDATE_H
