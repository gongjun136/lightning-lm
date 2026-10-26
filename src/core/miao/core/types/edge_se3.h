//
// Created by xiang on 24-6-4.
//

#ifndef MIAO_EDGE_SE3_H
#define MIAO_EDGE_SE3_H

#include "core/common/math.h"
#include "core/graph/base_binary_edge.h"
#include "core/types/vertex_se3.h"

namespace lightning::miao {

/**
 * 两条SE3边之间的约束
 * 观测 = Tw1.inv * Tw2 = T12 =
 *     = [R1^T R2   R1^T t2 - R1^T t1 ]
 *     = [0^T       1                 ]
 *     -> [obs_R, obs_t]
 *        [0^T  , 1    ]
 *
 * 误差平移在前（与g2o保持一致，否则pose graph的info就不对）
 * 令 delta = measurement.inverse() * T12，平移残差是 delta.translation()，
 * 旋转残差是 delta 的单位四元数 xyz 部分，小角度下约为半个旋转向量。
 * @see @ref pgo_edge_residual_reading "二元边与先验边的残差尺度"
 *
 */
class EdgeSE3 : public BaseBinaryEdge<6, SE3, VertexSE3, VertexSE3> {
   public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /// 为了读g2o文件，这里的算法必须与g2o一致。

    // [pgo-binary-edge-residual]
    void ComputeError() override {
        SE3 from = ((VertexSE3*)(vertices_[0]))->Estimate();
        SE3 to = ((VertexSE3*)(vertices_[1]))->Estimate();

        SE3 T12 = from.inverse() * to;

        SE3 delta = measurement_.inverse() * T12;
        error_.head<3>() = delta.translation();
        error_.tail<3>() = delta.so3().unit_quaternion().coeffs().head<3>();
    }
    // [pgo-binary-edge-residual]
};

}  // namespace lightning::miao

#endif  // MIAO_EDGE_SE3_H
