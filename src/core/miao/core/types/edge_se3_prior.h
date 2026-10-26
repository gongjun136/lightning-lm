//
// Created by xiang on 25-7-14.
//

#ifndef MIAO_EDGE_SE3_PRIOR_H
#define MIAO_EDGE_SE3_PRIOR_H

#include "core/common/math.h"
#include "core/graph/base_unary_edge.h"
#include "core/types/vertex_se3.h"

namespace lightning::miao {

/**
 * 先验的SE3观测
 *
 * 误差平移在前（与g2o保持一致，否则pose graph的info就不对）
 * 平移残差 = T.translation - measurement.translation。
 * 旋转残差 = Log(T.rotation.inverse() * measurement.rotation)，单位为 rad。
 * @see @ref pgo_edge_residual_reading "二元边与先验边的残差尺度"
 *
 */
class EdgeSE3Prior : public BaseUnaryEdge<6, SE3, VertexSE3> {
   public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // [pgo-prior-edge-residual]
    void ComputeError() override {
        SE3 pose = ((VertexSE3*)(vertices_[0]))->Estimate();
        error_.head<3>() = pose.translation() - measurement_.translation();
        error_.tail<3>() = (pose.so3().inverse() * measurement_.so3()).log();
    }
    // [pgo-prior-edge-residual]

    // TODO: jacobian计算
};

}  // namespace lightning::miao
#endif  // LIGHTNING_EDGE_SE3_PRIOR_H
