//
// Created by xiang on 24-4-23.
//

#ifndef MIAO_BASE_BINARY_EDGE_H
#define MIAO_BASE_BINARY_EDGE_H

#include "base_fixed_sized_edge.h"

namespace lightning::miao {

template <int D, typename E, typename VertexXi, typename VertexXj>
class BaseBinaryEdge : public BaseFixedSizedEdge<D, E, VertexXi, VertexXj> {
   public:
    using VertexXiType = VertexXi;
    using VertexXjType = VertexXj;
    BaseBinaryEdge() : BaseFixedSizedEdge<D, E, VertexXi, VertexXj>(){};

   protected:
    typename BaseFixedSizedEdge<D, E, VertexXi, VertexXj>::template JacobianType<D, VertexXi::Dimension>&
        jacobian_oplus_xi_ = std::get<0>(this->jacobian_oplus_);  ///< 边误差对第一个顶点局部增量的 Jacobian 引用。
    typename BaseFixedSizedEdge<D, E, VertexXi, VertexXj>::template JacobianType<D, VertexXj::Dimension>&
        jacobian_oplus_xj_ = std::get<1>(this->jacobian_oplus_);  ///< 边误差对第二个顶点局部增量的 Jacobian 引用。
};

}  // namespace lightning::miao
#endif  // MIAO_BASE_BINARY_EDGE_H
