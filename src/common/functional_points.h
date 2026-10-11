//
// Created by xiang on 23-10-12.
//
#pragma once

#include "common/eigen_types.h"

namespace lightning {

/// 功能点定义
struct FunctionalPoint {
    FunctionalPoint() {}
    FunctionalPoint(const std::string& name, const SE3& pose) : name_(name), pose_(pose) {}

    std::string name_;  ///< 功能点的名称，用于选择定位初始化位置。
    SE3 pose_;  ///< 该功能点在地图坐标系中的位姿。
};

}  // namespace lightning
