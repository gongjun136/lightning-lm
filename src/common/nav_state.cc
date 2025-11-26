//
// Created by xiang on 2022/6/21.
//

#include "common/nav_state.h"

namespace lightning {
/// 向量状态变量索引配置
const std::vector<NavState::MetaInfo> NavState::vect_states_{
    {0, 0, 3},    // pos (位置): 索引0，维度0，自由度3
    {9, 9, 3},    // offset_t (外参平移): 索引9，维度9，自由度3
    {12, 12, 3},  // vel (速度): 索引12，维度12，自由度3
    {15, 15, 3},  // bg (陀螺零偏): 索引15，维度15，自由度3
    {18, 18, 3},  // ba (加计零偏): 索引18，维度18，自由度3
};

/// SO3旋转状态变量索引配置
const std::vector<NavState::MetaInfo> NavState::SO3_states_{
    {3, 3, 3},  // rot (姿态旋转): 索引3，维度3，自由度3
    {6, 6, 3},  // offset_R_li (外参旋转): 索引6，维度6，自由度3
};

/// S2球面状态变量索引配置
const std::vector<NavState::MetaInfo> NavState::S2_states_{
    {21, 21, 3},  // grav (重力向量): 索引21，维度21，dof=3(源空间3维→误差空间2维)
};

}  // namespace lightning
