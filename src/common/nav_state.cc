#include "common/nav_state.h"

namespace lightning {

const std::vector<NavState::MetaInfo> NavState::vect_states_{
    {NavState::kPosIdx, NavState::kPosIdx, NavState::kBlockDim},
    {NavState::kVelIdx, NavState::kVelIdx, NavState::kBlockDim},
    {NavState::kBgIdx, NavState::kBgIdx, NavState::kBlockDim},
    {NavState::kBaIdx, NavState::kBaIdx, NavState::kBlockDim},
    {NavState::kGravIdx, NavState::kGravIdx, NavState::kBlockDim},
};

const std::vector<NavState::MetaInfo> NavState::SO3_states_{
    {NavState::kRotIdx, NavState::kRotIdx, NavState::kBlockDim},
};

}  // namespace lightning
