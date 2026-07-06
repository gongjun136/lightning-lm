#!/usr/bin/env python3
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read(path: str) -> str:
    return (ROOT / path).read_text(encoding='utf-8')


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f'velocity propagation check failed: {message}')


nav_state = read('src/common/nav_state.h')
eskf_hpp = read('src/core/lio/eskf.hpp')
eskf_cc = read('src/core/lio/eskf.cc')
laser_mapping_cc = read('src/core/lio/laser_mapping.cc')
laser_mapping_h = read('src/core/lio/laser_mapping.h')
config = read('config/sany_livox_114_20260701_velprop_no_lba.yaml')

require(
    'vel_ += vec.middleRows(kVelIdx, kBlockDim) * dt;' in nav_state,
    'NavState::oplus must propagate velocity during IMU prediction',
)
require(
    '// vel_ += vec.middleRows(kVelIdx, kBlockDim) * dt;' not in nav_state,
    'velocity propagation must not remain commented out',
)
require(
    'bool propagate_velocity_ = false;' in eskf_hpp,
    'ESKF options must expose a default-off velocity propagation switch',
)
require(
    'if (!options_.propagate_velocity_)' in eskf_cc
    and 'segment<NavState::kBlockDim>(NavState::kVelIdx).setZero()' in eskf_cc,
    'ESKF::Predict must preserve legacy behavior when velocity propagation is disabled',
)
require(
    'eskf_options.propagate_velocity_ = propagate_velocity_;' in laser_mapping_cc,
    'LaserMapping must pass the YAML velocity propagation option into ESKF',
)
require(
    'propagate_velocity' in laser_mapping_cc,
    'LaserMapping must load fasterlio.propagate_velocity when present',
)
require(
    'bool propagate_velocity_ = false;' in laser_mapping_h,
    'LaserMapping must default velocity propagation to disabled',
)
require('propagate_velocity: true' in config, 'SANY 114 config must enable velocity propagation')

print('velocity propagation check passed')
