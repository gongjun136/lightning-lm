#!/usr/bin/env bash
# Single operator-facing entry point for SANY online localization.
set +e
set +u

WS=/home/nvidia/project/gj_ws
REPO="${WS}/lightning-lm"
MAP_DIR="${WS}/maps/sany_4lidar_20260811_a4defa3_solid"

export LIGHTNING_LM_INSTALL_SETUP="${REPO}/install/setup.bash"
export SANY_MAP_PATH="${MAP_DIR}"
export LIGHTNING_LM_OUT_ROOT="${WS}/runs"
export SANY_LIDAR_LAYOUT=3
export SANY_ENABLE_CAN_OBSERVATION=1
export SANY_WHEEL_SPEED_TOPIC=/SpeThrCAN4_topic
export SANY_RECORD_BAG=0
export SANY_POSRES_TIMEOUT_SECONDS=5
export SANY_MIN_FREE_GB=20
export SANY_ENABLE_POSRES_WATCHDOG=0
export LIGHTNING_LM_RUN_MODE=diagnostic
export LIGHTNING_LM_COMPUTE_PROFILE=1
export LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD=0
unset LIGHTNING_LM_CONFIG

cd "${REPO}" || exit 1

exec bash scripts/run_sany_online_diagnostics.sh \
  "sany_4lidar_can_on_$(date +%Y%m%d_%H%M%S)"
