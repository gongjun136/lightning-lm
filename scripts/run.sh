#!/usr/bin/env bash
# Single operator-facing entry point for SANY online localization.
set +e
set +u

export ROS_DOMAIN_ID="${ROS_DOMAIN_ID:-42}"

WS="${LIGHTNING_LM_WS:-/home/nvidia/project/gj_ws}"
REPO="${LIGHTNING_LM_REPO_DIR:-${WS}/lightning-lm}"
MAP_DIR="${WS}/maps/sany_4lidar_20260811_a4defa3_solid"

export LIGHTNING_LM_INSTALL_SETUP="${LIGHTNING_LM_INSTALL_SETUP:-${REPO}/install/setup.bash}"
export SANY_MAP_PATH="${SANY_MAP_PATH:-${MAP_DIR}}"
export LIGHTNING_LM_OUT_ROOT="${LIGHTNING_LM_OUT_ROOT:-${WS}/runs}"
export SANY_LIDAR_LAYOUT="${SANY_LIDAR_LAYOUT:-3}"
export SANY_ENABLE_CAN_OBSERVATION="${SANY_ENABLE_CAN_OBSERVATION:-1}"
export SANY_WHEEL_SPEED_TOPIC="${SANY_WHEEL_SPEED_TOPIC:-/SpeThrCAN4_topic}"
export SANY_RECORD_BAG="${SANY_RECORD_BAG:-0}"
export SANY_POSE_VEL_TIMEOUT_SECONDS="${SANY_POSE_VEL_TIMEOUT_SECONDS:-5}"
export SANY_MIN_FREE_GB="${SANY_MIN_FREE_GB:-20}"
export SANY_ENABLE_POSE_VEL_WATCHDOG="${SANY_ENABLE_POSE_VEL_WATCHDOG:-0}"
export LIGHTNING_LM_RUN_MODE="${LIGHTNING_LM_RUN_MODE:-diagnostic}"
case "${LIGHTNING_LM_RUN_MODE}" in
  diagnostic)
    export LIGHTNING_LM_COMPUTE_PROFILE="${LIGHTNING_LM_COMPUTE_PROFILE:-1}"
    export LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD="${LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD:-0}"
    ;;
  production)
    export LIGHTNING_LM_COMPUTE_PROFILE="${LIGHTNING_LM_COMPUTE_PROFILE:-0}"
    export LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD="${LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD:-1}"
    ;;
  *)
    echo "ERROR: LIGHTNING_LM_RUN_MODE must be diagnostic or production." >&2
    exit 1
    ;;
esac

cd "${REPO}" || exit 1

exec bash scripts/run_sany_online_diagnostics.sh \
  "sany_${SANY_LIDAR_LAYOUT}lidar_$(date +%Y%m%d_%H%M%S)"
