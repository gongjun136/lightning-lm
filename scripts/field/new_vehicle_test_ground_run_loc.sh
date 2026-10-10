#!/usr/bin/env bash
# New vehicle / four-LiDAR test-ground entry; copy to ${SANY_WS}/run_loc.sh.
# Start all four LiDAR drivers separately before this entry.
set -euo pipefail
source "$HOME/.config/sany/workspace.bash"
export ROS_DOMAIN_ID="${ROS_DOMAIN_ID:-42}"

export LIGHTNING_LM_CONFIG="$LIGHTNING_LM_REPO_DIR/config/reproduction/multi_lidar/sany_4livox/sany_4lidar_new_vehicle_test_ground_localization_solid.yaml"
export SANY_MAP_PATH="$SANY_WS/maps/new_vehicle_4lidar_20260922_201743_calib_user_20261010_v1"
export SANY_LIDAR_LAYOUT=4
export LIGHTNING_LM_RUN_MODE=production
export SANY_RECORD_BAG=0
# This vehicle's CAN speed conversion has not been provided or validated.
export SANY_ENABLE_CAN_OBSERVATION=0
exec bash "$LIGHTNING_LM_REPO_DIR/scripts/run_sany_lidar_loc.sh" "$@"
