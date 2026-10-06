#!/usr/bin/env bash
set -euo pipefail
repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
if [[ "$(basename -- "$(dirname -- "$repo")")" == src ]]; then
  default_ws="$(dirname -- "$(dirname -- "$repo")")"
else
  default_ws="$(dirname -- "$repo")/lightning_lm_ws"
fi
ws="${LIGHTNING_LM_BUILD_WS:-$default_ws}"
jobs="${LIGHTNING_BUILD_JOBS:-2}"
[[ "$jobs" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid LIGHTNING_BUILD_JOBS' >&2; exit 2; }
cache_args=()
if [[ "${1:-}" == --cmake-clean-cache ]]; then cache_args+=(--cmake-clean-cache); shift; fi
[[ $# -eq 0 ]] || { echo 'Usage: build_workspace.sh [--cmake-clean-cache]' >&2; exit 2; }
for package in common_msgs/cgi430_interfaces cgi430_driver livox_ros_driver2 Livox-SDK2 lightning-lm; do
  [[ -f "$ws/src/$package/package.xml" ]] || { echo "Missing workspace package: $ws/src/$package" >&2; exit 2; }
done
export MAKEFLAGS="-j$jobs -l$jobs"
export CMAKE_BUILD_PARALLEL_LEVEL="$jobs"
set +u
source /opt/ros/humble/setup.bash
set -u
cd "$ws"
# Drop the old external CGI interface cache entry; colcon supplies this workspace's dependency prefix.
colcon build --executor sequential --event-handlers desktop_notification- "${cache_args[@]}" \
  --cmake-args -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -Ucgi430_interfaces_DIR
