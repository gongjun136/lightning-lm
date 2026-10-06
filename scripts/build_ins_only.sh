#!/usr/bin/env bash
set -euo pipefail
repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
ws="${LIGHTNING_LM_BUILD_WS:-$(dirname -- "$repo")/lightning_lm_ws}"
sdk="${CGI430_SDK_WS:-$(dirname -- "$repo")/sdk/cgi430_sdk}"
jobs="${LIGHTNING_BUILD_JOBS:-2}"
[[ "$jobs" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid LIGHTNING_BUILD_JOBS' >&2; exit 2; }
cache_args=()
if [[ "${1:-}" == --cmake-clean-cache ]]; then cache_args+=(--cmake-clean-cache); shift; fi
[[ $# -eq 0 ]] || { echo 'Usage: build_ins_only.sh [--cmake-clean-cache]' >&2; exit 2; }
[[ -d "$sdk/src/cgi430_interfaces" && -d "$ws/src" ]] || { echo 'Set CGI430_SDK_WS and LIGHTNING_LM_BUILD_WS to existing workspaces' >&2; exit 2; }
# colcon adds its own Make -j/-l unless MAKEFLAGS is set; CMAKE_BUILD_PARALLEL_LEVEL alone is insufficient.
export MAKEFLAGS="-j$jobs -l$jobs"
export CMAKE_BUILD_PARALLEL_LEVEL="$jobs"
set +u
source /opt/ros/humble/setup.bash
set -u
cd "$sdk"
colcon build --packages-select cgi430_interfaces --executor sequential --event-handlers desktop_notification- "${cache_args[@]}" --cmake-args -DCMAKE_BUILD_TYPE=Release
set +u
source "$sdk/install/local_setup.bash"
set -u
cd "$ws"
colcon build --packages-up-to lightning_lm --executor sequential --event-handlers desktop_notification- "${cache_args[@]}" --cmake-args -DBUILD_TESTING=ON "-Dcgi430_interfaces_DIR=$sdk/install/cgi430_interfaces/share/cgi430_interfaces/cmake"
