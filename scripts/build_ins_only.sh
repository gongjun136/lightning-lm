#!/usr/bin/env bash
set -euo pipefail
repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
ws="${LIGHTNING_LM_BUILD_WS:-$(dirname -- "$repo")/lightning_lm_ws}"
sdk="${CGI430_SDK_WS:-$(dirname -- "$repo")/sdk/cgi430_sdk}"
jobs="${LIGHTNING_BUILD_JOBS:-2}"
[[ "$jobs" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid LIGHTNING_BUILD_JOBS' >&2; exit 2; }
[[ -d "$sdk/src/cgi430_interfaces" && -d "$ws/src" ]] || { echo 'Set CGI430_SDK_WS and LIGHTNING_LM_BUILD_WS to existing workspaces' >&2; exit 2; }
set +u
source /opt/ros/humble/setup.bash
set -u
cd "$sdk"
colcon build --packages-select cgi430_interfaces --executor sequential --event-handlers desktop_notification- --cmake-args -DCMAKE_BUILD_TYPE=Release
set +u
source "$sdk/install/local_setup.bash"
set -u
cd "$ws"
CMAKE_BUILD_PARALLEL_LEVEL="$jobs" colcon build --packages-up-to lightning_lm --executor sequential --event-handlers desktop_notification- --cmake-args -DBUILD_TESTING=ON "-Dcgi430_interfaces_DIR=$sdk/install/cgi430_interfaces/share/cgi430_interfaces/cmake"
