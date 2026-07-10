#!/usr/bin/env bash
# ROS2 online/replay frontend template. The YAML is intentionally supplied from
# outside the source tree through LIGHTNING_LM_CONFIG.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="${LIGHTNING_LM_REPO_DIR:-$(cd "$script_dir/.." && pwd)}"
ros_setup="${LIGHTNING_LM_ROS_SETUP:-/opt/ros/humble/setup.bash}"
install_setup="${LIGHTNING_LM_INSTALL_SETUP:-$repo_dir/install/setup.bash}"
config_path="${LIGHTNING_LM_CONFIG:-}"

if [[ -z "$config_path" || ! -f "$config_path" ]]; then
  echo "LIGHTNING_LM_CONFIG is required and must exist: $config_path" >&2
  exit 2
fi

set +u
source "$ros_setup"
source "$install_setup"
set -u

exec ros2 run lightning run_frontend_online --config "$(realpath "$config_path")" "$@"
