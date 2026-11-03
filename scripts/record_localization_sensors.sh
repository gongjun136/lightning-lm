#!/usr/bin/env bash
# 主 Orin 的原始传感器与精度评价采集入口；不启动算法或修改 CAN。
set -eo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "$HOME/.config/sany/workspace.bash" ]]; then
  source "$HOME/.config/sany/workspace.bash"
fi
source "${LIGHTNING_LM_ROS_SETUP:-/opt/ros/humble/setup.bash}"
setup="${LIGHTNING_LM_INSTALL_SETUP:-${LIGHTNING_LM_WS:-$script_dir/../..}/install/setup.bash}"
if [[ ! -f "$setup" ]]; then
  printf '未找到工作空间环境：%s\n' "$setup" >&2
  exit 2
fi
source "$setup"
exec python3 "$script_dir/record_localization_sensors.py" "$@"
