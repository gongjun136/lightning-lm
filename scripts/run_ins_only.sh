#!/usr/bin/env bash
set -euo pipefail
repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
[[ $# -ge 1 ]] || { echo "Usage: $0 /absolute/site.yaml [run_loc_online flags]" >&2; exit 2; }
config="$(realpath -- "$1")"
shift
if [[ "$(basename -- "$(dirname -- "$repo")")" == src ]]; then
  default_ws="$(dirname -- "$(dirname -- "$repo")")"
else
  default_ws="$(dirname -- "$repo")/lightning_lm_ws"
fi
install="${LIGHTNING_LM_INSTALL_SETUP:-${LIGHTNING_LM_BUILD_WS:-$default_ws}/install/local_setup.bash}"
set +u
source /opt/ros/humble/setup.bash
source "$install"
set -u
export ROS_DOMAIN_ID="${ROS_DOMAIN_ID:-42}"
python3 - "$config" <<'PY'
import sys,yaml
with open(sys.argv[1]) as f: c=yaml.safe_load(f)
if c.get('system',{}).get('localization_mode')!='ins_only':
    raise SystemExit('run_ins_only.sh requires system.localization_mode: ins_only')
PY
cd "$repo"
exec ros2 run lightning_lm run_loc_online --config "$config" "$@"
