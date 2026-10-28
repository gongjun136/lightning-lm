#!/usr/bin/env bash
# Offline reads only. Sourcing the configured workspace does not start ROS nodes.
set -euo pipefail
script_dir="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]}")")"
# Apply scheduling before ROS setup invokes its Python helpers, too. exec keeps
# one short-lived analyzer process; this never changes localization's affinity.
if [[ "${LIGHTNING_DIAG_SCHEDULING_READY:-0}" != 1 ]]; then
  exec python3 -B - "${script_dir}" "$@" <<'PY'
import json, os, sys
sys.path.insert(0, sys.argv[1])
from diagnose_localization import low_impact
try:
    low_impact()
except RuntimeError as exc:
    print(json.dumps({'error': str(exc)}, ensure_ascii=False))
    sys.exit(2)
os.environ['LIGHTNING_DIAG_SCHEDULING_READY'] = '1'
os.execv('/bin/bash', ['bash', os.path.join(sys.argv[1], 'diagnose_localization.sh'), *sys.argv[2:]])
PY
fi
if [[ -r "${HOME}/.config/sany/workspace.bash" ]]; then
  source "${HOME}/.config/sany/workspace.bash"
fi
# Sidecars work without ROS; a bare MCAP needs its recorded message definitions.
set +u
if [[ -r /opt/ros/humble/setup.bash ]]; then source /opt/ros/humble/setup.bash; fi
if [[ -n "${LIGHTNING_LM_INSTALL_SETUP:-}" && -r "${LIGHTNING_LM_INSTALL_SETUP}" ]]; then
  source "${LIGHTNING_LM_INSTALL_SETUP}"
fi
set -u
export PYTHONDONTWRITEBYTECODE=1
exec python3 "${script_dir}/diagnose_localization.py" "$@"
