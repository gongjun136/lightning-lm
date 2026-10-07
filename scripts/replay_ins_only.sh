#!/usr/bin/env bash
set -euo pipefail
[[ $# -ge 2 ]] || { echo "Usage: $0 site_replay.yaml bag [ros2 bag play options]; start run_sany_ins_only.sh separately" >&2; exit 2; }
dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
config="$1"; bag="$2"; shift 2
python3 - "$config" <<'PY'
import sys,yaml
with open(sys.argv[1]) as f: c=yaml.safe_load(f)
if not c['ins_only'].get('use_sim_time',False):
    raise SystemExit('Replay requires ins_only.use_sim_time: true')
PY
topic_text="$(python3 "$dir/ins_topics.py" "$config")"
mapfile -t topics <<< "$topic_text"
# Never replay old business outputs, TF or an old /clock from the input bag.
exec ros2 bag play "$bag" "$@" --clock --topics "${topics[@]}"
