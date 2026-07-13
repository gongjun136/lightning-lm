#!/usr/bin/env bash
# Formal 1x online interface validation for the five public frontend topics.
set -euo pipefail

bag="${1:?usage: run_sany_online_contract.sh BAG CONFIG OUTPUT_DIR [FAULT_CONTRACT]}"
config="${2:?missing external YAML config}"
output_dir="${3:?missing output directory}"
fault_contract="${4:-}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
default_repo="$(cd "$script_dir/../../../.." && pwd)"
repo="${LIGHTNING_LM_REPO_DIR:-$default_repo}"
script_file="$(realpath "${BASH_SOURCE[0]}")"
bag="$(realpath "$bag")"
config="$(realpath "$config")"
repo="$(realpath "$repo")"
output_dir="$(realpath -m "$output_dir")"
git_pointer="$(sed -n 's/^gitdir: //p' "$repo/.git" 2>/dev/null || true)"
if [[ -n "$git_pointer" && "$git_pointer" =~ ^([A-Za-z]):/(.*)$ ]]; then
  git_drive="$(printf '%s' "${BASH_REMATCH[1]}" | tr '[:upper:]' '[:lower:]')"
  git_dir="/mnt/$git_drive/${BASH_REMATCH[2]}"
elif [[ -n "$git_pointer" ]]; then
  git_dir="$(realpath -m "$repo/$git_pointer")"
else
  git_dir="$repo/.git"
fi
if [[ ! -d "$git_dir" ]]; then
  echo "cannot resolve Git metadata directory: $git_dir" >&2
  exit 2
fi
cpu_set="${BENCH_CPUSET:-0-7}"
ros_domain_id="${SANY_ROS_DOMAIN_ID:-217}"
record_dir="$output_dir/topic_record"
logs="$output_dir/logs"
analyzer="$repo/scripts/analyze_frontend_topic_bag.py"
fault_analyzer="$script_dir/analyze_sany_fault_recovery.py"
selected_analyzer="$analyzer"
if [[ -n "$fault_contract" ]]; then
  fault_contract="$(realpath "$fault_contract")"
  selected_analyzer="$fault_analyzer"
fi
template="$repo/scripts/run_frontend_online.sh"
binary="$repo/install/lightning/lib/lightning/run_frontend_online"
source_file="$repo/src/app/run_frontend_online.cc"

if [[ ! -f "$bag/metadata.yaml" || ! -f "$config" || ! -x "$template" || ! -f "$selected_analyzer" || ! -x "$binary" ]]; then
  echo "missing bag, config, online template, analyzer, or installed binary" >&2
  exit 2
fi
if [[ -n "$fault_contract" && ! -f "$fault_contract" ]]; then
  echo "fault contract does not exist: $fault_contract" >&2
  exit 2
fi
if [[ ! "$ros_domain_id" =~ ^[0-9]+$ ]] || (( ros_domain_id > 232 )); then
  echo "SANY_ROS_DOMAIN_ID must be an integer in [0,232]" >&2
  exit 2
fi
if [[ "$source_file" -nt "$binary" ]]; then
  echo "installed online binary is older than run_frontend_online.cc; rebuild and install first" >&2
  exit 2
fi
if [[ -e "$output_dir" ]]; then
  echo "refusing to overwrite: $output_dir" >&2
  exit 2
fi
case "$output_dir/" in
  "$repo/"*|"$bag/"*)
    echo "output directory must be outside the repository and input bag" >&2
    exit 2
    ;;
esac
mkdir -p "$logs"

set +u
source /opt/ros/humble/setup.bash
source "$repo/install/setup.bash"
set -u
export ROS_DOMAIN_ID="$ros_domain_id"
export ROS_LOCALHOST_ONLY=1
if [[ "$(realpath "$(ros2 pkg prefix lightning)")" != "$(realpath "$repo/install/lightning")" ]]; then
  echo "ros2 package prefix does not resolve to the requested install tree" >&2
  exit 2
fi
if [[ -n "$(git --git-dir="$git_dir" --work-tree="$repo" status --porcelain)" ]]; then
  echo "formal online validation requires a clean Lightning-LM worktree" >&2
  exit 2
fi
existing_nodes="$(timeout 5s ros2 node list 2>/dev/null || true)"
if [[ -n "$existing_nodes" ]]; then
  echo "dedicated ROS_DOMAIN_ID=$ROS_DOMAIN_ID is not isolated; existing nodes:" >&2
  printf '%s\n' "$existing_nodes" >&2
  exit 2
fi
export ROS_LOG_DIR="$output_dir/ros_logs"
mkdir -p "$ROS_LOG_DIR"
cd "$output_dir"
{
  printf 'ros_distro=%s\n' "${ROS_DISTRO:-}"
  printf 'package_prefix=%s\n' "$(ros2 pkg prefix lightning)"
  printf 'ros_domain_id=%s\n' "$ROS_DOMAIN_ID"
  printf 'ros_localhost_only=%s\n' "$ROS_LOCALHOST_ONLY"
  printf 'kernel='; uname -a
  ldd "$binary" | sed -E 's/ \(0x[0-9a-f]+\)//g'
} > "$output_dir/environment_identity.txt"

bag_duration_s="$(python3 - "$bag/metadata.yaml" <<'PY'
import sys, yaml
data=yaml.safe_load(open(sys.argv[1], encoding="utf-8"))
print(data["rosbag2_bagfile_information"]["duration"]["nanoseconds"] / 1e9)
PY
)"
player_timeout_s="$(python3 -c 'import math,sys; print(math.ceil(float(sys.argv[1])+60.0))' "$bag_duration_s")"

python3 - "$output_dir/run_metadata.json" "$bag" "$config" "$cpu_set" "$binary" "$selected_analyzer" "$template" "$repo" "$script_file" "$bag_duration_s" "$output_dir/environment_identity.txt" "$fault_contract" "$ROS_DOMAIN_ID" "$git_dir" <<'PY'
import hashlib, json, os, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path

out, bag, config, cpus, binary, analyzer, template, repo, controller, duration, environment, fault_contract, ros_domain_id, git_dir = sys.argv[1:]
def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()
def identity(path):
    path=Path(path).resolve()
    return {"path":str(path),"size":path.stat().st_size,"sha256":digest(path)}
bag_path=Path(bag).resolve()
bag_files=sorted(p for p in bag_path.iterdir() if p.is_file() and (p.name=="metadata.yaml" or p.suffix in {".db3",".mcap"}))
git_args=["git",f"--git-dir={git_dir}",f"--work-tree={repo}"]
git_head=subprocess.run(git_args+["rev-parse","HEAD"],check=True,text=True,capture_output=True).stdout.strip()
git_status=subprocess.run(git_args+["status","--porcelain"],check=True,text=True,capture_output=True).stdout.splitlines()
fingerprint_payload={
    "bag":{"path":str(bag_path),"files":[identity(path) for path in bag_files]},
    "config":identity(config),"binary":identity(binary),"analyzer":identity(analyzer),
    "template":identity(template),"controller":identity(controller),"git_head":git_head,
    "environment":identity(environment),
    "git_status_porcelain":git_status,"cpu_set":cpus,"play_rate":1.0,
    "clock_publish_frequency_hz":100.0,"bag_duration_s":float(duration),
    "ros_domain_id":int(ros_domain_id),"ros_localhost_only":True,
}
if fault_contract:
    contract=json.loads(Path(fault_contract).read_text(encoding="utf-8"))
    if contract.get("complete_input_pass") is not True:
        raise SystemExit("fault contract is not from a complete input pass")
    if Path(contract.get("output_bag", "")).resolve() != bag_path:
        raise SystemExit("fault contract output_bag does not match the actual playback bag")
    actual_files=[identity(path) for path in bag_files]
    contracted_files=contract.get("bag_files")
    if not isinstance(contracted_files, list) or actual_files != contracted_files:
        raise SystemExit("fault contract bag file hashes do not match the actual playback bag")
    fingerprint_payload["fault_contract"]=identity(fault_contract)
fingerprint=hashlib.sha256(json.dumps(fingerprint_payload,sort_keys=True,ensure_ascii=False).encode()).hexdigest()
payload={
    "status":"running","started_at":datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
    "experiment_fingerprint":fingerprint,"fingerprint_payload":fingerprint_payload,
    "record_started_before_node":True,
    "recorded_topics":["/slamPoseRaw_topic","/final_points_topic","/slamSafety_topic","/slamState_topic","/SystemState","/clock"],
}
path=Path(out); temporary=path.with_suffix(path.suffix+".tmp")
temporary.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+"\n",encoding="utf-8"); os.replace(temporary,path)
PY

recorder_pid=""
node_pid=""
player_pid=""
stop_group() {
  local pid="${1:-}"
  [[ -z "$pid" ]] && return 0
  if kill -0 "$pid" 2>/dev/null; then
    kill -INT -- "-$pid" 2>/dev/null || true
    for _ in {1..50}; do
      kill -0 "$pid" 2>/dev/null || return 0
      sleep 0.1
    done
    kill -TERM -- "-$pid" 2>/dev/null || true
    for _ in {1..50}; do
      kill -0 "$pid" 2>/dev/null || return 0
      sleep 0.1
    done
    kill -KILL -- "-$pid" 2>/dev/null || true
    for _ in {1..20}; do
      kill -0 "$pid" 2>/dev/null || return 0
      sleep 0.1
    done
  fi
}
cleanup() {
  local exit_code=$?
  stop_group "$player_pid"
  stop_group "$node_pid"
  stop_group "$recorder_pid"
  if [[ "$exit_code" -ne 0 && -f "$output_dir/run_metadata.json" ]]; then
    python3 - "$output_dir/run_metadata.json" "$exit_code" <<'PY' || true
import json, os, sys
from datetime import datetime, timezone
from pathlib import Path
path=Path(sys.argv[1]); payload=json.loads(path.read_text(encoding="utf-8"))
if payload.get("status") == "running":
    payload.update({"status":"aborted","aborted_at":datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),"controller_exit_code":int(sys.argv[2])})
    temporary=path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    os.replace(temporary,path)
PY
  fi
  return "$exit_code"
}
trap cleanup EXIT

graph_count() {
  local topic="$1"
  local label="$2"
  ros2 topic info "$topic" 2>/dev/null | awk -F ': ' -v label="$label" '$1 == label {print $2; exit}'
}
at_least_one() {
  local value="${1:-0}"
  [[ "$value" =~ ^[0-9]+$ ]] && (( value >= 1 ))
}

setsid taskset -c "$cpu_set" ros2 bag record \
  --storage sqlite3 --output "$record_dir" \
  /slamPoseRaw_topic /final_points_topic /slamSafety_topic /slamState_topic /SystemState /clock \
  >"$logs/recorder.stdout.log" 2>"$logs/recorder.stderr.log" &
recorder_pid=$!
sleep 2
if ! kill -0 "$recorder_pid" 2>/dev/null; then
  echo "recorder exited before node startup" >&2
  wait "$recorder_pid" || true
  exit 3
fi

setsid taskset -c "$cpu_set" /usr/bin/time -v -o "$output_dir/node_time.txt" \
  env LIGHTNING_LM_REPO_DIR="$repo" \
      LIGHTNING_LM_INSTALL_SETUP="$repo/install/setup.bash" \
      LIGHTNING_LM_CONFIG="$config" \
      "$template" \
  >"$logs/node.stdout.log" 2>"$logs/node.stderr.log" &
node_pid=$!
ready=false
readiness_deadline=$((SECONDS + 30))
while (( SECONDS < readiness_deadline )); do
  if ! kill -0 "$node_pid" 2>/dev/null || ! kill -0 "$recorder_pid" 2>/dev/null; then
    break
  fi
  if at_least_one "$(graph_count /livox/lidar_192_168_1_114 'Subscription count')" \
    && at_least_one "$(graph_count /livox/lidar_192_168_1_127 'Subscription count')" \
    && at_least_one "$(graph_count /livox/lidar_192_168_1_187 'Subscription count')" \
    && at_least_one "$(graph_count /livox/lidar_192_168_1_195 'Subscription count')" \
    && at_least_one "$(graph_count /livox/imu_192_168_1_114 'Subscription count')" \
    && at_least_one "$(graph_count /slamPoseRaw_topic 'Publisher count')" \
    && at_least_one "$(graph_count /slamPoseRaw_topic 'Subscription count')" \
    && at_least_one "$(graph_count /final_points_topic 'Publisher count')" \
    && at_least_one "$(graph_count /final_points_topic 'Subscription count')" \
    && at_least_one "$(graph_count /slamSafety_topic 'Subscription count')" \
    && at_least_one "$(graph_count /slamState_topic 'Subscription count')" \
    && at_least_one "$(graph_count /SystemState 'Subscription count')"; then
    ready=true
    break
  fi
  sleep 0.1
done
if [[ "$ready" != true ]]; then
  echo "frontend/recorder ROS graph did not become ready before playback" >&2
  exit 4
fi

setsid taskset -c "$cpu_set" timeout --signal=INT --kill-after=10s "${player_timeout_s}s" \
  ros2 bag play "$bag" --rate 1.0 --clock 100.0 --read-ahead-queue-size 1000 \
  >"$logs/player.stdout.log" 2>"$logs/player.stderr.log" &
player_pid=$!
set +e
wait "$player_pid"
player_rc=$?
set -e
player_pid=""
sleep 3

stop_group "$node_pid"
set +e
wait "$node_pid"
node_rc=$?
set -e
node_pid=""

stop_group "$recorder_pid"
set +e
wait "$recorder_pid"
recorder_rc=$?
set -e
recorder_pid=""

set +e
if [[ -n "$fault_contract" ]]; then
  python3 "$fault_analyzer" "$record_dir" \
    --fault-contract "$fault_contract" \
    --input-bag "$bag" \
    --output "$output_dir/fault_recovery_contract.json"
  analyzer_rc=$?
else
  python3 "$analyzer" "$record_dir" \
    --output-json "$output_dir/topic_contract.json" \
    --output-csv "$output_dir/topic_timing.csv" \
    --enforce
  analyzer_rc=$?
fi
set -e

python3 - "$output_dir/run_metadata.json" "$output_dir/artifact_manifest.json" "$output_dir" \
  "$player_rc" "$node_rc" "$recorder_rc" "$analyzer_rc" <<'PY'
import hashlib, json, os, sys
from datetime import datetime, timezone
from pathlib import Path

metadata_path, manifest_path, root, player, node, recorder, analyzer = sys.argv[1:]
metadata_path=Path(metadata_path); manifest_path=Path(manifest_path); root=Path(root)
exit_codes={"player":int(player),"node":int(node),"recorder":int(recorder),"analyzer":int(analyzer)}
passed=all(value == 0 for value in exit_codes.values())
payload=json.loads(metadata_path.read_text(encoding="utf-8"))
payload.update({
    "status":"completed" if passed else "failed",
    "finished_at":datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
    "exit_codes":exit_codes,
})
def atomic(path,value):
    temporary=path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(value,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    os.replace(temporary,path)
atomic(metadata_path,payload)
def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()
artifacts=[]
for path in sorted(root.rglob("*")):
    if path.is_file() and path != manifest_path and not path.name.endswith(".tmp"):
        artifacts.append({"path":str(path.relative_to(root)),"size":path.stat().st_size,"sha256":digest(path)})
manifest={
    "created_at":datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
    "experiment_fingerprint":payload["experiment_fingerprint"],"contract_passed":passed,
    "exit_codes":exit_codes,"artifacts":artifacts,
}
atomic(manifest_path,manifest)
PY

if [[ "$player_rc" -ne 0 || "$node_rc" -ne 0 || "$recorder_rc" -ne 0 || "$analyzer_rc" -ne 0 ]]; then
  echo "online contract failed: player=$player_rc node=$node_rc recorder=$recorder_rc analyzer=$analyzer_rc" >&2
  exit 5
fi
echo "online contract passed: $output_dir"
