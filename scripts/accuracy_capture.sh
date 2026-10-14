#!/usr/bin/env bash
# Local-only data capture and isolated, input-allowlist ROS 2 replay.
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  accuracy_capture.sh record --config FILE --role reference|algorithm|combined|outputs --output NEW_DIR [--duration-seconds N] [--dry-run]
  accuracy_capture.sh playback --config FILE --bag INPUT_BAG --output NEW_DIR --domain-id N [--post-wait-seconds N] [--dry-run]

Run record locally on BOTH computers using the same source_recording_id in
confirmed config copies. reference records CGI topics plus receive-only candump
on the existing CAN interface; algorithm records listed sensor inputs and outputs.
combined is optional only after verifying cross-host DDS discovery and clock audit.
Copy each complete local session directory to the offline analysis host afterwards.

Source ROS 2 Humble and the workspace overlay first. This script never starts an
algorithm/CGI driver, changes CAN settings, transmits CAN, calls save services or
replays recorded output/TF/CGI topics. Keep original bags and can_raw.log immutable.

playback starts an output recorder and replays ONLY capture.input_topics at 1x,
with /clock, ROS_LOCALHOST_ONLY=1 and the explicitly chosen isolated ROS_DOMAIN_ID.
Start the experiment algorithm separately in that same isolated domain with
use_sim_time=true. Choose the input/config for the intended mapping/localization
phase. All capture output files are new; existing directories are refused.

--duration-seconds 0 (default): record until Ctrl-C; orderly rosbag finalization.
--post-wait-seconds 5 (default): record output tail after input playback finishes.
--dry-run: validate config and print topic/command plan; no ROS or CAN operations.
Exit codes: 0 clean stop, 2 contract/argument/preflight failure, 3 child/record failure.
USAGE
}

[[ $# -eq 0 ]] && { usage; exit 2; }
[[ "$1" == "--help" || "$1" == "-h" ]] && { usage; exit 0; }
mode="$1"
shift
config="" role="" output="" bag="" domain="" duration=0 post_wait=5 dry_run=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --help|-h) usage; exit 0 ;;
    --config|--role|--output|--bag|--domain-id|--duration-seconds|--post-wait-seconds)
      [[ $# -ge 2 ]] || { echo "missing value for $1" >&2; exit 2; }
      key="$1" value="$2"
      case "$key" in
        --config) config="$value" ;; --role) role="$value" ;;
        --output) output="$value" ;; --bag) bag="$value" ;;
        --domain-id) domain="$value" ;; --duration-seconds) duration="$value" ;;
        --post-wait-seconds) post_wait="$value" ;;
      esac
      shift 2 ;;
    --dry-run) dry_run=1; shift ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
[[ "$mode" == "record" || "$mode" == "playback" ]] || { usage; exit 2; }
[[ -n "$config" && -n "$output" ]] || { usage; exit 2; }
[[ "$duration" =~ ^[0-9]+$ && "$post_wait" =~ ^[0-9]+$ ]] || { echo "durations must be nonnegative integer seconds" >&2; exit 2; }
if [[ "$mode" == "playback" ]]; then
  [[ -n "$bag" && "$domain" =~ ^[0-9]+$ && "$domain" -ge 1 && "$domain" -le 232 ]] || { echo "playback requires --bag and an isolated --domain-id 1..232" >&2; exit 2; }
  role=playback
  export ROS_DOMAIN_ID="$domain" ROS_LOCALHOST_ONLY=1
else
  [[ "$role" =~ ^(reference|algorithm|combined|outputs)$ ]] || { echo "record requires --role" >&2; exit 2; }
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
helper="$script_dir/accuracy_export_bag.py"
python="${ACCURACY_PYTHON:-python3}"
config="$(realpath "$config")"
[[ ! -e "$output" ]] || { echo "output already exists: $output" >&2; exit 2; }
"$python" "$helper" check-config --config "$config" --purpose capture --role "$role"
if [[ "$dry_run" == 1 ]]; then
  "$python" - "$config" "$role" "$output" "$bag" <<'PY'
import sys, yaml
c = yaml.safe_load(open(sys.argv[1], encoding='utf-8'))
print('role=' + sys.argv[2], 'new_output=' + sys.argv[3], 'source_recording_id=' + c['source_recording_id'])
print('input replay allowlist:', *c['capture']['input_topics'], sep='\n  ')
print('output streams:', *c['capture']['output_topics'], sep='\n  ')
print('CAN preflight is read-only; reference/combined also require candump and current UP/bitrate match.')
if sys.argv[2] == 'playback':
    print('playback bag:', sys.argv[4], 'rate=1.0, /clock enabled; output/TF/CGI topics excluded')
PY
  exit 0
fi
command -v ros2 >/dev/null || { echo "source ROS 2 Humble and workspace overlay first" >&2; exit 2; }
[[ "${ROS_DISTRO:-}" == "humble" ]] || { echo "this capture contract is verified for ROS_DISTRO=humble" >&2; exit 2; }
"$python" "$helper" plan --config "$config" --role "$role" --output "$output"
output="$(realpath "$output")"
mkdir "$output/logs"
metadata="$output/capture.metadata.json"
if [[ "$mode" == "playback" ]]; then
  "$python" - "$metadata" "$bag" <<'PY'
import json, pathlib, sys, yaml
path = pathlib.Path(sys.argv[1])
meta = json.loads(path.read_text())
bag = pathlib.Path(sys.argv[2]).resolve()
info = yaml.safe_load((bag / 'metadata.yaml').read_text())['rosbag2_bagfile_information']
start = int(info['starting_time']['nanoseconds_since_epoch'])
meta['playback_source_bag'] = str(bag)
meta['playback_source_interval'] = {'start_ns': start, 'end_ns': start + int(info['duration']['nanoseconds'])}
path.write_text(json.dumps(meta, indent=2) + '\n')
PY
fi
mapfile -t topics < "$output/topics.txt"
mapfile -t required < "$output/required_topics.txt"
mapfile -t inputs < "$output/input_topics.txt"
storage="$("$python" -c 'import json,sys; print(json.load(open(sys.argv[1]))["storage_id"])' "$metadata")"

recorder_pid="" raw_pid="" clock_pid="" player_pid="" can_interface="" started=0 cleanup_done=0 forced=0 child_failed=0 requested_interrupt=0
stop_child() {
  local pid="$1" name="$2" signaled=0 status=0
  [[ -n "$pid" ]] || return 0
  if kill -0 "$pid" 2>/dev/null; then
    if kill -INT "$pid" 2>/dev/null; then signaled=1; fi
    for ((i=0; i<40; i++)); do
      kill -0 "$pid" 2>/dev/null || break
      sleep 0.25
    done
    if kill -0 "$pid" 2>/dev/null; then
      forced=1
      kill -TERM "$pid" 2>/dev/null || true
    fi
  fi
  if wait "$pid" 2>/dev/null; then status=0; else status=$?; fi
  printf '%s\t%s\n' "$name" "$status" >> "$output/logs/child_exit_status.tsv"
  if [[ "$status" != 0 ]]; then
    if [[ "$status" != 130 || ( "$signaled" == 0 && "$requested_interrupt" == 0 ) ]]; then
      child_failed=1
    fi
  fi
}
cleanup() {
  local code="$?"
  [[ "$cleanup_done" == 0 ]] || return
  cleanup_done=1
  trap - EXIT INT TERM
  stop_child "$player_pid" player
  stop_child "$recorder_pid" recorder
  stop_child "$raw_pid" candump
  # Python background jobs can inherit ignored SIGINT from bash; TERM is normal
  # for this flushed per-second clock log and does not affect bag finalization.
  if [[ -n "$clock_pid" ]]; then
    kill -TERM "$clock_pid" 2>/dev/null || true
    clock_status=0
    if wait "$clock_pid" 2>/dev/null; then clock_status=0; else clock_status=$?; fi
    printf 'clock\t%s\n' "$clock_status" >> "$output/logs/child_exit_status.tsv"
    [[ "$clock_status" == 0 || "$clock_status" == 143 ]] || child_failed=1
  fi
  [[ "$forced" == 0 && "$child_failed" == 0 ]] || code=3
  [[ "$started" == 1 ]] || code=2
  if [[ "$started" == 1 ]]; then
    [[ -f "$output/bag/metadata.yaml" ]] || code=3
    if [[ "$role" == "reference" || "$role" == "combined" ]]; then
      [[ -s "$output/can_raw.log" ]] || code=3
    fi
  fi
  if [[ -n "$can_interface" ]]; then
    ip -details -statistics -json link show dev "$can_interface" > "$output/can_link_end.json" 2> "$output/logs/can_link_end.stderr.log" || code=3
  fi
  "$python" "$helper" mark --metadata "$metadata" --event end --exit-code "$code" || code=3
  if [[ "$code" == 0 ]]; then
    echo "capture finalized: $output"
  else
    echo "capture requires review (exit $code): $output/logs" >&2
  fi
  exit "$code"
}
trap cleanup EXIT
trap 'requested_interrupt=1; exit 0' INT
trap 'exit 3' TERM

# Record the actual visible graph and publishers; no remote commands are run.
timeout 20 ros2 topic list > "$output/topics_observed.txt"
"$python" - "$output/required_topics.txt" "$output/topics_observed.txt" <<'PY'
import sys
required = set(open(sys.argv[1]).read().splitlines())
observed = set(open(sys.argv[2]).read().splitlines())
missing = sorted(required - observed)
if missing:
    raise SystemExit('required local topics missing: ' + ', '.join(missing))
PY
for topic in "${required[@]}"; do
  name="${topic//\//_}"
  timeout 15 ros2 topic info --verbose "$topic" > "$output/logs/topic${name}.txt"
  if ! grep -Eq '^Publisher count: [1-9][0-9]*$' "$output/logs/topic${name}.txt"; then
    echo "required topic has no observed publisher: $topic" >&2
    exit 2
  fi
done
if [[ "$role" == "reference" || "$role" == "combined" ]]; then
  command -v candump >/dev/null || { echo "candump (can-utils) is required to retain all CAN frames" >&2; exit 2; }
  can_interface="$("$python" -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"]["capture"]["can_interface"])' "$metadata")"
  ip -details -statistics -json link show dev "$can_interface" > "$output/can_link_start.json"
  "$python" - "$metadata" "$output/can_link_start.json" <<'PY'
import json, sys
c = json.load(open(sys.argv[1]))['config']['capture']
links = json.load(open(sys.argv[2]))
if len(links) != 1:
    raise SystemExit('expected one existing CAN interface')
link = links[0]
info = link.get('linkinfo', {})
data = info.get('info_data', {})
if info.get('info_kind') != 'can' or 'UP' not in link.get('flags', []) or data.get('state') in ('STOPPED', 'BUS-OFF'):
    raise SystemExit('CAN is not UP and usable; recording refuses to configure it')
actual = data.get('bittiming', {}).get('bitrate')
if actual != c['expected_can_bitrate']:
    raise SystemExit(f'CAN bitrate {actual} != verified expected {c["expected_can_bitrate"]}; no changes made')
PY
  # Explicit zero mask observes all CAN IDs; error mask also retains error frames.
  # -d reports socket drops. Receiving never changes the existing CAN settings.
  candump -L -d "$can_interface,0:0,#FFFFFFFF" > "$output/can_raw.log" 2> "$output/logs/candump.stderr.log" &
  raw_pid=$!
fi

"$python" "$helper" clock-log --output "$output/clock_samples.jsonl" > "$output/logs/clock.stdout.log" 2> "$output/logs/clock.stderr.log" &
clock_pid=$!
"$python" "$helper" mark --metadata "$metadata" --event start
started=1
ros2 bag record --storage "$storage" --output "$output/bag" --qos-profile-overrides-path "$output/qos.yaml" "${topics[@]}" > "$output/logs/recorder.stdout.log" 2> "$output/logs/recorder.stderr.log" &
recorder_pid=$!
sleep 2
kill -0 "$recorder_pid" 2>/dev/null || { echo "rosbag recorder exited during startup" >&2; exit 3; }
if [[ -n "$raw_pid" ]]; then
  kill -0 "$raw_pid" 2>/dev/null || { echo "candump exited during startup" >&2; exit 3; }
fi
if [[ "$mode" == "playback" ]]; then
  ros2 bag play "$bag" --rate 1.0 --clock --disable-keyboard-controls --topics "${inputs[@]}" > "$output/logs/playback.stdout.log" 2> "$output/logs/playback.stderr.log" &
  player_pid=$!
  player_status=0
  if wait "$player_pid"; then player_status=0; else player_status=$?; fi
  printf 'player\t%s\n' "$player_status" >> "$output/logs/child_exit_status.tsv"
  player_pid=""
  [[ "$player_status" == 0 ]] || { echo "input playback failed" >&2; exit 3; }
  sleep "$post_wait"
else
  begin=$SECONDS
  while kill -0 "$recorder_pid" 2>/dev/null; do
    kill -0 "$clock_pid" 2>/dev/null || { echo "clock audit stopped before recording ended" >&2; exit 3; }
    if [[ -n "$raw_pid" ]] && ! kill -0 "$raw_pid" 2>/dev/null; then
      echo "candump stopped before recording ended" >&2
      exit 3
    fi
    [[ "$duration" == 0 || $((SECONDS - begin)) -lt "$duration" ]] || break
    sleep 1
  done
fi
kill -0 "$recorder_pid" 2>/dev/null || { echo "recorder exited before requested stop" >&2; exit 3; }
exit 0
