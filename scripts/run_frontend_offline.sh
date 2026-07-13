#!/usr/bin/env bash
# Canonical monitored runner for every offline Lightning-LM frontend dataset.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_frontend_offline.sh --bag ROS2_BAG --config YAML --output-dir DIR [options]

Required:
  --bag PATH                    SQLite3 ROS2 bag directory
  --config PATH                 Lightning-LM YAML configuration
  --output-dir PATH             Exact run output directory

Options:
  --sequence NAME               Dataset sequence label (default: bag directory name)
  --repeat INDEX                Repeat label (default: 1)
  --inventory-json PATH         Optional frozen benchmark inventory
  --playback-rate RATE          Sensor-time pacing; 0 disables pacing (default: 0)
  --max-lidar-frames COUNT      Stop after COUNT fused frames (default: 0, all)
  --cpu-set LIST                Linux CPU list (default: 0-7)
  --cpu-count COUNT             Allocated logical CPUs (default: 8)
  --completion-tolerance SEC    Allowed trajectory tail difference (default: 0.25)
  --max-output-gap SEC          Maximum adjacent LiDAR pose gap (default: 0.20)
  --watchdog-margin SEC         Extra watchdog wall time (default: 300)
  --experiment-fingerprint ID   Optional frozen experiment identifier
  --development-exposed BOOL    Metadata field (default: false)
  --confirmatory BOOL           Metadata field (default: true)
  --output-imu-tum PATH         Override IMU TUM path
  --output-lidar-tum PATH       Override primary-LiDAR TUM path
  --output-rear-axle-tum PATH   Override rear-axle TUM path
  --output-map PATH             Override PCD map path
  --output-frame-stats PATH     Override fused-frame CSV path
  --wait-ui BOOL                Wait for UI close (default: false)
  -h, --help                    Show this help

Paths supplied to output options may be absolute or relative to --output-dir.
Additional run_frontend_offline flags may be passed after `--`.
EOF
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "$script_dir/.." && pwd)"
ros_setup="${LIGHTNING_LM_ROS_SETUP:-/opt/ros/humble/setup.bash}"

bag_dir=""
config_path=""
run_dir=""
sequence=""
repeat="1"
inventory_json=""
playback_rate="0"
max_lidar_frames="0"
cpu_set="0-7"
allocated_cpus="8"
completion_tolerance="0.25"
max_output_gap="0.20"
watchdog_margin="300"
experiment_fingerprint="unfrozen_manual_run"
development_exposed="false"
confirmatory="true"
wait_ui="false"
imu_tum=""
lidar_tum=""
rear_tum=""
map_pcd=""
frame_stats=""
extra_args=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bag) bag_dir="${2:?missing value for --bag}"; shift 2 ;;
    --config) config_path="${2:?missing value for --config}"; shift 2 ;;
    --output-dir) run_dir="${2:?missing value for --output-dir}"; shift 2 ;;
    --sequence) sequence="${2:?missing value for --sequence}"; shift 2 ;;
    --repeat) repeat="${2:?missing value for --repeat}"; shift 2 ;;
    --inventory-json) inventory_json="${2:?missing value for --inventory-json}"; shift 2 ;;
    --playback-rate) playback_rate="${2:?missing value for --playback-rate}"; shift 2 ;;
    --max-lidar-frames) max_lidar_frames="${2:?missing value for --max-lidar-frames}"; shift 2 ;;
    --cpu-set) cpu_set="${2:?missing value for --cpu-set}"; shift 2 ;;
    --cpu-count) allocated_cpus="${2:?missing value for --cpu-count}"; shift 2 ;;
    --completion-tolerance) completion_tolerance="${2:?missing value for --completion-tolerance}"; shift 2 ;;
    --max-output-gap) max_output_gap="${2:?missing value for --max-output-gap}"; shift 2 ;;
    --watchdog-margin) watchdog_margin="${2:?missing value for --watchdog-margin}"; shift 2 ;;
    --experiment-fingerprint) experiment_fingerprint="${2:?missing value for --experiment-fingerprint}"; shift 2 ;;
    --development-exposed) development_exposed="${2:?missing value for --development-exposed}"; shift 2 ;;
    --confirmatory) confirmatory="${2:?missing value for --confirmatory}"; shift 2 ;;
    --output-imu-tum) imu_tum="${2:?missing value for --output-imu-tum}"; shift 2 ;;
    --output-lidar-tum) lidar_tum="${2:?missing value for --output-lidar-tum}"; shift 2 ;;
    --output-rear-axle-tum) rear_tum="${2:?missing value for --output-rear-axle-tum}"; shift 2 ;;
    --output-map) map_pcd="${2:?missing value for --output-map}"; shift 2 ;;
    --output-frame-stats) frame_stats="${2:?missing value for --output-frame-stats}"; shift 2 ;;
    --wait-ui) wait_ui="${2:?missing value for --wait-ui}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    --) shift; extra_args=("$@"); break ;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$bag_dir" || -z "$config_path" || -z "$run_dir" ]]; then
  usage >&2
  exit 2
fi
if [[ ! -d "$bag_dir" || ! -f "$bag_dir/metadata.yaml" ]]; then
  echo "SQLite3 ROS2 bag not found: $bag_dir" >&2
  exit 2
fi
if [[ ! -f "$config_path" ]]; then
  echo "configuration not found: $config_path" >&2
  exit 2
fi
if [[ -n "$inventory_json" && ! -f "$inventory_json" ]]; then
  echo "inventory not found: $inventory_json" >&2
  exit 2
fi

bag_dir="$(realpath "$bag_dir")"
config_path="$(realpath "$config_path")"
run_dir="$(realpath -m "$run_dir")"
if [[ -z "$sequence" ]]; then sequence="$(basename "$bag_dir")"; fi

for owned in results logs run_metadata.txt bag_contract.json watchdog_status.json resource_samples.csv resource_summary.json monitor.stop; do
  if [[ -e "$run_dir/$owned" ]]; then
    echo "refusing to overwrite existing run artifact: $run_dir/$owned" >&2
    exit 2
  fi
done
mkdir -p "$run_dir/results" "$run_dir/logs"

resolve_output() {
  local configured="$1"
  local fallback="$2"
  if [[ -z "$configured" ]]; then
    printf '%s\n' "$fallback"
  elif [[ "${configured:0:1}" == "/" ]]; then
    printf '%s\n' "$configured"
  else
    printf '%s\n' "$run_dir/$configured"
  fi
}

imu_tum="$(resolve_output "$imu_tum" "$run_dir/results/trajectory_imu.tum")"
lidar_tum="$(resolve_output "$lidar_tum" "$run_dir/results/trajectory_lidar114.tum")"
rear_tum="$(resolve_output "$rear_tum" "$run_dir/results/trajectory_rear_axle.tum")"
map_pcd="$(resolve_output "$map_pcd" "$run_dir/results/map_lio.pcd")"
frame_stats="$(resolve_output "$frame_stats" "$run_dir/results/frame_stats.csv")"
mkdir -p "$(dirname "$imu_tum")" "$(dirname "$lidar_tum")" "$(dirname "$rear_tum")" \
  "$(dirname "$map_pcd")" "$(dirname "$frame_stats")"

set +u
source "$ros_setup"
source "$repo_dir/install/setup.bash"
set -u
prefix="$(realpath "$(ros2 pkg prefix lightning)")"
expected_prefix="$(realpath "$repo_dir/install/lightning")"
if [[ "$prefix" != "$expected_prefix" ]]; then
  echo "lightning resolves to unexpected install prefix: $prefix (expected $expected_prefix)" >&2
  exit 2
fi
binary="$prefix/lib/lightning/run_frontend_offline"
monitor="$script_dir/monitor_process_tree.py"
inspector="$script_dir/inspect_rosbag2_sqlite.py"
timing_extractor="$script_dir/extract_frontend_timing.py"
if [[ ! -x "$binary" || ! -f "$monitor" || ! -f "$inspector" || ! -f "$timing_extractor" ]]; then
  echo "missing standard-build binary or runner helper; run colcon build first" >&2
  exit 2
fi

inspect_args=(--bag "$bag_dir" --config "$config_path" --output "$run_dir/bag_contract.json")
if [[ -n "$inventory_json" ]]; then
  inspect_args+=(--inventory "$inventory_json" --sequence "$sequence")
fi
python3 "$inspector" "${inspect_args[@]}"
readarray -t contract < <(python3 -c 'import json,sys
p=json.load(open(sys.argv[1],encoding="utf-8"))
print(p["primary_lidar_topic"])
print(p["expected_last_lidar_s"])
print(p["sensor_duration_s"])
print(p["expected_end_source"])
print(p["payload_file_count"])
print(p["payload_sha256"])' "$run_dir/bag_contract.json")
primary_lidar_topic="${contract[0]}"
expected_end="${contract[1]}"
sensor_duration="${contract[2]}"
expected_end_source="${contract[3]}"
bag_payload_file_count="${contract[4]}"
bag_payload_sha256="${contract[5]}"

watchdog_timeout="$(python3 -c 'import math,sys
duration,rate,margin=map(float,sys.argv[1:])
if not all(math.isfinite(v) for v in (duration,rate,margin)) or duration <= 0 or rate < 0 or margin < 0: raise SystemExit(2)
effective_rate=rate if rate > 0 else 1.0
print(max(1, math.ceil(duration/effective_rate+margin)))' "$sensor_duration" "$playback_rate" "$watchdog_margin")"

write_watchdog() {
  local status="$1"
  python3 -c 'import json,os,sys,time
path,status,timeout,margin,duration,rate=sys.argv[1:]
payload={"status":status,"timeout_s":float(timeout),"margin_s":float(margin),"sensor_duration_s":float(duration),"playback_rate":float(rate),"updated_wall_ns":time.time_ns()}
temporary=path+".tmp"
with open(temporary,"w",encoding="utf-8") as stream: json.dump(payload,stream,indent=2,sort_keys=True); stream.write("\n")
os.replace(temporary,path)' "$run_dir/watchdog_status.json" "$status" "$watchdog_timeout" "$watchdog_margin" "$sensor_duration" "$playback_rate"
}

algorithm_pid=""
monitor_pid=""
cleanup_group() {
  local pid="${1:-}"
  [[ -z "$pid" ]] && return 0
  kill -INT -- "-$pid" 2>/dev/null || true
  sleep 1
  kill -TERM -- "-$pid" 2>/dev/null || true
  sleep 1
  kill -KILL -- "-$pid" 2>/dev/null || true
}
cleanup() {
  touch "$run_dir/monitor.stop" 2>/dev/null || true
  cleanup_group "$algorithm_pid"
}
trap cleanup EXIT

export OMP_NUM_THREADS="$allocated_cpus"
export OPENBLAS_NUM_THREADS="$allocated_cpus"
export MKL_NUM_THREADS="$allocated_cpus"
start_ns="$(date +%s%N)"
write_watchdog "armed"
setsid taskset -c "$cpu_set" "$binary" \
  --input_bag="$bag_dir" \
  --config="$config_path" \
  --output_tum="$imu_tum" \
  --output_lidar_tum="$lidar_tum" \
  --output_rear_axle_tum="$rear_tum" \
  --output_map="$map_pcd" \
  --output_frame_stats_csv="$frame_stats" \
  --playback_rate="$playback_rate" \
  --max_lidar_frames="$max_lidar_frames" \
  --wait_ui="$wait_ui" \
  "${extra_args[@]}" \
  >"$run_dir/logs/algorithm.stdout.log" 2>"$run_dir/logs/algorithm.stderr.log" &
algorithm_pid=$!
write_watchdog "running"

rm -f "$run_dir/monitor.stop"
setsid taskset -c "$cpu_set" python3 "$monitor" \
  --pid "$algorithm_pid" \
  --stop-file "$run_dir/monitor.stop" \
  --csv "$run_dir/resource_samples.csv" \
  --summary "$run_dir/resource_summary.json" \
  --allocated-cpus "$allocated_cpus" \
  --interval 0.2 \
  >"$run_dir/logs/resource_monitor.log" 2>&1 &
monitor_pid=$!

watchdog_status="completed"
watchdog_deadline=$(( $(date +%s) + watchdog_timeout ))
while true; do
  process_state="$(ps -o stat= -p "$algorithm_pid" 2>/dev/null | awk '{print substr($1,1,1)}' || true)"
  if [[ -z "$process_state" || "$process_state" == "Z" ]]; then break; fi
  if (( $(date +%s) >= watchdog_deadline )); then
    watchdog_status="timed_out"
    write_watchdog "$watchdog_status"
    cleanup_group "$algorithm_pid"
    break
  fi
  sleep 1
done
set +e
wait "$algorithm_pid"
algorithm_rc=$?
set -e
algorithm_pid=""
if [[ "$watchdog_status" == "completed" ]]; then write_watchdog "$watchdog_status"; fi
end_ns="$(date +%s%N)"
touch "$run_dir/monitor.stop"
wait "$monitor_pid" 2>/dev/null || true
monitor_pid=""
trap - EXIT

trajectory_lines=0
last_stamp=0
invalid_count=0
nonmonotonic_count=0
maximum_output_gap_s=0
excessive_output_gap_count=0
if [[ -f "$lidar_tum" ]]; then
  readarray -t checks < <(python3 -c 'import math,sys
last=-math.inf; count=invalid=nonmono=0; gaps=[]
for line in open(sys.argv[1],encoding="utf-8"):
    try: values=[float(v) for v in line.split()]
    except ValueError: invalid+=1; continue
    if len(values)!=8 or not all(math.isfinite(v) for v in values): invalid+=1; continue
    if values[0] <= last: nonmono+=1; continue
    qnorm=math.sqrt(sum(v*v for v in values[4:8]))
    if abs(qnorm-1.0)>1e-3: invalid+=1; continue
    if math.isfinite(last): gaps.append(values[0]-last)
    last=values[0]; count+=1
limit=float(sys.argv[2])
print(count); print(last if math.isfinite(last) else 0); print(invalid); print(nonmono)
print(max(gaps) if gaps else 0.0); print(sum(gap>limit for gap in gaps))' "$lidar_tum" "$max_output_gap")
  trajectory_lines="${checks[0]}"
  last_stamp="${checks[1]}"
  invalid_count="${checks[2]}"
  nonmonotonic_count="${checks[3]}"
  maximum_output_gap_s="${checks[4]}"
  excessive_output_gap_count="${checks[5]}"
fi

completion="incomplete"
if [[ "$max_lidar_frames" -gt 0 ]]; then
  completion="limited_frame_run"
elif awk -v actual="$last_stamp" -v expected="$expected_end" -v tol="$completion_tolerance" 'BEGIN { exit !(actual >= expected-tol) }'; then
  completion="reached_final_lidar"
fi
wall_ns=$((end_ns-start_ns))
wall_time_s="$(awk -v value="$wall_ns" 'BEGIN { printf "%.6f", value/1000000000.0 }')"
timing_csv="$run_dir/results/processing_timing.csv"
timing_summary="$run_dir/results/processing_timing_summary.json"
python3 "$timing_extractor" \
  --log "$run_dir/logs/algorithm.stderr.log" \
  --csv "$timing_csv" \
  --summary "$timing_summary" \
  --wall-time-s "$wall_time_s" \
  --sensor-duration-s "$sensor_duration" \
  --trajectory-frames "$trajectory_lines" \
  --completion "$completion" \
  --max-lidar-frames "$max_lidar_frames" \
  --algorithm-rc "$algorithm_rc" \
  --watchdog-status "$watchdog_status" \
  --playback-rate "$playback_rate" \
  --wait-ui "$wait_ui"
readarray -t timing_checks < <(python3 -c 'import json,sys
p=json.load(open(sys.argv[1],encoding="utf-8"))
e=p["end_to_end"]
print(p["status"]); print(p["stage_count"])
for key in ("trajectory_frames_per_wall_s","wall_ms_per_trajectory_frame","realtime_factor","processing_speed_x"):
    value=e[key]
    print("null" if value is None else format(value,".9g"))' "$timing_summary")
timing_status="${timing_checks[0]}"
timing_stage_count="${timing_checks[1]}"
trajectory_frames_per_wall_s="${timing_checks[2]}"
wall_ms_per_trajectory_frame="${timing_checks[3]}"
realtime_factor="${timing_checks[4]}"
processing_speed_x="${timing_checks[5]}"

{
  echo "method=lightning_lm"
  echo "sequence=$sequence"
  echo "repeat=$repeat"
  echo "bag=$bag_dir"
  echo "bag_metadata_sha256=$(sha256sum "$bag_dir/metadata.yaml" | awk '{print $1}')"
  echo "bag_payload_file_count=$bag_payload_file_count"
  echo "bag_payload_sha256=$bag_payload_sha256"
  echo "primary_lidar_topic=$primary_lidar_topic"
  echo "expected_end_source=$expected_end_source"
  echo "cpu_set=$cpu_set"
  echo "allocated_cpus=$allocated_cpus"
  echo "play_rate=$playback_rate"
  echo "experiment_fingerprint=$experiment_fingerprint"
  echo "development_exposed=$development_exposed"
  echo "confirmatory=$confirmatory"
  echo "config=$config_path"
  echo "config_sha256=$(sha256sum "$config_path" | awk '{print $1}')"
  echo "algorithm_binary=$binary"
  echo "algorithm_binary_sha256=$(sha256sum "$binary" | awk '{print $1}')"
  echo "runner_sha256=$(sha256sum "$0" | awk '{print $1}')"
  echo "bag_contract_sha256=$(sha256sum "$run_dir/bag_contract.json" | awk '{print $1}')"
  if [[ -n "$inventory_json" ]]; then echo "inventory_sha256=$(sha256sum "$inventory_json" | awk '{print $1}')"; fi
  echo "output_imu_tum=$imu_tum"
  echo "output_lidar_tum=$lidar_tum"
  echo "output_rear_axle_tum=$rear_tum"
  echo "output_map=$map_pcd"
  echo "output_frame_stats=$frame_stats"
  echo "output_processing_timing_csv=$timing_csv"
  echo "output_processing_timing_summary=$timing_summary"
  echo "max_lidar_frames=$max_lidar_frames"
  echo "watchdog_status=$watchdog_status"
  echo "watchdog_timeout_s=$watchdog_timeout"
  echo "watchdog_margin_s=$watchdog_margin"
  echo "sensor_duration_s=$sensor_duration"
  echo "expected_last_lidar_end_s=$expected_end"
  echo "completion_tolerance_s=$completion_tolerance"
  echo "completion=$completion"
  echo "algorithm_rc=$algorithm_rc"
  echo "wall_time_s=$wall_time_s"
  echo "processing_timing_status=$timing_status"
  echo "processing_timing_stage_count=$timing_stage_count"
  echo "trajectory_frames_per_wall_s=$trajectory_frames_per_wall_s"
  echo "wall_ms_per_trajectory_frame=$wall_ms_per_trajectory_frame"
  echo "realtime_factor=$realtime_factor"
  echo "processing_speed_x=$processing_speed_x"
  echo "trajectory_lines=$trajectory_lines"
  echo "last_stamp=$last_stamp"
  echo "invalid_count=$invalid_count"
  echo "nonmonotonic_count=$nonmonotonic_count"
  echo "maximum_allowed_output_gap_s=$max_output_gap"
  echo "maximum_output_gap_s=$maximum_output_gap_s"
  echo "excessive_output_gap_count=$excessive_output_gap_count"
  echo "completed_at=$(date --iso-8601=seconds)"
} >"$run_dir/run_metadata.txt"

expected_completion="reached_final_lidar"
if [[ "$max_lidar_frames" -gt 0 ]]; then expected_completion="limited_frame_run"; fi
missing_output=0
for output in "$imu_tum" "$lidar_tum" "$rear_tum" "$map_pcd" "$frame_stats" "$timing_csv" "$timing_summary" \
              "$run_dir/resource_samples.csv" "$run_dir/resource_summary.json"; do
  if [[ ! -s "$output" ]]; then echo "missing or empty output: $output" >&2; missing_output=1; fi
done
if [[ "$algorithm_rc" -ne 0 || "$watchdog_status" != "completed" || "$completion" != "$expected_completion" || \
      "$trajectory_lines" -lt 10 || "$invalid_count" -ne 0 || "$nonmonotonic_count" -ne 0 || \
      "$excessive_output_gap_count" -ne 0 || "$timing_status" != "ok" || "$timing_stage_count" -lt 1 || \
      "$missing_output" -ne 0 ]]; then
  echo "run failed contract: rc=$algorithm_rc watchdog=$watchdog_status completion=$completion lines=$trajectory_lines gaps=$excessive_output_gap_count" >&2
  exit 4
fi
echo "completed method=lightning_lm sequence=$sequence repeat=$repeat lines=$trajectory_lines output=$run_dir"
