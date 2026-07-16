#!/usr/bin/env bash
# Canonical monitored runner for offline Lightning-LM localization.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_loc_offline.sh --bag ROS2_BAG --config YAML --map MAP_DIR --output-dir DIR [options]

Required:
  --bag PATH                    SQLite3 ROS2 bag directory
  --config PATH                 Lightning-LM YAML configuration
  --map PATH                    Tiled map directory containing index.txt
  --output-dir PATH             Exact run output directory

Options:
  --sequence NAME               Dataset sequence label (default: bag directory name)
  --repeat INDEX                Repeat label (default: 1)
  --reference-tum PATH          Optional reference trajectory for timestamp-aligned error analysis
  --playback-rate RATE          Sensor-time pacing; 0 disables pacing (default: 0)
  --max-lidar-frames COUNT      Stop after COUNT fused frames (default: 0, all)
  --cpu-set LIST                Linux CPU list (default: 0-7)
  --cpu-count COUNT             Allocated logical CPUs (default: 8)
  --completion-tolerance SEC    Allowed trajectory tail difference (default: 0.25)
  --max-output-gap SEC          Maximum adjacent localization pose gap (default: 0.30)
  --max-trajectory-speed MPS   Diagnostic speed threshold (default: 4.0)
  --max-trajectory-z-range M   Diagnostic Z-range threshold (default: 1.0)
  --max-trajectory-axis-range M  Diagnostic X/Y-range threshold (default: 40.0)
  --watchdog-margin SEC         Extra watchdog wall time (default: 300)
  --wait-ui BOOL                Wait for UI close (default: false)
  --publish-topics BOOL         Initialize ROS 2 publishers (default: false)
  --use-config-initial-pose BOOL  Use YAML initial pose (default: true)
  -h, --help                    Show this help

Arguments after `--` are forwarded unchanged to run_loc_offline.

Legacy environment variables are still accepted:
  LIGHTNING_LM_INPUT_BAG, LIGHTNING_LM_CONFIG, LIGHTNING_LM_MAP_PATH,
  LIGHTNING_LM_OUT_ROOT, LIGHTNING_LM_RUN_NAME, LIGHTNING_LM_REFERENCE_TUM.
EOF
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="${LIGHTNING_LM_REPO_DIR:-$(cd "$script_dir/.." && pwd)}"
ros_setup="${LIGHTNING_LM_ROS_SETUP:-/opt/ros/humble/setup.bash}"
install_setup="${LIGHTNING_LM_INSTALL_SETUP:-$repo_dir/install/setup.bash}"

bag_dir="${LIGHTNING_LM_INPUT_BAG:-}"
config_path="${LIGHTNING_LM_CONFIG:-$repo_dir/config/default.yaml}"
map_path="${LIGHTNING_LM_MAP_PATH:-$repo_dir/data/new_map}"
run_dir=""
out_root="${LIGHTNING_LM_OUT_ROOT:-$repo_dir/runs}"
run_name="${LIGHTNING_LM_RUN_NAME:-run_loc_offline_$(date +%Y%m%d_%H%M%S)}"
sequence=""
repeat="1"
reference_tum="${LIGHTNING_LM_REFERENCE_TUM:-}"
playback_rate="0"
max_lidar_frames="0"
cpu_set="0-7"
allocated_cpus="8"
completion_tolerance="0.25"
max_output_gap="0.30"
max_trajectory_speed="4.0"
max_trajectory_z_range="1.0"
max_trajectory_axis_range="40.0"
watchdog_margin="300"
wait_ui="false"
publish_topics="false"
use_config_initial_pose="true"
extra_args=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bag) bag_dir="${2:?missing value for --bag}"; shift 2 ;;
    --config) config_path="${2:?missing value for --config}"; shift 2 ;;
    --map|--map-path) map_path="${2:?missing value for --map}"; shift 2 ;;
    --output-dir) run_dir="${2:?missing value for --output-dir}"; shift 2 ;;
    --sequence) sequence="${2:?missing value for --sequence}"; shift 2 ;;
    --repeat) repeat="${2:?missing value for --repeat}"; shift 2 ;;
    --reference-tum) reference_tum="${2:?missing value for --reference-tum}"; shift 2 ;;
    --playback-rate) playback_rate="${2:?missing value for --playback-rate}"; shift 2 ;;
    --max-lidar-frames) max_lidar_frames="${2:?missing value for --max-lidar-frames}"; shift 2 ;;
    --cpu-set) cpu_set="${2:?missing value for --cpu-set}"; shift 2 ;;
    --cpu-count) allocated_cpus="${2:?missing value for --cpu-count}"; shift 2 ;;
    --completion-tolerance) completion_tolerance="${2:?missing value for --completion-tolerance}"; shift 2 ;;
    --max-output-gap) max_output_gap="${2:?missing value for --max-output-gap}"; shift 2 ;;
    --max-trajectory-speed) max_trajectory_speed="${2:?missing value for --max-trajectory-speed}"; shift 2 ;;
    --max-trajectory-z-range) max_trajectory_z_range="${2:?missing value for --max-trajectory-z-range}"; shift 2 ;;
    --max-trajectory-axis-range) max_trajectory_axis_range="${2:?missing value for --max-trajectory-axis-range}"; shift 2 ;;
    --watchdog-margin) watchdog_margin="${2:?missing value for --watchdog-margin}"; shift 2 ;;
    --wait-ui) wait_ui="${2:?missing value for --wait-ui}"; shift 2 ;;
    --publish-topics) publish_topics="${2:?missing value for --publish-topics}"; shift 2 ;;
    --use-config-initial-pose) use_config_initial_pose="${2:?missing value for --use-config-initial-pose}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    --) shift; extra_args=("$@"); break ;;
    *)
      if [[ "${1:0:1}" != "-" && -z "${LIGHTNING_LM_RUN_NAME:-}" && -z "$run_dir" ]]; then
        run_name="$1"
        shift
      else
        echo "unknown argument: $1" >&2
        usage >&2
        exit 2
      fi
      ;;
  esac
done

if [[ -z "$bag_dir" || -z "$config_path" || -z "$map_path" ]]; then
  usage >&2
  exit 2
fi
if [[ -z "$run_dir" ]]; then
  out_root="$(mkdir -p "$out_root" && cd "$out_root" && pwd)"
  run_dir="$out_root/$run_name"
fi
if [[ ! -d "$bag_dir" || ! -f "$bag_dir/metadata.yaml" ]]; then
  echo "SQLite3 ROS2 bag not found: $bag_dir" >&2
  exit 2
fi
if [[ ! -f "$config_path" ]]; then
  echo "configuration not found: $config_path" >&2
  exit 2
fi
if [[ ! -f "$map_path/index.txt" ]]; then
  echo "tiled map index not found: $map_path/index.txt" >&2
  exit 2
fi
if [[ -n "$reference_tum" && ! -f "$reference_tum" ]]; then
  echo "reference trajectory not found: $reference_tum" >&2
  exit 2
fi

bag_dir="$(realpath "$bag_dir")"
config_path="$(realpath "$config_path")"
map_path="$(realpath "$map_path")"
run_dir="$(realpath -m "$run_dir")"
if [[ -n "$reference_tum" ]]; then reference_tum="$(realpath "$reference_tum")"; fi
if [[ -z "$sequence" ]]; then sequence="$(basename "$bag_dir")"; fi

for owned in results logs run_metadata.txt bag_contract.json watchdog_status.json resource_samples.csv resource_summary.json monitor.stop; do
  if [[ -e "$run_dir/$owned" ]]; then
    echo "refusing to overwrite existing run artifact: $run_dir/$owned" >&2
    exit 2
  fi
done
mkdir -p "$run_dir/results" "$run_dir/logs"

set +u
source "$ros_setup"
source "$install_setup"
set -u
prefix="$(realpath "$(ros2 pkg prefix lightning)")"
expected_prefix="$(realpath "$repo_dir/install/lightning")"
if [[ "$prefix" != "$expected_prefix" ]]; then
  echo "lightning resolves to unexpected install prefix: $prefix (expected $expected_prefix)" >&2
  exit 2
fi

binary="$prefix/lib/lightning/run_loc_offline"
monitor="$script_dir/monitor_process_tree.py"
inspector="$script_dir/inspect_rosbag2_sqlite.py"
timing_extractor="$script_dir/extract_frontend_timing.py"
analyzer="$script_dir/analyze_localization_run.py"
if [[ ! -x "$binary" || ! -f "$monitor" || ! -f "$inspector" || ! -f "$timing_extractor" || ! -f "$analyzer" ]]; then
  echo "missing standard-build binary or runner helper; run colcon build first" >&2
  exit 2
fi

python3 "$inspector" --bag "$bag_dir" --config "$config_path" --output "$run_dir/bag_contract.json"
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
with open(temporary,"w",encoding="utf-8") as stream:
    json.dump(payload,stream,indent=2,sort_keys=True)
    stream.write("\n")
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

trajectory_tum="$run_dir/results/trajectory_loc.tum"
lidar_loc_tum="$run_dir/results/trajectory_lidar_loc.tum"
localization_csv="$run_dir/results/localization_stats.csv"
frame_stats="$run_dir/results/frame_stats.csv"
timing_csv="$run_dir/results/processing_timing.csv"
timing_summary="$run_dir/results/processing_timing_summary.json"
analysis_summary="$run_dir/results/localization_summary.json"
error_csv="$run_dir/results/trajectory_reference_errors.csv"

export OMP_NUM_THREADS="$allocated_cpus"
export OPENBLAS_NUM_THREADS="$allocated_cpus"
export MKL_NUM_THREADS="$allocated_cpus"
start_ns="$(date +%s%N)"
write_watchdog "armed"
setsid taskset -c "$cpu_set" "$binary" \
  --input_bag="$bag_dir" \
  --config="$config_path" \
  --map_path="$map_path" \
  --output_tum="$trajectory_tum" \
  --output_lidar_loc_tum="$lidar_loc_tum" \
  --output_csv="$localization_csv" \
  --output_frame_stats_csv="$frame_stats" \
  --playback_rate="$playback_rate" \
  --max_lidar_frames="$max_lidar_frames" \
  --wait_ui="$wait_ui" \
  --publish_topics="$publish_topics" \
  --use_config_initial_pose="$use_config_initial_pose" \
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
if [[ -f "$trajectory_tum" ]]; then
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
print(max(gaps) if gaps else 0.0); print(sum(gap>limit for gap in gaps))' "$trajectory_tum" "$max_output_gap")
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

analyzer_args=(--localization-csv "$localization_csv" --trajectory-tum "$trajectory_tum" \
  --lidar-loc-tum "$lidar_loc_tum" --summary-json "$analysis_summary" --errors-csv "$error_csv" \
  --max-speed "$max_trajectory_speed" --max-z-range "$max_trajectory_z_range" \
  --max-axis-range "$max_trajectory_axis_range")
if [[ -n "$reference_tum" ]]; then analyzer_args+=(--reference-tum "$reference_tum"); fi
if [[ -s "$localization_csv" && -s "$trajectory_tum" ]]; then
  python3 "$analyzer" "${analyzer_args[@]}"
fi

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
  echo "method=lightning_lm_offline_localization"
  echo "sequence=$sequence"
  echo "repeat=$repeat"
  echo "bag=$bag_dir"
  echo "bag_metadata_sha256=$(sha256sum "$bag_dir/metadata.yaml" | awk '{print $1}')"
  echo "bag_payload_file_count=$bag_payload_file_count"
  echo "bag_payload_sha256=$bag_payload_sha256"
  echo "primary_lidar_topic=$primary_lidar_topic"
  echo "expected_end_source=$expected_end_source"
  echo "map_path=$map_path"
  echo "map_index_sha256=$(sha256sum "$map_path/index.txt" | awk '{print $1}')"
  echo "cpu_set=$cpu_set"
  echo "allocated_cpus=$allocated_cpus"
  echo "play_rate=$playback_rate"
  echo "publish_topics=$publish_topics"
  echo "use_config_initial_pose=$use_config_initial_pose"
  echo "config=$config_path"
  echo "config_sha256=$(sha256sum "$config_path" | awk '{print $1}')"
  echo "algorithm_binary=$binary"
  echo "algorithm_binary_sha256=$(sha256sum "$binary" | awk '{print $1}')"
  echo "runner_sha256=$(sha256sum "$0" | awk '{print $1}')"
  echo "bag_contract_sha256=$(sha256sum "$run_dir/bag_contract.json" | awk '{print $1}')"
  if [[ -n "$reference_tum" ]]; then echo "reference_tum=$reference_tum"; fi
  echo "output_tum=$trajectory_tum"
  echo "output_lidar_loc_tum=$lidar_loc_tum"
  echo "output_localization_csv=$localization_csv"
  echo "output_frame_stats=$frame_stats"
  echo "output_processing_timing_csv=$timing_csv"
  echo "output_processing_timing_summary=$timing_summary"
  echo "output_localization_summary=$analysis_summary"
  echo "output_reference_errors=$error_csv"
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
  echo "max_trajectory_speed_mps=$max_trajectory_speed"
  echo "max_trajectory_z_range_m=$max_trajectory_z_range"
  echo "max_trajectory_axis_range_m=$max_trajectory_axis_range"
  echo "completed_at=$(date --iso-8601=seconds)"
} >"$run_dir/run_metadata.txt"

expected_completion="reached_final_lidar"
if [[ "$max_lidar_frames" -gt 0 ]]; then expected_completion="limited_frame_run"; fi
missing_output=0
for output in "$trajectory_tum" "$lidar_loc_tum" "$localization_csv" "$frame_stats" "$timing_csv" "$timing_summary" \
              "$run_dir/resource_samples.csv" "$run_dir/resource_summary.json"; do
  if [[ ! -s "$output" ]]; then echo "missing or empty output: $output" >&2; missing_output=1; fi
done
physical_diagnostic_pass=0
if [[ -s "$analysis_summary" ]]; then
  physical_diagnostic_pass="$(python3 -c 'import json,sys
p=json.load(open(sys.argv[1],encoding="utf-8"))
print(int(bool(p["trajectory"]["motion"]["passed"]) and bool(p["lidar_loc_trajectory"]["motion"]["passed"])))' "$analysis_summary")"
fi
if [[ "$algorithm_rc" -ne 0 || "$watchdog_status" != "completed" || "$completion" != "$expected_completion" || \
      "$trajectory_lines" -lt 10 || "$invalid_count" -ne 0 || "$nonmonotonic_count" -ne 0 || \
      "$excessive_output_gap_count" -ne 0 || "$timing_status" != "ok" || "$timing_stage_count" -lt 1 || \
      "$missing_output" -ne 0 ]]; then
  echo "run failed contract: rc=$algorithm_rc watchdog=$watchdog_status completion=$completion lines=$trajectory_lines gaps=$excessive_output_gap_count physical_diagnostic=$physical_diagnostic_pass" >&2
  exit 4
fi
echo "completed method=lightning_lm_offline_localization sequence=$sequence repeat=$repeat lines=$trajectory_lines physical_diagnostic=$physical_diagnostic_pass output=$run_dir"
