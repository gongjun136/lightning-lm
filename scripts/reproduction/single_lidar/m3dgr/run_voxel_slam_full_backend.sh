#!/usr/bin/env bash
# Reproducible ROS1 full Voxel-SLAM baseline on normalized M3DGR bags.
set -euo pipefail

bag="${1:?usage: run_voxel_slam_full_backend.sh BAG SEQUENCE OUTPUT_DIR}"
sequence="${2:?missing sequence}"
run_dir="${3:?missing output directory}"
repeat="${4:-1}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
lightning_repo="$(cd "$script_dir/../../../.." && pwd)"
voxel_ws="${VOXEL_SLAM_WS:-/mnt/f/SLAM_AI_KnowledgeBase/code/WSL_Ubuntu_20.04/ros1_ws/ws_voxel_slam}"
bundle="$voxel_ws/src/Voxel-SLAM/reproduction/m3dgr"
config="$bundle/config/voxel_slam_mid360_frontend.yaml"
binary="$voxel_ws/devel/lib/voxel_slam/voxelslam"
monitor="$bundle/scripts/monitor_process_tree.py"
extractor="$script_dir/extract_voxel_slam_loop_events.py"
cpu_set="${BENCH_CPUSET:-0-7}"
cpu_count="${BENCH_CPU_COUNT:-8}"
play_rate="${BENCH_PLAY_RATE:-1.0}"
loop_icp_eigval="${VOXEL_LOOP_ICP_EIGVAL:-}"
ros_port="${BENCH_ROS_PORT:-11331}"
finish_timeout="${BENCH_FINISH_TIMEOUT_S:-1800}"
inventory="${BENCH_INVENTORY_JSON:-/mnt/f/SLAM_AI_KnowledgeBase/code/_m3dgr_work/bench/inventory/bag_inventory.json}"
completion_tolerance="${BENCH_COMPLETION_TOLERANCE_S:-0.25}"
maximum_allowed_output_gap="${BENCH_MAX_OUTPUT_GAP_S:-0.20}"

for required in "$bag" "$config" "$binary" "$monitor" "$extractor" "$inventory"; do
  [[ -e "$required" ]] || { echo "missing dependency: $required" >&2; exit 2; }
done
if [[ -d "$run_dir" ]] && [[ -n "$(find "$run_dir" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
  echo "refusing to overwrite non-empty output: $run_dir" >&2
  exit 2
fi

mkdir -p "$run_dir/logs" "$run_dir/results" \
  "$run_dir/data/new_map/backend_diagnostics" "$run_dir/voxel_output"
run_dir="$(realpath "$run_dir")"
trajectory=""
output_tum="$run_dir/results/trajectory_slam_opt.tum"

export ROS_MASTER_URI="http://127.0.0.1:$ros_port"
export ROS_IP=127.0.0.1
export ROS_HOME="$run_dir/ros_home"
export ROS_LOG_DIR="$run_dir/logs/ros"
export OMP_NUM_THREADS="$cpu_count"
export OPENBLAS_NUM_THREADS="$cpu_count"
export MKL_NUM_THREADS="$cpu_count"
mkdir -p "$ROS_HOME" "$ROS_LOG_DIR"

set +u
source /opt/ros/noetic/setup.bash
source "$voxel_ws/devel/setup.bash"
set -u

core_pid=""
algorithm_pid=""
monitor_pid=""
play_pid=""
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
  cleanup_group "$play_pid"
  cleanup_group "$algorithm_pid"
  cleanup_group "$core_pid"
}
trap cleanup EXIT

setsid taskset -c "$cpu_set" roscore -p "$ros_port" >"$run_dir/logs/roscore.log" 2>&1 &
core_pid=$!
master_deadline=$(( $(date +%s) + 60 ))
until rosparam list >/dev/null 2>&1; do
  kill -0 "$core_pid" 2>/dev/null || { echo "roscore exited during startup" >&2; exit 3; }
  (( $(date +%s) < master_deadline )) || { echo "roscore startup timeout" >&2; exit 3; }
  sleep 1
done
rosparam set /use_sim_time true
rosparam load "$config"
rosparam set /General/enable_backend true
rosparam set /General/is_save_map 1
rosparam set /General/save_path "$run_dir/voxel_output/"
rosparam set /General/bagname "$sequence"
if [[ -n "$loop_icp_eigval" ]]; then
  rosparam set /Loop/icp_eigval "$loop_icp_eigval"
fi
rosparam set /finish false

setsid taskset -c "$cpu_set" "$binary" >"$run_dir/logs/algorithm.log" 2>&1 &
algorithm_pid=$!
subscriber_deadline=$(( $(date +%s) + 90 ))
until python3 -c 'import rosgraph
s=dict(rosgraph.Master("/voxel_baseline_probe").getSystemState()[1])
raise SystemExit(0 if s.get("/livox/mid360/lidar") and s.get("/livox/mid360/imu") else 1)' \
    >/dev/null 2>&1; do
  kill -0 "$algorithm_pid" 2>/dev/null || { echo "Voxel-SLAM exited during startup" >&2; exit 3; }
  (( $(date +%s) < subscriber_deadline )) || { echo "sensor subscriber timeout" >&2; exit 3; }
  sleep 1
done

rm -f "$run_dir/monitor.stop"
setsid taskset -c "$cpu_set" python3 "$monitor" \
  --pid "$algorithm_pid" --stop-file "$run_dir/monitor.stop" \
  --csv "$run_dir/resource_samples.csv" --summary "$run_dir/resource_summary.json" \
  --allocated-cpus "$cpu_count" --interval 0.2 \
  >"$run_dir/logs/resource_monitor.log" 2>&1 &
monitor_pid=$!

start_ns="$(date +%s%N)"
setsid taskset -c "$cpu_set" rosbag play "$bag" --clock --quiet --wait-for-subscribers \
  --rate "$play_rate" --topics /livox/mid360/lidar /livox/mid360/imu \
  >"$run_dir/logs/rosbag_play.log" 2>&1 &
play_pid=$!
wait "$play_pid"
play_pid=""
rosparam set /finish true

deadline=$(( $(date +%s) + finish_timeout ))
while [[ -z "$trajectory" ]]; do
  trajectory="$(find "$run_dir/voxel_output" -mindepth 2 -maxdepth 2 \
    -name alidarState.txt -type f -size +0c -printf '%T@ %p\n' 2>/dev/null \
    | sort -nr | head -n1 | cut -d' ' -f2-)"
  [[ -n "$trajectory" ]] && break
  kill -0 "$algorithm_pid" 2>/dev/null || { echo "Voxel-SLAM exited before saving trajectory" >&2; exit 4; }
  (( $(date +%s) < deadline )) || { echo "Voxel-SLAM backend finish timeout" >&2; exit 4; }
  sleep 2
done
previous_size=-1
stable_count=0
while (( stable_count < 3 )); do
  current_size="$(stat -c %s "$trajectory")"
  if [[ "$current_size" == "$previous_size" ]]; then
    stable_count=$((stable_count + 1))
  else
    stable_count=0
    previous_size="$current_size"
  fi
  sleep 2
done
end_ns="$(date +%s%N)"

awk 'NF >= 8 {print $1,$2,$3,$4,$5,$6,$7,$8}' "$trajectory" >"$output_tum"
python3 "$extractor" --log "$run_dir/logs/algorithm.log" --trajectory "$output_tum" \
  --output "$run_dir/data/new_map/backend_diagnostics/btc_loop_candidates.csv"

readarray -t inventory_values < <(python3 -c 'import json,sys
rows=json.load(open(sys.argv[1],encoding="utf-8"))
row=next((item for item in rows if item.get("sequence")==sys.argv[2]),None)
if row is None: raise SystemExit(f"missing inventory row: {sys.argv[2]}")
print(row["sensor_duration_s"])
print(row["lidar_last_end_ns"]*1e-9)
print(row["lidar_count"])' "$inventory" "$sequence")
sensor_duration="${inventory_values[0]}"
expected_last_lidar_end="${inventory_values[1]}"
expected_lidar_frames="${inventory_values[2]}"

readarray -t trajectory_checks < <(python3 -c 'import math,sys
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
print(count); print(last if math.isfinite(last) else 0.0); print(invalid); print(nonmono)
print(max(gaps) if gaps else 0.0); print(sum(gap>limit for gap in gaps))' "$output_tum" "$maximum_allowed_output_gap")
trajectory_lines="${trajectory_checks[0]}"
last_stamp="${trajectory_checks[1]}"
invalid_count="${trajectory_checks[2]}"
nonmonotonic_count="${trajectory_checks[3]}"
maximum_output_gap_s="${trajectory_checks[4]}"
excessive_output_gap_count="${trajectory_checks[5]}"
completion="incomplete"
if awk -v actual="$last_stamp" -v expected="$expected_last_lidar_end" -v tol="$completion_tolerance" \
  'BEGIN { exit !(actual >= expected-tol) }'; then
  completion="reached_final_lidar"
fi
pcd_count="$(find "$(dirname "$trajectory")" -maxdepth 1 -type f -name '*.pcd' | wc -l)"

touch "$run_dir/monitor.stop"
wait "$monitor_pid" 2>/dev/null || true
monitor_pid=""
cleanup_group "$algorithm_pid"
wait "$algorithm_pid" 2>/dev/null || true
algorithm_pid=""
cleanup_group "$core_pid"
wait "$core_pid" 2>/dev/null || true
core_pid=""
trap - EXIT

wall_s="$(awk -v ns="$((end_ns-start_ns))" 'BEGIN {printf "%.6f", ns/1e9}')"
cat >"$run_dir/run_metadata.txt" <<EOF
method=voxel_slam_full_backend
sequence=$sequence
repeat=$repeat
bag=$(realpath "$bag")
bag_size_bytes=$(stat -c %s "$bag")
bag_sha256=$(sha256sum "$bag" | awk '{print $1}')
config=$(realpath "$config")
config_sha256=$(sha256sum "$config" | awk '{print $1}')
algorithm_binary=$(realpath "$binary")
algorithm_binary_sha256=$(sha256sum "$binary" | awk '{print $1}')
runner_sha256=$(sha256sum "$0" | awk '{print $1}')
inventory=$(realpath "$inventory")
inventory_sha256=$(sha256sum "$inventory" | awk '{print $1}')
cpu_set=$cpu_set
allocated_cpus=$cpu_count
play_rate=$play_rate
loop_icp_eigval=${loop_icp_eigval:-config_default}
wall_time_s=$wall_s
output_tum=$output_tum
sensor_duration_s=$sensor_duration
expected_lidar_frames=$expected_lidar_frames
expected_last_lidar_end_s=$expected_last_lidar_end
completion_tolerance_s=$completion_tolerance
completion=$completion
trajectory_lines=$trajectory_lines
last_stamp=$last_stamp
invalid_count=$invalid_count
nonmonotonic_count=$nonmonotonic_count
maximum_allowed_output_gap_s=$maximum_allowed_output_gap
maximum_output_gap_s=$maximum_output_gap_s
excessive_output_gap_count=$excessive_output_gap_count
pcd_count=$pcd_count
completed_at=$(date --iso-8601=seconds)
EOF
if [[ "$completion" != "reached_final_lidar" || "$trajectory_lines" -lt 10 || "$invalid_count" -ne 0 || \
      "$nonmonotonic_count" -ne 0 || "$excessive_output_gap_count" -ne 0 || "$pcd_count" -ne "$trajectory_lines" || \
      ! -s "$run_dir/resource_summary.json" ]]; then
  echo "Voxel-SLAM full backend failed contract: completion=$completion lines=$trajectory_lines pcds=$pcd_count invalid=$invalid_count nonmono=$nonmonotonic_count gaps=$excessive_output_gap_count" >&2
  exit 4
fi
echo "completed method=voxel_slam_full_backend sequence=$sequence output=$run_dir"
