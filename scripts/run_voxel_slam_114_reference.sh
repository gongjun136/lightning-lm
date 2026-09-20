#!/usr/bin/env bash
set -eo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_voxel_slam_114_reference.sh --bag ROS1_BAG --output-dir DIR --sequence NAME [options]

Required:
  --bag PATH             ROS1 bag containing Livox 114 lidar and IMU topics
  --output-dir PATH      Exact output directory; must not already contain files
  --sequence NAME        Output sequence name used by Voxel-SLAM

Options:
  --voxel-ws PATH        Read-only ws_voxel_slam workspace
  --lidar-topic TOPIC    ROS1 PointCloud2 topic (default: legacy 114 topic)
  --imu-topic TOPIC      ROS1 IMU topic (default: legacy 114 topic)
  --launch-file PATH     Voxel-SLAM launch file relative to its package
  --config-file PATH     Voxel-SLAM config file used for provenance
  --ros-port PORT        Private ROS master port (default: 11331)
  --shutdown-wait SEC    Maximum wait for final optimization (default: 900)
  --cpu-set LIST         Linux CPU list reserved for the run (default: 0-7)
  --cpu-count COUNT      Number of reserved logical CPUs (default: 8)
  --play-rate RATE       rosbag playback rate (default: 1.0)
  --repeat INDEX         Formal repeat identifier (default: 1)
  -h, --help             Show this help
EOF
}

bag=""
output_dir=""
sequence=""
voxel_ws="/mnt/f/SLAM_AI_KnowledgeBase/code/WSL_Ubuntu_20.04/ros1_ws/ws_voxel_slam"
lidar_topic="/livox/lidar_192_168_1_114"
imu_topic="/livox/imu_192_168_1_114"
launch_rel="launch/vxlm_sany_20260701_livox_pc2_114.launch"
config_rel="config/sany_20260701_livox_pc2_114.yaml"
ros_port=11331
shutdown_wait=900
cpu_set="0-7"
cpu_count=8
play_rate=1.0
repeat=1

while (($#)); do
  case "$1" in
    --bag) bag="${2:?missing value for --bag}"; shift 2 ;;
    --output-dir) output_dir="${2:?missing value for --output-dir}"; shift 2 ;;
    --sequence) sequence="${2:?missing value for --sequence}"; shift 2 ;;
    --voxel-ws) voxel_ws="${2:?missing value for --voxel-ws}"; shift 2 ;;
    --lidar-topic) lidar_topic="${2:?missing value for --lidar-topic}"; shift 2 ;;
    --imu-topic) imu_topic="${2:?missing value for --imu-topic}"; shift 2 ;;
    --launch-file) launch_rel="${2:?missing value for --launch-file}"; shift 2 ;;
    --config-file) config_rel="${2:?missing value for --config-file}"; shift 2 ;;
    --ros-port) ros_port="${2:?missing value for --ros-port}"; shift 2 ;;
    --shutdown-wait) shutdown_wait="${2:?missing value for --shutdown-wait}"; shift 2 ;;
    --cpu-set) cpu_set="${2:?missing value for --cpu-set}"; shift 2 ;;
    --cpu-count) cpu_count="${2:?missing value for --cpu-count}"; shift 2 ;;
    --play-rate) play_rate="${2:?missing value for --play-rate}"; shift 2 ;;
    --repeat) repeat="${2:?missing value for --repeat}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$bag" && -n "$output_dir" && -n "$sequence" ]] || {
  usage >&2
  exit 2
}
[[ -f "$bag" ]] || { echo "bag does not exist: $bag" >&2; exit 2; }
[[ "$sequence" =~ ^[A-Za-z0-9_.-]+$ ]] || {
  echo "sequence must contain only letters, digits, dot, underscore, or dash" >&2
  exit 2
}
[[ "$ros_port" =~ ^[0-9]+$ ]] && ((ros_port >= 1024 && ros_port <= 65535)) || {
  echo "invalid ROS port: $ros_port" >&2
  exit 2
}
[[ "$shutdown_wait" =~ ^[0-9]+$ ]] && ((shutdown_wait > 0)) || {
  echo "invalid shutdown wait: $shutdown_wait" >&2
  exit 2
}
[[ "$cpu_count" =~ ^[0-9]+$ ]] && ((cpu_count > 0)) || {
  echo "invalid CPU count: $cpu_count" >&2
  exit 2
}
[[ "$repeat" =~ ^[0-9]+$ ]] && ((repeat > 0)) || {
  echo "invalid repeat index: $repeat" >&2
  exit 2
}
awk -v value="$play_rate" 'BEGIN { exit !(value > 0) }' || {
  echo "invalid playback rate: $play_rate" >&2
  exit 2
}

setup="$voxel_ws/devel/setup.bash"
recorder="$voxel_ws/src/Voxel-SLAM/reproduction/m3dgr/scripts/trajectory_recorder.py"
binary="$voxel_ws/devel/lib/voxel_slam/voxelslam"
launch_file="$voxel_ws/src/Voxel-SLAM/VoxelSLAM/$launch_rel"
config_file="$voxel_ws/src/Voxel-SLAM/VoxelSLAM/$config_rel"
monitor="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/monitor_process_tree.py"
[[ -f "$setup" ]] || { echo "Voxel-SLAM setup does not exist: $setup" >&2; exit 2; }
[[ -f "$recorder" ]] || { echo "trajectory recorder does not exist: $recorder" >&2; exit 2; }
[[ -x "$binary" ]] || { echo "Voxel-SLAM binary does not exist: $binary" >&2; exit 2; }
[[ -f "$launch_file" ]] || { echo "Voxel-SLAM launch file does not exist: $launch_file" >&2; exit 2; }
[[ -f "$config_file" ]] || { echo "Voxel-SLAM config does not exist: $config_file" >&2; exit 2; }
[[ -f "$monitor" ]] || { echo "resource monitor does not exist: $monitor" >&2; exit 2; }

if [[ -d "$output_dir" ]] && find "$output_dir" -mindepth 1 -print -quit | grep -q .; then
  echo "output directory is not empty: $output_dir" >&2
  exit 2
fi
mkdir -p "$output_dir/logs" "$output_dir/results" "$output_dir/data"
output_dir="$(realpath "$output_dir")"
bag="$(realpath "$bag")"

source /opt/ros/noetic/setup.bash
rosbag info --yaml "$bag" >"$output_dir/logs/input_bag_info.yaml"
readarray -t bag_contract < <(python3 -c 'import sys,rosbag,yaml
info=yaml.safe_load(open(sys.argv[1],encoding="utf-8"))
lidar_topic=sys.argv[3]
topic=next((row for row in info.get("topics",[]) if row.get("topic")==lidar_topic),None)
if topic is None: raise SystemExit(f"missing lidar topic: {lidar_topic}")
with rosbag.Bag(sys.argv[2],"r") as bag:
    times=[msg.header.stamp.to_sec() for _,msg,_ in bag.read_messages(topics=[lidar_topic])]
if not times: raise SystemExit("empty Livox 114 lidar index")
print(topic["messages"])
print(info["duration"])
print(info["start"])
print(info["end"])
print(min(times))
print(max(times))' "$output_dir/logs/input_bag_info.yaml" "$bag" "$lidar_topic")
expected_lidar_frames="${bag_contract[0]}"
sensor_duration_s="${bag_contract[1]}"
bag_start_s="${bag_contract[2]}"
bag_end_s="${bag_contract[3]}"
lidar_first_s="${bag_contract[4]}"
lidar_last_s="${bag_contract[5]}"

# shellcheck disable=SC1091
source "$setup"
set -u
export ROS_MASTER_URI="http://127.0.0.1:${ros_port}"
export ROS_HOSTNAME=127.0.0.1
export ROS_HOME="$output_dir/ros_home"
export ROS_LOG_DIR="$output_dir/logs/ros"
export OMP_NUM_THREADS="$cpu_count"
export OPENBLAS_NUM_THREADS="$cpu_count"
export MKL_NUM_THREADS="$cpu_count"
mkdir -p "$ROS_HOME" "$ROS_LOG_DIR"

core_pid=""
launch_pid=""
recorder_pid=""
monitor_pid=""
play_pid=""
cleanup() {
  set +e
  touch "$output_dir/monitor.stop" 2>/dev/null || true
  for pid in "$play_pid" "$recorder_pid" "$launch_pid" "$core_pid"; do
    if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then
      kill -INT -- "-$pid" 2>/dev/null || kill -INT "$pid" 2>/dev/null || true
    fi
  done
  sleep 1
  for pid in "$play_pid" "$recorder_pid" "$launch_pid" "$core_pid"; do
    if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then
      kill -TERM -- "-$pid" 2>/dev/null || kill -TERM "$pid" 2>/dev/null || true
    fi
  done
  wait "$play_pid" "$recorder_pid" "$launch_pid" "$monitor_pid" "$core_pid" 2>/dev/null || true
}
trap cleanup EXIT INT TERM

setsid taskset -c "$cpu_set" roscore -p "$ros_port" >"$output_dir/logs/roscore.log" 2>&1 &
core_pid=$!
for _ in $(seq 1 100); do
  rosparam list >/dev/null 2>&1 && break
  sleep 0.1
done
rosparam list >/dev/null 2>&1 || { echo "ROS master did not start" >&2; exit 1; }
rosparam set /use_sim_time true

setsid taskset -c "$cpu_set" roslaunch "$launch_file" \
  rviz:=false save_path:="$output_dir/data/" bagname:="$sequence" \
  >"$output_dir/logs/voxel_slam.log" 2>&1 &
launch_pid=$!
for _ in $(seq 1 300); do
  rosnode list 2>/dev/null | grep -qx /voxelslam && break
  kill -0 "$launch_pid" 2>/dev/null || {
    echo "Voxel-SLAM exited during startup; see $output_dir/logs/voxel_slam.log" >&2
    exit 1
  }
  sleep 0.1
done
rosnode list 2>/dev/null | grep -qx /voxelslam || {
  echo "Voxel-SLAM node did not become ready" >&2
  exit 1
}

setsid taskset -c "$cpu_set" python3 "$recorder" \
  --mode tf --topic /tf --parent-frame camera_init --child-frame aft_mapped \
  --tx -0.011 --ty -0.02329 --tz 0.04412 \
  --output "$output_dir/results/trajectory_voxel_frontend.tum" \
  --status "$output_dir/results/trajectory_voxel_frontend.status.json" \
  >"$output_dir/logs/trajectory_recorder.log" 2>&1 &
recorder_pid=$!

rm -f "$output_dir/monitor.stop"
setsid taskset -c "$cpu_set" python3 "$monitor" \
  --pid "$launch_pid" --stop-file "$output_dir/monitor.stop" \
  --csv "$output_dir/resource_samples.csv" --summary "$output_dir/resource_summary.json" \
  --allocated-cpus "$cpu_count" --interval 0.2 \
  >"$output_dir/logs/resource_monitor.log" 2>&1 &
monitor_pid=$!

start_ns="$(date +%s%N)"
setsid taskset -c "$cpu_set" rosbag play "$bag" --clock --quiet --wait-for-subscribers \
  --rate "$play_rate" \
  --topics "$lidar_topic" "$imu_topic" \
  >"$output_dir/logs/rosbag_play.log" 2>&1 &
play_pid=$!
wait "$play_pid"
play_pid=""

# Let subscriber queues drain before requesting the final global optimization.
sleep 5
rosparam set /finish true

optimized_state="$output_dir/data/$sequence/alidarState.txt"
deadline=$((SECONDS + shutdown_wait))
stable_size=-1
stable_checks=0
while ((stable_checks < 5)); do
  if ((SECONDS >= deadline)); then
    echo "Voxel-SLAM did not write a stable optimized trajectory within ${shutdown_wait}s" >&2
    exit 1
  fi
  if ! kill -0 "$launch_pid" 2>/dev/null; then
    echo "Voxel-SLAM exited before the optimized trajectory was stable" >&2
    exit 1
  fi
  if [[ -s "$optimized_state" ]]; then
    current_size=$(stat -c %s "$optimized_state")
    if ((current_size == stable_size)); then
      stable_checks=$((stable_checks + 1))
    else
      stable_size=$current_size
      stable_checks=0
    fi
  fi
  sleep 1
done
end_ns="$(date +%s%N)"

# Voxel-SLAM writes the final trajectory and then enters ros::spin().  Shut
# down the node only after the result has remained stable for five seconds.
rosnode kill /voxelslam >/dev/null
shutdown_deadline=$((SECONDS + 30))
while kill -0 "$launch_pid" 2>/dev/null && ((SECONDS < shutdown_deadline)); do
  sleep 1
done
if kill -0 "$launch_pid" 2>/dev/null; then
  echo "Voxel-SLAM did not exit after ROS shutdown" >&2
  exit 1
fi
wait "$launch_pid" || true
launch_pid=""

if kill -0 "$recorder_pid" 2>/dev/null; then
  kill -INT "$recorder_pid" 2>/dev/null || true
  wait "$recorder_pid" 2>/dev/null || true
fi
recorder_pid=""

optimized_tum="$output_dir/results/trajectory_voxel_opt.tum"
[[ -s "$optimized_state" ]] || {
  echo "optimized Voxel-SLAM trajectory was not written: $optimized_state" >&2
  exit 1
}
awk 'NF >= 8 {print $1, $2, $3, $4, $5, $6, $7, $8}' \
  "$optimized_state" >"$optimized_tum"

opt_lines=$(wc -l <"$optimized_tum")
frontend_lines=$(wc -l <"$output_dir/results/trajectory_voxel_frontend.tum")
((opt_lines >= 3)) || { echo "optimized trajectory has only $opt_lines poses" >&2; exit 1; }
((frontend_lines >= 3)) || { echo "front-end trajectory has only $frontend_lines poses" >&2; exit 1; }

readarray -t trajectory_contract < <(python3 -c 'import math,sys
last=-math.inf; valid=invalid=nonmono=0
for line in open(sys.argv[1],encoding="utf-8"):
    try: values=[float(value) for value in line.split()]
    except ValueError: invalid+=1; continue
    if len(values)!=8 or not all(math.isfinite(value) for value in values): invalid+=1; continue
    if values[0] <= last: nonmono+=1; continue
    norm=math.sqrt(sum(value*value for value in values[4:8]))
    if abs(norm-1.0)>1e-3: invalid+=1; continue
    last=values[0]; valid+=1
print(valid); print(invalid); print(nonmono); print(last if math.isfinite(last) else 0.0)' "$optimized_tum")
valid_lines="${trajectory_contract[0]}"
invalid_lines="${trajectory_contract[1]}"
nonmonotonic_lines="${trajectory_contract[2]}"
last_stamp="${trajectory_contract[3]}"
output_ratio="$(awk -v actual="$valid_lines" -v expected="$expected_lidar_frames" 'BEGIN {printf "%.9f", actual/expected}')"
final_lidar_gap_s="$(awk -v expected="$lidar_last_s" -v actual="$last_stamp" 'BEGIN {printf "%.9f", expected-actual}')"
completion="incomplete"
if awk -v gap="$final_lidar_gap_s" 'BEGIN {exit !(gap >= -0.25 && gap <= 5.0)}'; then
  completion="reached_final_lidar"
fi

touch "$output_dir/monitor.stop"
wait "$monitor_pid" 2>/dev/null || true
monitor_pid=""
wall_time_s="$(awk -v ns="$((end_ns-start_ns))" 'BEGIN {printf "%.6f", ns/1e9}')"
voxel_repo="$voxel_ws/src/Voxel-SLAM"
voxel_commit="$(git -C "$voxel_repo" rev-parse HEAD 2>/dev/null || echo unknown)"
voxel_dirty="$(git -C "$voxel_repo" status --porcelain 2>/dev/null | wc -l)"
cat >"$output_dir/run_metadata.txt" <<EOF
method=voxel_slam_114_reference
sequence=$sequence
repeat=$repeat
bag=$bag
bag_size_bytes=$(stat -c %s "$bag")
bag_sha256=$(sha256sum "$bag" | awk '{print $1}')
bag_start_s=$bag_start_s
bag_end_s=$bag_end_s
sensor_duration_s=$sensor_duration_s
expected_lidar_frames=$expected_lidar_frames
lidar_first_s=$lidar_first_s
lidar_last_s=$lidar_last_s
lidar_timestamp_source=message_header_stamp
output_ratio=$output_ratio
final_lidar_gap_s=$final_lidar_gap_s
completion=$completion
optimized_lines=$opt_lines
frontend_lines=$frontend_lines
valid_lines=$valid_lines
invalid_lines=$invalid_lines
nonmonotonic_lines=$nonmonotonic_lines
last_stamp=$last_stamp
voxel_repo=$(realpath "$voxel_repo")
voxel_commit=$voxel_commit
voxel_dirty_paths=$voxel_dirty
algorithm_binary=$(realpath "$binary")
algorithm_binary_sha256=$(sha256sum "$binary" | awk '{print $1}')
launch_file=$(realpath "$launch_file")
launch_file_sha256=$(sha256sum "$launch_file" | awk '{print $1}')
config_file=$(realpath "$config_file")
config_file_sha256=$(sha256sum "$config_file" | awk '{print $1}')
runner_sha256=$(sha256sum "$0" | awk '{print $1}')
cpu_set=$cpu_set
allocated_cpus=$cpu_count
play_rate=$play_rate
lidar_topic=$lidar_topic
imu_topic=$imu_topic
wall_time_s=$wall_time_s
completed_at=$(date --iso-8601=seconds)
EOF

if ((valid_lines < 3 || invalid_lines != 0 || nonmonotonic_lines != 0)) || \
   [[ "$completion" != "reached_final_lidar" ]] || \
   ! awk -v ratio="$output_ratio" 'BEGIN {exit !(ratio >= 0.80 && ratio <= 1.01)}' || \
   [[ ! -s "$output_dir/resource_summary.json" ]]; then
  echo "Voxel-SLAM reference failed contract: completion=$completion final_gap=$final_lidar_gap_s valid=$valid_lines expected=$expected_lidar_frames ratio=$output_ratio invalid=$invalid_lines nonmono=$nonmonotonic_lines" >&2
  exit 1
fi

trap - EXIT INT TERM
cleanup
echo "completed method=ws_voxel_slam_114 sequence=$sequence optimized_lines=$opt_lines frontend_lines=$frontend_lines output=$output_dir"
