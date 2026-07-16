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
  --ros-port PORT        Private ROS master port (default: 11331)
  --shutdown-wait SEC    Maximum wait for final optimization (default: 900)
  -h, --help             Show this help
EOF
}

bag=""
output_dir=""
sequence=""
voxel_ws="/mnt/f/SLAM_AI_KnowledgeBase/code/WSL_Ubuntu_20.04/ros1_ws/ws_voxel_slam"
ros_port=11331
shutdown_wait=900

while (($#)); do
  case "$1" in
    --bag) bag="${2:?missing value for --bag}"; shift 2 ;;
    --output-dir) output_dir="${2:?missing value for --output-dir}"; shift 2 ;;
    --sequence) sequence="${2:?missing value for --sequence}"; shift 2 ;;
    --voxel-ws) voxel_ws="${2:?missing value for --voxel-ws}"; shift 2 ;;
    --ros-port) ros_port="${2:?missing value for --ros-port}"; shift 2 ;;
    --shutdown-wait) shutdown_wait="${2:?missing value for --shutdown-wait}"; shift 2 ;;
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

setup="$voxel_ws/devel/setup.bash"
recorder="$voxel_ws/src/Voxel-SLAM/reproduction/m3dgr/scripts/trajectory_recorder.py"
[[ -f "$setup" ]] || { echo "Voxel-SLAM setup does not exist: $setup" >&2; exit 2; }
[[ -f "$recorder" ]] || { echo "trajectory recorder does not exist: $recorder" >&2; exit 2; }

if [[ -d "$output_dir" ]] && find "$output_dir" -mindepth 1 -print -quit | grep -q .; then
  echo "output directory is not empty: $output_dir" >&2
  exit 2
fi
mkdir -p "$output_dir/logs" "$output_dir/results" "$output_dir/data"

# shellcheck disable=SC1091
source /opt/ros/noetic/setup.bash
# shellcheck disable=SC1090
source "$setup"
set -u
export ROS_MASTER_URI="http://127.0.0.1:${ros_port}"
export ROS_HOSTNAME=127.0.0.1

core_pid=""
launch_pid=""
recorder_pid=""
cleanup() {
  set +e
  for pid in "$recorder_pid" "$launch_pid" "$core_pid"; do
    if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then
      kill -INT "$pid" 2>/dev/null || true
    fi
  done
  sleep 1
  for pid in "$recorder_pid" "$launch_pid" "$core_pid"; do
    if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then
      kill -TERM "$pid" 2>/dev/null || true
    fi
  done
  wait "$recorder_pid" "$launch_pid" "$core_pid" 2>/dev/null || true
}
trap cleanup EXIT INT TERM

roscore -p "$ros_port" >"$output_dir/logs/roscore.log" 2>&1 &
core_pid=$!
for _ in $(seq 1 100); do
  rosparam list >/dev/null 2>&1 && break
  sleep 0.1
done
rosparam list >/dev/null 2>&1 || { echo "ROS master did not start" >&2; exit 1; }
rosparam set /use_sim_time true

roslaunch voxel_slam vxlm_sany_20260701_livox_pc2_114.launch \
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

python3 "$recorder" \
  --mode tf --topic /tf --parent-frame camera_init --child-frame aft_mapped \
  --tx -0.011 --ty -0.02329 --tz 0.04412 \
  --output "$output_dir/results/trajectory_voxel_frontend.tum" \
  --status "$output_dir/results/trajectory_voxel_frontend.status.json" \
  >"$output_dir/logs/trajectory_recorder.log" 2>&1 &
recorder_pid=$!

rosbag play "$bag" --clock --quiet --wait-for-subscribers \
  --topics /livox/lidar_192_168_1_114 /livox/imu_192_168_1_114 \
  >"$output_dir/logs/rosbag_play.log" 2>&1

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

trap - EXIT INT TERM
cleanup
echo "completed method=ws_voxel_slam_114 sequence=$sequence optimized_lines=$opt_lines frontend_lines=$frontend_lines output=$output_dir"
