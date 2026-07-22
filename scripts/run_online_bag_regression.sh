#!/usr/bin/env bash
# Run an online SLAM or localization node against a ROS 2 bag at 1x speed.
set -euo pipefail

if [[ $# -lt 4 || $# -gt 5 ]]; then
  echo "usage: $0 <slam|loc> <run-name> <config.yaml> <bag-dir> [map-dir]" >&2
  exit 2
fi

mode="$1"
run_name="$2"
config_path="$(realpath "$3")"
bag_path="$(realpath "$4")"
map_path="${5:-}"
if [[ "$mode" != "slam" && "$mode" != "loc" ]]; then
  echo "mode must be slam or loc" >&2
  exit 2
fi
if [[ "$mode" == "loc" && -z "$map_path" ]]; then
  echo "localization mode requires a map directory" >&2
  exit 2
fi
if [[ -n "$map_path" ]]; then
  map_path="$(realpath "$map_path")"
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="${LIGHTNING_LM_REPO_DIR:-$(cd "$script_dir/.." && pwd)}"
out_root="${LIGHTNING_LM_OUT_ROOT:-$repo_dir/runs/online_regression}"
run_dir="$out_root/$run_name"
post_wait="${LIGHTNING_LM_POST_WAIT_SECONDS:-15}"
discovery_delay="${LIGHTNING_LM_DISCOVERY_DELAY_SECONDS:-10}"
embedded_playback="${LIGHTNING_LM_EMBEDDED_PLAYBACK:-1}"
ros_setup="${LIGHTNING_LM_ROS_SETUP:-/opt/ros/humble/setup.bash}"
install_setup="${LIGHTNING_LM_INSTALL_SETUP:-$repo_dir/install/setup.bash}"

mkdir -p "$run_dir/logs"
run_dir="$(realpath "$run_dir")"
set +u
source "$ros_setup"
source "$install_setup"
set -u

export ROS_DOMAIN_ID="${LIGHTNING_LM_ROS_DOMAIN_ID:-83}"
node_pid=""
recorder_pid=""

stop_process() {
  local pid="$1"
  [[ -z "$pid" ]] && return
  if kill -0 "$pid" 2>/dev/null; then
    kill -INT "$pid" 2>/dev/null || true
    for _ in $(seq 1 20); do
      kill -0 "$pid" 2>/dev/null || return 0
      sleep 0.5
    done
    kill -TERM "$pid" 2>/dev/null || true
  fi
}

cleanup() {
  stop_process "$recorder_pid"
  stop_process "$node_pid"
}
trap cleanup EXIT INT TERM

{
  echo "mode=$mode"
  echo "config_path=$config_path"
  echo "bag_path=$bag_path"
  echo "map_path=$map_path"
  echo "playback_rate=1.0"
  echo "post_wait_seconds=$post_wait"
  echo "discovery_delay_seconds=$discovery_delay"
  echo "embedded_playback=$embedded_playback"
  echo "ros_domain_id=$ROS_DOMAIN_ID"
  date --iso-8601=seconds
} > "$run_dir/run_metadata.txt"

cd "$run_dir"
exe_prefix="$(ros2 pkg prefix lightning)/lib/lightning"
if [[ "$embedded_playback" == "1" ]]; then
  if [[ "$mode" == "slam" ]]; then
    /usr/bin/time -f 'elapsed_seconds=%e\nmax_rss_kb=%M\nuser_seconds=%U\nsystem_seconds=%S' \
      -o node_resource.txt \
      "$exe_prefix/run_slam_online" --config "$config_path" --bag "$bag_path" --playback_rate 1.0 \
      --post_wait_seconds "$post_wait" --save_map ./data/online_map \
      > logs/node.stdout.log 2> logs/node.stderr.log
  else
    /usr/bin/time -f 'elapsed_seconds=%e\nmax_rss_kb=%M\nuser_seconds=%U\nsystem_seconds=%S' \
      -o node_resource.txt \
      "$exe_prefix/run_loc_online" --config "$config_path" --map "$map_path" --bag "$bag_path" \
      --playback_rate 1.0 --post_wait_seconds "$post_wait" --output_tum ./trajectory_loc_online.tum \
      --output_high_frequency_tum ./trajectory_loc_online_high_frequency.tum \
      > logs/node.stdout.log 2> logs/node.stderr.log
  fi
  trap - EXIT INT TERM
  echo "online regression complete: $run_dir"
  exit 0
fi

if [[ "$mode" == "slam" ]]; then
  "$exe_prefix/run_slam_online" --config "$config_path" \
    > logs/node.stdout.log 2> logs/node.stderr.log &
else
  "$exe_prefix/run_loc_online" --config "$config_path" --map "$map_path" \
    > logs/node.stdout.log 2> logs/node.stderr.log &
fi
node_pid=$!

sleep 5
if ! kill -0 "$node_pid" 2>/dev/null; then
  wait "$node_pid" || true
  echo "online node exited before playback" >&2
  exit 3
fi

if [[ "$mode" == "loc" ]]; then
  ros2 bag record -o output_pose_bag /slamPoseRaw_topic /PosRes \
    > logs/recorder.stdout.log 2> logs/recorder.stderr.log &
  recorder_pid=$!
  sleep 2
fi

playback_start="$(date +%s.%N)"
/usr/bin/time -f 'elapsed_seconds=%e\nmax_rss_kb=%M\nuser_seconds=%U\nsystem_seconds=%S' \
  -o playback_resource.txt \
  ros2 bag play "$bag_path" --rate 1.0 --delay "$discovery_delay" \
  --read-ahead-queue-size 1000 --disable-keyboard-controls \
  > logs/playback.stdout.log 2> logs/playback.stderr.log
playback_end="$(date +%s.%N)"
{
  echo "playback_wall_start=$playback_start"
  echo "playback_wall_end=$playback_end"
} >> run_metadata.txt

sleep "$post_wait"
if [[ "$mode" == "slam" ]]; then
  timeout 240 ros2 service call /lightning/save_map lightning/srv/SaveMap "{map_id: online_map}" \
    > logs/save_map.stdout.log 2> logs/save_map.stderr.log
else
  stop_process "$recorder_pid"
  wait "$recorder_pid" 2>/dev/null || true
  recorder_pid=""
fi

stop_process "$node_pid"
wait "$node_pid" 2>/dev/null || true
node_pid=""
trap - EXIT INT TERM

if [[ "$mode" == "loc" ]]; then
  python3 "$script_dir/extract_pose_tum_from_rosbag2.py" \
    --bag output_pose_bag --topic /slamPoseRaw_topic --output trajectory_loc_online.tum \
    > logs/extract_pose.log 2>&1
fi

echo "online regression complete: $run_dir"
