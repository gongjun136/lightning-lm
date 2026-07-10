#!/usr/bin/env bash
# Reproducible ROS2 offline frontend template. Keep experiment YAML and outputs
# outside the repository by setting LIGHTNING_LM_CONFIG and LIGHTNING_LM_OUT_ROOT.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="${LIGHTNING_LM_REPO_DIR:-$(cd "$script_dir/.." && pwd)}"
ros_setup="${LIGHTNING_LM_ROS_SETUP:-/opt/ros/humble/setup.bash}"
install_setup="${LIGHTNING_LM_INSTALL_SETUP:-$repo_dir/install/setup.bash}"

set +u
source "$ros_setup"
source "$install_setup"
set -u

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  exec ros2 run lightning run_frontend_offline --help
fi

bag_dir="${LIGHTNING_LM_INPUT_BAG:-}"
config_path="${LIGHTNING_LM_CONFIG:-$repo_dir/config/default.yaml}"
out_root="${LIGHTNING_LM_OUT_ROOT:-$repo_dir/runs}"
run_name="${LIGHTNING_LM_RUN_NAME:-run_frontend_offline_$(date +%Y%m%d_%H%M%S)}"
wait_ui="${LIGHTNING_LM_WAIT_UI:-false}"
max_lidar_frames="${LIGHTNING_LM_MAX_LIDAR_FRAMES:-0}"
if [[ $# -gt 0 && "${1:0:1}" != "-" ]]; then
  run_name="$1"
  shift
fi

if [[ -z "$bag_dir" || ! -e "$bag_dir" ]]; then
  echo "LIGHTNING_LM_INPUT_BAG is required and must exist: $bag_dir" >&2
  exit 2
fi
if [[ ! -f "$config_path" ]]; then
  echo "configuration not found: $config_path" >&2
  exit 2
fi

bag_dir="$(realpath "$bag_dir")"
config_path="$(realpath "$config_path")"
out_root="$(mkdir -p "$out_root" && cd "$out_root" && pwd)"
run_dir="$out_root/$run_name"
result_dir="$run_dir/results"
log_dir="$run_dir/logs"
mkdir -p "$result_dir" "$log_dir"

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

imu_tum="$(resolve_output "${LIGHTNING_LM_OUTPUT_TUM:-}" "$result_dir/trajectory_imu.tum")"
lidar_tum="$(resolve_output "${LIGHTNING_LM_OUTPUT_LIDAR_TUM:-}" "$result_dir/trajectory_lidar114.tum")"
rear_tum="$(resolve_output "${LIGHTNING_LM_OUTPUT_REAR_AXLE_TUM:-}" "$result_dir/trajectory_rear_axle.tum")"
map_pcd="$(resolve_output "${LIGHTNING_LM_OUTPUT_MAP:-}" "$result_dir/map_lio.pcd")"
frame_stats="$(resolve_output "${LIGHTNING_LM_OUTPUT_FRAME_STATS:-}" "$result_dir/frame_stats.csv")"

{
  echo "executable=run_frontend_offline"
  echo "repo_dir=$repo_dir"
  echo "bag_dir=$bag_dir"
  echo "config_path=$config_path"
  echo "run_name=$run_name"
  echo "run_dir=$run_dir"
  echo "imu_tum=$imu_tum"
  echo "lidar_tum=$lidar_tum"
  echo "rear_tum=$rear_tum"
  echo "map_pcd=$map_pcd"
  echo "frame_stats=$frame_stats"
  echo "wait_ui=$wait_ui"
  echo "max_lidar_frames=$max_lidar_frames"
  echo "extra_args=$*"
  date --iso-8601=seconds
} > "$run_dir/run_metadata.txt"

(
  cd "$run_dir"
  ros2 run lightning run_frontend_offline \
    --input_bag "$bag_dir" \
    --config "$config_path" \
    --output_tum "$imu_tum" \
    --output_lidar_tum "$lidar_tum" \
    --output_rear_axle_tum "$rear_tum" \
    --output_map "$map_pcd" \
    --output_frame_stats_csv "$frame_stats" \
    --wait_ui="$wait_ui" \
    --max_lidar_frames="$max_lidar_frames" \
    "$@"
) > "$log_dir/run_frontend_offline.stdout.log" \
  2> "$log_dir/run_frontend_offline.stderr.log"

echo "completed: $run_dir"
