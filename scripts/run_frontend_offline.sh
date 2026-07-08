#!/usr/bin/env bash
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
  exe_path="$(ros2 pkg prefix lightning)/lib/lightning/run_frontend_offline"
  set +e
  "$exe_path" --help
  help_status=$?
  set -e
  [[ "$help_status" -eq 0 || "$help_status" -eq 1 ]] && exit 0
  exit "$help_status"
fi

bag_dir="${LIGHTNING_LM_INPUT_BAG:-}"
config_path="${LIGHTNING_LM_CONFIG:-$repo_dir/config/default.yaml}"
out_root="${LIGHTNING_LM_OUT_ROOT:-$repo_dir/runs}"
run_name="${LIGHTNING_LM_RUN_NAME:-run_frontend_offline_$(date +%Y%m%d_%H%M%S)}"
wait_ui="${LIGHTNING_LM_WAIT_UI:-true}"
max_lidar_frames="${LIGHTNING_LM_MAX_LIDAR_FRAMES:-0}"
if [[ $# -gt 0 && "${1:0:1}" != "-" ]]; then
  run_name="$1"
  shift
fi

if [[ -z "$bag_dir" ]]; then
  echo "LIGHTNING_LM_INPUT_BAG is required" >&2
  exit 2
fi
if [[ ! -f "$config_path" ]]; then
  echo "config not found: $config_path" >&2
  exit 2
fi

run_dir="$out_root/$run_name"
result_dir="$run_dir/results"
log_dir="$run_dir/logs"
tum_path="${LIGHTNING_LM_OUTPUT_TUM:-$result_dir/${run_name}.tum}"
mkdir -p "$result_dir" "$log_dir"

{
  echo "executable=run_frontend_offline"
  echo "repo_dir=$repo_dir"
  echo "bag_dir=$bag_dir"
  echo "config_path=$config_path"
  echo "run_name=$run_name"
  echo "run_dir=$run_dir"
  echo "tum_path=$tum_path"
  echo "wait_ui=$wait_ui"
  echo "max_lidar_frames=$max_lidar_frames"
  echo "extra_args=$*"
  date --iso-8601=seconds
} > "$run_dir/run_metadata.txt"

ros2 run lightning run_frontend_offline \
  --input_bag "$bag_dir" \
  --config "$config_path" \
  --output_tum "$tum_path" \
  --wait_ui="$wait_ui" \
  --max_lidar_frames="$max_lidar_frames" \
  "$@" \
  > "$log_dir/run_frontend_offline.stdout.log" \
  2> "$log_dir/run_frontend_offline.stderr.log"
