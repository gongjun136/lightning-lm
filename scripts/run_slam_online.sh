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
  exe_path="$(ros2 pkg prefix lightning)/lib/lightning/run_slam_online"
  set +e
  "$exe_path" --help
  help_status=$?
  set -e
  [[ "$help_status" -eq 0 || "$help_status" -eq 1 ]] && exit 0
  exit "$help_status"
fi

config_path="${LIGHTNING_LM_CONFIG:-$repo_dir/config/default.yaml}"
out_root="${LIGHTNING_LM_OUT_ROOT:-$repo_dir/runs}"
run_name="${LIGHTNING_LM_RUN_NAME:-run_slam_online_$(date +%Y%m%d_%H%M%S)}"
if [[ $# -gt 0 && "${1:0:1}" != "-" ]]; then
  run_name="$1"
  shift
fi

if [[ ! -f "$config_path" ]]; then
  echo "config not found: $config_path" >&2
  exit 2
fi

run_dir="$out_root/$run_name"
log_dir="$run_dir/logs"
mkdir -p "$log_dir"

{
  echo "executable=run_slam_online"
  echo "repo_dir=$repo_dir"
  echo "config_path=$config_path"
  echo "run_name=$run_name"
  echo "run_dir=$run_dir"
  echo "extra_args=$*"
  date --iso-8601=seconds
} > "$run_dir/run_metadata.txt"

ros2 run lightning run_slam_online \
  --config "$config_path" \
  "$@" \
  > "$log_dir/run_slam_online.stdout.log" \
  2> "$log_dir/run_slam_online.stderr.log"
