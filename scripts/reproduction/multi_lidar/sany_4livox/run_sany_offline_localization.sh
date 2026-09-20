#!/usr/bin/env bash
# Build a SANY 4-LiDAR tiled map and immediately validate offline localization on it.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_sany_offline_localization.sh [options]

Options:
  --bag PATH                    ROS2 bag directory
  --mapping-config PATH         YAML used by the SLAM/map-building stage
  --localization-config PATH    YAML used by the localization stage
  --config PATH                 Legacy: use one YAML for both stages
  --output-root PATH            Parent directory for both runs
  --run-prefix NAME             Prefix for run directories
  --sequence NAME               Sequence label (default: SANY_4lidar)
  --repeat INDEX                Repeat label (default: 1)
  --playback-rate RATE          Sensor-time pacing; 0 disables pacing (default: 0)
  --max-lidar-frames COUNT      Stop both stages after COUNT fused frames (default: 0, all)
  --completion-tolerance SEC    Allowed trajectory tail difference (default: 2.5 for this SANY bag)
  --cpu-set LIST                Linux CPU list (default: 0-7)
  --cpu-count COUNT             Allocated logical CPUs (default: 8)
  --watchdog-margin SEC         Extra watchdog wall time for each stage (default: 300)
  --wait-ui BOOL                Wait for UI close in each stage (default: false)
  -h, --help                    Show this help

Defaults target the local SANY dataset layout under /mnt/f/datasets/SANY/mid360/data_20260701.
EOF
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="${LIGHTNING_LM_REPO_DIR:-$(cd "$script_dir/../../../.." && pwd)}"
data_dir="${LIGHTNING_LM_SANY_DATA_DIR:-/mnt/f/datasets/SANY/mid360/data_20260701}"

bag_dir="${LIGHTNING_LM_INPUT_BAG:-$data_dir/rosbag2_2026_07_01-17_30_42_merged}"
legacy_config_path="${LIGHTNING_LM_CONFIG:-}"
mapping_config_path="${LIGHTNING_LM_MAPPING_CONFIG:-${legacy_config_path:-$repo_dir/config/reproduction/multi_lidar/sany_4livox/sany_4lidar_mapping.yaml}}"
localization_config_path="${LIGHTNING_LM_LOCALIZATION_CONFIG:-${legacy_config_path:-$repo_dir/config/reproduction/multi_lidar/sany_4livox/sany_4lidar_localization_solid.yaml}}"
output_root="${LIGHTNING_LM_OUT_ROOT:-$data_dir/lightning_lm_4lidar_offline_localization_20260713}"
run_prefix="${LIGHTNING_LM_RUN_NAME:-sany_4lidar_offline_loc_$(date +%Y%m%d_%H%M%S)}"
sequence="SANY_4lidar"
repeat="1"
playback_rate="0"
max_lidar_frames="0"
completion_tolerance="2.5"
cpu_set="0-7"
allocated_cpus="8"
watchdog_margin="300"
wait_ui="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bag) bag_dir="${2:?missing value for --bag}"; shift 2 ;;
    --mapping-config) mapping_config_path="${2:?missing value for --mapping-config}"; shift 2 ;;
    --localization-config) localization_config_path="${2:?missing value for --localization-config}"; shift 2 ;;
    --config)
      mapping_config_path="${2:?missing value for --config}"
      localization_config_path="$mapping_config_path"
      shift 2
      ;;
    --output-root) output_root="${2:?missing value for --output-root}"; shift 2 ;;
    --run-prefix) run_prefix="${2:?missing value for --run-prefix}"; shift 2 ;;
    --sequence) sequence="${2:?missing value for --sequence}"; shift 2 ;;
    --repeat) repeat="${2:?missing value for --repeat}"; shift 2 ;;
    --playback-rate) playback_rate="${2:?missing value for --playback-rate}"; shift 2 ;;
    --max-lidar-frames) max_lidar_frames="${2:?missing value for --max-lidar-frames}"; shift 2 ;;
    --completion-tolerance) completion_tolerance="${2:?missing value for --completion-tolerance}"; shift 2 ;;
    --cpu-set) cpu_set="${2:?missing value for --cpu-set}"; shift 2 ;;
    --cpu-count) allocated_cpus="${2:?missing value for --cpu-count}"; shift 2 ;;
    --watchdog-margin) watchdog_margin="${2:?missing value for --watchdog-margin}"; shift 2 ;;
    --wait-ui) wait_ui="${2:?missing value for --wait-ui}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ ! -d "$bag_dir" || ! -f "$bag_dir/metadata.yaml" ]]; then
  echo "SQLite3 ROS2 bag not found: $bag_dir" >&2
  exit 2
fi
if [[ ! -f "$mapping_config_path" ]]; then
  echo "mapping configuration not found: $mapping_config_path" >&2
  exit 2
fi
if [[ ! -f "$localization_config_path" ]]; then
  echo "localization configuration not found: $localization_config_path" >&2
  exit 2
fi

bag_dir="$(realpath "$bag_dir")"
mapping_config_path="$(realpath "$mapping_config_path")"
localization_config_path="$(realpath "$localization_config_path")"
output_root="$(mkdir -p "$output_root" && cd "$output_root" && pwd)"
slam_run_dir="$output_root/${run_prefix}_slam_map"
loc_run_dir="$output_root/${run_prefix}_localization"

if [[ -e "$slam_run_dir" || -e "$loc_run_dir" ]]; then
  echo "refusing to overwrite existing run directories:" >&2
  echo "  $slam_run_dir" >&2
  echo "  $loc_run_dir" >&2
  exit 2
fi

common_args=(
  --bag "$bag_dir"
  --sequence "$sequence"
  --repeat "$repeat"
  --playback-rate "$playback_rate"
  --max-lidar-frames "$max_lidar_frames"
  --completion-tolerance "$completion_tolerance"
  --cpu-set "$cpu_set"
  --cpu-count "$allocated_cpus"
  --watchdog-margin "$watchdog_margin"
  --wait-ui "$wait_ui"
)

bash "$repo_dir/scripts/run_slam_offline.sh" \
  "${common_args[@]}" \
  --config "$mapping_config_path" \
  --output-dir "$slam_run_dir"

map_path="$slam_run_dir/data/new_map"
reference_tum="$slam_run_dir/results/trajectory_slam.tum"
if [[ ! -f "$map_path/index.txt" || ! -s "$reference_tum" ]]; then
  echo "SLAM stage did not produce required map/reference trajectory" >&2
  exit 4
fi

bash "$repo_dir/scripts/run_loc_offline.sh" \
  "${common_args[@]}" \
  --config "$localization_config_path" \
  --map "$map_path" \
  --reference-tum "$reference_tum" \
  --output-dir "$loc_run_dir"

{
  echo "method=lightning_lm_sany_offline_slam_then_localization"
  echo "sequence=$sequence"
  echo "repeat=$repeat"
  echo "bag=$bag_dir"
  echo "mapping_config=$mapping_config_path"
  echo "localization_config=$localization_config_path"
  echo "slam_run_dir=$slam_run_dir"
  echo "localization_run_dir=$loc_run_dir"
  echo "map_path=$map_path"
  echo "reference_tum=$reference_tum"
  echo "localization_summary=$loc_run_dir/results/localization_summary.json"
  echo "reference_errors=$loc_run_dir/results/trajectory_reference_errors.csv"
  echo "completed_at=$(date --iso-8601=seconds)"
} > "$output_root/${run_prefix}_summary.txt"

echo "completed SANY offline localization pipeline"
echo "slam_run=$slam_run_dir"
echo "loc_run=$loc_run_dir"
