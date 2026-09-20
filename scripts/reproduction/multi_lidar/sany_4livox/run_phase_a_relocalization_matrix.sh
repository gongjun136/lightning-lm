#!/usr/bin/env bash
set -uo pipefail

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 BAG_ROOT CONFIG MAP OUTPUT_ROOT [DATASET ...]" >&2
  exit 2
fi

bag_root="$1"
config="$2"
map="$3"
output_root="$4"
shift 4
datasets=("$@")
if [[ ${#datasets[@]} -eq 0 ]]; then
  datasets=(data1 data2 data3 data4 data5)
fi

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
runner="$repo_dir/scripts/run_loc_offline.sh"
offsets=(0 3 5 7 10 15 20 25 30 35)
runner_failures=0

for dataset in "${datasets[@]}"; do
  bag="$bag_root/$dataset"
  if [[ ! -f "$bag/metadata.yaml" ]]; then
    echo "missing bag metadata: $bag/metadata.yaml" >&2
    exit 3
  fi
  first_sensor_time=""
  for offset in "${offsets[@]}"; do
    printf -v tag "%02d" "$offset"
    output="$output_root/$dataset/offset_${tag}s"
    if [[ -e "$output" ]]; then
      echo "existing output: $output" >&2
      exit 4
    fi

    start_sensor_time="0"
    if [[ "$offset" -gt 0 ]]; then
      start_sensor_time="$(awk -v start="$first_sensor_time" -v delta="$offset" \
        'BEGIN { printf "%.6f", start + delta }')"
    fi

    if bash "$runner" \
      --bag "$bag" \
      --config "$config" \
      --map "$map" \
      --output-dir "$output" \
      --sequence "SANY_4lidar_${dataset}_offset${tag}" \
      --repeat 1 \
      --playback-rate 0 \
      --start-sensor-time "$start_sensor_time" \
      --max-lidar-frames 90 \
      --cpu-set 0-15 \
      --cpu-count 16 \
      --completion-tolerance 2.5 \
      --max-output-gap 0.30 \
      --watchdog-margin 1200 \
      --wait-ui false \
      --use-config-initial-pose false; then
      echo "trial ${dataset}/${tag}s runner_ok"
    else
      rc=$?
      runner_failures=$((runner_failures + 1))
      echo "trial ${dataset}/${tag}s runner_rc=$rc"
    fi

    if [[ "$offset" -eq 0 ]]; then
      stats="$output/results/localization_stats.csv"
      first_sensor_time="$(awk -F, 'NR == 2 { print $2; exit }' "$stats")"
      if [[ -z "$first_sensor_time" ]]; then
        echo "cannot determine first localization timestamp: $stats" >&2
        exit 5
      fi
      echo "dataset ${dataset} first_localization_time=$first_sensor_time"
    fi
  done
done

echo "batch_runner_failures=$runner_failures"
