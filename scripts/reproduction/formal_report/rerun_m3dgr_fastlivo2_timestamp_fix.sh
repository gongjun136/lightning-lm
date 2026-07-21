#!/usr/bin/env bash
set -u

# Rerun only the FAST-LIVO2 cells invalidated by the odometry timestamp fix.
# The order is the FAST-LIVO2 projection of the frozen 20260717 schedule.
base="${1:-/mnt/f/SLAM_AI_KnowledgeBase/code/WSL_Ubuntu_22.04/lightning-lm/runs/formal_report_20260717/m3dgr_frontend_v2}"
runner=/mnt/f/SLAM_AI_KnowledgeBase/code/WSL_Ubuntu_20.04/ros1_ws/ws_fastlivo2/src/FAST-LIVO2/reproduction/m3dgr/scripts/run_frontend_offline.sh
inventory=/mnt/f/SLAM_AI_KnowledgeBase/code/_m3dgr_work/bench/inventory/bag_inventory.json
bag_root=/mnt/f/SLAM_AI_KnowledgeBase/code/_m3dgr_work/ros1_normalized
fingerprint=6365ed995ff4124ba7ad414066c1838beef778a7c4ca64bf9acb24f917fdb031
items=(
  Dark01:1 Dark01:3 Z-Rough-Road01:1 Z-Rough-Road01:2
  Dark01:2 Grass02:1 Grass02:2 Outdoor04:1
  Z-Rough-Road01:3 Grass02:3 Outdoor04:2 Outdoor04:3
)

metadata_value() {
  local path="$1"
  local key="$2"
  awk -F= -v key="$key" '$1 == key { print $2 }' "$path"
}

index=0
for item in "${items[@]}"; do
  index=$((index + 1))
  sequence="${item%%:*}"
  repeat="${item##*:}"
  port=$((15500 + index))
  run_dir=$(printf '%s/%s/fastlivo2_lio/repeat_%02d' "$base" "$sequence" "$repeat")
  bag="$bag_root/$sequence/$sequence.bag"
  if [[ -f "$run_dir/run_metadata.txt" ]]; then
    metadata="$run_dir/run_metadata.txt"
    completion=$(metadata_value "$metadata" completion)
    invalid=$(metadata_value "$metadata" invalid_count)
    nonmonotonic=$(metadata_value "$metadata" nonmonotonic_count)
    contract=$(metadata_value "$metadata" algorithm_contract)
    subscribers=$(metadata_value "$metadata" subscriber_contract)
    if [[ "$completion" == reached_final_lidar && "$invalid" == 0 && "$nonmonotonic" == 0 && \
          "$contract" == passed && "$subscribers" == passed ]]; then
      echo "SKIP [$index/${#items[@]}] $sequence repeat=$repeat: valid result already exists"
      continue
    fi
  fi
  if [[ -e "$run_dir" ]]; then
    echo "existing run directory: $run_dir" >&2
    exit 2
  fi

  echo "START [$index/${#items[@]}] $sequence repeat=$repeat port=$port"
  /usr/bin/env \
    BENCH_PLAY_RATE=1.0 \
    BENCH_CPUSET=0-7 \
    BENCH_CPU_COUNT=8 \
    BENCH_INVENTORY_JSON="$inventory" \
    BENCH_EXPERIMENT_FINGERPRINT="$fingerprint" \
    BENCH_DEVELOPMENT_EXPOSED=false \
    BENCH_CONFIRMATORY=true \
    BENCH_ROS_PORT="$port" \
    bash "$runner" "$bag" "$sequence" "$repeat" "$run_dir"
  runner_rc=$?

  metadata="$run_dir/run_metadata.txt"
  if [[ ! -f "$metadata" ]]; then
    echo "missing metadata: $metadata (runner rc=$runner_rc)" >&2
    exit 3
  fi
  completion=$(metadata_value "$metadata" completion)
  invalid=$(metadata_value "$metadata" invalid_count)
  nonmonotonic=$(metadata_value "$metadata" nonmonotonic_count)
  contract=$(metadata_value "$metadata" algorithm_contract)
  subscribers=$(metadata_value "$metadata" subscriber_contract)
  frames=$(metadata_value "$metadata" trajectory_lines)
  max_gap=$(metadata_value "$metadata" maximum_output_gap_s)
  gaps=$(metadata_value "$metadata" excessive_output_gap_count)
  echo "DONE [$index/${#items[@]}] $sequence repeat=$repeat rc=$runner_rc frames=$frames max_gap=$max_gap gaps=$gaps"

  if [[ "$completion" != reached_final_lidar || "$invalid" != 0 || "$nonmonotonic" != 0 || \
        "$contract" != passed || "$subscribers" != passed ]]; then
    echo "run contract failed: $sequence repeat=$repeat" >&2
    exit 4
  fi
done
