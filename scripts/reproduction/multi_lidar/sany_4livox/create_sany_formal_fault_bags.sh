#!/usr/bin/env bash
set -euo pipefail

set +u
source /opt/ros/humble/setup.bash
set -u

input="${1:?usage: create_sany_formal_fault_bags.sh INPUT OUTPUT_ROOT}"
root="${2:?missing output root}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
creator="$script_dir/create_sany_fault_bag.py"

run() {
  local name="$1"
  local topic="$2"
  local mode="$3"
  local start="$4"
  local duration="$5"
  printf 'creating %s\n' "$name"
  python3 "$creator" "$input" "$root/$name" \
    --scenario "$name" --target-topic "$topic" --mode "$mode" \
    --start-offset-s "$start" --duration-s "$duration"
}

mkdir -p "$root"
run missing_imu114_all /livox/imu_192_168_1_114 drop_all 40 20
run interrupt_imu114_40_50 /livox/imu_192_168_1_114 interrupt 40 10
run interrupt_lidar127_40_60 /livox/lidar_192_168_1_127 interrupt 40 20
run interrupt_lidar114_40_60 /livox/lidar_192_168_1_114 interrupt 40 20
