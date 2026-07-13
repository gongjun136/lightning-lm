#!/usr/bin/env bash
set -uo pipefail

fault_root="${1:?usage: run_sany_formal_fault_matrix.sh FAULT_ROOT CONFIG OUTPUT_ROOT}"
config="${2:?missing config}"
output_root="${3:?missing output root}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
runner="$script_dir/run_sany_online_contract.sh"
scenarios=(
  missing_imu114_all
  interrupt_imu114_40_50
  interrupt_lidar127_40_60
  interrupt_lidar114_40_60
)

mkdir -p "$output_root"
failures=0
for scenario in "${scenarios[@]}"; do
  bag="$fault_root/$scenario"
  contract="$bag/fault_contract.json"
  output="$output_root/$scenario"
  printf 'running %s\n' "$scenario"
  if bash "$runner" "$bag" "$config" "$output" "$contract"; then
    printf 'passed %s\n' "$scenario"
  else
    rc=$?
    printf 'failed %s rc=%s\n' "$scenario" "$rc" >&2
    failures=$((failures + 1))
  fi
done
printf 'fault scenarios=%s failures=%s\n' "${#scenarios[@]}" "$failures"
(( failures == 0 ))
