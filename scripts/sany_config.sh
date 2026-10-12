#!/usr/bin/env bash
# Resolve SANY templates before input checks; leave other datasets' configurations unchanged.
resolve_sany_config() {
  local input="$1" run_directory="$2"
  if grep -q '^sany_vehicle_config:' "$input"; then
    python3 "${script_dir}/resolve_sany_config.py" \
      --base-config "$input" \
      --vehicle-config "${SANY_VEHICLE_CONFIG:-/opt/conf/vehicle_config.yaml}" \
      --output "$run_directory/config.yaml"
  else
    printf '%s\n' "$input"
  fi
}
