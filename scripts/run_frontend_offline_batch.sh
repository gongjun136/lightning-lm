#!/usr/bin/env bash
# Run one SQLite3 ROS2 bag with multiple Lightning-LM frontend YAML files.
set -uo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_frontend_offline_batch.sh --bag ROS2_BAG --output-root DIR [options] [-- RUNNER_OPTIONS]

Required:
  --bag PATH              SQLite3 ROS2 bag directory
  --output-root PATH      Parent directory for all batch outputs

Options:
  --config-dir PATH       YAML directory (default: SANY four-LiDAR matrix configs)
  --config NAME_OR_PATH   Run only this YAML; may be specified multiple times
  --repeats COUNT         Runs per YAML (default: 1; use 3 for the formal 7x3 matrix)
  --dry-run               Print commands without creating outputs
  -h, --help              Show this help

Arguments after `--` are forwarded unchanged to run_frontend_offline.sh, for example:
  -- --playback-rate 1.0 --max-lidar-frames 160 --wait-ui false

Without --config, every *.yaml file in --config-dir is run in filename order.
Each run is written to OUTPUT_ROOT/<config-name>/repeat_NN/.
EOF
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "$script_dir/.." && pwd)"
runner="$script_dir/run_frontend_offline.sh"
config_dir="$repo_dir/config/reproduction/multi_lidar/sany_4livox/matrix"
bag_dir=""
output_root=""
repeats="1"
dry_run="false"
selected_configs=()
runner_args=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bag) bag_dir="${2:?missing value for --bag}"; shift 2 ;;
    --output-root) output_root="${2:?missing value for --output-root}"; shift 2 ;;
    --config-dir) config_dir="${2:?missing value for --config-dir}"; shift 2 ;;
    --config) selected_configs+=("${2:?missing value for --config}"); shift 2 ;;
    --repeats) repeats="${2:?missing value for --repeats}"; shift 2 ;;
    --dry-run) dry_run="true"; shift ;;
    -h|--help) usage; exit 0 ;;
    --) shift; runner_args=("$@"); break ;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$bag_dir" || -z "$output_root" ]]; then
  usage >&2
  exit 2
fi
if [[ ! "$repeats" =~ ^[1-9][0-9]*$ ]]; then
  echo "--repeats must be a positive integer: $repeats" >&2
  exit 2
fi
if [[ ! -x "$runner" ]]; then
  echo "offline runner is missing or not executable: $runner" >&2
  exit 2
fi
if [[ ! -d "$bag_dir" || ! -f "$bag_dir/metadata.yaml" ]]; then
  echo "SQLite3 ROS2 bag not found: $bag_dir" >&2
  exit 2
fi
if [[ ! -d "$config_dir" ]]; then
  echo "configuration directory not found: $config_dir" >&2
  exit 2
fi

bag_dir="$(realpath "$bag_dir")"
config_dir="$(realpath "$config_dir")"
output_root="$(realpath -m "$output_root")"

configs=()
if [[ ${#selected_configs[@]} -eq 0 ]]; then
  shopt -s nullglob
  configs=("$config_dir"/*.yaml)
  shopt -u nullglob
else
  for selected in "${selected_configs[@]}"; do
    candidate="$selected"
    if [[ ! -f "$candidate" ]]; then
      candidate="$config_dir/$selected"
    fi
    if [[ ! -f "$candidate" && "$candidate" != *.yaml ]]; then
      candidate="$candidate.yaml"
    fi
    if [[ ! -f "$candidate" ]]; then
      echo "configuration not found: $selected" >&2
      exit 2
    fi
    configs+=("$(realpath "$candidate")")
  done
fi
if [[ ${#configs[@]} -eq 0 ]]; then
  echo "no YAML configurations found in: $config_dir" >&2
  exit 2
fi

mapfile -t configs < <(printf '%s\n' "${configs[@]}" | sort -u)
declare -A seen_names=()
for config in "${configs[@]}"; do
  name="$(basename "$config" .yaml)"
  if [[ -n "${seen_names[$name]:-}" ]]; then
    echo "duplicate configuration name would overwrite outputs: $name" >&2
    exit 2
  fi
  seen_names[$name]=1
done

summary="$output_root/batch_summary.tsv"
if [[ "$dry_run" == "false" ]]; then
  if [[ -e "$summary" ]]; then
    echo "refusing to overwrite existing batch summary: $summary" >&2
    exit 2
  fi
  mkdir -p "$output_root"
  printf 'config\tconfig_sha256\trepeat\tstatus\treturncode\toutput_dir\n' >"$summary"
fi

total=$(( ${#configs[@]} * repeats ))
completed=0
failures=0
echo "batch plan: configs=${#configs[@]} repeats=$repeats runs=$total output=$output_root"

for ((repeat = 1; repeat <= repeats; repeat++)); do
  for config in "${configs[@]}"; do
    name="$(basename "$config" .yaml)"
    run_dir="$output_root/$name/$(printf 'repeat_%02d' "$repeat")"
    command=(
      bash "$runner"
      --bag "$bag_dir"
      --config "$config"
      --output-dir "$run_dir"
      --sequence "$name"
      --repeat "$repeat"
      "${runner_args[@]}"
    )
    completed=$((completed + 1))
    printf '[%d/%d] ' "$completed" "$total"
    printf '%q ' "${command[@]}"
    printf '\n'
    if [[ "$dry_run" == "true" ]]; then
      continue
    fi

    if "${command[@]}"; then
      returncode=0
      status="passed"
    else
      returncode=$?
      status="failed"
      failures=$((failures + 1))
    fi
    config_sha256="$(sha256sum "$config" | awk '{print $1}')"
    printf '%s\t%s\t%d\t%s\t%d\t%s\n' \
      "$name" "$config_sha256" "$repeat" "$status" "$returncode" "$run_dir" >>"$summary"
  done
done

if [[ "$dry_run" == "true" ]]; then
  echo "dry-run complete: $total commands"
  exit 0
fi
echo "batch complete: runs=$total failures=$failures summary=$summary"
if (( failures > 0 )); then
  exit 4
fi
