#!/usr/bin/env bash
# Production preset for the shared SANY online runner. Explicit caller values
# still win, which makes controlled A/B deployment possible without editing it.
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export LIGHTNING_LM_RUN_MODE="production"
export LIGHTNING_LM_COMPUTE_PROFILE="${LIGHTNING_LM_COMPUTE_PROFILE:-0}"
export LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD="${LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD:-1}"
export SANY_RECORD_BAG="${SANY_RECORD_BAG:-0}"

exec bash "${script_dir}/run_sany_online_diagnostics.sh" "$@"
