#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
binary="${repo_dir}/install/lightning_lm/lib/lightning_lm/align_maps_offline"
if [[ ! -x "${binary}" ]]; then
  binary="${repo_dir}/bin/align_maps_offline"
fi
if [[ ! -x "${binary}" ]]; then
  echo 'align_maps_offline is not built; run: MAKEFLAGS="-j4" colcon build --packages-up-to lightning_lm' >&2
  exit 2
fi
exec "${binary}" "$@"
