#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
binary="${repo_dir}/install/lightning/lib/lightning/align_maps_offline"
if [[ ! -x "${binary}" ]]; then
  binary="${repo_dir}/bin/align_maps_offline"
fi
if [[ ! -x "${binary}" ]]; then
  echo "align_maps_offline is not built; run: colcon build --packages-select lightning" >&2
  exit 2
fi
exec "${binary}" "$@"
