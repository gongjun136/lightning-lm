#!/usr/bin/env bash
set -euo pipefail
# [build-ins-delegate]
# Compatibility entry point: all packages now share one workspace.
exec bash "$(dirname -- "${BASH_SOURCE[0]}")/build_workspace.sh" "$@"
# [build-ins-delegate]
