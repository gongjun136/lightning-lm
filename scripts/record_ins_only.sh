#!/usr/bin/env bash
set -euo pipefail
[[ $# -eq 2 ]] || { echo "Usage: $0 site.yaml output_bag (source ROS/workspaces first)" >&2; exit 2; }
dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
topic_text="$(python3 "$dir/ins_topics.py" "$1")"
mapfile -t topics <<< "$topic_text"
exec ros2 bag record -s mcap -o "$2" "${topics[@]}" /PosRes /slamPoseRaw_topic \
  /localization/pose_vel /localization/loc_status /localization/fault_status \
  /localization/ins_diagnostics /diagnostics/heartbeat/lightning_slam
