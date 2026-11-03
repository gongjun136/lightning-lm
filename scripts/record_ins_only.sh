#!/usr/bin/env bash
set -euo pipefail
# [ins-record-inputs]
[[ $# -eq 2 ]] || { echo "Usage: $0 site.yaml output_bag (source ROS/workspaces first)" >&2; exit 2; }
dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
topic_text="$(python3 "$dir/ins_topics.py" "$1")"
mapfile -t topics <<< "$topic_text"
output_text="$(python3 "$dir/ins_topics.py" --outputs "$1")"
mapfile -t outputs <<< "$output_text"
# [ins-record-inputs]
# [ins-record-launch]
exec ros2 bag record -s mcap -o "$2" "${topics[@]}" "${outputs[@]}"
# [ins-record-launch]
