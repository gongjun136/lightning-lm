#!/usr/bin/env bash
# Run SANY multi-LiDAR localization together with a lossless compressed MCAP
# recorder, a /PosRes watchdog, and incident snapshots. Localization is never
# restarted by this script; the in-process global relocalizer owns recovery.
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="${LIGHTNING_LM_REPO_DIR:-$(cd -- "${script_dir}/.." && pwd)}"
ros_setup="${LIGHTNING_LM_ROS_SETUP:-/opt/ros/humble/setup.bash}"
install_setup="${LIGHTNING_LM_INSTALL_SETUP:-${repo_dir}/install/setup.bash}"
config_path="${LIGHTNING_LM_CONFIG:-${repo_dir}/config/reproduction/multi_lidar/sany_3livox/sany_3lidar_localization_blind5.yaml}"
map_path="${SANY_MAP_PATH:-}"
out_root="${LIGHTNING_LM_OUT_ROOT:-/home/nvidia/project/gj_ws/runs}"
run_name="${1:-sany_loc_diag_$(date +%Y%m%d_%H%M%S)}"
qos_file="${SANY_RECORD_QOS_FILE:-${script_dir}/config/sany_localization_record_qos.yaml}"
topic_wait_seconds="${SANY_TOPIC_WAIT_SECONDS:-60}"
posres_timeout_seconds="${SANY_POSRES_TIMEOUT_SECONDS:-2}"
min_free_gb="${SANY_MIN_FREE_GB:-20}"

usage() {
  cat <<'EOF'
Usage:
  LIGHTNING_LM_CONFIG=/absolute/config.yaml \
  SANY_MAP_PATH=/absolute/map/path \
  scripts/run_sany_online_diagnostics.sh [run_name]

Required before launch:
  1. Start every Livox driver configured by LIGHTNING_LM_CONFIG and the primary IMU.
  2. Start: ros2 run livox_ros_driver2 pointcloud_zstd_compressor

Useful environment variables:
  LIGHTNING_LM_INSTALL_SETUP  Built workspace setup.bash
  LIGHTNING_LM_OUT_ROOT       Run root (default: ~/project/gj_ws/runs)
  SANY_TOPIC_WAIT_SECONDS     Compressed-input discovery timeout (default: 60)
  SANY_POSRES_TIMEOUT_SECONDS Declare loss after this silence (default: 2)
  SANY_MIN_FREE_GB            Refuse to start below this free space (default: 20)
  SANY_RECORD_QOS_FILE        QoS override YAML
  SANY_COMPRESSED_LIDAR_TOPICS
                              Optional space-separated recorder topic override
  SANY_IMU_TOPIC              Optional primary IMU topic override

The run continues until localization exits or Ctrl-C. A /PosRes loss only
records a snapshot; it does not stop or restart localization.
EOF
}

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi
[[ "${run_name}" != */* ]] || fail "run_name must not contain '/'."
[[ -r "${ros_setup}" ]] || fail "ROS setup not found: ${ros_setup}"
[[ -r "${install_setup}" ]] || fail "workspace setup not found: ${install_setup}"
[[ -r "${config_path}" ]] || fail "config not found: ${config_path}"
[[ -r "${qos_file}" ]] || fail "QoS file not found: ${qos_file}"
[[ "${out_root}" == /* ]] || fail "LIGHTNING_LM_OUT_ROOT must be absolute."
[[ "${config_path}" == /* ]] || fail "LIGHTNING_LM_CONFIG must be absolute."
if [[ -n "${map_path}" && "${map_path}" != /* ]]; then
  fail "SANY_MAP_PATH must be absolute when set."
fi
if [[ -n "${map_path}" && ! -e "${map_path}" ]]; then
  fail "SANY_MAP_PATH does not exist: ${map_path}"
fi
[[ "${topic_wait_seconds}" =~ ^[1-9][0-9]*$ ]] || fail "SANY_TOPIC_WAIT_SECONDS must be positive."
[[ "${posres_timeout_seconds}" =~ ^[1-9][0-9]*$ ]] || fail "SANY_POSRES_TIMEOUT_SECONDS must be positive."
[[ "${min_free_gb}" =~ ^[1-9][0-9]*$ ]] || fail "SANY_MIN_FREE_GB must be positive."

set +u
source "${ros_setup}"
source "${install_setup}"
set -u

command -v ros2 >/dev/null 2>&1 || fail "ros2 is unavailable after sourcing the workspace."
command -v timeout >/dev/null 2>&1 || fail "timeout is unavailable."
command -v python3 >/dev/null 2>&1 || fail "python3 is unavailable."
grep -Fqx mcap <<<"$(ros2 bag list storage)" || fail "MCAP storage plugin is not installed."

sensor_topics="$({ python3 - "${config_path}" <<'PY'
import re
import sys

import yaml

root = yaml.safe_load(open(sys.argv[1], encoding="utf-8")) or {}
common = root.get("common") or {}
multi = root.get("multi_lidar") or {}
configured = multi.get("topics") or {}

if multi.get("enabled", False):
    lidar_topics = []
    for key, topic in configured.items():
        match = re.fullmatch(r"lidar_(\d+)", str(key))
        if match and topic:
            lidar_topics.append((int(match.group(1)), str(topic)))
    lidar_topics.sort()
    primary_id = int(multi.get("primary_lidar_id", 0))
    imu_topic = configured.get(f"imu_{primary_id}") or common.get("imu_topic")
else:
    topic = common.get("lidar_topic")
    lidar_topics = [(0, str(topic))] if topic else []
    imu_topic = common.get("imu_topic")

if not lidar_topics:
    raise SystemExit("no LiDAR topics found in the localization config")
if not imu_topic:
    raise SystemExit("primary IMU topic is missing from the localization config")

print(str(imu_topic))
for _, topic in lidar_topics:
    print(topic if topic.endswith("/zstd") else f"{topic}/zstd")
PY
} 2>&1)" || fail "failed to read sensor topics from ${config_path}: ${sensor_topics}"
mapfile -t configured_sensor_topics <<<"${sensor_topics}"
(( ${#configured_sensor_topics[@]} >= 2 )) || fail "config must provide at least one LiDAR and one IMU topic."

imu_topic="${SANY_IMU_TOPIC:-${configured_sensor_topics[0]}}"
if [[ -n "${SANY_COMPRESSED_LIDAR_TOPICS:-}" ]]; then
  read -r -a lidar_topics <<<"${SANY_COMPRESSED_LIDAR_TOPICS}"
else
  lidar_topics=("${configured_sensor_topics[@]:1}")
fi
(( ${#lidar_topics[@]} > 0 )) || fail "at least one compressed LiDAR topic is required."

record_topics=(
  "${lidar_topics[@]}"
  "${imu_topic}"
  /PosRes
  /slamPoseRaw_topic
  /localization/fault_status
  /localization/loc_status
  /localization/pipeline_diagnostics
  /tf
  /tf_static
  /rosout
)

mkdir -p "${out_root}"
out_root="$(cd -- "${out_root}" && pwd)"
run_dir="${out_root}/${run_name}"
[[ ! -e "${run_dir}" ]] || fail "run directory already exists: ${run_dir}"
available_kb="$(df -Pk "${out_root}" | awk 'NR==2 {print $4}')"
required_kb=$((min_free_gb * 1024 * 1024))
((available_kb >= required_kb)) || fail "less than ${min_free_gb} GiB free under ${out_root}"

mkdir -p "${run_dir}/bag" "${run_dir}/logs" "${run_dir}/snapshots" "${run_dir}/results"
config_path="$(realpath "${config_path}")"
cp -- "${config_path}" "${run_dir}/config.yaml"

wait_for_inputs() {
  local deadline=$((SECONDS + topic_wait_seconds))
  local listed topic missing
  while ((SECONDS < deadline)); do
    listed="$(ros2 topic list)"
    missing=0
    for topic in "${lidar_topics[@]}" "${imu_topic}"; do
      grep -Fqx "${topic}" <<<"${listed}" || missing=$((missing + 1))
    done
    ((missing == 0)) && return 0
    sleep 1
  done
  echo "Missing required topics:" >&2
  listed="$(ros2 topic list)"
  for topic in "${lidar_topics[@]}" "${imu_topic}"; do
    grep -Fqx "${topic}" <<<"${listed}" || echo "  ${topic}" >&2
  done
  return 1
}

snapshot_incident() {
  local index="$1"
  local snapshot="${run_dir}/snapshots/loss_$(printf '%03d' "${index}")_$(date +%Y%m%d_%H%M%S).txt"
  {
    echo "captured_at=$(date --iso-8601=ns)"
    echo "reason=PosRes silent for at least ${posres_timeout_seconds}s"
    echo
    echo "[pipeline_diagnostics]"
    timeout 3 ros2 topic echo --once /localization/pipeline_diagnostics 2>&1 || true
    echo
    echo "[fault_status]"
    timeout 3 ros2 topic echo --once /localization/fault_status 2>&1 || true
    echo
    echo "[nodes]"
    ros2 node list 2>&1 || true
    echo
    echo "[topics]"
    ros2 topic list -t 2>&1 || true
    echo
    echo "[input_topic_info]"
    for topic in "${lidar_topics[@]}" "${imu_topic}"; do
      echo "--- ${topic}"
      ros2 topic info --verbose "${topic}" 2>&1 || true
    done
    echo
    echo "[processes]"
    ps -eo pid,ppid,stat,psr,pcpu,pmem,rss,etimes,comm,args --sort=-pcpu 2>&1 || true
    echo
    echo "[memory]"
    free -h 2>&1 || true
    echo
    echo "[disk]"
    df -h "${out_root}" 2>&1 || true
    echo
    echo "[network]"
    ip -s link 2>&1 || true
    echo
    echo "[clock]"
    timedatectl status 2>&1 || true
    command -v chronyc >/dev/null 2>&1 && chronyc tracking 2>&1 || true
  } >"${snapshot}"
}

watch_posres() {
  local seen=false
  local lost=false
  local incident=0
  while true; do
    if timeout "${posres_timeout_seconds}" ros2 topic echo --once /PosRes >/dev/null 2>&1; then
      if [[ "${lost}" == true ]]; then
        echo "$(date --iso-8601=ns),RECOVERED,${incident}" >>"${run_dir}/logs/posres_watchdog.csv"
      elif [[ "${seen}" == false ]]; then
        echo "$(date --iso-8601=ns),FIRST_POSE,0" >>"${run_dir}/logs/posres_watchdog.csv"
      fi
      seen=true
      lost=false
    elif [[ "${seen}" == true && "${lost}" == false ]]; then
      incident=$((incident + 1))
      lost=true
      echo "$(date --iso-8601=ns),LOST,${incident}" >>"${run_dir}/logs/posres_watchdog.csv"
      snapshot_incident "${incident}"
    fi
  done
}

child_pids=()
recorder_pid=""
watchdog_pid=""
tegrastats_pid=""
algorithm_pid=""
stopping=false
stop_children() {
  [[ "${stopping}" == false ]] || return 0
  stopping=true
  trap - INT TERM
  local pid
  # SIGINT lets rosbag2 and localization flush MCAP/trajectories. The shell
  # watchdog and tegrastats do not own buffered artifacts, so TERM is safer
  # and avoids background shells ignoring SIGINT.
  [[ -z "${recorder_pid}" ]] || kill -INT "${recorder_pid}" 2>/dev/null || true
  [[ -z "${algorithm_pid}" ]] || kill -INT "${algorithm_pid}" 2>/dev/null || true
  [[ -z "${watchdog_pid}" ]] || kill -TERM "${watchdog_pid}" 2>/dev/null || true
  [[ -z "${tegrastats_pid}" ]] || kill -TERM "${tegrastats_pid}" 2>/dev/null || true
  for pid in "${child_pids[@]}"; do
    wait "${pid}" 2>/dev/null || true
  done
}
trap 'stop_children; exit 130' INT TERM
trap stop_children EXIT

echo "Waiting for ${#lidar_topics[@]} Zstd LiDAR topics and the primary IMU..."
wait_for_inputs || fail "required compressed sensor inputs did not appear."
for topic in "${lidar_topics[@]}"; do
  actual_type="$(ros2 topic type "${topic}")"
  [[ "${actual_type}" == "livox_ros_driver2/msg/CompressedPointCloud2" ]] ||
    fail "${topic} has type ${actual_type}, expected livox_ros_driver2/msg/CompressedPointCloud2"
done
actual_type="$(ros2 topic type "${imu_topic}")"
[[ "${actual_type}" == "sensor_msgs/msg/Imu" ]] ||
  fail "${imu_topic} has type ${actual_type}, expected sensor_msgs/msg/Imu"

{
  echo "started_at=$(date --iso-8601=ns)"
  echo "repo_dir=${repo_dir}"
  echo "config_path=${config_path}"
  echo "map_path=${map_path:-<from-config>}"
  echo "run_dir=${run_dir}"
  echo "compressed_lidar_topics=${lidar_topics[*]}"
  echo "primary_imu_topic=${imu_topic}"
  echo "git_commit=$(git -C "${repo_dir}" rev-parse HEAD 2>/dev/null || echo unavailable)"
  echo "config_sha256=$(sha256sum "${config_path}" | awk '{print $1}')"
  echo "kernel=$(uname -a)"
  echo "ros_distro=${ROS_DISTRO:-unknown}"
} >"${run_dir}/run_metadata.txt"
git -C "${repo_dir}" status --short >"${run_dir}/git_status.txt" 2>&1 || true
ros2 topic list -t >"${run_dir}/topics_at_start.txt" 2>&1 || true

ros2 bag record \
  --storage mcap \
  --storage-preset-profile fastwrite \
  --max-cache-size 1073741824 \
  --qos-profile-overrides-path "${qos_file}" \
  --output "${run_dir}/bag/localization_incident" \
  "${record_topics[@]}" \
  >"${run_dir}/logs/rosbag.stdout.log" \
  2>"${run_dir}/logs/rosbag.stderr.log" &
recorder_pid=$!
child_pids+=("${recorder_pid}")

watch_posres &
watchdog_pid=$!
child_pids+=("${watchdog_pid}")

if command -v tegrastats >/dev/null 2>&1; then
  tegrastats --interval 1000 >"${run_dir}/logs/tegrastats.log" 2>&1 &
  tegrastats_pid=$!
  child_pids+=("${tegrastats_pid}")
fi

algorithm_args=(
  ros2 run lightning run_loc_online
  --config="${config_path}"
  --output_tum="${run_dir}/results/trajectory_global.tum"
  --output_high_frequency_tum="${run_dir}/results/trajectory_high_frequency.tum"
)
if [[ -n "${map_path}" ]]; then
  algorithm_args+=(--map="${map_path}")
fi

echo "Recording to ${run_dir}; Ctrl-C stops the run cleanly."
set +e
stdbuf -oL -eL "${algorithm_args[@]}" \
  >"${run_dir}/logs/run_loc_online.stdout.log" \
  2>"${run_dir}/logs/run_loc_online.stderr.log" &
algorithm_pid=$!
child_pids+=("${algorithm_pid}")
wait "${algorithm_pid}"
algorithm_status=$?
set -e

echo "algorithm_exit_code=${algorithm_status}" >>"${run_dir}/run_metadata.txt"
echo "finished_at=$(date --iso-8601=ns)" >>"${run_dir}/run_metadata.txt"
stop_children
trap - EXIT
echo "Diagnostics saved to ${run_dir}"
exit "${algorithm_status}"
