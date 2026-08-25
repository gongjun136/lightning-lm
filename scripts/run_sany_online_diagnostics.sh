#!/usr/bin/env bash
# Run SANY multi-LiDAR localization with an optional raw-sensor MCAP
# recorder and optional /PosRes watchdog incident snapshots. Localization is
# never restarted by this script; the in-process global relocalizer owns recovery.
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="${LIGHTNING_LM_REPO_DIR:-$(cd -- "${script_dir}/.." && pwd)}"
ros_setup="${LIGHTNING_LM_ROS_SETUP:-/opt/ros/humble/setup.bash}"
install_setup="${LIGHTNING_LM_INSTALL_SETUP:-${repo_dir}/install/setup.bash}"
lidar_layout="${SANY_LIDAR_LAYOUT:-4}"
case "${lidar_layout}" in
  3)
    default_config_path="${repo_dir}/config/reproduction/multi_lidar/sany_3livox/sany_3lidar_localization_solid.yaml"
    ;;
  4)
    default_config_path="${repo_dir}/config/reproduction/multi_lidar/sany_4livox/sany_4lidar_localization_solid.yaml"
    ;;
  *)
    echo "ERROR: SANY_LIDAR_LAYOUT must be 3 or 4." >&2
    exit 1
    ;;
esac
config_path="${LIGHTNING_LM_CONFIG:-${default_config_path}}"
map_path="${SANY_MAP_PATH:-}"
out_root="${LIGHTNING_LM_OUT_ROOT:-/home/nvidia/project/gj_ws/runs}"
run_name="${1:-sany_loc_diag_$(date +%Y%m%d_%H%M%S)}"
qos_file="${SANY_RECORD_QOS_FILE:-${script_dir}/config/sany_localization_record_qos.yaml}"
record_bag="${SANY_RECORD_BAG:-1}"
topic_wait_seconds="${SANY_TOPIC_WAIT_SECONDS:-60}"
posres_timeout_seconds="${SANY_POSRES_TIMEOUT_SECONDS:-2}"
enable_posres_watchdog="${SANY_ENABLE_POSRES_WATCHDOG:-0}"
min_free_gb="${SANY_MIN_FREE_GB:-20}"
wheel_speed_topic="${SANY_WHEEL_SPEED_TOPIC:-/SpeThrCAN4_topic}"
enable_can_observation="${SANY_ENABLE_CAN_OBSERVATION:-1}"
run_mode="${LIGHTNING_LM_RUN_MODE:-diagnostic}"
case "${run_mode}" in
  diagnostic)
    default_compute_profile=1
    default_reduce_nonessential_overhead=0
    ;;
  production)
    default_compute_profile=0
    default_reduce_nonessential_overhead=1
    ;;
  *)
    echo "ERROR: LIGHTNING_LM_RUN_MODE must be diagnostic or production." >&2
    exit 1
    ;;
esac
compute_profile="${LIGHTNING_LM_COMPUTE_PROFILE:-${default_compute_profile}}"
reduce_nonessential_overhead="${LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD:-${default_reduce_nonessential_overhead}}"
cpu_affinity="${LIGHTNING_LM_CPU_AFFINITY:-}"

usage() {
  cat <<'EOF'
Usage:
  LIGHTNING_LM_CONFIG=/absolute/config.yaml \
  SANY_MAP_PATH=/absolute/map/path \
  scripts/run_sany_online_diagnostics.sh [run_name]

Required before launch:
  1. Start every Livox driver configured by LIGHTNING_LM_CONFIG.
     No-bag mode requires the primary IMU; bag mode requires all configured IMUs.
  2. When CAN observation is enabled, start the CAN bridge that publishes
     /SpeThrCAN4_topic.

Useful environment variables:
  SANY_LIDAR_LAYOUT           Select bundled localization YAML: 3 or 4 (default: 4)
  LIGHTNING_LM_INSTALL_SETUP  Built workspace setup.bash
  LIGHTNING_LM_OUT_ROOT       Run root (default: ~/project/gj_ws/runs)
  SANY_RECORD_BAG             Record a background MCAP: 1=yes, 0=no (default: 1)
  SANY_TOPIC_WAIT_SECONDS     Sensor-input discovery timeout (default: 60)
  SANY_ENABLE_POSRES_WATCHDOG Watch /PosRes and capture loss snapshots: 1=yes, 0=no (default: 0)
  SANY_POSRES_TIMEOUT_SECONDS Declare loss after this silence (default: 2)
  SANY_MIN_FREE_GB            Refuse to start below this free space (default: 20)
  SANY_RECORD_QOS_FILE        QoS override YAML
  SANY_ENABLE_CAN_OBSERVATION Fuse CAN wheel speed: 1=yes, 0=no (default: 1)
  SANY_WHEEL_SPEED_TOPIC      Motor-speed topic (default: /SpeThrCAN4_topic)
  SANY_IMU_TOPIC              Optional primary IMU topic override
  LIGHTNING_LM_RUN_MODE       diagnostic or production (this entry defaults to diagnostic;
                              scripts/run.sh fixes the normal field preset)
  LIGHTNING_LM_COMPUTE_PROFILE
                              Emit machine-readable compute timing: 1=yes, 0=no
                              (mode default: diagnostic=1, production=0)
  LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD
                              Skip high-rate diagnostic I/O: 1=yes, 0=no
                              (mode default: diagnostic=0, production=1)
  LIGHTNING_LM_LIO_THREADS    Override compute_budget.lio_threads
  LIGHTNING_LM_NDT_THREADS    Override compute_budget.ndt_threads
  LIGHTNING_LM_SOLID_ICP_WORKERS
                              Override compute_budget.solid_icp_workers
  LIGHTNING_LM_SOLID_WORKER_NICE
                              Override SOLiD worker nice floor: 0=unchanged, 1-19=yield
  LIGHTNING_LM_SOLID_CPU_AFFINITY
                              Optional per-SOLiD-worker CPU list, e.g. 1-5,7;
                              must be a subset of the process-wide CPU affinity
  LIGHTNING_LM_CPU_AFFINITY   Optional process-wide taskset CPU list, e.g. 0-9;
                              it does not isolate LIO/NDT/SOLiD within the process

The run continues until localization exits or Ctrl-C. When the watchdog is
enabled, a /PosRes loss only records a snapshot; it does not stop or restart
localization.
When recording, every raw PointCloud2, every configured IMU, and the CAN topic
are added to the bag even when CAN observation is disabled.
EOF
}

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

validate_optional_thread_count() {
  local name="$1"
  local value="${!name:-}"
  [[ -z "${value}" ]] && return 0
  [[ "${value}" =~ ^[1-9][0-9]*$ ]] || fail "${name} must be an integer in [1, 128]."
  ((value <= 128)) || fail "${name} must be an integer in [1, 128]."
}

validate_optional_nice() {
  local name="$1"
  local value="${!name:-}"
  [[ -z "${value}" ]] && return 0
  [[ "${value}" =~ ^[0-9]+$ ]] || fail "${name} must be an integer in [0, 19]."
  ((value <= 19)) || fail "${name} must be an integer in [0, 19]."
}

cpu_list_is_subset() {
  local subset="$1"
  local superset="$2"
  python3 - "${subset}" "${superset}" <<'PY'
import sys


def expand_cpu_list(text):
    cpus = set()
    for item in text.split(","):
        bounds = [part.strip() for part in item.strip().split("-", 1)]
        if not bounds[0] or len(bounds) > 2 or (len(bounds) == 2 and not bounds[1]):
            raise ValueError("invalid CPU list")
        first = int(bounds[0])
        last = int(bounds[-1])
        if first < 0 or last < first:
            raise ValueError("invalid CPU range")
        cpus.update(range(first, last + 1))
    return cpus


try:
    requested = expand_cpu_list(sys.argv[1])
    allowed = expand_cpu_list(sys.argv[2])
except ValueError:
    raise SystemExit(2)
raise SystemExit(0 if requested.issubset(allowed) else 1)
PY
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi
[[ "${run_name}" != */* ]] || fail "run_name must not contain '/'."
[[ -r "${ros_setup}" ]] || fail "ROS setup not found: ${ros_setup}"
[[ -r "${install_setup}" ]] || fail "workspace setup not found: ${install_setup}"
[[ -r "${config_path}" ]] || fail "config not found: ${config_path}"
[[ "${record_bag}" == "0" || "${record_bag}" == "1" ]] || fail "SANY_RECORD_BAG must be 0 or 1."
[[ "${enable_posres_watchdog}" == "0" || "${enable_posres_watchdog}" == "1" ]] ||
  fail "SANY_ENABLE_POSRES_WATCHDOG must be 0 or 1."
[[ "${enable_can_observation}" == "0" || "${enable_can_observation}" == "1" ]] ||
  fail "SANY_ENABLE_CAN_OBSERVATION must be 0 or 1."
[[ "${compute_profile}" == "0" || "${compute_profile}" == "1" ]] ||
  fail "LIGHTNING_LM_COMPUTE_PROFILE must be 0 or 1."
[[ "${reduce_nonessential_overhead}" == "0" || "${reduce_nonessential_overhead}" == "1" ]] ||
  fail "LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD must be 0 or 1."
validate_optional_thread_count LIGHTNING_LM_LIO_THREADS
validate_optional_thread_count LIGHTNING_LM_NDT_THREADS
validate_optional_thread_count LIGHTNING_LM_SOLID_ICP_WORKERS
validate_optional_nice LIGHTNING_LM_SOLID_WORKER_NICE
export SANY_ENABLE_CAN_OBSERVATION="${enable_can_observation}"
export LIGHTNING_LM_COMPUTE_PROFILE="${compute_profile}"
export LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD="${reduce_nonessential_overhead}"
if [[ "${record_bag}" == "1" ]]; then
  [[ -r "${qos_file}" ]] || fail "QoS file not found: ${qos_file}"
fi
[[ "${wheel_speed_topic}" == /* ]] || fail "SANY_WHEEL_SPEED_TOPIC must be an absolute topic name."
[[ "${out_root}" == /* ]] || fail "LIGHTNING_LM_OUT_ROOT must be absolute."
[[ "${config_path}" == /* ]] || fail "LIGHTNING_LM_CONFIG must be absolute."
if [[ -n "${map_path}" && "${map_path}" != /* ]]; then
  fail "SANY_MAP_PATH must be absolute when set."
fi
if [[ -n "${map_path}" && ! -e "${map_path}" ]]; then
  fail "SANY_MAP_PATH does not exist: ${map_path}"
fi
[[ "${topic_wait_seconds}" =~ ^[1-9][0-9]*$ ]] || fail "SANY_TOPIC_WAIT_SECONDS must be positive."
if [[ "${enable_posres_watchdog}" == "1" ]]; then
  [[ "${posres_timeout_seconds}" =~ ^[1-9][0-9]*$ ]] ||
    fail "SANY_POSRES_TIMEOUT_SECONDS must be positive."
fi
if [[ "${record_bag}" == "1" ]]; then
  [[ "${min_free_gb}" =~ ^[1-9][0-9]*$ ]] || fail "SANY_MIN_FREE_GB must be positive."
fi

set +u
source "${ros_setup}"
source "${install_setup}"
set -u

command -v ros2 >/dev/null 2>&1 || fail "ros2 is unavailable after sourcing the workspace."
command -v timeout >/dev/null 2>&1 || fail "timeout is unavailable."
command -v python3 >/dev/null 2>&1 || fail "python3 is unavailable."
command -v setsid >/dev/null 2>&1 || fail "setsid is unavailable."
if [[ -n "${cpu_affinity}" ]]; then
  command -v taskset >/dev/null 2>&1 || fail "taskset is required by LIGHTNING_LM_CPU_AFFINITY."
  taskset --cpu-list "${cpu_affinity}" true >/dev/null 2>&1 ||
    fail "LIGHTNING_LM_CPU_AFFINITY is invalid or unavailable in this cpuset: ${cpu_affinity}"
fi
if [[ -n "${LIGHTNING_LM_SOLID_CPU_AFFINITY:-}" ]]; then
  command -v taskset >/dev/null 2>&1 ||
    fail "taskset is required to validate LIGHTNING_LM_SOLID_CPU_AFFINITY."
  taskset --cpu-list "${LIGHTNING_LM_SOLID_CPU_AFFINITY}" true >/dev/null 2>&1 ||
    fail "LIGHTNING_LM_SOLID_CPU_AFFINITY is invalid or unavailable in this cpuset: ${LIGHTNING_LM_SOLID_CPU_AFFINITY}"
  effective_process_affinity="${cpu_affinity}"
  if [[ -z "${effective_process_affinity}" ]]; then
    effective_process_affinity="$(taskset --cpu-list --pid $$ 2>/dev/null)" ||
      fail "unable to read the launcher's effective CPU affinity."
    effective_process_affinity="${effective_process_affinity##*: }"
  fi
  cpu_list_is_subset "${LIGHTNING_LM_SOLID_CPU_AFFINITY}" "${effective_process_affinity}" ||
    fail "LIGHTNING_LM_SOLID_CPU_AFFINITY (${LIGHTNING_LM_SOLID_CPU_AFFINITY}) must be a subset of the effective process CPU affinity (${effective_process_affinity})."
fi
if [[ "${record_bag}" == "1" ]]; then
  grep -Fqx mcap <<<"$(ros2 bag list storage)" || fail "MCAP storage plugin is not installed."
fi

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
    imu_topics = []
    for key, topic in configured.items():
        match = re.fullmatch(r"lidar_(\d+)", str(key))
        if match and topic:
            lidar_topics.append((int(match.group(1)), str(topic)))
        match = re.fullmatch(r"imu_(\d+)", str(key))
        if match and topic:
            imu_topics.append((int(match.group(1)), str(topic)))
    lidar_topics.sort()
    imu_topics.sort()
    primary_id = int(multi.get("primary_lidar_id", 0))
    imu_topic = configured.get(f"imu_{primary_id}") or common.get("imu_topic")
else:
    topic = common.get("lidar_topic")
    lidar_topics = [(0, str(topic))] if topic else []
    imu_topic = common.get("imu_topic")
    imu_topics = [(0, str(imu_topic))] if imu_topic else []

if not lidar_topics:
    raise SystemExit("no LiDAR topics found in the localization config")
if not imu_topic:
    raise SystemExit("primary IMU topic is missing from the localization config")

print(str(imu_topic))
print(" ".join(topic for _, topic in imu_topics))
for _, topic in lidar_topics:
    print(topic)
PY
} 2>&1)" || fail "failed to read sensor topics from ${config_path}: ${sensor_topics}"
mapfile -t configured_sensor_topics <<<"${sensor_topics}"
(( ${#configured_sensor_topics[@]} >= 3 )) || fail "config must provide at least one LiDAR and one IMU topic."

imu_topic="${SANY_IMU_TOPIC:-${configured_sensor_topics[0]}}"
read -r -a configured_imu_topics <<<"${configured_sensor_topics[1]}"
lidar_topics=("${configured_sensor_topics[@]:2}")
record_imu_topics=()
required_imu_topics=("${imu_topic}")
record_topics=()
if [[ "${record_bag}" == "1" ]]; then
  record_imu_topics=("${configured_imu_topics[@]}")
  if [[ " ${record_imu_topics[*]} " != *" ${imu_topic} "* ]]; then
    record_imu_topics+=("${imu_topic}")
  fi
  required_imu_topics=("${record_imu_topics[@]}")
  record_topics=(
    "${lidar_topics[@]}"
    "${record_imu_topics[@]}"
    "${wheel_speed_topic}"
    /PosRes
    /localization/pose_vel
    /slamPoseRaw_topic
    /localization/fault_status
    /localization/loc_status
    /localization/pipeline_diagnostics
    /tf
    /tf_static
    /rosout
  )
fi
required_topics=("${lidar_topics[@]}" "${required_imu_topics[@]}")
if [[ "${enable_can_observation}" == "1" ]]; then
  required_topics+=("${wheel_speed_topic}")
fi

mkdir -p "${out_root}"
out_root="$(cd -- "${out_root}" && pwd)"
run_dir="${out_root}/${run_name}"
[[ ! -e "${run_dir}" ]] || fail "run directory already exists: ${run_dir}"
if [[ "${record_bag}" == "1" ]]; then
  available_kb="$(df -Pk "${out_root}" | awk 'NR==2 {print $4}')"
  required_kb=$((min_free_gb * 1024 * 1024))
  ((available_kb >= required_kb)) || fail "less than ${min_free_gb} GiB free under ${out_root}"
fi

mkdir -p "${run_dir}/logs" "${run_dir}/snapshots" "${run_dir}/results"
if [[ "${record_bag}" == "1" ]]; then
  mkdir -p "${run_dir}/bag"
fi
config_path="$(realpath "${config_path}")"
cp -- "${config_path}" "${run_dir}/config.yaml"

wait_for_inputs() {
  local deadline=$((SECONDS + topic_wait_seconds))
  local listed topic missing
  while ((SECONDS < deadline)); do
    listed="$(ros2 topic list)"
    missing=0
    for topic in "${required_topics[@]}"; do
      grep -Fqx "${topic}" <<<"${listed}" || missing=$((missing + 1))
    done
    ((missing == 0)) && return 0
    sleep 1
  done
  echo "Missing required topics:" >&2
  listed="$(ros2 topic list)"
  for topic in "${required_topics[@]}"; do
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
    timeout 3 ros2 topic echo --once --qos-reliability best_effort \
      /localization/pipeline_diagnostics 2>&1 || true
    echo
    echo "[fault_status]"
    timeout 3 ros2 topic echo --once --qos-reliability best_effort \
      /localization/fault_status 2>&1 || true
    echo
    echo "[nodes]"
    ros2 node list 2>&1 || true
    echo
    echo "[topics]"
    ros2 topic list -t 2>&1 || true
    echo
    echo "[input_topic_info]"
    for topic in "${required_topics[@]}"; do
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
  local incident=0
  local event silence_sec
  local monitor_pid=""
  local event_fifo="${run_dir}/logs/posres_watchdog.events"
  cleanup_posres_monitor() {
    if [[ -n "${monitor_pid}" ]]; then
      kill -TERM "${monitor_pid}" 2>/dev/null || true
      wait "${monitor_pid}" 2>/dev/null || true
      monitor_pid=""
    fi
    rm -f -- "${event_fifo}"
  }
  trap 'cleanup_posres_monitor; exit 0' INT TERM
  trap cleanup_posres_monitor EXIT

  rm -f -- "${event_fifo}"
  mkfifo "${event_fifo}"
  python3 -u "${script_dir}/monitor_posres_silence.py" \
    --timeout-sec "${posres_timeout_seconds}" \
    >"${event_fifo}" \
    2>"${run_dir}/logs/posres_watchdog.stderr.log" &
  monitor_pid=$!
  while IFS=, read -r event silence_sec; do
    case "${event}" in
      FIRST_POSE)
        echo "$(date --iso-8601=ns),FIRST_POSE,0,${silence_sec}" \
          >>"${run_dir}/logs/posres_watchdog.csv"
        ;;
      LOST)
        incident=$((incident + 1))
        echo "$(date --iso-8601=ns),LOST,${incident},${silence_sec}" \
          >>"${run_dir}/logs/posres_watchdog.csv"
        snapshot_incident "${incident}"
        ;;
      RECOVERED)
        echo "$(date --iso-8601=ns),RECOVERED,${incident},${silence_sec}" \
          >>"${run_dir}/logs/posres_watchdog.csv"
        ;;
    esac
  done <"${event_fifo}"
  set +e
  wait "${monitor_pid}"
  local monitor_status=$?
  set -e
  monitor_pid=""
  if ((monitor_status != 0)); then
    echo "$(date --iso-8601=ns),MONITOR_EXIT,${monitor_status},0" \
      >>"${run_dir}/logs/posres_watchdog.csv"
  fi
}

child_pids=()
recorder_pid=""
watchdog_pid=""
tegrastats_pid=""
algorithm_pid=""
stopping=false
profile_extracted=false
extract_compute_profile() {
  [[ "${profile_extracted}" == false ]] || return 0
  local algorithm_log="${run_dir}/logs/run_loc_online.stderr.log"
  [[ -f "${algorithm_log}" ]] || return 0
  grep -F "COMPUTE_BENCH_" "${algorithm_log}" \
    >"${run_dir}/results/compute_profile.log" || true
  local profile_record_count
  profile_record_count="$(wc -l <"${run_dir}/results/compute_profile.log")"
  if [[ -f "${run_dir}/run_metadata.txt" ]]; then
    echo "compute_profile_record_count=${profile_record_count}" \
      >>"${run_dir}/run_metadata.txt"
  fi
  profile_extracted=true
}
stop_children() {
  [[ "${stopping}" == false ]] || return 0
  stopping=true
  trap - INT TERM
  local pid
  # SIGINT lets rosbag2 and localization flush MCAP/trajectories. The shell
  # starts them in dedicated process groups so the Python ros2 launcher and
  # its child both receive the signal. The watchdog and tegrastats do not own
  # buffered artifacts, so TERM is safer for those helpers.
  [[ -z "${recorder_pid}" ]] || kill -INT -- "-${recorder_pid}" 2>/dev/null || true
  [[ -z "${algorithm_pid}" ]] || kill -INT -- "-${algorithm_pid}" 2>/dev/null || true
  [[ -z "${watchdog_pid}" ]] || kill -TERM "${watchdog_pid}" 2>/dev/null || true
  [[ -z "${tegrastats_pid}" ]] || kill -TERM "${tegrastats_pid}" 2>/dev/null || true
  for pid in "${child_pids[@]}"; do
    wait "${pid}" 2>/dev/null || true
  done
  extract_compute_profile
}
trap 'stop_children; exit 130' INT TERM
trap stop_children EXIT

if [[ "${record_bag}" == "1" ]]; then
  if [[ "${enable_can_observation}" == "1" ]]; then
    echo "Waiting for ${#lidar_topics[@]} raw LiDAR topics, ${#record_imu_topics[@]} IMU topics, and CAN wheel speed..."
  else
    echo "Waiting for ${#lidar_topics[@]} raw LiDAR topics and ${#record_imu_topics[@]} IMU topics; CAN observation is disabled."
  fi
else
  if [[ "${enable_can_observation}" == "1" ]]; then
    echo "Waiting for ${#lidar_topics[@]} raw LiDAR topics, the primary IMU, and CAN wheel speed..."
  else
    echo "Waiting for ${#lidar_topics[@]} raw LiDAR topics and the primary IMU; CAN observation is disabled."
  fi
fi
wait_for_inputs || fail "required sensor inputs did not appear."
for topic in "${lidar_topics[@]}"; do
  actual_type="$(ros2 topic type "${topic}")"
  [[ "${actual_type}" == "sensor_msgs/msg/PointCloud2" ]] ||
    fail "${topic} has type ${actual_type}, expected sensor_msgs/msg/PointCloud2"
done
for topic in "${required_imu_topics[@]}"; do
  actual_type="$(ros2 topic type "${topic}")"
  [[ "${actual_type}" == "sensor_msgs/msg/Imu" ]] ||
    fail "${topic} has type ${actual_type}, expected sensor_msgs/msg/Imu"
done
if [[ "${enable_can_observation}" == "1" ]]; then
  wheel_speed_type="$(ros2 topic type "${wheel_speed_topic}")"
  [[ "${wheel_speed_type}" == "geosun_msgs/msg/SpeThrCAN4" ]] ||
    fail "${wheel_speed_topic} has type ${wheel_speed_type}, expected geosun_msgs/msg/SpeThrCAN4"
fi

{
  echo "started_at=$(date --iso-8601=ns)"
  echo "repo_dir=${repo_dir}"
  echo "config_path=${config_path}"
  echo "map_path=${map_path:-<from-config>}"
  echo "run_dir=${run_dir}"
  echo "record_bag=${record_bag}"
  echo "enable_posres_watchdog=${enable_posres_watchdog}"
  echo "raw_lidar_topics=${lidar_topics[*]}"
  echo "enable_can_observation=${enable_can_observation}"
  echo "run_mode=${run_mode}"
  echo "compute_profile=${compute_profile}"
  echo "reduce_nonessential_overhead=${reduce_nonessential_overhead}"
  echo "lio_threads_override=${LIGHTNING_LM_LIO_THREADS:-<from-config>}"
  echo "ndt_threads_override=${LIGHTNING_LM_NDT_THREADS:-<from-config>}"
  echo "solid_icp_workers_override=${LIGHTNING_LM_SOLID_ICP_WORKERS:-<from-config>}"
  echo "solid_worker_nice_override=${LIGHTNING_LM_SOLID_WORKER_NICE:-<from-config>}"
  echo "solid_cpu_affinity_override=${LIGHTNING_LM_SOLID_CPU_AFFINITY:-<from-config>}"
  echo "cpu_affinity=${cpu_affinity:-<unrestricted>}"
  echo "wheel_speed_topic=${wheel_speed_topic}"
  echo "primary_imu_topic=${imu_topic}"
  echo "recorded_imu_topics=${record_imu_topics[*]:-<disabled>}"
  echo "git_commit=$(git -C "${repo_dir}" rev-parse HEAD 2>/dev/null || echo unavailable)"
  echo "config_sha256=$(sha256sum "${config_path}" | awk '{print $1}')"
  echo "kernel=$(uname -a)"
  echo "ros_distro=${ROS_DISTRO:-unknown}"
} >"${run_dir}/run_metadata.txt"
git -C "${repo_dir}" status --short >"${run_dir}/git_status.txt" 2>&1 || true
ros2 topic list -t >"${run_dir}/topics_at_start.txt" 2>&1 || true

if [[ "${record_bag}" == "1" ]]; then
  setsid ros2 bag record \
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
  sleep 1
  if ! kill -0 "${recorder_pid}" 2>/dev/null; then
    set +e
    wait "${recorder_pid}"
    recorder_status=$?
    set -e
    tail -n 40 "${run_dir}/logs/rosbag.stderr.log" >&2 || true
    fail "rosbag recorder exited during startup with status ${recorder_status}."
  fi
fi

if [[ "${enable_posres_watchdog}" == "1" ]]; then
  echo "wall_time,event,incident,silence_sec" >"${run_dir}/logs/posres_watchdog.csv"
  python3 -c 'import rclpy; from geosun_msgs.msg import PosRes'
  python3 "${script_dir}/monitor_posres_silence.py" --help >/dev/null
  watch_posres &
  watchdog_pid=$!
  child_pids+=("${watchdog_pid}")
fi

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
  --output_published_tum="${run_dir}/results/trajectory_published_rear_axle.tum"
)
if [[ -n "${map_path}" ]]; then
  algorithm_args+=(--map="${map_path}")
fi

if [[ "${record_bag}" == "1" ]]; then
  echo "Running ${run_mode} mode with background bag recording in ${run_dir}; Ctrl-C stops the run cleanly."
else
  echo "Running ${run_mode} mode without bag recording in ${run_dir}; Ctrl-C stops the run cleanly."
fi
algorithm_launcher=()
if [[ -n "${cpu_affinity}" ]]; then
  algorithm_launcher+=(taskset --cpu-list "${cpu_affinity}")
fi
algorithm_launcher+=(stdbuf -oL -eL)
set +e
setsid "${algorithm_launcher[@]}" "${algorithm_args[@]}" \
  >"${run_dir}/logs/run_loc_online.stdout.log" \
  2>"${run_dir}/logs/run_loc_online.stderr.log" &
algorithm_pid=$!
child_pids+=("${algorithm_pid}")
wait "${algorithm_pid}"
algorithm_status=$?
set -e

echo "algorithm_exit_code=${algorithm_status}" >>"${run_dir}/run_metadata.txt"
echo "finished_at=$(date --iso-8601=ns)" >>"${run_dir}/run_metadata.txt"
extract_compute_profile
stop_children
trap - EXIT
echo "${run_mode^} run artifacts saved to ${run_dir}"
exit "${algorithm_status}"
