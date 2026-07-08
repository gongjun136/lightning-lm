#!/usr/bin/env bash
# 通用离线 SLAM 运行模板，对应 src/app/run_slam_offline.cc。
#
# 适用场景：
#   完整离线建图流程，包括前端 LIO、关键帧、后端/地图保存；用于正式复现实验和地图产物导出。
#
# 必需环境变量：
#   LIGHTNING_LM_INPUT_BAG=/path/to/rosbag2
#       输入 ROS2 bag 目录，例如 .../data/Outdoor01_mid360_ros2。
#
# 常用可选环境变量：
#   LIGHTNING_LM_CONFIG=/path/to/config.yaml
#       配置文件；默认使用仓库 config/default.yaml。
#   LIGHTNING_LM_OUT_ROOT=/path/to/output_root
#       实验输出根目录；默认使用仓库 runs/。
#   LIGHTNING_LM_WAIT_UI=true|false
#       离线处理完是否等待 Pangolin UI 关闭；需要查看地图时设为 true。
#   LIGHTNING_LM_OUTPUT_TUM=/path/to/output.tum
#       前端逐帧 LIO TUM 轨迹；默认写入 results/<run_name>.tum。
#       程序还会额外导出同名前缀的 _keyframes_lio.tum 和 _keyframes_opt.tum。
#   LIGHTNING_LM_REPO_DIR=/path/to/lightning-lm
#       从数据集目录调用脚本时可显式指定仓库路径。
#
# 位置参数：
#   第一个非 -- 开头参数作为 run_name；省略时自动生成带时间戳的名称。
#   其余参数会原样追加给 run_slam_offline。
#
# 输出组织：
#   $LIGHTNING_LM_OUT_ROOT/<run_name>/run_metadata.txt
#   $LIGHTNING_LM_OUT_ROOT/<run_name>/logs/run_slam_offline.{stdout,stderr}.log
#   $LIGHTNING_LM_OUT_ROOT/<run_name>/results/<run_name>.tum
#   $LIGHTNING_LM_OUT_ROOT/<run_name>/results/<run_name>_keyframes_lio.tum
#   $LIGHTNING_LM_OUT_ROOT/<run_name>/results/<run_name>_keyframes_opt.tum
#   $LIGHTNING_LM_OUT_ROOT/<run_name>/data/new_map/global.pcd
#   程序会在 $LIGHTNING_LM_OUT_ROOT/<run_name>/ 内执行，便于归档地图和相对路径产物。
#
# 示例：
#   LIGHTNING_LM_INPUT_BAG=/path/to/rosbag2 \
#   LIGHTNING_LM_CONFIG=/path/to/config.yaml \
#   LIGHTNING_LM_OUT_ROOT=/path/to/repro/runs \
#   LIGHTNING_LM_WAIT_UI=false \
#   scripts/run_slam_offline.sh outdoor01_slam_test
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="${LIGHTNING_LM_REPO_DIR:-$(cd "$script_dir/.." && pwd)}"
ros_setup="${LIGHTNING_LM_ROS_SETUP:-/opt/ros/humble/setup.bash}"
install_setup="${LIGHTNING_LM_INSTALL_SETUP:-$repo_dir/install/setup.bash}"

set +u
source "$ros_setup"
source "$install_setup"
set -u

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  exe_path="$(ros2 pkg prefix lightning)/lib/lightning/run_slam_offline"
  set +e
  "$exe_path" --help
  help_status=$?
  set -e
  [[ "$help_status" -eq 0 || "$help_status" -eq 1 ]] && exit 0
  exit "$help_status"
fi

bag_dir="${LIGHTNING_LM_INPUT_BAG:-}"
config_path="${LIGHTNING_LM_CONFIG:-$repo_dir/config/default.yaml}"
out_root="${LIGHTNING_LM_OUT_ROOT:-$repo_dir/runs}"
run_name="${LIGHTNING_LM_RUN_NAME:-run_slam_offline_$(date +%Y%m%d_%H%M%S)}"
wait_ui="${LIGHTNING_LM_WAIT_UI:-true}"
if [[ $# -gt 0 && "${1:0:1}" != "-" ]]; then
  run_name="$1"
  shift
fi

if [[ -z "$bag_dir" ]]; then
  echo "LIGHTNING_LM_INPUT_BAG is required" >&2
  exit 2
fi
if [[ ! -e "$bag_dir" ]]; then
  echo "input bag not found: $bag_dir" >&2
  exit 2
fi
if [[ ! -f "$config_path" ]]; then
  echo "config not found: $config_path" >&2
  exit 2
fi

bag_dir="$(realpath "$bag_dir")"
config_path="$(realpath "$config_path")"
out_root="$(mkdir -p "$out_root" && cd "$out_root" && pwd)"
run_dir="$out_root/$run_name"
result_dir="$run_dir/results"
log_dir="$run_dir/logs"
mkdir -p "$result_dir" "$log_dir"
if [[ -n "${LIGHTNING_LM_OUTPUT_TUM:-}" ]]; then
  tum_path="$LIGHTNING_LM_OUTPUT_TUM"
  [[ "${tum_path:0:1}" != "/" ]] && tum_path="$run_dir/$tum_path"
else
  tum_path="$result_dir/${run_name}.tum"
fi

{
  echo "executable=run_slam_offline"
  echo "repo_dir=$repo_dir"
  echo "bag_dir=$bag_dir"
  echo "config_path=$config_path"
  echo "run_name=$run_name"
  echo "run_dir=$run_dir"
  echo "tum_path=$tum_path"
  echo "wait_ui=$wait_ui"
  echo "extra_args=$*"
  date --iso-8601=seconds
} > "$run_dir/run_metadata.txt"

(
  cd "$run_dir"
  ros2 run lightning run_slam_offline \
    --input_bag "$bag_dir" \
    --config "$config_path" \
    --output_tum "$tum_path" \
    --wait_ui="$wait_ui" \
    "$@"
) > "$log_dir/run_slam_offline.stdout.log" \
  2> "$log_dir/run_slam_offline.stderr.log"
