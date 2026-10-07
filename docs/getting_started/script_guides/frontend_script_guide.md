@page frontend_script_guide 前端脚本导读

这组三个入口都围绕 LIO 前端。在线脚本接 ROS 话题，单次离线脚本交给 C++ 读取 bag，批量脚本再编排多个单次任务。

| 入口 | 先看哪里 | 输入由谁送入算法 |
|---|---|---|
| `run_frontend_online.sh` | @ref frontend_online_guide "在线前端" | 驱动或外部 `ros2 bag play` 发布 ROS 话题 |
| `run_frontend_offline.sh` | @ref frontend_offline_guide "单次离线前端" | C++ 离线程序读取指定 bag |
| `run_frontend_offline_batch.sh` | @ref frontend_batch_guide "批量离线前端" | 逐个调用单次离线入口 |

@anchor frontend_online_guide
# run_frontend_online.sh：配置交给实时前端

## 先看核心调用

@snippet{lineno} run_frontend_online.sh frontend-online-launch

`LIGHTNING_LM_CONFIG` 指定必需的外部 YAML，脚本将它转为绝对路径；`"$@"` 透传其余程序 flags。`exec` 用 ROS 启动命令替换当前 Shell，没有后台等待循环。脚本不改变 cwd，也不建立统一日志目录。

```text
外部 YAML → 加载 ROS/overlay → C++ run_frontend_online
驱动或外部回放 → ROS 订阅 → LaserMapping → 前端输出
```

## 按执行顺序阅读

| 阶段 | 目的 | 阅读入口 |
|---|---|---|
| 1. 确定路径和配置 | 找到仓库、安装环境与真实 YAML | @ref frontend_online_settings "路径与配置" |
| 2. 加载环境 | 让 ROS 找到程序和消息包 | @ref frontend_online_environment "环境" |
| 3. 进入程序 | 启动实时前端，等待 ROS 输入 | 上面的核心调用 |

@anchor frontend_online_settings
### 1. 路径与配置

@snippet{lineno} run_frontend_online.sh frontend-online-settings

这里的配置缺省为空；未设置或文件不存在就退出。拿到 `config_path` 后才加载运行环境。

@anchor frontend_online_environment
### 2. 加载环境

@snippet{lineno} run_frontend_online.sh frontend-online-environment

两次 `source` 后，核心调用启动安装空间中的可执行文件。接着读 @ref run_frontend_online.cc "run_frontend_online.cc 的 main()"：构造 `FrontendNode`，初始化 `LaserMapping`，进入 `rclcpp::spin`。算法主线见 @ref guide_lio "传感器与 LIO"。

<details>
<summary>运行设置与输出</summary>

`LIGHTNING_LM_REPO_DIR`、`LIGHTNING_LM_ROS_SETUP` 和 `LIGHTNING_LM_INSTALL_SETUP` 可覆盖仓库和环境路径。驱动/回放需另行启动，并与前端处于同一 ROS domain。输出 topic、UI 与诊断由配置和程序 flags 决定。

完整源码：@ref run_frontend_online.sh "run_frontend_online.sh"；精确边界：@ref script_contracts "脚本契约"。

</details>

@anchor frontend_offline_guide
# run_frontend_offline.sh：处理一次离线数据并归档

## 先看核心调用

@snippet{lineno} run_frontend_offline.sh frontend-offline-launch

`binary` 是本仓库安装空间中的 `run_frontend_offline`。三个核心输入是 `bag_dir`、`config_path` 和各输出路径；C++ 自己读取 bag。`setsid` 建立进程组，`taskset` 指定 CPU，`&` 让 Shell 继续执行监控，`algorithm_pid=$!` 为等待与采样提供对象。

| 参数 | 与算法的交接 |
|---|---|
| `--input_bag / --config` | 数据与前端 YAML；单/多雷达由配置决定 |
| `--output_tum / --output_lidar_tum / --output_rear_axle_tum` | IMU、主 LiDAR、后轴三种轨迹 |
| `--output_map / --output_frame_stats_csv` | LIO 点云与融合帧统计 |
| `--playback_rate / --max_lidar_frames` | 处理节奏与停止条件 |

## 按执行顺序阅读

| 阶段 | 目的 | 阅读入口 |
|---|---|---|
| 1. 解析输入和输出 | 确认 bag/YAML，保留本次产物 | @ref frontend_offline_prepare "运行准备" |
| 2. 加载环境、读取输入合同 | 确认二进制来源，取得末帧时间和数据时长 | @ref frontend_offline_prepare "环境与输入合同" |
| 3. 注册退出处理、启动算法 | 允许异常退出时清理算法进程组 | 核心调用；此前的 `trap cleanup EXIT` |
| 4. 采样、等待和超时处理 | 记录资源，避免无限等待 | @ref frontend_offline_wait "监控与等待" |
| 5. 汇总并判定运行结果 | 记录完成程度、计时和产物完整性 | @ref frontend_offline_finish "结果与退出" |

@anchor frontend_offline_prepare
### 1–2. 运行准备与输入合同

@snippet{lineno} run_frontend_offline.sh frontend-offline-inputs

bag 必须有 `metadata.yaml`，配置必须存在。目录中若已有脚本拥有的 logs/results/metadata 等产物，会拒绝覆盖。输出参数的相对路径由 `resolve_output()` 放到本次目录下。

@snippet{lineno} run_frontend_offline.sh frontend-offline-environment

加载环境后核对包 prefix，避免跑到其他安装树；随后检查算法与采样/计时工具。

@snippet{lineno} run_frontend_offline.sh frontend-offline-bag-contract

`bag_contract.json` 提供主 LiDAR、预期末帧、传感器时长与输入指纹。脚本据时长、回放倍率和 margin 计算 watchdog 超时，然后进入核心调用。

@anchor frontend_offline_wait
### 4. 启动资源采样，再等待算法

@snippet{lineno} run_frontend_offline.sh frontend-offline-monitor

资源采样发生在算法启动之后，因为需要其 PID。主脚本随后循环查看算法是否结束或超过 deadline；超时时清理算法进程组，再取回退出码。

@snippet{lineno} run_frontend_offline.sh frontend-offline-watchdog

`set +e` 允许失败的算法返回后继续归档。停止采样、取消 EXIT trap 后，进入统计阶段。首次阅读先跟踪 `algorithm_rc`、`watchdog_status` 和 `completion`，再看具体的 Python/awk 统计。

@anchor frontend_offline_finish
### 5. 结果与退出

轨迹检查、计时提取和 metadata 写入后，脚本用下面的条件决定是否返回失败：

@snippet{lineno} run_frontend_offline.sh frontend-offline-finish

因此退出码同时反映算法运行与 runner 的结果合同。限帧任务允许 `limited_frame_run`；完整任务要求到达最后 LiDAR。进入 C++ 时读 @ref run_frontend_offline.cc "main()"，沿 `LaserMapping::Init → rosbag.Go → 最后一帧 flush → 地图导出` 继续。

<details>
<summary>参数、计时与产物：需要调整实验时展开</summary>

@snippet{lineno} run_frontend_offline.sh frontend-offline-usage

`--` 后的 flags 进入 `extra_args`，再原样传给 C++。默认不按传感器时间节拍限速，使用 CPU `0-7`、8 个逻辑核；正式实验可提供 inventory 与指纹。主要结果在 results/，日志在 logs/，资源采样、bag 合同和运行信息在运行根目录。

计时提取只在算法等待结束后执行：

@snippet{lineno} run_frontend_offline.sh frontend-offline-timing

完整源码：@ref run_frontend_offline.sh "run_frontend_offline.sh"。

</details>

@anchor frontend_batch_guide
# run_frontend_offline_batch.sh：重复调度单次前端

## 先看核心调用

@snippet{lineno} run_frontend_offline_batch.sh frontend-batch-command

外层先遍历 repeat，内层遍历 YAML。`command` 每次构造一个完整的 @ref frontend_offline_guide "单次离线任务"，输出放到 `<配置名>/repeat_NN/`。此脚本本身没有 C++ 算法；进入 C++ 的位置在子脚本中。

## 按执行顺序阅读

1. 解析 bag、输出根和可选配置；检查输入、runner 与 repeats。
2. 得到 `configs`：未指定 `--config` 时按文件名顺序取目录下全部 YAML，否则解析所选名称/路径；检查名称冲突。
3. 拒绝覆盖旧 `batch_summary.tsv`，建立批次汇总；按上面的双层循环执行任务。
4. 每个任务记录状态、退出码和 YAML 指纹；继续后续任务，最后根据失败总数退出。

@anchor frontend_batch_results
### 任务失败后为什么还能继续？

@snippet{lineno} run_frontend_offline_batch.sh frontend-batch-results

子脚本放在 `if` 的条件中，失败会进入 `else`，保留退出码并增加计数。循环完成后，有失败就返回 4：

@snippet{lineno} run_frontend_offline_batch.sh frontend-batch-finish

<details>
<summary>配置选择与运行参数</summary>

@snippet{lineno} run_frontend_offline_batch.sh frontend-batch-configs

@snippet{lineno} run_frontend_offline_batch.sh frontend-batch-usage

`--dry-run` 打印命令后跳过执行与产物写入。`--` 后的参数交给单次 runner，批量层只负责调度与汇总。

完整源码：@ref run_frontend_offline_batch.sh "run_frontend_offline_batch.sh"。

</details>
