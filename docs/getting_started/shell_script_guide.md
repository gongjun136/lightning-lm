@page shell_script_guide 脚本目录与在线定位入口

# 从任务找到脚本

Shell 脚本负责环境、输入、输出和进程编排；定位与建图算法由它启动的 C++ 程序实现。先按下表选择导读，再阅读执行顺序；精确参数与副作用查 @ref script_contracts "脚本契约"。每个脚本都提供主线导读与完整源码。导读中的片段从真实脚本提取，按语义区段定位，行号由工具自动显示。完整分组见 @ref guide_scripts "脚本启动主线"。

现场操作从仓库外的 `${SANY_WS}/run.sh` 进入，它加载集中路径配置并固定 `lidar / ins_only` 分支。开发联调时可直接阅读或调用仓库脚本。`fusion` 仅为规划，尚无可运行入口。`production / diagnostic` 是 LiDAR 分支的开销与采集设置，与定位方式分别管理。

## 定位、前端与建图

| 要做的事 | 脚本 | 职责与输入 | 主要结果 |
|---|---|---|---|
| SANY LiDAR 地图定位 | @ref sany_lidar_launcher_guide "run_sany_lidar_loc.sh 导读" · @ref run_sany_lidar_loc.sh "源码" | 现场参数、LiDAR/IMU/CAN 预检、运行证据、可选诊断和进程管理；详见下文 | 三类 TUM、日志、配置快照、运行信息；可选 MCAP/资源记录/静默快照 |
| SANY CGI-430 纯组合导航 | @ref sany_ins_launcher_guide "run_sany_ins_only.sh 导读" · @ref run_sany_ins_only.sh "源码" | 必需现场 YAML；确认 `ins_only` 后启动 `InsLocSystem`；详见下文 | 导航审核记录、发布位姿和状态；驱动、录包另行启动 |
| 通用在线定位联调 | @ref generic_online_launcher_guide "run_loc_online.sh 导读" · @ref run_loc_online.sh "源码" | ROS 话题输入；建立日志目录、透传算法 flags；详见下文 | 基本运行信息、stdout/stderr；其他产物由 flags 决定 |
| 离线地图定位 | @ref loc_offline_guide "run_loc_offline.sh 导读" · @ref run_loc_offline.sh "源码" | bag、配置、地图；启动定位、资源监控和结果分析 | 轨迹、定位统计与误差报告、日志和结果 bag |
| 在线前端 | @ref frontend_online_guide "run_frontend_online.sh 导读" · @ref run_frontend_online.sh "源码" | 加载环境、启动实时前端；传感器或外部 `rosbag play` 提供输入 | 前端 ROS 输出、点云与里程计 |
| 单次离线前端 | @ref frontend_offline_guide "run_frontend_offline.sh 导读" · @ref run_frontend_offline.sh "源码" | bag、配置、输出目录；管理处理、超时和资源采样 | TUM、PCD、帧统计、耗时/资源记录与运行信息 |
| 批量离线前端 | @ref frontend_batch_guide "run_frontend_offline_batch.sh 导读" · @ref run_frontend_offline_batch.sh "源码" | 同一 bag、多配置和重复次数；依次调用单次入口 | 各配置独立运行目录和批次汇总 |
| 在线建图 | @ref slam_online_guide "run_slam_online.sh 导读" · @ref run_slam_online.sh "源码" | 建立运行目录、启动完整在线 SLAM；输入驱动另行启动 | 日志、在线 SLAM 输出；地图由配置/服务控制保存 |
| 离线建图 | @ref slam_offline_guide "run_slam_offline.sh 导读" · @ref run_slam_offline.sh "源码" | bag、配置、输出目录；编排前端、后端和地图导出 | 轨迹、分块/全局地图、后端信息和重定位数据库 |

## 构建、地图与导航数据工具

| 脚本 | 具体作用 | 调用前需要知道 |
|---|---|---|
| @ref build_workspace_guide "build_workspace.sh 导读" · @ref build_workspace.sh "源码" | 检查所需包，加载 ROS，顺序构建统一工作区 | `LIGHTNING_LM_BUILD_WS` 指定工作区；`LIGHTNING_BUILD_JOBS` 默认 2；可用 `--cmake-clean-cache` |
| @ref build_ins_guide "build_ins_only.sh 导读" · @ref build_ins_only.sh "源码" | 将参数交给统一构建入口的兼容包装 | 不再单独构建一套 INS 工作区 |
| @ref install_dep_guide "install_dep.sh 导读" · @ref install_dep.sh "源码" | 安装基础开发依赖 | 执行 `sudo apt`，会修改系统软件包 |
| @ref align_maps_guide "align_maps_offline.sh 导读" · @ref align_maps_offline.sh "源码" | 将参数交给地图刚体对齐程序 | 输入/输出由程序 flags 指定；不启动定位节点 |
| @ref save_map_guide "save_default_map.sh 导读" · @ref save_default_map.sh "源码" | 请求 `/lightning/save_map` 服务保存 `new_map` | 在线建图节点必须已运行；此脚本不负责建图 |
| @ref ins_record_guide "record_ins_only.sh 导读" · @ref record_ins_only.sh "源码" | 从现场 YAML 取 CGI/雷达输入列表，并录制业务输出、状态和心跳 | `site.yaml output_bag`；先加载 ROS/消息环境并设置一致的 domain；使用 MCAP |
| @ref ins_replay_guide "replay_ins_only.sh 导读" · @ref replay_ins_only.sh "源码" | 只回放配置中的输入话题，生成新 `/clock` | `site_replay.yaml bag [play flags]`；YAML 开启 `use_sim_time`，定位节点在另一终端先启动 |

## 回归与复现实验

这些脚本用于回放、对比或冻结实验矩阵。脚本中的数据目录和实验条件属于相应场景，具体参数见 @ref script_contracts "复现与回归契约"。

| 脚本 | 具体作用与产物 |
|---|---|
| @ref online_regression_guide "run_online_bag_regression.sh 导读" · @ref run_online_bag_regression.sh "源码" | 编排在线节点、bag 输入回放和结果录制；可请求保存地图，结束时清理子进程 |
| @ref voxel_reference_guide "run_voxel_slam_114_reference.sh 导读" · @ref run_voxel_slam_114_reference.sh "源码" | 运行 ROS 1 Voxel-SLAM 参考轨迹实验，管理独立 ROS master、回放、录制及资源/轨迹提取 |
| @ref sany_offline_pipeline_guide "run_sany_offline_localization.sh 导读" · @ref reproduction/multi_lidar/sany_4livox/run_sany_offline_localization.sh "源码" | 先生成 SANY 地图，再用该地图运行离线定位，保存两阶段产物 |
| @ref phase_a_matrix_guide "run_phase_a_relocalization_matrix.sh 导读" · @ref reproduction/multi_lidar/sany_4livox/run_phase_a_relocalization_matrix.sh "源码" | 按数据集和传感器起始时间偏移运行重定位任务，保存各任务结果并打印 runner 失败总数 |
| @ref fault_bags_guide "create_sany_formal_fault_bags.sh 导读" · @ref reproduction/multi_lidar/sany_4livox/create_sany_formal_fault_bags.sh "源码" | 批量生成已定义的丢包/中断故障输入和场景合同 |
| @ref fault_matrix_guide "run_sany_formal_fault_matrix.sh 导读" · @ref reproduction/multi_lidar/sany_4livox/run_sany_formal_fault_matrix.sh "源码" | 逐个运行四个正式故障场景，调用在线合同脚本并汇总结果 |
| @ref sany_online_contract_guide "run_sany_online_contract.sh 导读" · @ref reproduction/multi_lidar/sany_4livox/run_sany_online_contract.sh "源码" | 管理在线前端、输入回放和输出录制，验收五个公开 Topic 及故障恢复 |
| @ref fastlivo_rerun_guide "rerun_m3dgr_fastlivo2_timestamp_fix.sh 导读" · @ref reproduction/formal_report/rerun_m3dgr_fastlivo2_timestamp_fix.sh "源码" | 按冻结的 M3DGR/FAST-LIVO2 矩阵重跑序列，保存指纹、运行清单和汇总；依赖报告对应的数据与工具布局 |
| @ref voxel_backend_guide "run_voxel_slam_full_backend.sh 导读" · @ref reproduction/single_lidar/m3dgr/run_voxel_slam_full_backend.sh "源码" | 运行 ROS 1 Voxel-SLAM 完整后端对比实验；管理 roscore、回放、结果和资源记录 |

@anchor sany_lidar_launcher_guide
# run_sany_lidar_loc.sh：启动主线导读

这份脚本完成 SANY LiDAR 定位的启动准备、程序运行和退出收尾。先看下面的核心启动代码，再沿阶段表回看它需要的参数、环境和输入。完整脚本见 @ref run_sany_lidar_loc.sh "源码页"；参数、诊断分支和输出文件放在本节末尾的 @ref sany_lidar_reference "按需参考"。

```text
${SANY_WS}/run.sh：加载路径配置，选择 lidar，固定发布预设
  → run_sany_lidar_loc.sh：准备环境与输入，管理定位和辅助进程
      → C++ run_loc_online：main → LocSystem → ROS 事件循环
```

## 先抓住核心：启动、等待、退出

@anchor sany_lidar_launch
### 把哪些参数交给 C++ 程序？

`algorithm_args` 是完整的程序命令与参数数组。`ros2 run` 的包名是 `lightning_lm`，可执行文件名是 `run_loc_online`。本段确定实际配置、地图和三份轨迹的保存位置：

@snippet{lineno} run_sany_lidar_loc.sh sany-lidar-launch-arguments

| 变量 | 从哪里来 | 交给程序做什么 |
|---|---|---|
| `config_path` | 显式 `LIGHTNING_LM_CONFIG`，或按雷达布局选择的内置 YAML | `--config`：定位系统的配置 |
| `map_path` | `SANY_MAP_PATH`，由根入口固定现场地图 | `--map`：覆盖配置中的地图路径 |
| `run_dir` | 输出根目录与本次 run name | 保存全局、高频和实际发布的后轴轨迹，以及运行日志 |

接下来，`algorithm_launcher` 补上可选的 CPU 核心限制和日志缓冲设置，`setsid` 真正发起程序启动：

@snippet{lineno} run_sany_lidar_loc.sh sany-lidar-launch-localization

读这段时抓住四处：

- `"${algorithm_args[@]}"`：把数组中的每一项作为独立参数传给程序。
- `setsid`：建立独立进程组，停止时可以同时通知 ROS 启动器和它启动的算法进程。
- 末尾的 `&`：后台启动，主脚本继续管理本次运行。
- `algorithm_pid=$!`：记住启动进程的 PID，后续用于等待、资源采样和收尾。

`stdout` 和 `stderr` 分别写入本次运行目录。这里直接启动 C++ 程序；通用 `run_loc_online.sh` 是另一条联调入口。

@anchor sany_lidar_wait
### 程序已经启动，主脚本为什么还没有结束？

算法启动后，按开关启动资源采样，然后主脚本在 `wait` 处等待算法进程结束：

@snippet{lineno} run_sany_lidar_loc.sh sany-lidar-wait-localization

`algorithm_status=$?` 保存算法退出状态；随后写入运行信息、提取计时并停止辅助进程，最后把该状态作为脚本退出码：

@snippet{lineno} run_sany_lidar_loc.sh sany-lidar-finish-run

收到 `Ctrl-C` 或 TERM 时走预先注册的 `trap`，调用同一收尾函数 `stop_children` 并返回 130。`stop_children` 向定位/录包进程组发送 INT，让轨迹和 MCAP 完成写入，并停止 watchdog、资源采样等辅助进程。

## 再按执行顺序理解准备工作

表中入口对应下面的稳定阅读区段。源码片段随文档构建从真实脚本提取；行号由工具显示，正文无需跟着代码位置变化修改。

| 执行阶段 | 这一步解决什么问题 | 阅读入口 |
|---|---|---|
| 1. 确定运行参数 | 本次用哪个配置、地图和运行名称？ | @ref sany_lidar_runtime "运行参数" |
| 2. 检查并加载环境 | 文件、开关、ROS 和安装环境是否可用？ | @ref sany_lidar_environment "环境准备" |
| 3. 确定输入话题 | 当前配置要求哪些传感器输入？ | @ref sany_lidar_inputs "输入清单" |
| 4. 准备运行目录 | 日志放哪里，怎样保存实际配置？ | @ref sany_lidar_artifacts "运行目录" |
| 5. 注册收尾并等待输入 | 怎样停止本次运行，输入是否已出现且类型正确？ | @ref sany_lidar_input_gate "启动前的输入门槛" |
| 6. 记录信息、启动可选采集 | 保留哪些运行信息，需要哪些辅助进程？ | @ref sany_lidar_collectors "证据与采集" |
| 7. 启动定位程序 | 怎样把配置、地图和输出路径交给 C++？ | @ref sany_lidar_launch "核心启动代码" |
| 8. 等待与收尾 | 怎样等待算法退出、停止辅助进程？ | @ref sany_lidar_wait "等待与退出" |

@anchor sany_lidar_runtime
### 1. 确定运行参数

输入是根入口导出的环境变量，以及可选的 run name。脚本先设置可覆盖的默认值，按 `SANY_LIDAR_LAYOUT` 选择内置三/四雷达配置，再解析本次运行的路径：

@snippet{lineno} run_sany_lidar_loc.sh sany-lidar-runtime-settings

`${变量:-默认值}` 表示外部有效值优先。因此，现场根入口固定的地图、三雷达、production 和关闭录包预设会传入本次运行。`LIGHTNING_LM_RUN_MODE` 同样可覆盖，决定计时、资源采样和诊断 I/O 的默认开关。此阶段得到 `config_path`、`map_path`、`out_root`、`run_name` 和各功能开关。

@anchor sany_lidar_environment
### 2. 检查并加载环境

先检查环境脚本和定位配置是否可读：

@snippet{lineno} run_sany_lidar_loc.sh sany-lidar-environment-files

再加载 ROS 与已构建的工作区安装环境，使 `ros2 run` 找到当前程序和消息包：

@snippet{lineno} run_sany_lidar_loc.sh sany-lidar-load-environment

这一阶段还检查地图和输出路径、开关值、线程数以及 CPU 核心限制；开启录包时检查 QoS 和 MCAP。检查失败会提前结束，此时还没有启动定位程序。

@anchor sany_lidar_inputs
### 3. 从配置确定必须等待的输入

脚本中的 Python 段解析 YAML，得到 `lidar_topics` 和 IMU 清单。关闭录包时需要主 IMU，开启录包时需要全部配置 IMU。最终拼成启动前的输入清单：

@snippet{lineno} run_sany_lidar_loc.sh sany-lidar-required-inputs

CAN 观测启用时，将轮速话题加入必需输入。`required_topics` 是稍后等待输入的依据，`record_topics` 是可选录包的采集清单，两者用途不同。首次阅读先理解这两个结果，再按需展开 YAML 解析细节。

@anchor sany_lidar_artifacts
### 4. 建立本次运行目录并保存配置

先把输出根与 run name 拼成 `run_dir`。同名目录存在时拒绝启动，以免覆盖旧证据：

@snippet{lineno} run_sany_lidar_loc.sh sany-lidar-run-directory

创建 logs/results/snapshots 等目录后，保存实际 YAML 和计算预算摘要：

@snippet{lineno} run_sany_lidar_loc.sh sany-lidar-config-snapshot

此阶段把“本次用了什么配置”和后续产物绑定到同一个目录。算法工作目录仍是仓库目录，三份轨迹通过明确的参数写入 `run_dir`。

@anchor sany_lidar_input_gate
### 5. 注册收尾，等待并核对传感器输入

先注册停止信号和退出处理，后续启动的辅助进程由统一收尾函数管理：

@snippet{lineno} run_sany_lidar_loc.sh sany-lidar-register-cleanup

然后调用 `wait_for_inputs` 等待必需话题，并检查 LiDAR、IMU 和启用的 CAN 输入类型：

@snippet{lineno} run_sany_lidar_loc.sh sany-lidar-check-inputs

这一关检查话题存在和类型，不代表已经完成频率、时延或定位质量验收。通过后才进入采集与算法启动阶段。

**函数定义与执行顺序：** `wait_for_inputs()`、`watch_pose_vel()` 和 `stop_children()` 出现在前面时只是定义函数；定义处不会执行函数体。上面这段 `wait_for_inputs || fail ...` 才是等待输入的调用点。首次阅读沿阶段表看调用点，遇到需要理解的细节再回到函数定义。

@anchor sany_lidar_collectors
### 6. 保存运行信息并按条件启动辅助采集

输入检查通过后，将版本、配置路径与哈希、最终开关写入 `run_metadata.txt`，将未提交改动清单和启动时话题名称/类型分别写入 `git_status.txt`、`topics_at_start.txt`，用于确认本次运行使用了什么代码、配置和输入。随后按下表启动辅助采集；主线继续进入 @ref sany_lidar_launch "定位程序启动"。以下输出路径均相对于本次 `run_dir`。

| 辅助采集 | 采集对象与用途 | 启动条件与顺序 | 主要产物 |
|---|---|---|---|
| MCAP 录包 | 保存传感器原始输入、定位输出与状态，供回放和事后对照 | `SANY_RECORD_BAG=1`；算法之前启动并检查录包进程仍在运行 | `bag/localization_incident/`、`logs/rosbag.stdout.log` / `.stderr.log` |
| 位姿静默 watchdog | 检测 `/localization/pose_vel` 消息接收是否中断，记录中断/恢复事件并保存中断现场 | `SANY_ENABLE_POSE_VEL_WATCHDOG=1`；算法之前启动 | `logs/pose_vel_watchdog.csv`、`snapshots/loss_*.txt` |
| `tegrastats` | 记录 Jetson/Orin 整机 CPU/GPU、内存、温度和功耗等，辅助排查机器资源压力 | 当前环境能找到该命令；算法之前启动，每 1000 ms 采样 | `logs/tegrastats.log` |
| 定位进程与线程资源采样 | 在本次算法进程组中找到真实 `run_loc_online`，记录该进程及线程的 CPU、调度等待、内存和 I/O | `LIGHTNING_LM_RESOURCE_PROFILE=1`；算法之后启动，因为需要 `algorithm_pid` | `results/process_resources.jsonl`、`logs/resource_monitor.stderr.log` |

录包与 watchdog 缺省关闭；切到 diagnostic 不会自动开启它们。进程资源采样缺省在 diagnostic 开启、production 关闭；显式子开关优先。`tegrastats` 只按命令是否可用决定启动，当前脚本没有独立开关，因此安装了该工具时，production 也会采集。

**watchdog 具体监测什么？** `watch_pose_vel()` 调用 `scripts/monitor_pose_vel_silence.py`，持续订阅 `/localization/pose_vel`（`lightning/msg/VehiclePose`）。`PoseVelSilenceMonitor.check_silence()` 根据本机单调时钟计算距上次消息接收的时间，不检查消息内的坐标、速度或 Header 时间戳。收到首条消息前不触发静默告警；收到首条后，达到 `SANY_POSE_VEL_TIMEOUT_SECONDS`（缺省 5 秒）且再经过缺省 0.25 秒确认窗口仍未收到消息，才报告 `LOST`。实际报告时刻还受定时器和进程调度影响。

`FIRST_POSE`、`LOST`、`RECOVERED` 写入 CSV，字段为墙钟时间、事件、事件编号和静默秒数。一次持续中断只记录一次 `LOST` 并抓取一次快照；再次收到消息时记录 `RECOVERED`，后续新的中断再增加编号。`snapshot_incident()` 抓取流水线诊断、故障状态、ROS 节点/话题与输入发布订阅信息，以及进程、内存、磁盘、网卡和时钟同步状态，便于回查中断时的现场。诊断/故障话题各等待最多 3 秒，快照各项依次采集，不是同一瞬间的原子快照。

`LOST` 表示监测端没有收到输出，可能涉及输入中断、输出门控、DDS 通信或调度延迟，不能直接等同于算法已定位失败；`RECOVERED` 也只表示消息恢复接收。输出持续到达但位姿错误、速度异常或时间戳陈旧，不会触发此静默检测。watchdog 只留存证据，不停止/重启定位，也不调用重定位；定位恢复由 C++ 内部逻辑负责。

**tegrastats 有什么作用？** 它是 NVIDIA Jetson 的整机资源统计工具。脚本执行 `tegrastats --interval 1000` 并保存原始输出，用来观察 CPU 各核负载/频率、RAM/SWAP、GPU 活跃度/频率（`GR3D_FREQ`）、内存控制器带宽使用率/频率（`EMC`），以及温度和各供电轨功耗；实际字段随设备和 Jetson Linux 版本变化。字段口径见 [NVIDIA tegrastats 官方说明](https://docs.nvidia.com/jetson/archives/r36.4.3/DeveloperGuide/AT/JetsonLinuxDevelopmentTools/TegrastatsUtility.html)。例如，输出变慢时可检查是否同时出现整机高负载、温度升高或频率下降，再结合进程采样与算法日志验证原因。它统计整台设备，包含驱动和其他业务进程，不能将全部 CPU/GPU 占用归给定位，也不能用它衡量定位精度或单个函数耗时。

**进程资源采样与 `tegrastats` 怎么配合？** `scripts/monitor_process_resources.py` 的 `discover_target()` 通过进程组定位真实 `run_loc_online`，而不是汇总进程组中所有进程；之后从 Linux `/proc` 读取该进程及各线程的 CPU、内存、读写增量、上下文切换和调度等待，并附带整机各核忙碌率及区间 CPU 消耗最高的 20 个进程作为背景。默认每 1 秒采样，间隔由 `LIGHTNING_LM_RESOURCE_INTERVAL_SEC` 控制（0.5–60 秒）。`tegrastats` 回答“整机是否繁忙”，进程采样帮助判断“定位自身及哪些线程消耗了 CPU、是否等待调度”。

采样中的 `cpu_core_equivalents=3.2` 表示区间平均消耗约 3.2 个逻辑核的计算时间，对应单核口径 `320%`；线程 `sched_wait_ns_delta` 是等待 CPU 调度的时间，不是互斥锁或算法队列等待，内核未启用相应统计时零值也不能证明没有等待。逐帧/模块耗时查 `results/compute_profile.log`，异步速度、队列与锁事件查 `results/causal_trace.csv`，分别受 `LIGHTNING_LM_COMPUTE_PROFILE` 和 `LIGHTNING_LM_CAUSAL_TRACE` 控制。详细口径见 @ref resource_diagnostics "资源与端到端延迟"、@ref causal_diagnostics "因果链与时间关联"。

以上启动与快照逻辑可对照 @ref run_sany_lidar_loc.sh "完整脚本" 中的 `watch_pose_vel()`、`snapshot_incident()` 和辅助采集调用点；具体启动片段见下方折叠参考。

@anchor sany_lidar_cpp_entry
## 到这里进入 C++，接下来读哪里？

`ros2 run lightning_lm run_loc_online` 解析安装空间中的可执行文件并启动它。后续源码阅读按下面的顺序进行：

1. @ref run_loc_online.cc "run_loc_online.cc 的 main()"：读取 flags，按 YAML 的 `system.localization_mode` 选择分支；lidar 使用 `LocSystem`。
2. @ref lightning::LocSystem::Init() "LocSystem::Init()"：使用刚传入的配置与地图，初始化定位系统和 ROS 输入输出。
3. @ref lightning::LocSystem::SetInitPose() "LocSystem::SetInitPose()"：设置初始位姿。
4. @ref lightning::LocSystem::Spin() "LocSystem::Spin()"：进入 ROS 事件循环，开始沿传感器回调处理在线数据。

接着转入 @ref online_localization_flow "在线定位端到端流程"，跟踪 `LIO → 地图匹配 → 定位 PGO → 业务输出`。Shell 继续负责运行管理，估计与地图匹配由 C++ 实现。

@anchor sany_lidar_reference
## 按需参考：参数、诊断分支与输出

@htmlonly[block]
<details>
<summary>参数与默认值：需要调整启动设置时展开</summary>
@endhtmlonly

可选第一个参数是 run name，省略为 `sany_<布局>lidar_YYYYmmdd_HHMMSS`。`-h / --help` 显示脚本用法。此入口按环境变量组装算法 flags，不提供通用脚本那样的剩余 flags 透传。

| 设置 | 脚本缺省值 | 用途与覆盖方式 |
|---|---|---|
| `ROS_DOMAIN_ID` | `42` | ROS 通信域；驱动与下游须使用同一 domain |
| `LIGHTNING_LM_WS` | `/home/nvidia/project/gj_ws` | 旧布局的 fallback 基准；当前部署应加载集中路径配置 |
| `LIGHTNING_LM_REPO_DIR` | `$LIGHTNING_LM_WS/lightning-lm` | 仓库路径，算法 cwd；统一 `src/` 布局应显式覆盖 |
| `LIGHTNING_LM_ROS_SETUP` | `/opt/ros/humble/setup.bash` | ROS 环境 |
| `LIGHTNING_LM_INSTALL_SETUP` | `$REPO/install/setup.bash` | 定位工作区安装环境；统一工作区应显式覆盖 |
| `LIGHTNING_LM_OUT_ROOT` | `$LIGHTNING_LM_WS/runs` | 本次日志和证据输出的根目录 |
| `SANY_LIDAR_LAYOUT` | `3` | 选择内置三/四雷达 YAML；不自动探测设备 |
| `LIGHTNING_LM_CONFIG` | 所选布局的正式定位 YAML | 可显式覆盖为绝对路径；输入话题从该配置解析 |
| `SANY_MAP_PATH` | `$LIGHTNING_LM_WS/maps/sany_4lidar_20260811_a4defa3_solid` | 保留旧入口的历史 fallback；现场根入口必须固定实际发布地图 |
| `SANY_ENABLE_CAN_OBSERVATION` / `SANY_WHEEL_SPEED_TOPIC` | `1` / `/SpeThrCAN4_topic` | CAN 轮速观测及输入话题；启用时预检该输入 |
| `SANY_RECORD_BAG` | `0` | 可选后台 MCAP；启用后采集全部配置雷达/IMU和 CAN 等话题 |
| `SANY_TOPIC_WAIT_SECONDS` | `60` | 等待必需传感器话题出现的时间 |
| `SANY_IMU_TOPIC` | 配置中的主 IMU | 可覆盖主 IMU 输入；录包时还检查全部配置 IMU |
| `SANY_RECORD_QOS_FILE` | `scripts/config/sany_localization_record_qos.yaml` | 录包 QoS 覆盖 |
| `SANY_MIN_FREE_GB` | `20` | 仅在录包时要求输出目录所在磁盘保有此空间 |
| `SANY_ENABLE_POSE_VEL_WATCHDOG` / `SANY_POSE_VEL_TIMEOUT_SECONDS` | `0` / `5` | 可选定位输出静默检测和现场快照；不自动重启 |

所有默认值都可由调用环境覆盖。项目根入口固定的三雷达、地图、`production` 和关闭录包等设置优先于仓库默认值，最终有效设置记录到 `run_metadata.txt`。

仓库入口的现场默认设置取自真实脚本：

@snippet{lineno} run_sany_lidar_loc.sh sany-lidar-field-defaults

@htmlonly[block]
</details>
@endhtmlonly

@htmlonly[block]
<details>
<summary>production / diagnostic：需要调整采集开销时展开</summary>
@endhtmlonly

| 开关 | diagnostic 缺省 | production 缺省 | 作用 |
|---|---|---|---|
| `LIGHTNING_LM_COMPUTE_PROFILE` | `1` | `0` | 算法热点计时，退出后提取 `COMPUTE_BENCH_*` |
| `LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD` | `0` | `1` | 裁剪非必要的高频诊断 I/O |
| `LIGHTNING_LM_RESOURCE_PROFILE` | `1` | `0` | 在本次进程组中找到定位进程，对其进程/线程进行 CPU、内存、I/O 等采样 |
| `LIGHTNING_LM_CAUSAL_TRACE` | `1` | `0` | 采集有限的异步速度、队列与锁事件 |

`LIGHTNING_LM_RUN_MODE` 缺省 `diagnostic`；显式子开关优先于模式默认值。录包与 watchdog 有独立开关，切到 diagnostic 不会自动开启它们。安装了 `tegrastats` 时，脚本还会记录 Orin 资源信息。

计算资源覆盖包括 `LIGHTNING_LM_LIO_THREADS`、`LIGHTNING_LM_NDT_THREADS`、`LIGHTNING_LM_NDT_MAX_POINTS`、`LIGHTNING_LM_SOLID_ICP_WORKERS`、`LIGHTNING_LM_SOLID_WORKER_NICE`。`LIGHTNING_LM_CPU_AFFINITY` 限定算法进程的核心，`LIGHTNING_LM_SOLID_CPU_AFFINITY` 限定 SOLiD worker 且须为前者有效核心的子集。`LIGHTNING_LM_RESOURCE_INTERVAL_SEC` 缺省 1 秒，范围 0.5–60 秒。实际预算同时受 YAML `compute_budget` 控制，见 @ref configuration "配置入口"。

两种模式的默认开关在脚本中集中选择：

@snippet{lineno} run_sany_lidar_loc.sh sany-lidar-mode-settings

@htmlonly[block]
</details>
@endhtmlonly

@htmlonly[block]
<details>
<summary>录包、watchdog 与资源监测：需要诊断现场时展开</summary>
@endhtmlonly

**后台 MCAP：** 在定位启动前录制原始输入和指定输出，先等待一秒并确认录包进程仍在运行，失败时结束启动。

@snippet{lineno} run_sany_lidar_loc.sh sany-lidar-start-recording

**位姿静默 watchdog：** 启用时后台执行 `watch_pose_vel`，记录它的 PID，退出时由 `stop_children` 管理；消息静默的判定、事件与快照内容见 @ref sany_lidar_collectors "第 6 步辅助采集说明"。

@snippet{lineno} run_sany_lidar_loc.sh sany-lidar-start-watchdog

**算法资源采样：** 算法已启动，使用 `algorithm_pid` 作为进程组标识，找到并采样其中真实的 `run_loc_online` 进程及线程，并将采样进程加入清理列表。

@snippet{lineno} run_sany_lidar_loc.sh sany-lidar-resource-monitor

@htmlonly[block]
</details>
@endhtmlonly

@htmlonly[block]
<details>
<summary>运行目录与输出文件：需要查看结果或排障时展开</summary>
@endhtmlonly

| 产物 | 用途 |
|---|---|
| `config.yaml`、`launch_compute_config.json` | 实际配置快照和计算预算摘要 |
| `run_metadata.txt`、`git_status.txt`、`topics_at_start.txt` | 本次版本、开关、输入图和退出状态 |
| `logs/run_loc_online.stdout.log` / `.stderr.log` | 算法原始输出；glog 和热点计时主要在 stderr |
| `results/trajectory_global.tum` | PGO 全局校正轨迹 |
| `results/trajectory_high_frequency.tum` | 定位器内部高频轨迹 |
| `results/trajectory_published_rear_axle.tum` | 实际发布边界后的后轴轨迹 |
| `results/compute_profile.log` | 从 stderr 提取的计时记录；采集开关决定内容 |
| `results/process_resources.jsonl` | 定位进程/线程资源采样，以及整机 CPU 和其他进程的背景统计 |
| `results/causal_trace.csv` | 开启追踪时保存的异步速度、队列与锁事件 |
| `logs/tegrastats.log` | 工具可用时每秒记录的 Jetson/Orin 整机资源、温度和功耗 |
| `bag/` | 开启录包时保存的 MCAP |
| `snapshots/loss_*.txt`、`logs/pose_vel_watchdog.csv` | 开启 watchdog 时的输出中断现场与首次接收/中断/恢复时间线 |

@htmlonly[block]
</details>
@endhtmlonly

@anchor sany_ins_launcher_guide
# run_sany_ins_only.sh：SANY 纯组合导航入口

## 先看核心调用

@snippet{lineno} run_sany_ins_only.sh ins-launch

仍然启动 C++ `run_loc_online`，通过现场 YAML 选择 `InsLocSystem` 分支。`config` 是第一个位置参数，剩余 flags 由 `"$@"` 透传；`exec` 替换 Shell，脚本没有后续的监控与等待阶段。

```text
现场 YAML → 环境与模式核对 → C++ run_loc_online / InsLocSystem
CGI-430 已解算的位置、姿态、速度 → 坐标与质量处理 → 业务输出
```

## 按执行顺序阅读

| 阶段 | 目的 | 阅读入口 |
|---|---|---|
| 1. 解析 YAML 和工作区 | 得到绝对配置路径与 overlay | @ref sany_ins_settings "配置与工作区" |
| 2. 加载环境和 domain | 让节点找到程序、消息包与输入 | @ref sany_ins_environment "环境" |
| 3. 核对定位模式 | 避免把 LiDAR YAML 交给 INS 入口 | @ref sany_ins_mode "模式门槛" |
| 4. 启动 C++ | 使用 CGI 导航解进行处理与发布 | 核心调用 |

@anchor sany_ins_settings
### 1. 配置与工作区

@snippet{lineno} run_sany_ins_only.sh ins-config

第一个参数必需；`realpath` 后 `shift` 消费它。仓库位于工作区 `src/` 时向上推导工作区，否则使用兼容布局。显式 `LIGHTNING_LM_INSTALL_SETUP` 优先，其次是 `LIGHTNING_LM_BUILD_WS`。

@anchor sany_ins_environment
### 2. 加载环境和通信域

@snippet{lineno} run_sany_ins_only.sh ins-environment

ROS 环境固定为 Humble，domain 缺省 42；驱动和下游必须使用一致的通信域。

@anchor sany_ins_mode
### 3. 核对模式再启动

@snippet{lineno} run_sany_ins_only.sh ins-mode

Python 段只核对 YAML 的模式，成功后进入核心调用。接着读 @ref run_loc_online.cc "main() 的 ins_only 分支"，沿 `InsLocSystem::Init → Spin` 继续；坐标转换、质量门控、点云去畸变和输出审核见 @ref ins_only_operation "纯组合导航操作契约"。

@htmlonly[block]
<details>
<summary>驱动、录包、回放与输出</summary>
@endhtmlonly

```bash
bash scripts/run_sany_ins_only.sh /absolute/site.yaml [run_loc_online flags]
```

CGI 与所需雷达驱动分别启动。节点按 YAML 的 `audit_root` 保存审核记录，支持 `--output_published_tum`，INS 分支拒绝 LIO/PGO 轨迹参数和内嵌 `--bag`。录制见 @ref ins_record_guide "导航录包导读"，回放见 @ref ins_replay_guide "输入回放导读"；这两个工具不会替你启动定位节点。

现场配置需固定 ENU 原点、高程类型、后轴外参和质量条件。完整源码：@ref run_sany_ins_only.sh "run_sany_ins_only.sh"。

@htmlonly[block]
</details>
@endhtmlonly

@anchor generic_online_launcher_guide
# run_loc_online.sh：通用在线定位联调入口

## 先看核心调用

@snippet{lineno} run_loc_online.sh loc-online-launch

`--config` 由脚本给出，`"$@"` 透传地图、bag、轨迹等程序 flags。外面的圆括号创建子 Shell 并切换到 `run_dir`，所以透传的相对地图、bag 和输出路径都以该目录为基准。没有 `&`，调用在前台等待，stdout/stderr 分别写入 logs/。

此入口可用于实时话题、外部回放，或传入 C++ 支持的 `--bag`。它是独立的通用联调入口；SANY 完整运行管理从 @ref sany_lidar_launcher_guide "SANY LiDAR 导读" 进入。

## 按执行顺序阅读

| 阶段 | 目的 | 阅读入口 |
|---|---|---|
| 1. 加载 ROS 与 overlay | 找到当前安装程序 | @ref generic_online_environment "环境" |
| 2. 确定配置和名称 | 区分脚本参数与程序 flags | @ref generic_online_settings "运行设置" |
| 3. 建立日志目录 | 保存本次输入路径与命令信息 | @ref generic_online_directory "运行目录" |
| 4. 启动并等待 C++ | 在 run dir 内运行定位 | 核心调用 |

@anchor generic_online_environment
### 1. 加载环境

@snippet{lineno} run_loc_online.sh loc-online-environment

此后遇到 `-h / --help` 会调用真实安装二进制的帮助并退出，不创建运行目录。

@anchor generic_online_settings
### 2. 确定配置与运行名称

@snippet{lineno} run_loc_online.sh loc-online-settings

第一个非 `-` 开头参数作为 run name，`shift` 后才将其余参数交给算法。配置必须存在，随后转为绝对路径。

@anchor generic_online_directory
### 3. 建立目录与基本运行信息

@snippet{lineno} run_loc_online.sh loc-online-directory

`mkdir -p` 允许复用同名目录；重复 run name 会截断基本日志和 metadata。准备完成后进入核心调用，脚本以算法调用的退出状态结束，没有独立的录包/采样收尾逻辑。

## 进入 C++ 后继续读哪里

从 @ref run_loc_online.cc "run_loc_online.cc 的 main()" 开始，LiDAR 模式沿 `LocSystem::Init → SetInitPose → Spin`，后续见 @ref online_localization_flow "在线定位数据流"。Shell 只处理环境和运行目录，算法与可选内嵌 bag 回放在 C++ 中执行。

@htmlonly[block]
<details>
<summary>调用示例、默认值与产物</summary>
@endhtmlonly

```bash
LIGHTNING_LM_CONFIG=/absolute/localization.yaml \
LIGHTNING_LM_OUT_ROOT=/absolute/runs \
bash scripts/run_loc_online.sh online_test --map /absolute/map
```

默认配置是 `$repo_dir/config/default.yaml`，输出根为 `$repo_dir/runs`，run name 可由 `LIGHTNING_LM_RUN_NAME` 设置，位置参数优先。ROS 和 overlay 路径可分别通过 `LIGHTNING_LM_ROS_SETUP`、`LIGHTNING_LM_INSTALL_SETUP` 覆盖。

基本产物为 `run_metadata.txt` 与 `logs/run_loc_online.{stdout,stderr}.log`；轨迹等产物由算法 flags 决定。完整参数边界见 @ref script_contracts "通用运行契约"，完整源码见 @ref run_loc_online.sh "run_loc_online.sh"。

@htmlonly[block]
</details>
@endhtmlonly
