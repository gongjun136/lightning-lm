@page shell_script_guide Shell 脚本用途与定位入口详解

# 从任务找到脚本

Shell 脚本负责环境、输入、输出和进程编排；定位与建图算法由它启动的 C++ 程序实现。先按下表选择入口，再阅读执行顺序；精确参数与副作用查 @ref script_contracts "脚本契约"。点击脚本名进入完整源码与行号，每个源码页都有用途摘要和说明页链接。

现场操作从仓库外的 `${SANY_WS}/run.sh` 进入，它加载集中路径配置并固定 `lidar / ins_only` 分支。开发联调时可直接阅读或调用仓库脚本。`fusion` 仅为规划，尚无可运行入口。`production / diagnostic` 是 LiDAR 分支的开销与采集设置，与定位方式分别管理。

## 定位、前端与建图

| 要做的事 | 脚本 | 职责与输入 | 主要结果 |
|---|---|---|---|
| SANY LiDAR 地图定位 | @ref run_sany_lidar_loc.sh "run_sany_lidar_loc.sh" | 现场参数、LiDAR/IMU/CAN 预检、运行证据、可选诊断和进程管理；详见下文 | 三类 TUM、日志、配置快照、运行信息；可选 MCAP/资源记录/静默快照 |
| SANY CGI-430 纯组合导航 | @ref run_sany_ins_only.sh "run_sany_ins_only.sh" | 必需现场 YAML；确认 `ins_only` 后启动 `InsLocSystem`；详见下文 | 导航审核记录、发布位姿和状态；驱动、录包另行启动 |
| 通用在线定位联调 | @ref run_loc_online.sh "run_loc_online.sh" | ROS 话题输入；建立日志目录、透传算法 flags；详见下文 | 基本运行信息、stdout/stderr；其他产物由 flags 决定 |
| 离线地图定位 | @ref run_loc_offline.sh "run_loc_offline.sh" | bag、配置、地图；启动定位、资源监控和结果分析 | 轨迹、定位统计与误差报告、日志和结果 bag |
| 在线前端 | @ref run_frontend_online.sh "run_frontend_online.sh" | 加载环境、启动实时前端；传感器或外部 `rosbag play` 提供输入 | 前端 ROS 输出、点云与里程计 |
| 单次离线前端 | @ref run_frontend_offline.sh "run_frontend_offline.sh" | bag、配置、输出目录；管理处理、超时和资源采样 | TUM、PCD、帧统计、耗时/资源记录与运行信息 |
| 批量离线前端 | @ref run_frontend_offline_batch.sh "run_frontend_offline_batch.sh" | 同一 bag、多配置和重复次数；依次调用单次入口 | 各配置独立运行目录和批次汇总 |
| 在线建图 | @ref run_slam_online.sh "run_slam_online.sh" | 建立运行目录、启动完整在线 SLAM；输入驱动另行启动 | 日志、在线 SLAM 输出；地图由配置/服务控制保存 |
| 离线建图 | @ref run_slam_offline.sh "run_slam_offline.sh" | bag、配置、输出目录；编排前端、后端和地图导出 | 轨迹、分块/全局地图、后端信息和重定位数据库 |

## 构建、地图与导航数据工具

| 脚本 | 具体作用 | 调用前需要知道 |
|---|---|---|
| @ref build_workspace.sh "build_workspace.sh" | 检查所需包，加载 ROS，顺序构建统一工作区 | `LIGHTNING_LM_BUILD_WS` 指定工作区；`LIGHTNING_BUILD_JOBS` 默认 2；可用 `--cmake-clean-cache` |
| @ref build_ins_only.sh "build_ins_only.sh" | 将参数交给统一构建入口的兼容包装 | 不再单独构建一套 INS 工作区 |
| @ref install_dep.sh "install_dep.sh" | 安装基础开发依赖 | 执行 `sudo apt`，会修改系统软件包 |
| @ref align_maps_offline.sh "align_maps_offline.sh" | 将参数交给地图刚体对齐程序 | 输入/输出由程序 flags 指定；不启动定位节点 |
| @ref save_default_map.sh "save_default_map.sh" | 请求 `/lightning/save_map` 服务保存 `new_map` | 在线建图节点必须已运行；此脚本不负责建图 |
| @ref record_ins_only.sh "record_ins_only.sh" | 从现场 YAML 取 CGI/雷达输入列表，并录制业务输出、状态和心跳 | `site.yaml output_bag`；先加载 ROS/消息环境并设置一致的 domain；使用 MCAP |
| @ref replay_ins_only.sh "replay_ins_only.sh" | 只回放配置中的输入话题，生成新 `/clock` | `site_replay.yaml bag [play flags]`；YAML 开启 `use_sim_time`，定位节点在另一终端先启动 |

## 回归与复现实验

这些脚本用于回放、对比或冻结实验矩阵。脚本中的数据目录和实验条件属于相应场景，具体参数见 @ref script_contracts "复现与回归契约"。

| 脚本 | 具体作用与产物 |
|---|---|
| @ref run_online_bag_regression.sh "run_online_bag_regression.sh" | 编排在线节点、bag 输入回放和结果录制；可请求保存地图，结束时清理子进程 |
| @ref run_voxel_slam_114_reference.sh "run_voxel_slam_114_reference.sh" | 运行 ROS 1 Voxel-SLAM 参考轨迹实验，管理独立 ROS master、回放、录制及资源/轨迹提取 |
| @ref reproduction/multi_lidar/sany_4livox/run_sany_offline_localization.sh "run_sany_offline_localization.sh" | 先生成或选定 SANY 地图，再运行离线地图定位，保存两阶段产物 |
| @ref reproduction/multi_lidar/sany_4livox/run_phase_a_relocalization_matrix.sh "run_phase_a_relocalization_matrix.sh" | 批量运行场景和初始偏移，汇总重定位接受、恢复及轨迹结果 |
| @ref reproduction/multi_lidar/sany_4livox/create_sany_formal_fault_bags.sh "create_sany_formal_fault_bags.sh" | 批量生成已定义的丢包/中断故障输入和场景合同 |
| @ref reproduction/multi_lidar/sany_4livox/run_sany_formal_fault_matrix.sh "run_sany_formal_fault_matrix.sh" | 逐个运行四个正式故障场景，调用在线合同脚本并汇总结果 |
| @ref reproduction/multi_lidar/sany_4livox/run_sany_online_contract.sh "run_sany_online_contract.sh" | 管理在线前端、输入回放和输出录制，验收五个公开 Topic 及故障恢复 |
| @ref reproduction/formal_report/rerun_m3dgr_fastlivo2_timestamp_fix.sh "rerun_m3dgr_fastlivo2_timestamp_fix.sh" | 按冻结的 M3DGR/FAST-LIVO2 矩阵重跑序列，保存指纹、运行清单和汇总；依赖报告对应的数据与工具布局 |
| @ref reproduction/single_lidar/m3dgr/run_voxel_slam_full_backend.sh "run_voxel_slam_full_backend.sh" | 运行 ROS 1 Voxel-SLAM 完整后端对比实验；管理 roscore、回放、结果和资源记录 |

@anchor sany_lidar_launcher_guide
# run_sany_lidar_loc.sh：SANY LiDAR 完整入口

原 SANY 诊断执行器的职责已合入此脚本。现在的调用链为：

```text
${SANY_WS}/run.sh：集中路径配置、固定定位方式和发布预设
  → scripts/run_sany_lidar_loc.sh：默认参数、预检、采集、进程管理
      → ros2 run lightning_lm run_loc_online
          → LocSystem → LIO → 地图定位 → PGO → 业务输出
```

脚本直接启动 C++ `run_loc_online`；通用 `run_loc_online.sh` 是另一条联调路径。录包、计时、资源采样和静默快照只影响运行管理与证据采集。watchdog 检测到 `/localization/pose_vel` 静默时保存快照，定位恢复由进程内重定位器负责。

## 参数和默认值

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

## production / diagnostic

| 开关 | diagnostic 缺省 | production 缺省 | 作用 |
|---|---|---|---|
| `LIGHTNING_LM_COMPUTE_PROFILE` | `1` | `0` | 算法热点计时，退出后提取 `COMPUTE_BENCH_*` |
| `LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD` | `0` | `1` | 裁剪非必要的高频诊断 I/O |
| `LIGHTNING_LM_RESOURCE_PROFILE` | `1` | `0` | 对算法进程组进行 CPU、内存、I/O 等采样 |
| `LIGHTNING_LM_CAUSAL_TRACE` | `1` | `0` | 采集有限的异步速度、队列与锁事件 |

`LIGHTNING_LM_RUN_MODE` 缺省 `diagnostic`；显式子开关优先于模式默认值。录包与 watchdog 有独立开关，切到 diagnostic 不会自动开启它们。安装了 `tegrastats` 时，脚本还会记录 Orin 资源信息。

计算资源覆盖包括 `LIGHTNING_LM_LIO_THREADS`、`LIGHTNING_LM_NDT_THREADS`、`LIGHTNING_LM_NDT_MAX_POINTS`、`LIGHTNING_LM_SOLID_ICP_WORKERS`、`LIGHTNING_LM_SOLID_WORKER_NICE`。`LIGHTNING_LM_CPU_AFFINITY` 限定算法进程的核心，`LIGHTNING_LM_SOLID_CPU_AFFINITY` 限定 SOLiD worker 且须为前者有效核心的子集。`LIGHTNING_LM_RESOURCE_INTERVAL_SEC` 缺省 1 秒，范围 0.5–60 秒。实际预算同时受 YAML `compute_budget` 控制，见 @ref configuration "配置入口"。

## 实际执行顺序

1. 导出 SANY 默认设置，读取模式和布局；进入仓库目录，检查配置、地图路径、环境脚本和开关。
2. 加载 ROS 与定位工作区环境，检查所需命令和 CPU 设置；录包时检查 MCAP/QoS。
3. 从 YAML 解析 LiDAR/IMU 输入；建立本次独立运行目录，保存配置和计算预算摘要。
4. 等待必需输入出现并检查消息类型；记录 Git 身份、运行开关和启动时 ROS 图。
5. 按开关启动 MCAP、静默 watchdog；安装了 `tegrastats` 时启动其采样。
6. 使用 `setsid` 为 C++ 定位程序建立独立进程组，重定向 stdout/stderr，传入配置、地图和三份轨迹路径。
7. 资源采样开启时，在算法启动后开始监测其进程组；主脚本等待算法退出。
8. 正常退出或收到停止信号后，停止辅助进程、等待产物刷新并提取计时记录；正常结束沿用算法退出状态，INT/TERM 路径返回 130。

运行目录预先存在时拒绝启动，以免覆盖证据。算法 cwd 是仓库目录，配置/地图/运行根要求绝对路径；轨迹路径由脚本明确指向本次运行目录。

## 输出文件应该怎样读

| 产物 | 用途 |
|---|---|
| `config.yaml`、`launch_compute_config.json` | 实际配置快照和计算预算摘要 |
| `run_metadata.txt`、`git_status.txt`、`topics_at_start.txt` | 本次版本、开关、输入图和退出状态 |
| `logs/run_loc_online.stdout.log` / `.stderr.log` | 算法原始输出；glog 和热点计时主要在 stderr |
| `results/trajectory_global.tum` | PGO 全局校正轨迹 |
| `results/trajectory_high_frequency.tum` | 定位器内部高频轨迹 |
| `results/trajectory_published_rear_axle.tum` | 实际发布边界后的后轴轨迹 |
| `results/compute_profile.log` | 从 stderr 提取的计时记录；采集开关决定内容 |
| `results/process_resources.jsonl` / `causal_trace.csv` | 开启对应采样时保存的资源与事件记录 |
| `bag/` | 开启录包时保存的 MCAP |
| `snapshots/`、`logs/pose_vel_watchdog.csv` | 开启 watchdog 时的静默现场与事件时间线 |

从 @ref run_sany_lidar_loc.sh "脚本源码" 的 `algorithm_args` 和 `setsid` 进入 @ref online_localization_flow "在线定位端到端流程"，随后阅读 @ref run_loc_online.cc "C++ main" 和 `LocSystem::Init / Spin`。

@anchor sany_ins_launcher_guide
# run_sany_ins_only.sh：SANY 纯组合导航入口

此脚本使用 CGI-430 已解算的位置、姿态和速度。它确认 YAML 的 `system.localization_mode: ins_only` 后启动同一个 C++ `run_loc_online`，程序在 `main` 中选择 `InsLocSystem` 分支。坐标转换、质量门控、点云去畸变和输出审核由 C++ 实现，见 @ref ins_only_operation "纯组合导航操作契约"。

```bash
bash scripts/run_sany_ins_only.sh /absolute/site.yaml [run_loc_online flags]
```

第一个参数是必需的现场 YAML，先通过 `realpath` 转为绝对路径，消费后其余参数直接传给程序。`LIGHTNING_LM_INSTALL_SETUP` 可指定安装环境；未设置时根据仓库位置推导工作区，`LIGHTNING_LM_BUILD_WS` 可覆盖推导值。ROS 环境固定加载 `/opt/ros/humble/setup.bash`；`ROS_DOMAIN_ID` 缺省 42。

执行顺序是：解析配置与工作区 → 加载 ROS/overlay → 设置 domain → 核对 YAML 模式 → `cd` 仓库 → `exec ros2 run lightning_lm run_loc_online --config ...`。此脚本不创建 LiDAR 入口的整套运行目录，也不管理录包/资源监测子进程；节点按 YAML 的 `audit_root` 保存自己的审核记录。

CGI 与所需雷达驱动分别启动。配置固定 ENU 原点、高程类型、后轴外参和质量条件后，才可用于现场。支持 `--output_published_tum`；INS 分支拒绝 LIO/PGO 轨迹参数和内嵌 `--bag`。回放时另用 @ref replay_ins_only.sh "输入回放脚本"，录包另用 @ref record_ins_only.sh "导航录包脚本"。

源码见 @ref run_sany_ins_only.sh "run_sany_ins_only.sh"。

@anchor generic_online_launcher_guide
# run_loc_online.sh：通用在线定位联调入口

此脚本提供少量通用运行管理：加载环境、选择配置、建立日志目录，然后把算法 flags 传给 C++。适用于实时话题联调、外部 `rosbag play`，也可传入程序支持的 `--bag`。它不带 SANY 的传感器预检、MCAP 录制和辅助进程管理。

```bash
LIGHTNING_LM_CONFIG=/absolute/localization.yaml \
LIGHTNING_LM_OUT_ROOT=/absolute/runs \
bash scripts/run_loc_online.sh online_test --map /absolute/map
```

可选第一个非 `-` 开头参数作为 run name，剩余参数原样透传。配置缺省 `$repo_dir/config/default.yaml`，输出根缺省 `$repo_dir/runs`，run name 可由 `LIGHTNING_LM_RUN_NAME` 设置；位置参数优先。ROS/安装环境由 `LIGHTNING_LM_ROS_SETUP` 和 `LIGHTNING_LM_INSTALL_SETUP` 指定。`-h / --help` 先加载环境，再调用真实安装二进制的帮助。

执行顺序是：加载环境 → 解析 run name → 检查并规范配置路径 → 创建日志目录和基本运行信息 → 在子 shell 中 `cd $run_dir` → 启动 C++。因此，透传的相对地图、bag 和输出路径相对于 run dir 解析。脚本不会拒绝已有同名目录，重复 run name 会截断基本日志和运行信息。

基本产物为 `run_metadata.txt`、`logs/run_loc_online.stdout.log` 和 `.stderr.log`；轨迹等额外产物由算法参数决定。脚本没有额外的辅助进程收尾逻辑，沿用算法调用的退出状态。完整参数表见 @ref script_contracts "通用运行契约"，源码见 @ref run_loc_online.sh "run_loc_online.sh"。
