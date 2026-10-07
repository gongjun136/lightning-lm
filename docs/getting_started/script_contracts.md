@page script_contracts Shell 启动与实验脚本契约

# 读取规则

CGI-430 的 @ref build_ins_only.sh "build_ins_only.sh"、@ref run_sany_ins_only.sh "run_sany_ins_only.sh"、@ref record_ins_only.sh "record_ins_only.sh"、@ref replay_ins_only.sh "replay_ins_only.sh"，以及坐标迁移/测试 Python 工具的参数、环境和文件副作用见 @ref ins_only_operation "纯组合导航操作契约"。

表中的“参数”指脚本自身消费的参数，不包含透传给二进制/子脚本的全部 flags；“环境”只列影响控制流、路径或资源的主要变量；cwd 指脚本是否主动改变工作目录；“副作用”包括文件写入、进程/ROS 图变化和包安装。精确默认值以脚本 `usage()` 与赋值语句为准。

点击表中的脚本名可查看完整源码与行号；顶部“脚本源码”收录仓库内全部 Shell 文件。源码页每次构建直接读取真实脚本。Bash 没有 C++ 式的类/函数 API 索引，阅读其参数、分支和调用链时结合本页契约。

## 脚本用途与阅读说明

@ref guide_scripts "脚本启动主线"：按定位、前端、建图、构建与导航数据、回归及复现实验分组阅读。每份导读先给实际启动命令，再说明执行顺序、输入输出与退出行为；源码页可直接返回对应导读。

## 两层定位入口

现场操作员使用项目根目录的 `${SANY_WS}/run.sh`，该文件负责加载集中路径配置并选择定位模式。它位于 Lightning-LM 仓库之外，具体内容由部署文档维护；本手册中的 @ref run_sany_lidar_loc.sh "scripts/run_sany_lidar_loc.sh" 是仓库内的 LiDAR 分支入口。

```text
${SANY_WS}/run.sh
  ├─ lidar    → scripts/run_sany_lidar_loc.sh → C++ run_loc_online / LocSystem
  └─ ins_only → scripts/run_sany_ins_only.sh → C++ run_loc_online / InsLocSystem
```

`lidar` 使用 LiDAR 地图定位及现有 IMU/CAN 辅助；`ins_only` 直接使用 CGI-430 导航结果。后续预留 `fusion`，但当前程序只接受前两种模式，融合方案和执行器尚未实现。`LIGHTNING_LM_RUN_MODE=production/diagnostic` 控制开销与诊断设置，与定位模式分别管理。

## 运维与通用入口

| 脚本 | 参数 / 主要环境 | cwd 与调用链 | 副作用 |
|---|---|---|---|
| @ref run_sany_lidar_loc.sh "run_sany_lidar_loc.sh" | 可选 run name；`LIGHTNING_LM_WS`, `LIGHTNING_LM_REPO_DIR`, `LIGHTNING_LM_RUN_MODE`, `ROS_DOMAIN_ID`, `SANY_*`, profile 变量 | `cd $REPO`；预检、配置快照 → `setsid ros2 run ... run_loc_online` → 采集与进程组收尾 | 写 logs/results/snapshots/config/metadata；可选录 bag、资源监控和静默快照；拒绝复用已有 run dir |
| @ref run_sany_ins_only.sh "run_sany_ins_only.sh" | 必需现场 YAML；其余算法 flags 透传；`LIGHTNING_LM_INSTALL_SETUP`, `LIGHTNING_LM_BUILD_WS`, `ROS_DOMAIN_ID` | source ROS/overlay；确认 `system.localization_mode: ins_only`；`cd $repo` → `exec ros2 run ... run_loc_online` | 由 `InsLocSystem` 写审核记录和发布结果；驱动、录包与回放分别启动 |
| @ref run_frontend_online.sh "run_frontend_online.sh" | `--config` 透传；`LIGHTNING_LM_{REPO_DIR,ROS_SETUP,INSTALL_SETUP,CONFIG}` | 不改 cwd；source 两个 setup → `ros2 run ... run_frontend_online` | 取决于节点配置：发布 ROS topic/UI/log |
| @ref run_slam_online.sh "run_slam_online.sh" | `--config` + 透传；同上，另有 OUT_ROOT/RUN_NAME | 创建并 `cd` run dir（同名目录可复用）→ `ros2 run ... run_slam_online` | 写 logs/run metadata；在线建图服务可能写地图 |
| @ref run_loc_online.sh "run_loc_online.sh" | `--config` + 透传；同上 | 创建并 `cd` run dir → `ros2 run ... run_loc_online` | 写 logs；加载地图并发布定位结果 |
| @ref align_maps_offline.sh "align_maps_offline.sh" | 全部参数透传 | 不改 cwd；解析安装 binary → `exec align_maps_offline` | 由二进制读写所给地图/结果路径 |
| @ref save_default_map.sh "save_default_map.sh" | 无 | 当前 cwd → `ros2 service call /lightning/save_map` | 请求在线节点写 `new_map` |
| @ref install_dep.sh "install_dep.sh" | 无 | 当前 cwd → `sudo apt install ...` | **系统级**安装 OpenCV/PCL/YAML/glog/gflags/ROS PCL 包 |

**代码依据：** 上述各脚本（`source/cd/exec/ros2/sudo` 调用与变量默认值）

## run_loc_online.sh 详细契约

这条脚本是通用联调的运行边界；SANY 现场完整启动边界见 @ref sany_lidar_launcher_guide "SANY LiDAR 入口详解"。它不实现定位算法，也不自行启动 rosbag；它负责选择环境和配置、建立可归档工作目录，然后把剩余 flags 交给 `run_loc_online`。

### 参数与环境

| 输入 | 默认值/解析 | 影响 |
|---|---|---|
| 第一个非 `-` 开头参数 | 时间戳形式 `run_loc_online_YYYYmmdd_HHMMSS` | 作为 run name；消费后不再传给二进制 |
| 其余参数 | 无 | 原样追加到 `run_loc_online`；可包含 `--map`、`--bag`、轨迹输出等程序 flags |
| `LIGHTNING_LM_REPO_DIR` | 脚本目录的上一级 | 配置、install setup、默认输出根的基准 |
| `LIGHTNING_LM_ROS_SETUP` | `/opt/ros/humble/setup.bash` | 第一个 source 的 ROS 环境 |
| `LIGHTNING_LM_INSTALL_SETUP` | `$repo_dir/install/setup.bash` | 第二个 source 的当前工作区 overlay |
| `LIGHTNING_LM_CONFIG` | `$repo_dir/config/default.yaml` | 传给程序的配置；运行前必须存在并转为绝对路径 |
| `LIGHTNING_LM_OUT_ROOT` | `$repo_dir/runs` | 持久运行目录根；创建后转为绝对路径 |
| `LIGHTNING_LM_RUN_NAME` | 时间戳名称 | 无位置 run name 时生效 |

`-h/--help` 仍会先 source 两个 setup，再通过 `ros2 pkg prefix lightning_lm` 找到安装空间中的真实二进制并调用其 help；help 分支不创建 run dir。

### 工作目录与调用链

```text
caller cwd
  → source $LIGHTNING_LM_ROS_SETUP
  → source $LIGHTNING_LM_INSTALL_SETUP
  → validate + realpath config
  → mkdir $LIGHTNING_LM_OUT_ROOT/$run_name/logs
  → subshell: cd $run_dir
  → ros2 run lightning_lm run_loc_online --config $config_path "$@"
```

主 shell 的 cwd 不变，因为 `cd` 位于子 shell；算法进程的 cwd 是 run dir。脚本只规范 config 和 output root，自行透传的 `--map`、`--bag`、`--output-*` 若使用相对路径，仍由算法进程相对于 run dir 解析。为避免换机器后含义改变，建议这些路径也传绝对路径。

### 文件与副作用

```text
$LIGHTNING_LM_OUT_ROOT/<run_name>/
├── run_metadata.txt
└── logs/
    ├── run_loc_online.stdout.log
    └── run_loc_online.stderr.log
```

`run_metadata.txt` 记录 executable、repo/config/run 路径、额外参数与 ISO 时间；标准输出和错误输出分开重定向。节点还会建立 ROS 订阅/发布、timer 与心跳；若程序 flags 请求 TUM 输出或地图配置允许持久化动态地图，还会产生额外文件。脚本没有“目标已存在则拒绝”的保护，同一个 run name 会截断上述 metadata 和日志。

`set -euo pipefail` 使 setup 不存在、配置不存在、目录创建失败或算法返回非零时脚本失败；配置缺失明确返回 2。算法级初始化失败和运行失败沿用 `ros2 run` 的退出状态。

在线定位内部调用与线程路径见 @ref online_localization_flow "在线定位端到端流程"。

**代码依据：** @ref run_loc_online.sh "scripts/run_loc_online.sh"（全部控制流）；@ref run_loc_online.cc "run_loc_online.cc"（可透传 flags 与算法副作用）

## 标准离线运行器

| 脚本 | 必需参数；重要可选参数 | 主要环境 | cwd / 调用链 / 副作用 |
|---|---|---|---|
| @ref run_frontend_offline.sh "run_frontend_offline.sh" | `--bag --config --output-dir`; sequence/repeat/playback/max frames/CPU/watchdog/UI/各输出 | ROS/INSTALL setup，repo，部分 legacy 输入/输出变量 | 不改主 shell cwd；检查 bag → inspector → `setsid taskset run_frontend_offline` + resource monitor → timing extractor；创建 results/logs、轨迹/PCD/CSV/JSON/metadata；拒绝覆盖受管产物 |
| @ref run_slam_offline.sh "run_slam_offline.sh" | `--bag --config --output-dir`; 另有 map/global/TUM/frame stats/backend-only | `LIGHTNING_LM_{INPUT_BAG,CONFIG,OUT_ROOT,RUN_NAME,WAIT_UI,OUTPUT_TUM,...}` | 同上，调用 `run_slam_offline`；额外写 tiled/global map、backend diagnostics/relocalization DB，并验证 map contract |
| @ref run_loc_offline.sh "run_loc_offline.sh" | `--bag --config --map --output-dir`; initial pose/reference/physical bounds 等 | 同类 setup/repo/bag/config/map/output 变量 | inspector → `run_loc_offline` + monitor → timing + localization analyzer；写轨迹、误差、结果 bag、CSV/JSON；拒绝覆盖受管产物 |
| @ref run_frontend_offline_batch.sh "run_frontend_offline_batch.sh" | `--bag --output-root`；`--config-dir` 或 `--config`；repeats/dry-run | repo | 枚举配置并逐个调用 `run_frontend_offline.sh`；为每组建立 run dir 和汇总，已有同名/失败计数影响退出码 |

这些脚本会用 `setsid` 建立进程组，trap/看门狗可能向整个算法进程组发送 INT、TERM、KILL；会导出 BLAS/OpenMP 线程数。输出目录不是临时目录，属于实验证据。

**代码依据：** `scripts/run_*_offline*.sh`（usage、resolve_output、cleanup_group、monitor 与最终契约判断）

## 在线诊断与回归

SANY 的诊断采集已合入 @ref run_sany_lidar_loc.sh "LiDAR 定位入口"，通过环境开关启用；不再经过独立的 SANY 诊断 Shell 执行器。

| 脚本 | 参数 / 环境 | cwd 与调用链 | 副作用 |
|---|---|---|---|
| @ref run_online_bag_regression.sh "run_online_bag_regression.sh" | bag/config/map/mode 与 playback 时序 | `cd run_dir`；启动在线 node + bag recorder/player；可调用 save_map service 与轨迹提取器 | 写 ROS bag、地图、轨迹和日志；结束时向子进程发信号 |

**代码依据：** 对应脚本（进程启动、trap/cleanup、目录创建与 ROS 命令）

## 复现与对比实验

| 脚本 | 参数/主要环境 | 调用链 | 主要副作用 |
|---|---|---|---|
| @ref reproduction/multi_lidar/sany_4livox/run_sany_offline_localization.sh "reproduction/.../run_sany_offline_localization.sh" | bag、mapping/localization config、map、输出、CPU、watchdog | 先 `run_slam_offline.sh` 生成/选地图，再 `run_loc_offline.sh` | 创建 mapping + localization 两阶段产物 |
| @ref reproduction/multi_lidar/sany_4livox/run_phase_a_relocalization_matrix.sh ".../run_phase_a_relocalization_matrix.sh" | bag/config/map/output matrix | 多数据集/offset 调用 `run_loc_offline.sh` | 批量 run dir、stats、失败汇总 |
| @ref reproduction/multi_lidar/sany_4livox/create_sany_formal_fault_bags.sh ".../create_sany_formal_fault_bags.sh" | 输入/root；scenario/mode/topic/start/duration | source ROS → Python fault-bag creator | 创建裁剪/删改 topic 的 ROS 2 故障 bag |
| @ref reproduction/multi_lidar/sany_4livox/run_sany_formal_fault_matrix.sh ".../run_sany_formal_fault_matrix.sh" | fault bags/config/output | 遍历场景调用在线 contract runner | 批量日志/结果，聚合失败码 |
| @ref reproduction/multi_lidar/sany_4livox/run_sany_online_contract.sh ".../run_sany_online_contract.sh" | bag/config/output/cpuset/fault contract | source ROS+workspace；record → online node → bag play → analyzers | 修改 ROS 图、录制结果 bag、写环境/commit identity、JSON/CSV；清理进程组 |
| @ref reproduction/formal_report/rerun_m3dgr_fastlivo2_timestamp_fix.sh "formal_report/rerun_m3dgr_fastlivo2_timestamp_fix.sh" | 脚本内实验矩阵及已有 runner | 对多个序列/repeat 调用外部/仓库 runner并核验指纹/轨迹 | 大批量报告实验目录与 inventory；依赖本机数据布局 |
| @ref reproduction/single_lidar/m3dgr/run_voxel_slam_full_backend.sh "single_lidar/m3dgr/run_voxel_slam_full_backend.sh" | ROS1 数据/资源环境 | 启动 ROS1 roscore + Voxel-SLAM + rosbag + monitor + extractor | 使用独立 ROS master 端口；写轨迹/PCD/资源与日志；信号清理 |
| @ref run_voxel_slam_114_reference.sh "run_voxel_slam_114_reference.sh" | bag/output、ROS1 workspace、topic/frame/CPU | ROS1 roscore → roslaunch → recorder/monitor → rosbag play | 建 reference 轨迹与 contract；管理 ROS1 进程组 |

**代码依据：** `scripts/reproduction/**/*.sh`、`scripts/run_voxel_slam_114_reference.sh`（循环矩阵、子运行器和输出/清理逻辑）

# 脚本修改检查清单

修改任一 `.sh` 时至少检查：

1. 参数：usage、解析 case、默认值、透传边界是否一致。
2. 环境：新变量是否有默认值、是否 export、是否进入 metadata/config snapshot。
3. cwd：相对路径在哪个 `cd` 前后求值；安装 binary 与 repo 是否可能混用。
4. 调用链：source 顺序、ROS domain/master、子进程组与信号是否正确。
5. 副作用：覆盖保护、目录所有权、临时/持久产物、ROS service/topic、系统安装。
6. 退出：trap 是否保留真实算法退出码；超时/验证失败是否区分。

文档同步判定见 @ref documentation_rules "基于代码 diff 的文档更新规则"。
