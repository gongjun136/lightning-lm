@page online_localization_flow 在线定位端到端流程

# 阅读目标与边界

此页描述缺省 `lidar` 分支。启动配置 `system.localization_mode: ins_only` 时，入口改为独立 `InsLocSystem`，流程见 @ref ins_only_operation "CGI-430 导航解到业务输出"。

本页沿生产入口 `run_loc_online` 跟踪一帧 IMU/LiDAR 数据，直到地图匹配、定位 PGO 和 ROS 2 输出。优先理解输入怎样形成 LIO 状态、地图匹配如何接受、PGO 怎样融合以及最终输出为何被门控。算法公式和类成员细节留在 @ref laser_mapping_module "LaserMapping" 与 @ref localization_module "定位模块"，关键计算的入口见下表。启动对象图与完整关闭顺序放在本页末尾的 @ref online_loc_runtime_reference "按需参考"。

> 范围说明：本文描述 `src/app/run_loc_online.cc` 的在线模式。`--bag` 只是在进程内增加一个回放生产者，不会把算法切换为离线模式。标为 **推断** 的结论来自调用关系，尚未由运行时故障注入验证。

# 一眼看懂主路径

@dot
digraph online_loc_overview {
  rankdir=LR;
  node [shape=box, fontsize=10];
  sany [label="run_sany_lidar_loc.sh\n默认设置、预检、采集、进程管理"];
  generic [label="run_loc_online.sh\n通用环境、cwd、日志"];
  main [label="main\nLocSystem"];
  ros [label="ROS callbacks\nIMU / LiDAR"];
  sensor [label="sensor_proc_\n有序传感器队列"];
  lio [label="LaserMapping\nLIO"];
  locq [label="lidar_loc_proc_cloud_\n最新定位帧"];
  match [label="LidarLoc\n地图匹配/重定位"];
  pgo [label="localization PGO\n融合/高频外推"];
  output [label="LocSystem\nTF / pose / cloud / health"];
  sany -> main;
  generic -> main [label="通用启动路径", style=dashed];
  main -> ros -> sensor -> lio -> locq -> match -> pgo -> output;
  wheel [label="wheel-speed callback"];
  wheel -> lio [label="direct observation", style=dashed];
}
@enddot

最短阅读链如下：

1. @ref run_sany_lidar_loc.sh "run_sany_lidar_loc.sh"：确定 SANY 参数、输入预检、运行目录和启动 flags。@ref run_loc_online.sh "run_loc_online.sh" 是通用联调入口，不在这条现场调用链中。
2. @ref run_loc_online.cc "run_loc_online.cc" 的 `main()`：创建 `LocSystem`、设置初始位姿、进入 ROS spin、最后保存轨迹并关闭。
3. @ref lightning::LocSystem::Init() "LocSystem::Init()"：加载地图/输出约束，建立订阅、发布和回调桥接。
4. @ref lightning::loc::Localization::ProcessIMUMsg() "ProcessIMUMsg()" 与 @ref lightning::loc::Localization::ProcessLidarMsg() "ProcessLidarMsg()"：把输入统一送入传感器队列。
5. @ref lightning::loc::Localization "Localization" 的 `ProcessSensorInput()`、`LidarOdomProcCloud()` 与 `DrainLioOutputs()`：保持时序并驱动 LIO。
6. @ref lightning::loc::LidarLoc::ProcessCloud() "LidarLoc::ProcessCloud()"：在活动地图上跟踪，必要时转入全局重定位。
7. `Localization::LidarLocProcCloud()`：将定位结果交给 PGO，并通过回调回到 `LocSystem` 发布。

# 学习重点：关键计算与代码入口

| 单帧中的计算 | 应理解的关系 | 实现与公式专题 |
|---|---|---|
| 组帧与去畸变 | 多雷达统一到主雷达，点云统一到扫描末时刻 | `LaserMapping::SyncPackages`、`ImuProcess::UndistortPcl`；@ref imu_deskew "逐点变换公式" |
| IMU 预测与 LiDAR 更新 | 18 维状态预测，点面残差构造六维信息，再迭代注入 | `ESKF::Predict/Update`、`LaserMapping::ObsModel`；@ref eskf_theory "ESKF"、@ref lidar_residuals "残差与参数化" |
| 地图匹配与重定位 | LIO 提供初值，局部配准或检索候选提供地图约束 | `LidarLoc::Align/TryGlobalRelocalization`；@ref relocalization_theory "检索、粗配准与确认" |
| 定位 PGO 与外推 | 地图绝对约束与相对运动进入融合结果 | `PGO::ProcessLidarLoc/PubResult`；@ref pose_graph_theory "图约束与平滑" |
| 业务坐标与速度 | 主传感器到后轴的位姿与杆臂补偿，再判断发布许可 | `LocSystem::PublishLocalizationResult`；@ref output_contracts "参考点与输出门控" |

先沿下面的输入、LIO、匹配、PGO 与输出顺序读一遍，再进入相关模块推导。队列只重点看顺序、积压与新鲜度；逐个对象的销毁细节在修改线程或退出逻辑时回查。

# 1. 启动入口

## SANY 现场启动链

```text
${SANY_WS}/run.sh：加载集中路径配置，选择 lidar 分支
  → scripts/run_sany_lidar_loc.sh：默认设置、cd 仓库、预检与可选采集
  → ros2 run lightning_lm run_loc_online <config/map/trajectory flags>
  → src/app/run_loc_online.cc 的 main()
```

@ref run_sany_lidar_loc.sh "SANY LiDAR 入口" 内部直接组装 `algorithm_args=(ros2 run lightning_lm run_loc_online ...)`，不调用 `run_loc_online.sh`。`run_loc_online` 是由 `src/app/run_loc_online.cc` 编译并安装的 C++ 可执行程序，和同名 `.sh` 文件需要区分。

执行器先加载 ROS/工作区环境、校验配置并解析所需输入，建立运行目录及配置快照，再等待传感器输入并核对类型。可选录包与 watchdog 在算法启动前建立；它通过 `setsid` 在后台启动算法，再按设置挂接资源监控并 `wait` 算法进程。算法退出后记录退出状态、提取计算计时并清理子进程。因此预检、启动、等待、收尾都发生在执行器的一次运行中。

这条路径继承仓库入口设置的 cwd；执行器不将算法 cwd 切到 run dir，而是显式传入配置、地图和轨迹输出路径。录包、watchdog 和资源采样按各自开关管理，具体默认值见 @ref script_contracts "脚本契约"。

**代码依据：** @ref run_sany_lidar_loc.sh "SANY LiDAR 入口"（默认设置、`cd`、`algorithm_args`、`setsid`、`wait` 和收尾）；`src/app/CMakeLists.txt` 的 `add_executable(run_loc_online ...)` 与安装目标。

## 通用联调启动路径

@ref run_loc_online.sh "run_loc_online.sh" 提供另一条可复现运行边界，最终启动同一个 C++ 程序：

```text
调用者
  └─ source ROS 2 setup
  └─ source 工作区 install/setup.bash
  └─ 创建 runs/<run_name>/logs
  └─ cd runs/<run_name>
  └─ ros2 run lightning_lm run_loc_online --config <absolute-path> <extra-args>
```

- `LIGHTNING_LM_REPO_DIR` 决定仓库默认位置；`LIGHTNING_LM_ROS_SETUP` 与 `LIGHTNING_LM_INSTALL_SETUP` 决定两个被 source 的环境。
- `LIGHTNING_LM_CONFIG` 默认指向仓库 `config/default.yaml`，脚本在 `cd` 前将其转成绝对路径。
- `LIGHTNING_LM_OUT_ROOT` 与 `LIGHTNING_LM_RUN_NAME` 决定持久运行目录；第一个不以 `-` 开头的位置参数可覆盖 run name，剩余参数原样传给程序。
- 程序在 run dir 内执行，所以二进制 flags 中仍为相对路径的产物会相对于该目录生成。stdout/stderr 分别重定向到两个日志文件。
- 脚本创建目录和元数据，但不播放 bag、不等待传感器，也不删除旧 run dir；同名运行可能覆盖日志和 metadata。

完整变量表、失败条件与副作用见 @ref script_contracts "Shell 脚本契约"。

**代码依据：** @ref run_loc_online.sh "scripts/run_loc_online.sh"（`source`、`realpath`、`mkdir`、`cd`、`ros2 run` 与重定向语句）

## 进入算法程序

`main()` 初始化系统、地图和初始位姿后才进入正常输入处理；`LocSystem::SetInitPose()` 会打开 `loc_started_`。完整初始化和对象关系见 @ref online_loc_runtime_reference "按需参考"。

# 2. IMU 与 LiDAR 输入路径

## IMU

`LocSystem` 将 `sensor_msgs::msg::Imu` 转为内部 @ref lightning::IMU "IMU"，记录输入统计；只有 `loc_started_` 为 true 才调用 `Localization::ProcessIMUMsg()`。在线模式把 IMU 放进统一 sensor queue。worker 中 `ProcessIMUData()` 拒绝空值、非有限值、非正时间戳、重复时间戳与时间回退，再在 `processing_mutex_` 内交给 LIO，并驱动高频状态输出。

## LiDAR

单雷达模式使用 `common.lidar_topic`/`common.livox_lidar_topic`，多雷达模式按 @ref lightning::MultiLidarConfig "MultiLidarConfig" 为每个传感器建立订阅。`LocSystem` 先做输入时间戳门控；`Localization::ProcessLidarMsg()` 完成消息转换和点云预处理，再把点云及 `lidar_id` 放入 sensor queue。

worker 在 `LidarOdomProcCloud()` 中调用 `LaserMapping::ProcessPointCloud2()`，随后 `DrainLioOutputs()` 循环调用 @ref lightning::LaserMapping::RunDetailed() "RunDetailed()"，直到返回 `kNoData`。一次 ROS 点云输入不保证恰好产生一个 LIO 输出；内部同步/多雷达融合可能暂存或一次排出多个结果。

若配置了 sensor lag 阈值，输入积压会触发 LiDAR admission throttling；统一队列本身也会在容量溢出时丢最旧消息。两者的计数与当前 lag 会出现在 pipeline diagnostics 中。

## 轮速

轮速不进入统一 sensor queue。`LocSystem::ObserveWheelSpeedInput()` 完成 rpm 到 m/s 换算和输入统计后，直接调用 `Localization::ProcessWheelSpeed()`；后者在 lifecycle 共享锁下把观测交给 LIO，并在启用 IMU static hold 时更新独立的短时轮速历史。它因此是 LIO/静止判定的辅助观测，不是触发 `LidarLoc::ProcessCloud()` 的主帧。

**代码依据：** `src/core/system/loc_system.cc` 的订阅 lambda、`ProcessIMU/ProcessLidar()` 与 `ObserveWheelSpeedInput()`；`src/core/localization/localization.cpp` 的 `ProcessLidarMsg()`、`ProcessIMUMsg()`、`ProcessIMUData()`、`ProcessWheelSpeed()`、`LidarOdomProcCloud()`

# 3. LIO 输出到地图匹配

每次 `RunDetailed()` 产生 `kOutput` 时，`DrainLioOutputs()` 执行以下顺序：

1. 读取 LIO @ref lightning::NavState "NavState"，送入 `LidarLoc::ProcessLO()`，同时送入定位 PGO 的 `ProcessLidarOdom()`。
2. 取得定位用投影点云、发布点云、同帧多雷达统计和到达时间，组装成一个 `LidarLocInput`；这些字段作为一体传递，避免跨帧拼接。
3. 若 `lidar_loc.loc_on_kf=true`，只有 LIO 关键帧变化才提交地图匹配；否则每个 LIO 输出都提交。
4. 在线模式入 `lidar_loc_proc_cloud_`；队列默认容量为 1，重定位较慢时新帧会淘汰旧的待处理帧，以恢复到最新观测。

`LidarLoc::ProcessLO()` 自己维护有界 LO pose queue，并根据 LIO 状态更新可靠性滞回；`Align()` 按点云结束时刻从该队列分配对应 LO pose。时间对不上时会记录警告，而不是把任意最新姿态静默当成同帧先验。

**代码依据：** `src/core/localization/localization.cpp` 的 `DrainLioOutputs()`；`src/core/localization/lidar_loc/lidar_loc.cc` 的 `ProcessLO()`、`Align()`、`AssignLOPose()`

# 4. 跟踪、失锁与全局重定位

@dot
digraph localization_state {
  rankdir=LR;
  node [shape=ellipse];
  init [label="INITIALIZING\nloc_inited_=false"];
  track [label="GOOD\nlocal match accepted"];
  dr [label="FOLLOWING_DR\nmatch rejected"];
  relocal [label="global relocalization\nquery/refine/confirm"];
  init -> track [label="external/FP/global init accepted"];
  init -> relocal [label="no local initialization"];
  relocal -> track [label="candidate accepted"];
  track -> dr [label="NDT/consistency rejected"];
  dr -> track [label="next local match accepted"];
  dr -> relocal [label="consecutive failures >= threshold\nand backend ready"];
}
@enddot

`LidarLoc::Align()` 的主要分支是：

- 静止且已初始化：复用上一绝对位姿；按 `enable_parking_static` 决定有效状态。
- 未初始化：先尝试外部初始位姿，再按配置尝试功能点，最后轮询全局重定位；仍失败则保持 `INITIALIZING` 并返回。
- 已初始化：用相邻 LO pose 的增量外推本帧初值，加载该位姿附近 tile，执行局部 NDT；阈值区分 `min_init_confidence` 与 `min_tracking_confidence`。
- 局部结果还需通过 odometry delta 检查，以及可选全扫描地图一致性检查。拒绝结果不作为下一帧先验，而是输出 LO/DR 外推并标记 `FOLLOWING_DR`。
- 连续拒绝达到 `relocalization.lost_frame_threshold` 且全局 backend ready 时，将 `loc_inited_` 复位并重置 query，下一阶段转入全局重定位。

全局重定位通过 @ref lightning::loc::GlobalRelocalizer "GlobalRelocalizer" 统一接口选择 BTC、SOLiD 或 SOLiD-KISS 后端。候选不是直接接受：代码会执行候选预检、NDT/plane-ICP refinement、地图一致性和连续确认。接受后更新 `map→odom` 关系并在 `Localization::LidarLocProcCloud()` 中重置 PGO，避免沿用旧图优化状态。

`LocalizationStatus` 的源码含义为：`IDLE` 未定位、`INITIALIZING` 位姿无效、`GOOD` 正常定位、`FOLLOWING_DR` 异常匹配时跟随 DR、`FAIL` 定位失败。`lidar_loc_valid_` 属于 LiDAR 匹配结果，最终对外 `valid_` 由融合/DR 输出路径负责；不能仅凭一个 bool 判断整个系统健康。

**代码依据：** `src/core/localization/lidar_loc/lidar_loc.cc` 的 `Align()`、`Localize()` 与 `TryGlobalRelocalization()`；`src/core/localization/global_relocalizer.h`；`src/core/localization/localization_result.h`

# 5. 地图加载与更新策略

@ref lightning::TiledMap "TiledMap" 根据地图目录索引管理活动 tile。初始化/重定位可用 `LoadOnPose()` 载入候选附近地图；正常跟踪也按 LO guess 载入附近 tile。`maps.load_map_size` 与 `maps.unload_map_size` 控制活动范围，地图策略决定动态点云是否只在内存、短期/长期维护或持久化。

`LidarLoc::Init()` 启动 `update_map_thread_`。它与匹配共享地图/NDT 资源，因此调用者不应绕过现有接口直接替换活动地图。`Finish()` 先要求线程退出并 join；只有策略为 `PERSISTENT`、允许退出保存且没有外部 set-pose 条件时，才调用 `SaveToBin(true)`。因此“在线定位一定只读地图”不是普遍事实，是否落盘取决于地图策略配置。

**代码依据：** `src/core/localization/lidar_loc/lidar_loc.cc` 的 `Init()`、`UpdateGlobalMap()`、`UpdateMapThread()`、`Finish()`；`src/core/maps/tiled_map.cc` 的 `LoadMapIndex()` 与 `LoadOnPose()`

# 6. PGO、发布与输出门控

`LidarLocProcCloud()` 取得地图匹配结果后先更新运行统计。全局重定位成功会重置 PGO；静止地图冲突会把结果改为 `FAIL`、重置 PGO 并请求重新定位；正常结果进入 `PGO::ProcessLidarLoc()`。LIO 状态早已通过 `PGO::ProcessLidarOdom()` 进入同一融合器。

PGO 的高频结果进入容量 1 的输出队列。live output callback 用独立时间戳门拒绝非单调结果，再更新 `loc_result_`、TF/UI 和定位结果回调。全局优化结果走独立 callback，避免较旧全局结果让实时 ROS 输出时间倒退。

`LocSystem` 的发布边界还包含：

- 位姿：`/PosRes`、`/slamPoseRaw_topic`、`/localization/pose_vel`，以及可选 TF。
- 点云：`/LidarDataInv`、`/LidarDataInL`。
- 健康与可观测性：`/localization/fault_status`、`/localization/loc_status`、`/localization/pipeline_diagnostics`、`/localization/debug_message`、`/localization/path`。
- 100 ms 健康 timer 与 2 s path timer。

publication gate 综合“是否曾有 GOOD 结果”、连续失配、LiDAR 匹配新鲜度和时间戳等条件，决定 map frame 输出是否继续发布。因此看到 LIO/DR 内部仍运行，并不意味着外部 map 定位输出一定被放行。

固定地图变换在输出边界应用；`output.fixed_map_transform.target_frame` 必须与 `output.map_frame` 一致。不要在算法内部再叠加一次同类变换。

**代码依据：** `src/core/localization/localization.cpp` 的 `LidarLocProcCloud()` 与 `Init()` 输出回调；`src/core/system/loc_system.cc` 的 `Init()`、`PublishLocalizationResult()`、`PublishProcessedCloud()`、`PublishHealthStatus()`

# 7. 队列、时序与新鲜度

| 执行上下文 | 创建者 | 处理内容 | 容量/退出语义 |
|---|---|---|---|
| ROS executor | `LocSystem::Spin()` | subscription、timer 与发布回调入口 | 并发度取决于 `rclcpp::spin(node_)` 的 executor 实现 |
| `sensor_proc_` worker | `Localization::Init()` | IMU 与 LiDAR 的统一 `SensorInput` | 默认 10000；最旧待处理消息先丢；退出时 drain 后 join |
| `lidar_loc_proc_cloud_` worker | 同上 | LIO 投影点云的地图匹配 | 默认 1；慢定位时保留最新帧；可按配置跳帧；退出时 drain 后 join |
| `high_frequency_output_proc_` worker | 同上 | PGO 高频结果到 live output callback | 固定 1；只保留最新待发布结果；退出时 drain 后 join |
| `update_map_thread_` | `LidarLoc::Init()` | 根据当前位姿更新活动 tile/NDT 地图 | `LidarLoc::Finish()` 置退出标志并 join |
| 可选 bag thread | `main()` | 进程内 ROS 2 bag 回放 | 播放后等待并 shutdown；主线程 join |

@ref lightning::sys::AsyncMessageProcess "AsyncMessageProcess" 每次只从有界 deque 取一条消息到 worker 外执行。入队超过 `max_size_` 时从队头删除，因此丢弃的是最旧的尚未执行消息；`Quit()` 设置退出标志，但 `ProcLoop()` 仅在“退出且队列为空”时结束，故会 drain 已接收队列。正在执行的消息计入 pending，丢弃数和完成数可进入运行诊断。

统一 sensor worker 是正确性的组成部分：源码明确要求 LIO 按回调顺序观察 IMU 与 LiDAR；若把二者放到两个 worker，后到 IMU 可能越过点云，在线结果会偏离离线顺序。`processing_mutex_` 再串行保护实际 LIO 更新；`lifecycle_mutex_` 防止 `Finish()` 在算法仍被访问时释放成员。

**推断：** `rclcpp::spin(node_)` 当前通常使用单线程 executor，但源码没有把“所有 ROS 回调永久串行”声明为公共契约；本文只依赖 `Localization` 内部队列与锁提供的顺序保证。

**代码依据：** `src/core/system/async_message_process.h` 的 `Start()`、`AddMessage()`、`ProcLoop()`、`Quit()`；`src/core/localization/localization.cpp` 的 `Init()`、`ProcessSensorInput()` 与 `LidarOdomProcCloud()`

# 8. 配置索引：从现象回到代码

| 关注点 | 主要配置/环境 | 首个读取位置 |
|---|---|---|
| 地图与初始位姿 | `--map`, `system.map_path`, `online_localization.initial_pose` | `main()`、`LocSystem::Init()` |
| 输入 topic/QoS | `common.*_topic`, `system.lidar_qos`, multi-lidar 配置 | `LocSystem::Init()`、`LaserMapping::Init()` |
| 统一队列积压 | `system.online_sensor_queue_size`, `online_sensor_*lag_sec` | `Localization::Init()` |
| 定位队列新鲜度 | `system.online_lidar_loc_queue_size`, `enable_lidar_loc_skip`, `lidar_loc.loc_on_kf` | `Localization::Init()` |
| 局部匹配 | `lidar_loc.min_init_confidence`, `min_tracking_confidence`, `force_2d` | `LidarLoc::Init()` |
| 失锁/全局重定位 | `relocalization.*` 及 backend 选择 | `LidarLoc::Init()` |
| 地图活动范围/写回 | `maps.load_map_size`, `unload_map_size`, `dyn_cloud_policy` | `LidarLoc::Init()` / `TiledMap` |
| 发布门控 | `relocalization.lost_frame_threshold`, `system.localization_output_max_lidar_age_sec` | `LocSystem::Init()` |
| 输出坐标 | `output.map_frame`, `output.fixed_map_transform`, `primary_lidar_position_in_body` | `LocSystem::Init()` |
| 轮速观测 | `system.enable_wheel_speed_observation`, `SANY_ENABLE_CAN_OBSERVATION`, `SANY_WHEEL_SPEED_TOPIC` | `LocSystem::Init()` |

表只提供定位入口；默认值、单位和组合约束以相应读取代码与 `config/default.yaml` 为准。

# 9. 建议的源码阅读顺序

按下面的重点建立公式、算法步骤与代码的对应关系：

1. **单帧数据流：** `ProcessLidarMsg` → `ProcessSensorInput` → `LidarOdomProcCloud` → `DrainLioOutputs` → `LidarLocProcCloud`。
2. **前端计算：** `ImuProcess::UndistortPcl` → `ESKF::Predict/Update` ↔ `LaserMapping::ObsModel`；结合去畸变、预测和点面残差公式阅读。
3. **匹配与恢复：** `LidarLoc::Align` → `Localize` → `TryGlobalRelocalization` → `LocalizationResult`；区分跟踪失败、候选生成和最终接受。
4. **融合与输出：** PGO output callback → `PublishLocalizationResult` → publication gate → health/diagnostics；确认参考点、时间和有效性。

适合下断点或日志对照的符号：

- `LocSystem::ProcessLidar()`：确认 topic、lidar id 与输入时间戳是否进入系统。
- `Localization::ProcessSensorInput()`：确认统一队列顺序和等待时间。
- `Localization::DrainLioOutputs()`：确认 LIO 是否真正形成输出/关键帧。
- `LidarLoc::Align()`：观察 `loc_inited_`、`match_fail_count_` 和 initial/tracking 分支。
- `Localization::LidarLocProcCloud()`：观察匹配结果、重定位接受与 PGO 交接。
- `LocSystem::PublishLocalizationResult()`：确认最终门控与 ROS 输出。

下一步按 @ref guide_lio "LIO 前端"、@ref guide_localization "地图定位"、@ref guide_map_output "地图坐标与输出" 深入。需要修改启动或并发时，再查下面的初始化、所有权与关闭约束。

@anchor online_loc_runtime_reference

# 按需参考：启动、对象关系与退出

以下内容保留初始化、线程停止与对象管理的实现约束。第一次学习先完成上面的算法与数据主线；修改重初始化、并发或关闭行为时再核对本节。

## `main()` 生命周期

@dot
digraph online_loc_lifecycle {
  rankdir=TB;
  node [shape=box];
  init [label="rclcpp::init"];
  object [label="stack LocSystem loc"];
  system_init [label="loc.Init(config, map)"];
  pose [label="ReadInitialPose\nloc.SetInitPose"];
  bag [label="optional OnlineBagPlayer thread"];
  spin [label="loc.Spin"];
  finish [label="join bag thread\nloc.Finish"];
  save [label="optional TUM export\nrclcpp::shutdown"];
  init -> object -> system_init -> pose -> bag -> spin -> finish -> save;
}
@enddot

`main()` 在栈上持有 `LocSystem`。`Init()` 成功后才读取初始位姿；配置优先使用 `online_localization.initial_pose`，不存在时回退到 `offline_localization.initial_pose`，若均不可用则使用单位位姿。`SetInitPose()` 不只是保存位姿：它还把 `loc_started_` 置为 true，因此此前到达的传感器回调只做输入观测，不进入定位算法。

若传入 `--bag`，@ref lightning::OnlineBagPlayer "OnlineBagPlayer" 在额外线程中按 sensor time 发布消息，播放完成并等待 `--post_wait_seconds` 后调用 `rclcpp::shutdown()`，从而让主线程的 `Spin()` 返回。未传 `--bag` 时，进程持续 spin，直到外部 shutdown/信号结束 ROS 上下文。

`LocSystem` 析构函数还会调用 `Finish()`；`finished_` 让显式调用与析构调用保持幂等。

**代码依据：** @ref run_loc_online.cc "run_loc_online.cc" 的 `ReadInitialPose()` 与 `main()`；@ref lightning::LocSystem::SetInitPose() "LocSystem::SetInitPose()"、@ref lightning::LocSystem::Finish() "LocSystem::Finish()"

## 初始化与对象所有权

@dot
digraph online_loc_ownership {
  rankdir=LR;
  node [shape=box];
  main [label="main\nstack"];
  system [label="LocSystem"];
  ros_node [label="rclcpp::Node\npubs/subs/timers"];
  loc [label="Localization"];
  lio [label="LaserMapping"];
  lidar_loc [label="LidarLoc"];
  pgo [label="PGO"];
  map [label="TiledMap"];
  relocalizer [label="GlobalRelocalizer\nBTC / SOLiD"];
  main -> system;
  system -> ros_node;
  system -> loc;
  loc -> lio;
  loc -> lidar_loc;
  loc -> pgo;
  lidar_loc -> map;
  lidar_loc -> relocalizer;
}
@enddot

`LocSystem::Init()` 先创建在线模式的 `Localization`，解析 YAML，再创建 ROS 节点。地图路径采用 `--map` 非空值，否则采用 `system.map_path`；二者都没有时初始化失败。输出 map frame、固定地图变换、发布门控、LiDAR QoS、轮速输入和健康状态也在这一层验证。

`Localization::Init()` 在生命周期独占锁内建立三个算法对象：

- @ref lightning::LaserMapping "LaserMapping"：`is_in_slam_mode_=false`，只提供 LIO、去畸变点云与关键帧，不在此路径负责建图后端。
- @ref lightning::loc::LidarLoc "LidarLoc"：接收全局地图目录，加载 tile 地图与全局重定位数据库。
- @ref lightning::loc::PGO "PGO"：接收 LIO 与 LiDAR 定位约束，产生融合/外推结果。

随后 `LocSystem` 建立 IMU、单/多 LiDAR、可选 Livox custom message 和轮速订阅；建立位姿、点云、状态、诊断、debug message 与 path 发布器；最后把 `Localization` 的 TF、定位结果、全局结果和已处理点云回调连接回系统层。

**代码依据：** @ref lightning::LocSystem::Init() "LocSystem::Init()"；@ref lightning::loc::Localization::Init() "Localization::Init()"；`src/core/localization/lidar_loc/lidar_loc.cc` 的 `LidarLoc::Init()`

## 关闭顺序与对象寿命

正常退出顺序是：

1. ROS spin 返回；可选 bag thread 已 join。
2. `LocSystem::Finish()` 清除 debug sink，调用 `Localization::Finish()`。
3. `Localization::Finish()` 依次 `Quit()` sensor、LiDAR 定位、高频输出队列；每个 worker drain 并 join。
4. 取得 lifecycle 独占锁后，`LidarLoc::Finish()` 停止并 join 地图更新线程；随后退出 UI。
5. `LocSystem` 停止 heartbeat、关闭在线 TUM 文件；`main()` 可导出 PGO/高频轨迹并关闭 ROS。

这个顺序保证 worker 不再访问算法对象后才进入模块级清理。`Localization::Init()` 若在已运行对象上重新初始化，会先释放生命周期锁再调用 `Finish()`，避免持有独占锁等待正等待共享锁的 worker。

**关键约束：** 不要把 `LidarLoc::Finish()` 提前到队列 join 之前；否则定位 worker 可能仍在使用地图线程所维护的对象。不要移除 `LocSystem::Finish()` 的幂等保护；显式结束和析构都会进入它。

**代码依据：** @ref lightning::LocSystem::Finish() "LocSystem::Finish()"；@ref lightning::loc::Localization::Finish() "Localization::Finish()"；@ref lightning::loc::LidarLoc::Finish() "LidarLoc::Finish()"
