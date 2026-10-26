@page localization_module 定位与重定位

# 模块边界

@ref lightning::loc::Localization "Localization" 是定位数据面：它初始化 LIO、@ref lightning::loc::LidarLoc "LidarLoc" 和定位 PGO，接收 IMU/LiDAR/轮速并输出 @ref lightning::loc::LocalizationResult "LocalizationResult"。@ref lightning::LocSystem "LocSystem" 是 ROS 2 编排层，负责订阅、发布、TF、轨迹记录和功能安全心跳；两者不要混为一个对象。

如果目标是从入口系统阅读在线定位，请先看 @ref online_localization_flow "在线定位端到端流程"。本页先交代类边界，再用四个问题核对地图匹配和重定位的实际分支。每段片段均从真实源码生成；读完例子和边界即可返回主线，需要细查时再进入完整 API。

<b>代码依据：</b> `src/core/localization/localization.h`、`src/core/system/loc_system.h`（公开接口、成员所有权）

# 按主链路阅读 API

这些接口页先说明职责与调用关系，复杂入口再列主要步骤；关键代码按需展开。

| 阅读阶段 | 直接入口 |
|---|---|
| ROS 与算法装配 | @ref lightning::LocSystem::Init() "LocSystem::Init()" → @ref lightning::loc::Localization::Init() "Localization::Init()" |
| 接收 IMU 与点云 | @ref lightning::loc::Localization::ProcessIMUMsg() "ProcessIMUMsg()"、@ref lightning::loc::Localization::ProcessLidarMsg() "ProcessLidarMsg()" |
| 驱动前端 | @ref lightning::loc::Localization::ProcessSensorInput() "ProcessSensorInput()" → @ref lightning::loc::Localization::DrainLioOutputs() "DrainLioOutputs()" |
| 连接地图匹配与融合 | @ref lightning::loc::Localization::LidarLocProcCloud() "LidarLocProcCloud()" |
| 业务结果与退出 | @ref lightning::LocSystem::PublishLocalizationResult() "PublishLocalizationResult()"、@ref lightning::loc::Localization::Finish() "Localization::Finish()" |
| 理解积压与丢弃 | @ref lightning::sys::AsyncMessageProcess "AsyncMessageProcess" 的入队、消费与统计接口 |

# 数据路径与有序性

@dot
digraph localization_flow {
  rankdir=LR;
  node [shape=box];
  ros [label="ROS callbacks\nIMU/LiDAR"];
  wheel [label="wheel observation"];
  order [label="sensor_proc_\n按接收顺序串行"];
  lio [label="LaserMapping\nodometry/undistorted cloud"];
  loc [label="LidarLoc\nmap matching/relocalization"];
  pgo [label="localization PGO"];
  result [label="LocalizationResult\ncloud/pose callbacks"];
  ros -> order -> lio -> loc -> pgo -> result;
  wheel -> lio [style=dashed];
}
@enddot

- 在线模式把 IMU 和 LiDAR 包装成统一 `SensorInput`，交给 `sensor_proc_`；处理函数根据变体类型调用 `ProcessIMUData()` 或点云路径。
- 点云路径在 `processing_mutex_` 内把数据交给 LIO，并循环调用 @ref lightning::LaserMapping::RunDetailed() "RunDetailed()"，直到无更多数据。
- 每个 LIO 输出先把 LO 状态送入 `LidarLoc::ProcessLO()` 与定位 PGO；通过关键帧筛选后，定位用投影点云、发布点云和同帧状态组成 `LidarLocInput`，在线模式提交定位处理队列，离线模式直接调用定位处理函数。
- `LidarLoc` 以地图匹配/重定位结果形成定位状态；结果再送入 PGO、点云回调和高频输出路径。

输入回调的 `input_mutex_` 只保护接收/预处理临界区；真正保证 LIO 观察到 IMU 与 LiDAR 回调顺序的是统一 sensor 队列和 `processing_mutex_`。输出另有 `live_output_dispatch_mutex_`，防止两个生产线程逆序发布。

<b>代码依据：</b> `src/core/localization/localization.cpp`（`Init()`、`ProcessSensorInput()`、点云处理循环、`ProcessIMUMsg()`），`src/core/localization/localization.h`（队列和锁）

# 初始化、重定位和地图更新

1. `Localization::Init()` 在 `lifecycle_mutex_` 的独占锁下构造并初始化 LIO、UI、`LidarLoc`、PGO 和回调。
2. `LidarLoc::Init()` 加载配置及 tile 地图；其全局重定位器由配置选择 SOLiD/BTC 实现。
3. 初始化位姿可由外部设置；未初始化时点云路径可进入全局重定位，成功后转入连续地图匹配。
4. `LidarLoc` 自有 `update_map_thread_`，根据当前位姿更新活动 tile；析构负责结束并 join 该线程。

<b>代码依据：</b> @ref lightning::loc::Localization::Init() "Localization::Init()"，`src/core/localization/lidar_loc/lidar_loc.cc`（`Init()`、构造/析构、`UpdateMapThread()`），`src/core/localization/global_relocalizer.h` 及 SOLiD/BTC 实现

@anchor loc_read_prior
# 1. LIO 的位姿怎样成为地图匹配初值？

从 @ref lightning::loc::LidarLoc::Align() "Align()" 的跟踪分支读起。它在扫描末时刻调用 `AssignLOPose()`，随后用相邻两帧的 LO <b>相对运动</b>更新上一帧地图位姿：

@snippet{lineno} core/localization/lidar_loc/lidar_loc.cc loc-lo-prior

只需抓住 `delta = last_lo_pose_.inverse() * current_lo_pose_` 和 `last_abs_pose_ * delta`。`last_lo_pose_`、`current_lo_pose_` 在 odom 系，`last_abs_pose_` 在 map 系；这里把相对运动接到地图位姿后面，没有把 LO 的绝对坐标直接当成地图坐标。

接下来跟踪 `guess_from_lo` 的三个用途：初始化 `current_pose_esti`、选择附近地图、作为 `Localize()` 的输入位姿。`Localize()` 的 `pose` 参数是引用，因此调用后 `current_pose_esti` 已可能变成配准结果；`res_of_lo` 保存这个结果供后续验证。

@snippet{lineno} core/localization/lidar_loc/lidar_loc.cc loc-map-match-dispatch

<b>代入例子：</b> 只考虑同方向的一维平移，上一帧 map 位置为 100 m，LO 从 10 m 走到 11 m，则相对运动为 1 m，地图匹配初值为 101 m。不能直接拿 11 m 去加载地图。

<b>读到这里即可确认初值的来源与去向。</b> 这段只适用于已初始化且未走静止提前返回的跟踪路径；缺少任一 LO 有效标志时，初值保留 `last_abs_pose_`。初值如何变成最终结果，继续看下一节；点云时间匹配细节再查 @ref lightning::loc::LidarLoc::AssignLOPose() "AssignLOPose()"，完整实现见 @ref lidar_loc.cc "lidar_loc.cc"。

@anchor loc_read_acceptance
# 2. NDT 分数够高，为什么仍可能被拒绝？

先在 @ref lightning::loc::LidarLoc::Localize() "Localize()" 中核对基础门槛：

@snippet{lineno} core/localization/lidar_loc/lidar_loc.cc loc-ndt-acceptance

`&&` 要求收敛、分数有限、矩阵有限及分数达标同时成立。`loc_inited_` 决定使用跟踪阈值还是初始化阈值。返回前还会检查转换所得四元数；所以“分数高”本身不保证 `Localize()` 返回 true。

回到调用者 `Align()`，还有两道判断：

@snippet{lineno} core/localization/lidar_loc/lidar_loc.cc loc-tracking-acceptance

- `CheckLidarOdomValid(res_of_lo, ...)` 检查<b>原始配准结果</b>与 LO 的相对运动是否相符。前面虽然把修正按 0.1 权重平滑，检查仍使用保存的 `res_of_lo`，避免平滑掩盖大跳变。
- 地图一致性检查由 `enable_relocalization_map_consistency_` 控制；启用时，`consistency.passed` 可以把原本成功的结果改为失败。这里检查的是 `current_pose_esti`，不要把它与上一项的原始结果混淆。

<b>代入例子：</b> 假设 NDT 已收敛且分数超过阈值，但 LO 指示本帧前进 0.1 m，原始配准指示前进 3 m。在上一帧 LO 有效且平面运动差超过 `lidar_loc_odom_th_` 的条件下，结果仍会被拒绝；无需先读懂 NDT 求解器就能验证这条规则。

<b>此处可停止追踪配准内部。</b> 已确认“基础配准通过 → 里程计一致性 → 可选地图一致性”的接受条件；拒绝后的位姿和状态见下一节。需要核对平面运动差的公式时再进入 @ref lightning::loc::LidarLoc::CheckLidarOdomValid() "CheckLidarOdomValid()"；最终业务发布还要经过 @ref output_contracts "融合与输出门控"。

@anchor loc_read_recovery
# 3. 匹配失败后，本帧用什么位姿，何时重新找回？

在 `Align()` 的失败分支里，沿赋值和计数器阅读：

@snippet{lineno} core/localization/lidar_loc/lidar_loc.cc loc-rejected-match-fallback

`current_pose_esti = guess_from_lo` 将本帧退回 LO 外推。接下来的条件同时要求全局重定位器存在、ready 和失败计数达到阈值，才把 `loc_inited_` 改为 false，并清空旧查询及待确认候选。成功匹配分支会把 `match_fail_count_` 清零，因此该阈值针对连续拒绝，而非运行以来累计失败数。

再看同一函数稍后的状态赋值，就能把行为与日志中的状态名对应起来：

@snippet{lineno} core/localization/lidar_loc/lidar_loc.cc loc-tracking-status

<b>代入例子：</b> 假设失败阈值设为 3，且重定位器 ready。前两次连续拒绝仍保持已初始化，本帧位姿使用 LO 外推，状态为 `FOLLOWING_DR`；第三次拒绝令 `loc_inited_ = false`，本帧状态改为 `INITIALIZING`。后续进入未初始化分支，按外部初值、功能点与全局重定位的条件尝试恢复；这不表示第三次拒绝已完成全局搜索。

<b>读到状态赋值就可以返回主线。</b> `FOLLOWING_DR` 在这里是状态名，本段实际使用的是 `guess_from_lo`，不能据此推断调用了独立 DR 传感器。全局后端未 ready 时，达到计数阈值也不会触发上述复位；这些状态描述的是 LiDAR 匹配结果，不能替代融合 `valid_` 或对外发布许可。完整上下文见 @ref lightning::loc::LidarLoc::Align() "Align()"。

@anchor loc_read_confirmation
# 4. 为什么全局检索找到候选，还不能立即宣布定位恢复？

@ref lightning::loc::LidarLoc::TryGlobalRelocalization() "TryGlobalRelocalization()" 先检查后端与 LO 是否可用，再处理候选、细化和地图验证。这里聚焦最后的连续确认：只有通过前面检查的候选才会到达这段。

@snippet{lineno} core/localization/lidar_loc/lidar_loc.cc loc-confirmation-consistency

`refined_pose` 是候选对应的 map 系位姿，`current_lo_pose_` 是同帧 odom 系位姿。两者组合得到 `T_map_odom`，比较的是两套坐标系之间的关系。确认计数只有在待确认候选有效、时间戳严格变新且平移/旋转差在阈值内时才递增；否则从 1 重新开始。

随后检查计数。下面片段中的日志可略读，只看 `if`、`return false` 和之后的 `loc_inited_ = true`：

@snippet{lineno} core/localization/lidar_loc/lidar_loc.cc loc-confirmation-commit

<b>代入例子：</b> 假设要求确认 2 次，两次均通过细化和地图检查。只考虑同方向平移：第一帧 map=100 m、odom=10 m，第二帧 map=101 m、odom=11 m，两次 `T_map_odom` 都是 90 m。车辆移动了 1 m，候选关系仍一致，第二次可通过确认；若第二次 map=105 m、odom=11 m 且平移差超过门槛，计数退回 1，继续等待。

<b>在这里可以确认接受时机。</b> 当前判断要求正时间间隔，没有使用 `relocalization_confirmation_max_interval_` 作为最大间隔门控；同一时间戳不会递增计数。地图一致性开关关闭时相关检查会放行，确认次数设为 1 时也可能首帧通过，因此“必须两帧”不是固定规则。提交成功后，调用者 @ref lightning::loc::Localization::LidarLocProcCloud() "LidarLocProcCloud()" 根据 `relocalization_accepted` 重置 PGO，再继续融合与输出。检索和几何细化原理见 @ref relocalization_theory "全局重定位算法"，完整源码见 @ref lidar_loc.cc "lidar_loc.cc"。

# 继续跳转

- @ref lightning::loc::Localization::ProcessLidarMsg() "ProcessLidarMsg()"
- @ref lightning::loc::Localization::ProcessIMUMsg() "ProcessIMUMsg()"
- @ref lightning::loc::LidarLoc::ProcessCloud() "LidarLoc::ProcessCloud()"
- @ref lightning::loc::GlobalRelocalizer "GlobalRelocalizer"
- @ref lightning::loc::SolidRelocalizer "SolidRelocalizer"

# 深入阅读

@ref relocalization_theory "全局重定位：检索、几何验证与时间确认"、@ref lidar_residuals "LiDAR 残差、信息矩阵与配准参数化"、@ref pose_graph_theory "位姿图、增量求解与高频平滑"、@ref output_contracts "定位状态、车体参考点与 ROS 输出"。

# 按需参考：对象管理与并发

修改对象初始化、线程或退出逻辑时核对本节；首次阅读优先掌握上面的主路径与算法对应。

## 生命周期和并发检查表

| 对象/资源 | 拥有者 | 并发约束 |
|---|---|---|
| LIO、LidarLoc、PGO、UI | `Localization` 的 `shared_ptr` | 初始化/销毁与输入之间由 `lifecycle_mutex_` 协调 |
| sensor 队列 | `Localization` | 在线模式启动；串行保持跨传感器到达次序 |
| 定位点云队列 | `Localization` | 将 LIO 与地图匹配解耦；同帧状态和点云一起传递 |
| 高频输出队列 | `Localization` | 仅在线模式启动 |
| tile 更新线程 | `LidarLoc` | `match_mutex_` 防止替换/使用 NDT 对象冲突 |
| 定位结果 | `LidarLoc` | `result_mutex_` 保护读写 |

`Localization::~Localization() = default`，而 `AsyncMessageProcess` 没有用析构函数自动 join worker；因此调用者必须在对象销毁前走 `Localization::Finish()`。生产路径由 `LocSystem::Finish()` 保证这一点。`LidarLoc` 析构只为地图更新线程提供最后一道 join 保护，不能替代三个消息队列的显式停止。

<b>代码依据：</b> `src/core/localization/localization.h/.cpp`、`src/core/system/loc_system.cc`、`src/core/system/async_message_process.h`、`src/core/localization/lidar_loc/lidar_loc.cc`

## 关键约束

- 生命周期状态读写必须取得 `lifecycle_mutex_`；输入路径使用共享锁，重新初始化/停止路径使用独占锁。
- 多生产者不能绕过统一 sensor 队列直接并发调用 LIO；源码注释明确要求保持 IMU/LiDAR 回调顺序。
- 地图匹配结果、LIO 里程计和 DR 状态有独立互斥量，读取组合状态时要使用已有 getter，避免自行跨锁拼接。
- 运行统计、时间戳回退检测和静止检测各有专用锁；它们不等同于定位主状态锁。
