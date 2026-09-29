@page localization_module 定位与重定位

# 模块边界

@ref lightning::loc::Localization "Localization" 是定位数据面：它初始化 LIO、@ref lightning::loc::LidarLoc "LidarLoc" 和定位 PGO，接收 IMU/LiDAR/轮速并输出 @ref lightning::loc::LocalizationResult "LocalizationResult"。@ref lightning::LocSystem "LocSystem" 是 ROS 2 编排层，负责订阅、发布、TF、轨迹记录和功能安全心跳；两者不要混为一个对象。

如果目标是从入口系统阅读在线定位，请先看 @ref online_localization_flow "在线定位端到端流程"。本页用于按类边界回查实现，不重复脚本、ROS topic 和整条调用链。

**代码依据：** `src/core/localization/localization.h`、`src/core/system/loc_system.h`（公开接口、成员所有权）

# 数据路径与有序性

@dot
digraph localization_flow {
  rankdir=LR;
  node [shape=box];
  ros [label="ROS callbacks\nIMU/LiDAR/wheel"];
  order [label="sensor_proc_\n按接收顺序串行"];
  lio [label="LaserMapping\nodometry/undistorted cloud"];
  loc [label="LidarLoc\nmap matching/relocalization"];
  pgo [label="localization PGO"];
  result [label="LocalizationResult\ncloud/pose callbacks"];
  ros -> order -> lio -> loc -> pgo -> result;
}
@enddot

- 在线模式把 IMU 和 LiDAR 包装成统一 `SensorInput`，交给 `sensor_proc_`；处理函数根据变体类型调用 `ProcessIMUData()` 或点云路径。
- 点云路径在 `processing_mutex_` 内把数据交给 LIO，并循环调用 @ref lightning::LaserMapping::RunDetailed() "RunDetailed()"，直到无更多数据。
- 每个 LIO 输出同时送入 `LidarLoc::ProcessLO()` 与定位 PGO；完整去畸变点云和同帧状态进入定位处理队列。
- `LidarLoc` 以地图匹配/重定位结果形成定位状态；结果再送入 PGO、点云回调和高频输出路径。

输入回调的 `input_mutex_` 只保护接收/预处理临界区；真正保证 LIO 观察到 IMU 与 LiDAR 回调顺序的是统一 sensor 队列和 `processing_mutex_`。输出另有 `live_output_dispatch_mutex_`，防止两个生产线程逆序发布。

**代码依据：** `src/core/localization/localization.cpp`（`Init()`、`ProcessSensorInput()`、点云处理循环、`ProcessIMUMsg()`），`src/core/localization/localization.h`（队列和锁）

# 初始化、重定位和地图更新

1. `Localization::Init()` 在 `lifecycle_mutex_` 的独占锁下构造并初始化 LIO、UI、`LidarLoc`、PGO 和回调。
2. `LidarLoc::Init()` 加载配置及 tile 地图；其全局重定位器由配置选择 SOLiD/BTC 实现。
3. 初始化位姿可由外部设置；未初始化时点云路径可进入全局重定位，成功后转入连续地图匹配。
4. `LidarLoc` 自有 `update_map_thread_`，根据当前位姿更新活动 tile；析构负责结束并 join 该线程。

**代码依据：** `src/core/localization/localization.cpp:28` 起的 `Init()`，`src/core/localization/lidar_loc/lidar_loc.cc`（`Init()`、构造/析构、`UpdateMapThread()`），`src/core/localization/global_relocalizer.h` 及 SOLiD/BTC 实现

# 生命周期和并发检查表

| 对象/资源 | 拥有者 | 并发约束 |
|---|---|---|
| LIO、LidarLoc、PGO、UI | `Localization` 的 `shared_ptr` | 初始化/销毁与输入之间由 `lifecycle_mutex_` 协调 |
| sensor 队列 | `Localization` | 在线模式启动；串行保持跨传感器到达次序 |
| 定位点云队列 | `Localization` | 将 LIO 与地图匹配解耦；同帧状态和点云一起传递 |
| 高频输出队列 | `Localization` | 仅在线模式启动 |
| tile 更新线程 | `LidarLoc` | `match_mutex_` 防止替换/使用 NDT 对象冲突 |
| 定位结果 | `LidarLoc` | `result_mutex_` 保护读写 |

`Localization::~Localization() = default`，而 `AsyncMessageProcess` 没有用析构函数自动 join worker；因此调用者必须在对象销毁前走 `Localization::Finish()`。生产路径由 `LocSystem::Finish()` 保证这一点。`LidarLoc` 析构只为地图更新线程提供最后一道 join 保护，不能替代三个消息队列的显式停止。

**代码依据：** `src/core/localization/localization.h/.cpp`、`src/core/system/loc_system.cc`、`src/core/system/async_message_process.h`、`src/core/localization/lidar_loc/lidar_loc.cc`

# 关键约束

- 生命周期状态读写必须取得 `lifecycle_mutex_`；输入路径使用共享锁，重新初始化/停止路径使用独占锁。
- 多生产者不能绕过统一 sensor 队列直接并发调用 LIO；源码注释明确要求保持 IMU/LiDAR 回调顺序。
- 地图匹配结果、LIO 里程计和 DR 状态有独立互斥量，读取组合状态时要使用已有 getter，避免自行跨锁拼接。
- 运行统计、时间戳回退检测和静止检测各有专用锁；它们不等同于定位主状态锁。

# 继续跳转

- @ref lightning::loc::Localization::ProcessLidarMsg() "ProcessLidarMsg()"
- @ref lightning::loc::Localization::ProcessIMUMsg() "ProcessIMUMsg()"
- @ref lightning::loc::LidarLoc::ProcessCloud() "LidarLoc::ProcessCloud()"
- @ref lightning::loc::GlobalRelocalizer "GlobalRelocalizer"
- @ref lightning::loc::SolidRelocalizer "SolidRelocalizer"

# 深入阅读

@ref relocalization_theory "全局重定位：检索、几何验证与时间确认"、@ref lidar_residuals "LiDAR 残差、信息矩阵与配准参数化"、@ref pose_graph_theory "位姿图、增量求解与高频平滑"、@ref output_contracts "定位状态、车体参考点与 ROS 输出"。
