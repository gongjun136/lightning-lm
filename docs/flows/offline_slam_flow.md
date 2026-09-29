@page offline_slam_flow 端到端：离线 SLAM 到地图包

# 入口契约

推荐调用 `scripts/run_slam_offline.sh`。它接收 bag、配置和唯一输出目录，source ROS/安装空间，验证 bag 与当前包前缀，然后以独立进程组运行安装后的 `run_slam_offline`。同时启动资源监控，最后验证轨迹、时序、地图 tile 和 global PCD。

**代码依据：** `scripts/run_slam_offline.sh`（usage、prefix 检查、setsid/taskset 调用和末尾 required_outputs 契约）

# 数据路径

@dot
digraph offline_slam {
  rankdir=LR;
  node [shape=box];
  script [label="run_slam_offline.sh"];
  app [label="run_slam_offline main"];
  bag [label="RosbagIO::Go"];
  ingest [label="ProcessIMU / ProcessPointCloud2"];
  sync [label="MultiLidarFrameAssembler + SyncPackages"];
  lio [label="LaserMapping::RunDetailed\nundistort → ESKF → map"];
  kf [label="Keyframe"];
  be [label="BackendPipeline\nBA → BTC → PGO/HBA"];
  export [label="map frame → TiledMap/PCD/PGM\nTUM + diagnostics + DB"];
  script -> app -> bag -> ingest -> sync -> lio -> kf -> be -> export;
}
@enddot

## 1. 初始化

`main()` 校验 flags，在栈上构造 @ref lightning::LaserMapping "LaserMapping" 并用 YAML 初始化。随后从 `system.with_loop_closing` 和 backend mode 选择禁用、新 `ba_btc_hba` 后端或旧 `LoopClosing`，最后按配置可选创建 UI。

**代码依据：** `src/app/run_slam_offline.cc:224-302`（flag 校验、`LaserMapping::Init`、backend 分支和 UI 注入）

## 2. 输入注册与同步

入口从配置取得 topic；多雷达时以 primary sensor 的 IMU topic 为准，并按 bag 中真实 topic type 为每个雷达注册 `PointCloud2` 或 Livox 回调。IMU 回调调用 @ref lightning::LaserMapping "LaserMapping" 的 `ProcessIMU()`，雷达回调调用 `ProcessPointCloud2()`，每次回调后都执行局部 `drain()`。

多雷达输入先进入 @ref lightning::MultiLidarFrameAssembler "MultiLidarFrameAssembler"；组帧结果写入与 LiDAR 队列平行的统计队列。@ref lightning::LaserMapping::SyncPackages() "SyncPackages()" 只有在 LiDAR 帧和覆盖其结束时刻的 IMU 数据满足条件时才构造 `MeasureGroup`。

**代码依据：** `src/app/run_slam_offline.cc:281-291、352-421`、`src/core/lio/laser_mapping.cc:1330-1462`、`src/core/lio/multi_lidar_fusion.h`（topic 选择、回调注册、组帧入队和同步 API）

## 3. 前端状态更新

`drain()` 循环调用 @ref lightning::LaserMapping::RunDetailed() "RunDetailed()"：

1. 无同步数据时返回 `kNoData`，本次回调结束。
2. IMU 处理器完成初始化、预测和点云去畸变。
3. 当前帧降采样后，用 IVox 邻域建立点面/点点约束，由 ESKF 迭代更新状态。
4. 跟踪健康时维护增量地图；达到位移/转角阈值时 @ref lightning::LaserMapping::MakeKF() "MakeKF()" 创建关键帧。
5. 有可发布状态时返回 `kOutput`，入口写 TUM，并记录帧位姿。

首帧、IMU 未初始化、点数过少或跟踪不健康等分支可能“消费但不输出”；因此 `consumed_frames` 与 `output_frames` 不保证相等。

**代码依据：** `src/core/lio/laser_mapping.cc:641-1155、1266-1329`、`src/app/run_slam_offline.cc:352-371`（`RunDetailed` 分支、地图/关键帧更新以及计数逻辑）

## 4. 后端

每轮 drain 读取全部关键帧，只把尚未分发的部分送到后端。新后端离线模式不会创建 keyframe worker，`AddKeyframe()` 直接调用处理；若 HBA 启用，HBA 线程仍在初始化时创建。bag 结束后入口执行 `FlushMultiLidar()`、再次 drain，并用 `WaitUntilIdle(true)` 请求最终全局优化。

**代码依据：** `src/app/run_slam_offline.cc:366-370、421-431`、`src/core/backend/backend_pipeline.cc:209-244、540-570`（关键帧游标、离线同步分支和最终等待）

## 5. 地图与诊断导出

至少一个关键帧是硬前置条件。入口选择 LIO 或优化位姿生成全局点云，估计导出坐标系，原位转换点云，然后写 tiled map、压缩 PCD，以及可选 PGM/YAML。轨迹、后端诊断、BTC 重定位库和 map metadata 随后写出。

导出目录若已存在，C++ 入口会递归删除后重建；标准 shell 运行器在更早阶段拒绝覆盖自己的受管产物。直接运行二进制不会获得脚本的这一层保护。

**代码依据：** `src/app/run_slam_offline.cc:433-573`（关键帧检查、`remove_all`、地图/轨迹/诊断导出）

# 对象生命周期

@dot
digraph lifetime {
  rankdir=TB;
  node [shape=box];
  main [label="main scope"];
  lio [label="LaserMapping (stack)"];
  bag [label="RosbagIO (stack)"];
  backend [label="BackendPipeline (shared_ptr)"];
  ui [label="PangolinWindow (shared_ptr, optional)"];
  main -> lio [label="construct/init"];
  main -> backend [label="construct/init"];
  main -> ui [label="Init success only"];
  main -> bag [label="register callbacks/Go"];
  bag -> lio [label="callbacks borrow by reference"];
  main -> backend [label="WaitUntilIdle"];
  main -> ui [label="Quit"];
}
@enddot

lambda 以引用捕获 `lio`、后端和统计变量，但它们都活到 `main()` 末尾，覆盖 `RosbagIO::Go()` 的同步回放期。后端析构调用 `Shutdown()` 并 join 可连接线程。

**代码依据：** `src/app/run_slam_offline.cc:245-421、575-594`、`src/core/backend/backend_pipeline.cc:207、553-573`（作用域、捕获与析构实现）

# 必须保持的约束

- 时间：输入 pacer 只改变壁钟节奏，不改变传感器时间；轨迹时间戳必须有限、严格递增。
- 同步：LiDAR 点云、起始时间、统计和预处理耗时队列必须同进同出；IMU 必须覆盖扫描结束时间。
- 坐标：地图点云、关键帧位姿、轨迹和诊断必须使用同一 `map_frame::Metadata` 转换。
- 后端：地图导出前必须完成最终优化，否则 `GetOptPose()` 可能仍变化。
- 输出：生产实验优先经 shell 运行器，让 watchdog、资源采样、哈希和完整性检查共同闭环。

**代码依据：** `scripts/run_slam_offline.sh`、`src/core/lio/laser_mapping.cc:1462-1551`、`src/app/run_slam_offline.cc:421-573`（输出校验、同步条件和导出顺序）
