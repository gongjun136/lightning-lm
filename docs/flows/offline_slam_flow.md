@page offline_slam_flow 端到端：离线 SLAM 到地图包

# 入口契约

推荐调用 `scripts/run_slam_offline.sh`。它接收 bag、配置和唯一输出目录，source ROS/安装空间，验证 bag 与当前包前缀，然后以独立进程组运行安装后的 `run_slam_offline`。同时启动资源监控，最后验证轨迹、时序、地图 tile 和 global PCD。

<b>代码依据：</b> @ref run_slam_offline.sh "运行器完整源码"（usage、prefix 检查、setsid/taskset 调用和末尾 required_outputs 契约）。脚本调用层见 @ref mapping_script_guide "建图脚本导读"。

@anchor offline_slam_code_reading
## 先验证三个容易混淆的问题

| 问题 | 短验证路径 |
|---|---|
| 消费帧、输出帧、关键帧是不是同一个计数？ | @ref offline_slam_drain_reading "排出结果与关键帧游标" |
| bag 结束后，为什么还要 Flush 和等待后端？ | @ref offline_slam_finish_reading "尾帧与最终优化" |
| 看到 global.pcd 就说明完整导出成功了吗？ | @ref offline_slam_export_reading "目录替换、轨迹重写与清单收尾" |

以下例子是对分支的手工推演，不表示已回放 bag 或实际覆盖地图。

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

@ref run_slam_offline.cc "离线 SLAM 入口" 校验 flags，在栈上构造 @ref lightning::LaserMapping "LaserMapping" 并用 YAML 初始化。随后从 `system.with_loop_closing` 和 backend mode 选择禁用、新 `ba_btc_hba` 后端或旧 `LoopClosing`，最后按配置可选创建 UI。

<b>代码依据：</b> @ref run_slam_offline.cc "入口中的 flags 校验、LaserMapping::Init、backend 分支和 UI 注入"。

## 2. 输入注册与同步

入口从配置取得 topic；多雷达时以 primary sensor 的 IMU topic 为准，并按 bag 中真实 topic type 为每个雷达注册 `PointCloud2` 或 Livox 回调。IMU 回调调用 @ref lightning::LaserMapping "LaserMapping" 的 `ProcessIMU()`，雷达回调调用 `ProcessPointCloud2()`，每次回调后都执行局部 `drain()`。

多雷达输入先进入 @ref lightning::MultiLidarFrameAssembler "MultiLidarFrameAssembler"；组帧结果写入与 LiDAR 队列平行的统计队列。@ref lightning::LaserMapping::SyncPackages() "SyncPackages()" 只有在 LiDAR 帧和覆盖其结束时刻的 IMU 数据满足条件时才构造 `MeasureGroup`。

<b>验证入口：</b> @ref run_slam_offline.cc "topic 选择与回调注册"；组帧等待条件见 @ref sensor_frame_ready_reading "何时排出融合帧"，IMU 队列消费见 @ref sensor_imu_sync_reading "同步边界"。

@anchor offline_slam_drain_reading
## 3. 前端状态更新

`drain()` 循环调用 @ref lightning::LaserMapping::RunDetailed() "RunDetailed()"：

1. 无同步数据时返回 `kNoData`，本次回调结束。
2. IMU 处理器完成初始化、预测和点云去畸变。
3. 当前帧降采样后，用 IVox 邻域建立点面/点点约束，由 ESKF 迭代更新状态。
4. 跟踪健康时维护增量地图；达到位移/转角阈值时 @ref lightning::LaserMapping::MakeKF() "MakeKF()" 创建关键帧。
5. 有可发布状态时返回 `kOutput`，入口写 TUM，并记录帧位姿。

首帧、IMU 未初始化、点数过少或跟踪不健康等分支可能“消费但不输出”；因此 `consumed_frames` 与 `output_frames` 不保证相等。

先看入口如何消费三种返回值，而不是先展开完整 LIO 算法：

@snippet{lineno} app/run_slam_offline.cc offline-slam-drain-status

`kNoData` 在计数前 `break`；其他返回值先增加 `consumed_frames`，只有 `kOutput` 增加 `output_frames`。其内部还要检查 `pose_is_ok_` 才保存帧位姿。TUM 写入函数另有文件是否打开和时间顺序检查，因此这三个数也不能直接等同于 TUM 行数。

接着看关键帧分发。`dispatched_keyframes` 在 lambda 外初始化为 0，跨多次 `drain()` 保留：

@snippet{lineno} app/run_slam_offline.cc offline-slam-keyframe-cursor

<b>手算例子：</b>若一次排出依次返回 `kConsumed → kOutput → kNoData`，则消费数增加 2、输出数增加 1。若上次已分发第 0 个关键帧，本次列表变成 `[KF0, KF1, KF2]`，游标从 1 开始，只分发 KF1、KF2。分发循环位于 `kOutput` 条件外，不能把“本次不是输出帧”推成“完全不查看关键帧列表”。

<b>读到这里：</b>能区分返回状态、计数和关键帧游标即可。`max_lidar_frames` 比较的是消费数，包含 `kConsumed`，达到上限后设退出标志。前端计算继续 @ref imu_deskew "去畸变"、@ref eskf_theory "ESKF"；新后端如何处理关键帧见 @ref backend_queue_ba_reading "入队与局部 BA"。完整代码在 @ref run_slam_offline.cc "drain lambda 所在入口"。

@anchor offline_slam_finish_reading
## 4. 后端

每轮 drain 读取全部关键帧，只把尚未分发的部分送到后端。新后端离线模式不会创建 keyframe worker，`AddKeyframe()` 直接调用处理；若 HBA 启用，HBA 线程仍在初始化时创建。bag 结束后入口执行 `FlushMultiLidar()`、再次 drain，并用 `WaitUntilIdle(true)` 请求最终全局优化。

@snippet{lineno} app/run_slam_offline.cc offline-slam-final-drain

按调用顺序读：`rosbag.Go()` 结束输入回放；`FlushMultiLidar()` 尝试处理 assembler 中尚未排出的尾部桶；`drain()` 仍沿前述 `RunDetailed()` 路径；最后新后端的 `WaitUntilIdle(true)` 等待并请求最终优化。Flush 不能凭空补出缺失 IMU，等待后端也不能把前端未输出的帧变成有效状态。

<b>时序例子：</b>假设正常回放结束、退出标志仍为 false、尾桶满足强制排出条件且 IMU 已覆盖扫描末，Flush 后这帧可以进入 drain；若 IMU 覆盖不足，它仍可能返回 `kNoData`。若此前 `max_lidar_frames` 已把 debug 命名空间中的 `flg_exit` 设为 true，末尾仍调用 Flush，但 `drain()` 中检查该标志的 while 循环不再执行，不能宣称“结束时无条件处理所有尾帧”。

<b>停止位置：</b>能连起回放、尾桶、前端排出和后端等待四步即可。`WaitUntilIdle(true)` 的 true 会请求全局优化，并非纯粹查询空闲；只在使用新后端时执行本片段的等待。实现见 @ref lightning::backend::BackendPipeline::WaitUntilIdle() "WaitUntilIdle()"，组帧前提见 @ref sensor_frame_ready_reading "融合帧准备条件"。

@anchor offline_slam_export_reading
## 5. 地图与诊断导出

至少一个关键帧是硬前置条件。入口选择 LIO 或优化位姿生成全局点云，估计导出坐标系，原位转换点云，然后写 tiled map、压缩 PCD，以及可选 PGM/YAML。轨迹、后端诊断、BTC 重定位库和 map metadata 随后写出。

最后写 `map_manifest.yaml`，绑定几何 SHA-256 和业务输出坐标声明；保存失败返回 4。新导出默认 `unknown`，不会复制共用定位 YAML 的历史新→老变换；只有明确配置 `map_export.alignment_status: not_required` 才声明原生坐标可直接使用。详见 @ref map_contracts "地图清单与启动契约"。

导出目录若已存在，C++ 入口会递归删除后重建；标准 shell 运行器在更早阶段拒绝覆盖自己的受管产物。直接运行二进制不会获得脚本的这一层保护。

先定位破坏旧目录的确切位置。它位于 `!FLAGS_backend_evaluation_only` 的地图导出分支，地面估计成功以后：

@snippet{lineno} app/run_slam_offline.cc offline-slam-export-replace

估计失败时，这个 lambda 先返回，尚未删除地图目录；成功后则先变换内存点云，再删除并重建目录。此后任何 tile、PCD、PGM 或附加产物失败，都没有恢复旧目录的代码。`backend_evaluation_only` 会跳过此地图分支，但仍可能写轨迹和后端诊断，不是只读模式。

再核对轨迹。播放期间写过 `output_tum`，最后会关闭流，并用统一 metadata 重写帧轨迹及其他轨迹文件：

@snippet{lineno} app/run_slam_offline.cc offline-slam-rewrite-trajectories

`output_tum` 对应保存的逐帧 LIO 位姿，另有关键帧 LIO、关键帧优化与修正后的逐帧轨迹；不能把四个文件都解释成同一条优化轨迹。最后的地图清单是收尾的一步：

@snippet{lineno} app/run_slam_offline.cc offline-slam-finalize-map

<b>状态例子：</b>假设 `global.pcd` 已写出，随后必要的 BTC 数据库或 manifest 保存失败，程序仍返回 4，目录中可能残留部分新文件。看到 PCD 存在只证明它曾写出，不能证明整个地图包已完成。反过来，地面估计失败保护的是此处的地图目录；早先打开的轨迹输出流可能已创建或截断文件，不能扩大成“失败前没有任何文件副作用”。

<b>读到这里：</b>能指出目录替换、轨迹重写与最终清单各自的位置即可。坐标与哈希契约继续 @ref map_frame_normalization_reading "同一导出变换"、@ref map_manifest_validation_reading "地图启动校验"；完整实现见 @ref run_slam_offline.cc "地图导出入口及 SaveFrameTrajectoryTum 等辅助函数"。

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

lambda 以引用捕获 `lio`、后端和统计变量，但它们都活到入口作用域末尾，覆盖 `RosbagIO::Go()` 的同步回放期。后端析构调用 `Shutdown()` 并 join 可连接线程。

<b>代码依据：</b> @ref run_slam_offline.cc "入口的局部对象与回调捕获"、@ref backend_pipeline.cc "后端析构和 Shutdown"。

# 必须保持的约束

- 时间：输入 pacer 只改变壁钟节奏，不改变传感器时间。轨迹应使用有限、严格递增时间；此入口的 `WriteTumPose()` 只显式检查正数和递增，不能单凭它证明非有限输入已排除。
- 同步：LiDAR 点云、起始时间、统计和预处理耗时队列必须同进同出；IMU 必须覆盖扫描结束时间。
- 坐标：地图点云、关键帧位姿、轨迹和诊断必须使用同一 `map_frame::Metadata` 转换。
- 后端：地图导出前必须完成最终优化，否则 `GetOptPose()` 可能仍变化。
- 输出：生产实验优先经 shell 运行器，让 watchdog、资源采样、哈希和完整性检查共同闭环。

<b>代码依据：</b> @ref run_slam_offline.sh "脚本输出校验"、@ref lightning::LaserMapping::SyncPackages() "同步条件"，以及本页的 @ref offline_slam_finish_reading "尾帧处理" 和 @ref offline_slam_export_reading "导出顺序"。
