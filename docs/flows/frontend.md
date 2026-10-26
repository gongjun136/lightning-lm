@page frontend_flow 纯前端：传感器到 LIO 状态

# 用途与入口

纯前端用于隔离时间同步、初始化、去畸变和滤波问题。离线入口 @ref run_frontend_offline.cc "run_frontend_offline.cc"、在线入口 @ref run_frontend_online.cc "run_frontend_online.cc" 都直接持有 `LaserMapping`，不创建地图定位器或完整 SLAM 后端。前端仍维护 IVox 局部地图并生成关键帧，“纯前端”不等于不建图；它导出的 LIO 点云也没有回环优化的效果。

@dot
digraph frontend_path { rankdir=LR; node [shape=box]; input [label="bag / ROS callbacks"]; sync [label="preprocess / multi-lidar / IMU sync"]; init [label="IMU init + deskew"]; update [label="ESKF + IVox"]; output [label="state / keyframe / trajectory"]; input -> sync -> init -> update -> output; }
@enddot

本页只核对入口的交接和结束条件：@ref frontend_offline_output_reading "哪些帧进入统计或轨迹" → @ref frontend_offline_finish_reading "尾帧、PCD 与退出成功" → @ref frontend_online_delivery_reading "在线排队、发布与退出"。完整 SLAM 的关键帧分发、最终优化和地图包导出另见 @ref offline_slam_flow "离线 SLAM 流程"；以下例子是源码手工推演。

@anchor frontend_offline_output_reading
# 1. 消费一帧，就会给三份轨迹各写一行吗？

离线 IMU 和 LiDAR 回调输入数据后都执行局部 `drain()`。先看消费计数的边界：

@snippet{lineno} app/run_frontend_offline.cc frontend-offline-consumption

`kNoData` 在计数之前退出；`kConsumed` 与 `kOutput` 都增加 `consumed_frames`。紧接着的多雷达统计分支位于 `kOutput` 判断之前，配置了 CSV 输出时也会记录已消费但未输出的帧：

@snippet{lineno} app/run_frontend_offline.cc frontend-offline-frame-statistics

状态和轨迹的条件更窄：

@snippet{lineno} app/run_frontend_offline.cc frontend-offline-trajectories

`output_frames` 只要求 `kOutput`；写轨迹还要求 `pose_is_ok_`。三次 `WriteTumPose()` 使用同一状态时间，但分别传入 IMU 位姿、主雷达位姿和重设原点后的后轴位姿。它们不是同一参考点下的重复文件，变换实现见 @ref lightning::RearAxlePoseTransformer "RearAxlePoseTransformer"。

每份流还单独检查是否打开及时间顺序：

@snippet{lineno} app/run_frontend_offline.cc frontend-offline-tum-guard

<b>手算例子：</b>假设一次 drain 返回 `kConsumed → kOutput → kNoData`，多雷达模式启用、CSV 与三条轨迹流均已打开，输出状态有效且时间大于 0、比各文件上一行都新。那么消费数增加 2，CSV 增加 2 行，输出数增加 1，三份轨迹各增加 1 行。若本次时间等于上一行，消费和输出数照增，TUM 不追加；某个路径为空则只跳过对应文件。

<b>停止位置与边界：</b>能区分消费、输出状态与真正写行即可。时间 helper 没有显式有限性检查或写入后流状态检查，不能从计数反推出文件完整性。`max_lidar_frames` 比较消费数，达到后设置退出标志。三种 LIO 返回状态见 @ref lio_drain_reading "排空与无输出条件"，完整调用顺序见 @ref run_frontend_offline.cc "离线前端入口"。

@anchor frontend_offline_finish_reading
# 2. bag 结束、PCD 存在或退出码 0，各自证明什么？

回放结束后先尝试排出尾部组帧，再沿同一 `drain()` 处理：

@snippet{lineno} app/run_frontend_offline.cc frontend-offline-final-flush

Flush 仍受组帧最少数量限制，也不能补出缺失 IMU。如果此前帧数上限把退出标志设为 true，上节的 `while (!flg_exit)` 不再运行，末尾调用 drain 也不会无条件处理所有尾帧。组帧和 IMU 边界分别见 @ref sensor_frame_ready_reading "时间桶排出条件"、@ref sensor_imu_sync_reading "IMU 覆盖与消费"。

随后导出的是可选单个前端 PCD：

@snippet{lineno} app/run_frontend_offline.cc frontend-offline-map-output

先读外层两个 `&&` 条件：既要指定 `output_map`，也要已经有关键帧，才尝试保存。`GetGlobalMap(true)` 使用 LIO 位姿；这里没有向 `BackendPipeline` 分发关键帧、等待最终回环优化、写 tiled map 和 manifest 的步骤，因此不能把该 PCD 当作完整定位地图包。

如果确实尝试保存且保存函数返回失败，这里显式返回 3。没有关键帧时整个分支跳过，之后即使消费数和输出数都是 0，仍可到达入口末尾：

@snippet{lineno} app/run_frontend_offline.cc frontend-offline-exit-success

<b>结束例子：</b>假设初始化、bag 读取和输出流打开均未提前失败，回放结束时仍没有关键帧。即便配置了 `output_map`，本次不会调用 PCD 保存函数，仍返回 0；若该路径已有旧文件，它也不能证明本次导出过地图。反过来，保存失败返回 3 时，早先的 TUM/CSV 可能已经创建或写入，整个流程没有文件回滚。

显式参数、LIO 初始化或已捕获的配置/输出打开错误返回 2；这些分支不能概括所有可能的异常。`wait_ui=true` 且 UI 存在时，算法处理完还会等待窗口关闭。`playback_rate` 只调节回放节奏，不能用吞吐速度推断在线性能。

<b>读到这里：</b>确认尾帧前提、PCD 保存条件和末尾无计数门槛的 `return 0` 即可。完整地图包要求继续 @ref offline_slam_export_reading "完整 SLAM 导出顺序"；脚本参数见 @ref script_contracts "Shell 启动与实验脚本契约"。

@anchor frontend_online_delivery_reading
# 3. 在线产生 kOutput，是否已经发布，退出前会排完队列吗？

`FrontendNode` 的 ROS 回调直接向 LIO 输入数据并执行 `Drain()`，不采用在线定位的三个 `AsyncMessageProcess`。但输出另有两个线程：一个发布位姿/点云，另一个周期发布状态。在线入口当前订阅 `PointCloud2` 与 IMU，离线入口则可按 bag 类型注册 Livox `CustomMsg`；不能只换配置就假定两者输入协议完全相同。

先看 `Drain()` 在持有 `state_mutex_` 时如何生成快照：

@snippet{lineno} app/run_frontend_online.cc frontend-online-enqueue

`output` 仅表示 `RunDetailed()==kOutput`，还需有效位姿且发布线程未报告失败才入队。队列上限 50 控制的是<b>等待发布的输出快照</b>。`tracking_normal_` 会被队列过载置 false，但这里没有用它阻止新快照入队。快照使用后轴位姿与 `GetScanUndist()` 的云；不要套用在线定位“发布完整云”的取数路径。

发布线程先把队首移动到局部变量、移出队列，释放锁后才调用两个 publisher：

@snippet{lineno} app/run_frontend_online.cc frontend-online-publish

<b>排队例子：</b>假设等待队列已有 F1…F50，新的 F51 到来，代码删除 F1 后入队 F51，等待队列变为 F2…F51。如果另一个快照已经被线程取走，它不在容量计数中，也不会被这次 `pop_front()` 取消。入队成功仅说明等待发布；两个 `publish()` 依次执行，也不等于订阅端已经收到。

若发布线程在运行期间捕获异常，会设置失败标志并清空等待队列；后续 `Drain()` 的 `!output_publish_failed_` 条件阻止再入队：

@snippet{lineno} app/run_frontend_online.cc frontend-online-publish-failure

这条路径没有重新启动输出线程。析构也明确清队列后再 join：

@snippet{lineno} app/run_frontend_online.cc frontend-online-shutdown

因此退出时不承诺把所有等待快照都发完。发布异常在工作线程内部处理，也不能直接推成 main 必然返回 2；main 的 2 对应其 try/catch 捕获到的异常，正常 spin 返回后是 0。

<b>停止位置：</b>能连起“LIO 输出 → 有界快照队列 → 锁外 publish → 失败/退出清空”即可。完整线程、状态计数与消息时间换算见 @ref run_frontend_online.cc "FrontendNode、OutputPublishLoop 和 StatePublishLoop"。这里的周期心跳与状态也不能独立证明轨迹持续成功。

# 算法调试入口

先核对 @ref sensor_pipeline "组帧与 IMU 同步"，再读 @ref imu_deskew "初始化和去畸变" → @ref eskf_theory "ESKF 预测与更新" → @ref lidar_residuals "残差与有效点"。缓存、关键帧和返回状态集中在 @ref laser_mapping_module "LIO 前端模块"；无需在两个 main 中重新推导算法。
