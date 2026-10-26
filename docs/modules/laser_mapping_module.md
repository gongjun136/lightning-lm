@page laser_mapping_module LIO 前端：LaserMapping

# 模块职责

@ref lightning::LaserMapping "LaserMapping" 是 LIO 前端边界：接收 ROS/PCL 点云与 IMU，完成多雷达组帧、时间同步、去畸变、ESKF 预测/更新、局部 IVox 地图维护和关键帧生成。它不负责 ROS bag 遍历、全局回环优化或最终地图包写出。

<b>代码依据：</b> @ref lightning::LaserMapping "LaserMapping API"、@ref run_slam_offline.cc "离线建图入口"（调用者承担后端与导出职责）。

# 输入、状态与输出

| 类别 | 内容 | 代码锚点 |
|---|---|---|
| 输入 | `PointCloud2`、Livox `CustomMsg` 或内部 `CloudPtr`；`IMUPtr`；可选轮速 | @ref lightning::LaserMapping "LaserMapping API" 的 `ProcessPointCloud2()`/`ProcessIMU()` |
| 配置 | `Options` + YAML 的传感器、体素、噪声、多雷达、滤波与运行预算 | @ref lightning::LaserMapping::Init() "Init()" |
| 长期状态 | IVox、ESKF、关键帧列表、IMU/LiDAR 缓冲、多雷达 assembler | @ref lightning::LaserMapping "LaserMapping 成员" |
| 单帧状态 | `MeasureGroup`、去畸变/降采样点云、对应关系、`NavState` | @ref lightning::LaserMapping::RunDetailed() "RunDetailed()" |
| 输出 | 当前状态、健康度、关键帧、全局点云、统计 | getter 与 @ref lightning::LaserMapping "LaserMapping API" |

# 主执行路径

@dot
digraph laser_mapping {
  rankdir=TB;
  node [shape=box];
  cloud [label="ProcessPointCloud2"];
  assemble [label="preprocess / multi-lidar assemble"];
  queue [label="lidar/time/stats queues"];
  imu [label="ProcessIMU → imu queue"];
  sync [label="SyncPackages"];
  predict [label="ImuProcess: init/predict/undistort"];
  down [label="voxel downsample"];
  obs [label="ObsModel: IVox neighbors/residuals"];
  eskf [label="ESKF iterated update"];
  health [label="tracking health + load controller"];
  map [label="MapIncremental"];
  kf [label="MakeKF"];
  cloud -> assemble -> queue -> sync;
  imu -> sync;
  sync -> predict -> down -> obs -> eskf -> health;
  health -> map;
  health -> kf;
}
@enddot

<b>代码依据：</b> @ref lightning::LaserMapping::RunDetailed() "RunDetailed()"、@ref lightning::LaserMapping::SyncPackages() "SyncPackages()"、@ref lightning::LaserMapping::MakeKF() "MakeKF()"。下节先核对在线调用者怎样驱动这些步骤，再进入算法专题。

# 沿在线定位核对四个问题

本节对应 `LocSystem + Localization` 的 LiDAR 在线分支。`Localization` 负责排队和派发，`LaserMapping` 负责同步与估计；不要把调用者的队列当作前端内部的 IMU/LiDAR 缓冲。每段先提出一个问题，阅读片段时只追踪指出的语句；其他分支可按需从 API 页进入完整源码。

@anchor lio_sensor_queue_reading
## 1. IMU 与点云为什么由同一个 worker 处理？

先对比两种消息的入队位置，再看 worker 的分派条件。下面的点云片段来自 `PointCloud2` 重载；执行到这里前已经完成输入检查、预处理和可选的 LIO 跳帧判断。

@htmlonly[block]
<details open><summary>展开源码：两种载荷进入 sensor_proc_，再按 is_imu 分派</summary>
@endhtmlonly

@snippet{lineno} core/localization/localization.cpp localization-lidar-enqueue

@snippet{lineno} core/localization/localization.cpp localization-imu-enqueue

@snippet{lineno} core/localization/localization.cpp localization-sensor-dispatch

@htmlonly[block]
</details>
@endhtmlonly

两次 `AddMessage()` 的接收者都是 `sensor_proc_`。聚合初始化中的 `false/true` 对应 @ref lightning::loc::Localization::SensorInput "SensorInput" 的 `is_imu`：点云载荷放在 `cloud`，IMU 载荷放在 `imu`。`ProcessSensorInput()` 只据此选择处理函数；它没有再创建两个工作线程。队列处理函数在 @ref lightning::loc::Localization::Init() "Init()" 中绑定，单 worker 的实现见 @ref lightning::sys::AsyncMessageProcess "AsyncMessageProcess"。

<b>用例核对：</b> 若 IMU A、点云 B、IMU C 依次进入同一个队列，且没有容量淘汰，worker 就按 A → B → C 调用处理函数。这里保留的是<b>入队顺序</b>，没有按消息时间戳重新排序；重复/回退 IMU 仍会在 @ref lightning::loc::Localization::ProcessIMUData() "ProcessIMUData()" 中被拒绝。在线模式之外，两种输入走直接调用分支。

找到“相同队列对象”和“一个类型分派”后，就能确认这层编排。只有需要核对某一帧为何被拒绝时，才继续读输入校验或 @ref sensor_pipeline "传感器预处理"。

@anchor lio_drain_reading
## 2. 点云到了，为什么可能暂时没有 LIO 输出？

在 @ref lightning::LaserMapping::SyncPackages() "SyncPackages()" 中，下面这个条件检查 IMU 是否已经覆盖当前扫描的结束时刻。两个时间戳的单位都是秒；此前还会检查缓冲是否为空并确定本帧起止时间。

@htmlonly[block]
<details open><summary>展开源码：等待 IMU 覆盖扫描末时刻，再区分无数据、已消费和有效输出</summary>
@endhtmlonly

@snippet{lineno} core/lio/laser_mapping.cc lio-wait-for-imu

@snippet{lineno} core/localization/localization.cpp localization-drain-status

@htmlonly[block]
</details>
@endhtmlonly

`SyncPackages()` 返回 false 时，@ref lightning::LaserMapping::RunDetailed() "RunDetailed()" 返回 `kNoData`。调用者 @ref lightning::loc::Localization::DrainLioOutputs() "DrainLioOutputs()" 的 `break` 结束本次排空；`kConsumed` 则经过 `continue` 继续尝试下一帧，只有 `kOutput` 才读取状态并调用 `ProcessLO()`、`ProcessLidarOdom()`。因此 `kConsumed` 不能当成一帧有效定位结果，也不能把它当成“以后没有数据了”。

<b>用例核对：</b> 扫描结束时刻为 12.10 s，最新 IMU 只有 12.09 s，其他条件正常时仍需等待。接受新的 12.11 s IMU 后，`ProcessIMUData()` 会再次驱动排空：

@htmlonly[block]
<details><summary>展开源码：IMU 更新后也会调用 DrainLioOutputs()</summary>
@endhtmlonly

@snippet{lineno} core/localization/localization.cpp localization-imu-drain

@htmlonly[block]
</details>
@endhtmlonly

此时前述时间覆盖条件已满足，但 IMU 初始化、点数与其他前端检查仍可能使该帧只得到 `kConsumed`。一次输入也可能让多帧缓存依次可处理，所以不能用 ROS 回调次数推断 LIO 输出次数。

能区分 `break`、`continue` 和状态交接后，本问题已核对完。进一步追踪具体失败条件时，进入 `RunDetailed()` 中对应的返回分支；去畸变计算见 @ref imu_deskew "IMU 与逐点去畸变"。

@anchor lio_keyframe_reading
## 3. 只在关键帧做地图定位，会不会停止其他帧的 LIO 更新？

先看 @ref lio_drain_reading "上节的 LIO 排空循环"：`ProcessLO(lo_state)` 与 `ProcessLidarOdom(lo_state)` 已经执行。随后 `DrainLioOutputs()` 才到达下面的关键帧选择；`kf == lio_kf_` 比较的是当前与上次提交的关键帧对象。

@htmlonly[block]
<details open><summary>展开源码：关键帧未变化时跳过地图匹配，新关键帧在线入队或离线直调</summary>
@endhtmlonly

@snippet{lineno} core/localization/localization.cpp localization-keyframe-selection

@snippet{lineno} core/localization/localization.cpp localization-map-enqueue

@htmlonly[block]
</details>
@endhtmlonly

第一个片段的 `continue` 跳过本次循环剩余的地图匹配提交，下一次循环仍会尝试 `RunDetailed()`。关键帧发生变化时，代码先保存 `lio_kf_ = kf`，再进入第二个片段；在线模式调用 `AddMessage()`，离线模式直接调用 `LidarLocProcCloud()`。`loc_on_kf_` 关闭时则不检查关键帧是否变化，每个 `kOutput` 都走地图匹配提交路径。

<b>用例核对：</b> 两帧 S1、S2 都产生 `kOutput`，但 `GetKeyframe()` 返回同一个 K；假设 S1 已提交 K，开启 `lidar_loc.loc_on_kf` 后，S2 仍会将 LO 状态交给 LidarLoc 和 PGO，但不会再提交地图匹配。这里确认的是调用边界；进入匹配队列后仍受跳帧、容量淘汰策略影响，不能把“提交”理解成“完成匹配”。

读到这里就能解释关键帧开关影响哪一步。需要理解何时产生新 K，再读 @ref lightning::LaserMapping::MakeKF() "MakeKF()"；需要跟踪提交后的处理，转到 @ref localization_module "地图定位模块"。

@anchor lio_localization_queue_reading
## 4. 定位队列容量为 1，究竟会丢哪一帧？

先核对容量来自哪里，再看入队与执行。下面限制的是 `lidar_loc_proc_cloud_` 的等待队列；它与上面的传感器队列是两个对象。

@htmlonly[block]
<details open><summary>展开源码：可覆盖的默认容量、队尾入队/队头淘汰、取出后执行</summary>
@endhtmlonly

@snippet{lineno} core/localization/localization.cpp localization-map-queue-limit

@snippet{lineno} core/system/async_message_process.h async-message-capacity

@snippet{lineno} core/system/async_message_process.h async-message-execute-one

@htmlonly[block]
</details>
@endhtmlonly

`default_loc_queue_size = 1` 可由 `system.online_lidar_loc_queue_size` 覆盖，零值会使初始化失败。`push_back()` 把新消息加到队尾，超限后 `pop_front()` 淘汰队头，因此删掉的是最旧的<b>待处理</b>帧。worker 先把一帧移入局部变量 `msg`，从队列删除并解锁，再执行 `custom_func_(msg)`；入队淘汰不会取消这个正在执行的调用。

假设容量为 1，A 正在计算，B 和 C 均通过入队前的跳帧判断，且都在 A 完成前到达：

| 时刻 | 正在执行 | 等待队列 | 可从代码确认的结果 |
|---|---|---|---|
| A 已取出 | A | 空 | A 已经不占等待队列容量 |
| B 入队 | A | B | 队列未超限 |
| C 入队 | A | B、C → C | 队头 B 被淘汰，丢弃计数加一 |
| A 完成 | C | 空 | worker 下一次取出 C |

所以“容量为 1”允许同时有一帧执行、一帧等待。@ref lightning::sys::AsyncMessageProcess::PendingCount() "PendingCount()" 包含这两部分，出现 2 不一定表示容量失效；@ref lightning::sys::AsyncMessageProcess::DroppedCount() "DroppedCount()" 不统计周期跳帧或退出后拒收的消息。

确认默认值、队头/队尾和执行帧所在位置后即可停止。若要修改容量或结束顺序，再进入 @ref lightning::sys::AsyncMessageProcess "队列完整 API 与源码"；在线各队列如何接力见 @ref online_localization_flow "在线定位流程"。

# 关键不变量

1. `time_buffer_[i]` 与 `lidar_buffer_[i]` 表示同一帧；多雷达统计/预处理耗时也必须按相同顺序出队。
2. 输入时间回退会触发缓冲清理或拒绝；轨迹消费者不能假定每个输入都产生输出。
3. ESKF 更新前必须完成 IMU 初始化和点云去畸变；点数低于阈值时不得构造不可靠观测。
4. `MakeKF()` 生成的 ID 单调递增；后端以关键帧 shared pointer 原地更新优化位姿。
5. `GetGlobalMap(use_lio_pose, ...)` 的位姿选择必须与 backend mode 一致：无后端用 LIO pose，有后端用 optimized pose。
6. 多雷达外参与 primary frame 的定义同时影响融合点云、地图射线原点和导出结果，不能只改一处。

<b>代码依据：</b> @ref lightning::LaserMapping::SyncPackages() "SyncPackages()"、@ref lightning::LaserMapping::MakeKF() "MakeKF()"、@ref lightning::LaserMapping::GetGlobalMap() "GetGlobalMap()"、@ref run_slam_offline.cc "离线建图入口"（全局图姿态选择和传感器原点）。

# 失败与观测

| 现象 | 代码层含义 | 先查 |
|---|---|---|
| consumed 增长、output 不增长 | 同步帧被消费但初始化/健康门限未通过 | IMU 初始化、点数、tracking health |
| partial frame 增多 | 多雷达组帧超时或缺传感器 | `MultiLidarFrameStats`、late/tolerance/invalid drops |
| 后端没有新 keyframe | 位移/角度阈值未达到或跟踪不健康 | `Options::kf_*`、`MakeKF()` 分支 |
| 轨迹跳变 | 时间同步、外参、观测退化或轮速融合异常 | 最大 IMU gap、残差/有效特征、wheel stats |
| 地图与轨迹不一致 | 导出时 pose 类型或 map-frame metadata 不一致 | `use_lio_pose`、metadata 传递 |

# 深入算法

@ref sensor_pipeline "传感器预处理、多雷达组帧与算力预算" → @ref imu_deskew "IMU 初始化、时间积分与点云去畸变" → @ref eskf_theory "18 维迭代 ESKF：从预测到观测注入" → @ref lidar_residuals "LiDAR 残差、信息矩阵与配准参数化"。

# 按需参考：对象管理与并发

修改对象初始化、线程或退出逻辑时核对本节；首次阅读优先掌握上面的主路径与算法对应。

## 生命周期与所有权

- 构造时仅保存 `Options` 并建立成员默认值；`Init(yaml)` 才创建预处理器、IMU 处理器、IVox 并装配 ESKF 观测函数。
- `LaserMapping` 独占自己的值状态（ESKF、队列、计数器），但通过 shared pointer 持有大对象与 UI；关键帧 shared pointer 会被入口和后端继续持有。
- 析构只释放当前帧点云指针并记录日志；类本身没有后台线程需要 join。并发来自调用者、OpenMP 以及可选 UI，而不是 `LaserMapping` 自己创建的 `std::thread`。

<b>代码依据：</b> @ref lightning::LaserMapping "LaserMapping 的构造、析构与成员"、@ref lightning::LaserMapping::Init() "Init()"。

## 线程与锁

`mtx_buffer_` 保护 LiDAR、时间、IMU、帧统计和预处理耗时队列。`wheel_speed_mutex_` 单独保护轮速缓冲与统计。若在线调用者让 IMU 与 LiDAR 回调并发，这两个锁是输入一致性的边界；`RunDetailed()` 的算法状态不是为多个并发消费者设计的。

点匹配可使用 OpenMP。观测 Hessian/gradient 采用固定块存储，使求和结果不依赖工作线程数；这既是性能设计也是数值可复现约束。

<b>推断：</b> 公开 API 没有声明 `RunDetailed()` 可重入，且它修改大量无锁单帧成员，因此应由单一消费线程串行调用。在线定位以统一 sensor 队列串行消费；在线建图由当前 ROS executor 回调直接驱动，见 @ref online_slam_flow "在线建图：ROS 输入到地图服务"。

<b>代码依据：</b> @ref lightning::LaserMapping "LaserMapping 的锁与单帧成员"、@ref lightning::SlamSystem "SlamSystem"、@ref lightning::sys::AsyncMessageProcess "AsyncMessageProcess"（在线调用者与单线程处理队列）。
