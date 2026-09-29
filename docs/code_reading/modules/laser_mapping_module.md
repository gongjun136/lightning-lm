@page laser_mapping_module 核心模块样板：LaserMapping

# 模块职责

@ref lightning::LaserMapping "LaserMapping" 是 LIO 前端边界：接收 ROS/PCL 点云与 IMU，完成多雷达组帧、时间同步、去畸变、ESKF 预测/更新、局部 IVox 地图维护和关键帧生成。它不负责 ROS bag 遍历、全局回环优化或最终地图包写出。

**代码依据：** `src/core/lio/laser_mapping.h:84-330`、`src/app/run_slam_offline.cc`（公开 API 与调用者承担的后端/导出职责）

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

**代码依据：** `src/core/lio/laser_mapping.cc:566-1155、1266-1655`（输入、`RunDetailed`、关键帧、同步、地图与观测模型实现）

# 生命周期与所有权

- 构造时仅保存 `Options` 并建立成员默认值；`Init(yaml)` 才创建预处理器、IMU 处理器、IVox 并装配 ESKF 观测函数。
- `LaserMapping` 独占自己的值状态（ESKF、队列、计数器），但通过 shared pointer 持有大对象与 UI；关键帧 shared pointer 会被入口和后端继续持有。
- 析构只释放当前帧点云指针并记录日志；类本身没有后台线程需要 join。并发来自调用者、OpenMP 以及可选 UI，而不是 `LaserMapping` 自己创建的 `std::thread`。

**代码依据：** `src/core/lio/laser_mapping.cc:29-69`、`src/core/lio/laser_mapping.h:119-126、335-478`（初始化、成员类型和析构体）

# 线程与锁

`mtx_buffer_` 保护 LiDAR、时间、IMU、帧统计和预处理耗时队列。`wheel_speed_mutex_` 单独保护轮速缓冲与统计。若在线调用者让 IMU 与 LiDAR 回调并发，这两个锁是输入一致性的边界；`RunDetailed()` 的算法状态不是为多个并发消费者设计的。

点匹配可使用 OpenMP。观测 Hessian/gradient 采用固定块存储，使求和结果不依赖工作线程数；这既是性能设计也是数值可复现约束。

**推断：** 公开 API 没有声明 `RunDetailed()` 可重入，且它修改大量无锁单帧成员，因此应由单一消费线程串行调用。在线系统当前也以异步消息处理队列满足这一条件。

**代码依据：** `src/core/lio/laser_mapping.h:374-478`、`src/core/system/slam.cc`、`src/core/system/async_message_process.h`（锁、共享临时状态和单线程处理队列）

# 关键不变量

1. `time_buffer_[i]` 与 `lidar_buffer_[i]` 表示同一帧；多雷达统计/预处理耗时也必须按相同顺序出队。
2. 输入时间回退会触发缓冲清理或拒绝；轨迹消费者不能假定每个输入都产生输出。
3. ESKF 更新前必须完成 IMU 初始化和点云去畸变；点数低于阈值时不得构造不可靠观测。
4. `MakeKF()` 生成的 ID 单调递增；后端以关键帧 shared pointer 原地更新优化位姿。
5. `GetGlobalMap(use_lio_pose, ...)` 的位姿选择必须与 backend mode 一致：无后端用 LIO pose，有后端用 optimized pose。
6. 多雷达外参与 primary frame 的定义同时影响融合点云、地图射线原点和导出结果，不能只改一处。

**代码依据：** `src/core/lio/laser_mapping.cc`、`src/app/run_slam_offline.cc:448-501`（队列同步、关键帧、全局图姿态选择和传感器原点）

# 失败与观测

| 现象 | 代码层含义 | 先查 |
|---|---|---|
| consumed 增长、output 不增长 | 同步帧被消费但初始化/健康门限未通过 | IMU 初始化、点数、tracking health |
| partial frame 增多 | 多雷达组帧超时或缺传感器 | `MultiLidarFrameStats`、late/tolerance/invalid drops |
| 后端没有新 keyframe | 位移/角度阈值未达到或跟踪不健康 | `Options::kf_*`、`MakeKF()` 分支 |
| 轨迹跳变 | 时间同步、外参、观测退化或轮速融合异常 | 最大 IMU gap、残差/有效特征、wheel stats |
| 地图与轨迹不一致 | 导出时 pose 类型或 map-frame metadata 不一致 | `use_lio_pose`、metadata 传递 |

# 复制此样板扩展模块

后续模块页保持相同七段：职责边界；输入/状态/输出；主执行路径；生命周期/所有权；线程/锁；关键不变量；失败与观测。每个结论至少附一个文件和一个符号/行段；只有调用关系不能证明语义时标“推断”。
