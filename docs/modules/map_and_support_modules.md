@page map_and_support_modules 地图与支撑模块

# 地图存储与导出

@ref lightning::TiledMap "TiledMap" 将静态/动态点云按空间块组织，索引与块文件支持按位姿加载；@ref lightning::MapChunk "MapChunk" 在需要时把块点云载入内存。NDT voxel 统计由 `UpdateVoxel()` 增量更新。`map_frame` 负责地图坐标归一化和 metadata，`navigation_map_export` 把点云及可选射线帧导出为 PGM/YAML。

@dot
digraph map_products {
  rankdir=LR;
  node [shape=box];
  kf [label="optimized keyframes"];
  cloud [label="global cloud"];
  tile [label="TiledMap\nindex + chunks"];
  frame [label="map_frame\nnormalized metadata"];
  nav [label="navigation export\nPGM + YAML"];
  kf -> cloud;
  kf -> tile;
  cloud -> frame -> nav;
  tile -> frame;
}
@enddot

关键约束：tile 索引和实际块文件必须成套；坐标归一化后的 offset/原点信息必须与导出的点云、栅格和 metadata 一致；`LoadOnPose()` 会改变内存中的活动块集合，不能视为只读查询。

**代码依据：** `src/core/maps/tiled_map.h/.cc`、`tiled_map_chunk.h/.cc`、`map_frame.h/.cc`、`navigation_map_export.h/.cc`，`src/app/run_slam_offline.cc`（最终导出顺序）

# 2.5D/2D 地图

@ref lightning::g2p5::G2P5 "G2P5" 接收前端关键帧，并维护可供显示/发布的栅格表示。前端绘制使用 @ref lightning::sys::AsyncMessageProcess "AsyncMessageProcess"；后端位姿更新会唤醒独立重绘线程。析构和重新初始化必须停止两条路径后再释放地图对象。

**代码依据：** `src/core/g2p5/g2p5.h/.cc`（关键帧回调、`draw_frontend_map_thread_`、`draw_backend_map_thread_` 和 `RenderBack()`）

# IVox、多雷达与预处理

- `PointCloudPreprocess` 把不同厂商 ROS 点类型转换为内部 `PointXYZIT`，并做盲区/线数/时间等处理。
- `MultiLidarFusion` 依据 lidar id、外参和时间窗口聚合输入；flush 是离线输入结束时不可省略的边界动作。
- @ref lightning::IVox "IVox" 是 LIO 的增量近邻地图；节点实现负责体素内采样/近邻候选，LIO 的匹配与地图更新依赖其参数一致性。

**代码依据：** `src/core/lio/pointcloud_preprocess.h/.cc`、`multi_lidar_fusion.h/.cc`、`src/core/ivox3d/ivox3d.h`、`src/core/lio/laser_mapping.cc`

# ROS bag、文件与 YAML 适配

@ref lightning::RosbagIO "RosbagIO" 通过 `AddHandle(topic, callback)` 建立 topic 到回调的映射，在 `Go()` 或 `GoRealtime()` 的调用线程中同步派发。未注册 topic 不进入回调；两个读取函数都忽略回调的 bool 返回值。`YAML_IO` 和 `file_io` 提供同步文件操作，不拥有后台线程。

@anchor support_bag_dispatch_reading
## 回调返回 false，或者重复注册话题，会怎样？

先看注册使用的是 `emplace`：

@snippet{lineno} wrapper/bag_io.h bag-topic-registration

同一个键已存在时不会覆盖旧回调。再看 `Go()` 的整个循环；`iter->second(msg)` 没有接收返回值，也没有用它控制 break：

@snippet{lineno} wrapper/bag_io.cc bag-dispatch-stop

<b>手工推演：</b>先给 `/imu` 注册 A，再给同一话题注册 B，随后读取 `/imu`、一个未注册话题、`/imu` 三条消息。假设没有异常或退出请求，即使 A 每次返回 false，仍调用 A 两次、B 零次。若 `sleep_usec=1000`，三条消息之后各等待一次，固定等待共 3 ms，此外还有读取和回调耗时；未注册话题也占这一等待。

<b>边界与停止点：</b>循环依靠 bag 结束或 debug 命名空间中的 `flg_exit` 返回；异常向外传播。退出标志在派发与固定等待之后检查，不能抢占正在运行的回调，进入函数前置位也不等同于在读取第一条消息前检查。构造器的 `dataset_type` 当前不改变这套分发。确认注册语句、未使用的返回值和退出检查位置即可；完整 API 见 @ref lightning::RosbagIO::AddHandle "AddHandle()"、@ref lightning::RosbagIO::Go "Go()"，离线算法如何设置停止条件见 @ref offline_slam_drain_reading "离线帧消费与上限"。

@anchor support_bag_clock_reading
## 两倍速回放改变哪一种时间，回调慢了会怎样？

`GoRealtime()` 先要求倍率有限且严格为正：

@snippet{lineno} wrapper/bag_io.cc bag-realtime-rate

节拍使用 bag 存储时间 `msg->time_stamp`，起点对应 `steady_clock`，不是消息内的 header 时间：

@snippet{lineno} wrapper/bag_io.cc bag-realtime-clock

第一条读出的消息就建立时钟锚点，即使它没有注册回调；每条消息先等目标时刻，再查话题和执行回调。函数本身不重写序列化消息的 header。`Go()` 没有这段存储时间节拍；两种接口不能仅按名称推断出相同的回放时序。

<b>手算：</b>两条消息的 bag 时间为 100.0 s、100.2 s，`playback_rate=2` 时第二条目标为实际播放起点之后 0.1 s。它们原有的 header 即使是 87.0 s、87.2 s，也不会因此变成 100.0 s、100.1 s。若第一条回调已耗时 0.15 s，第二条目标时刻已经过去，`sleep_until` 不会再补等 0.1 s；处理仍是同步串行，不会因设置二倍速而获得额外处理能力。

<b>边界与停止点：</b>这里不保证实时期限，也没有对每条 bag 时间额外做单调性检查；目标落在过去时不会等待。退出标志在等待与派发后才检查，不能中断 `sleep_until`。在线内置回放由 `OnlineBagPlayer::Run()` 调用此接口，直接离线 SLAM/定位/前端入口调用的是 `Go()`。读到目标时间公式与派发顺序即可；完整源码见 @ref lightning::RosbagIO::GoRealtime "GoRealtime()"、@ref online_bag_player.h "OnlineBagPlayer"，在线系统的等待范围见 @ref online_slam_flow "在线建图入口"。

# miao 优化器

`src/core/miao` 被编译为 `miao.core` 与 `miao.utils` 两个共享库。调用方通过配置选择算法和线性求解器，向图加入顶点/边，再执行初始化与迭代；SLAM 后端和定位 PGO 都复用这套图优化抽象。顶点、边和鲁棒核主要用 `shared_ptr` 交给图持有。

**代码依据：** `src/core/miao/CMakeLists.txt`、`src/core/miao/core/graph/optimizer.h`、`graph.h`、`src/core/backend/backend_pipeline.cc`（优化器装配和使用）

# UI、公共模型与运行时工具

| 模块 | 作用 | 生命周期/副作用 |
|---|---|---|
| `src/ui` | Pangolin 窗口、车体/点云/轨迹对象 | `Init()` 创建可视资源，`Quit()` 结束；无 UI 模式不得假设窗口存在 |
| `src/common` | `NavState`、`Keyframe`、参数与点类型 | `Keyframe::Ptr` 跨前后端共享；优化位姿可在后端线程更新 |
| `src/utils` | timing、profiling、诊断与线程辅助 | 部分工具直接参与 Safety 与线程调度，不能一概当作只读统计；按实际调用者核对影响 |
| `src/io` | YAML/文件序列化 | 同步磁盘副作用；错误需沿返回值或异常边界处理 |

**代码依据：** `src/ui`、`src/common/keyframe.h`、`src/common/nav_state.h`、`src/utils`、`src/io` 的公开头文件及 `src/CMakeLists.txt` 源清单

`src/utils/radar_dtc_state.h` 的 `RadarDtcState::Evaluate()` 直接参与 Safety 判定：通信恢复和同步历史 DTC
保留诊断可见性，但新鲜的历史 FAIL 不单独触发当前故障；当前 DTC、证据过期和四台设备可用性仍参与判断。
完整的故障优先级与主从版本要求见 @ref dtc_delivery 。

`src/utils/thread_scheduling.cc` 的在线角色策略在创建工作线程前初始化，
验证各掩码属于启动允许集合，应用后读回实际亲和性；不匹配时拒绝继续。
角色赋值日志可与 `/proc/<pid>/task/*/status` 对照，但不能单独证明调度时限。
没有启用角色环境变量时不会改变调用线程的亲和性。

## 计算预算与诊断工具

| API 入口 | 职责与阅读边界 |
|---|---|
| @ref lightning::compute::ComputeBudget "ComputeBudget"、@ref lightning::compute::LoadComputeBudget "LoadComputeBudget()" | 以调用方默认值为起点，读取 YAML 与环境覆盖，校验后整体提交；加载函数不创建线程或实际应用亲和性 |
| @ref lightning::profiling::Stopwatch "Stopwatch"、@ref lightning::profiling::MultiStageTimingWindow "MultiStageTimingWindow" | 区分单调经过时间、调用线程 CPU 与含并发线程的进程 CPU；窗口由 Add() 驱动，输入停止后不自动输出未到期窗口 |
| @ref lightning::profiling::CausalEvent "CausalEvent"、@ref lightning::profiling::CausalTraceEnabled "CausalTraceEnabled()"、@ref lightning::profiling::RecordCausalEvent "RecordCausalEvent()" | 首次调用按输出路径初始化共享追踪器；事件经有界队列异步写入，允许因容量或锁竞争丢弃诊断记录 |

接口页说明单位、无效值、配置覆盖和生命周期约束；预算使用场景见 @ref sensor_pipeline "传感器与算力预算"，
计时日志解释见 @ref resource_diagnostics "资源与端到端延迟"，事件文件与完整性检查见 @ref causal_diagnostics "因果链与时间关联"。

# 阅读时的跨模块约束

- 坐标系名称必须和变换方向一起读；不能仅从变量含 `map`/`odom` 猜乘法方向。
- `Keyframe` 同时携带 LIO pose 与 optimized pose；地图导出应使用哪一个由流程阶段决定。
- PCL shared pointer 表示共享所有权，不等于对象内容线程安全；跨线程写点云仍需调用侧同步。
- bag 回放是同步派发，但在线 ROS 输入不是；复用代码时必须保留两种时序差异。
- 文件写出、地图块换入换出、UI 和 ROS 发布都是可观察副作用，文档变更检查不能只看算法类。

**推断：** `src/utils` 中某些诊断工具按设计应只观测而不影响结果；是否完全无时序影响需要性能测试，源码静态阅读只能确认其调用边界。


# 深入阅读

@ref map_contracts "地图包、坐标归一化与栅格导出"、@ref sensor_pipeline "传感器预处理、多雷达组帧与算力预算"、@ref output_contracts "定位状态、车体参考点与 ROS 输出"、@ref pose_graph_theory "位姿图、增量求解与高频平滑"。
