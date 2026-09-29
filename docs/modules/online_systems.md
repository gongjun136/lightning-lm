@page online_systems 在线系统编排

# 两个系统对象

@ref lightning::SlamSystem "SlamSystem" 封装在线建图，@ref lightning::LocSystem "LocSystem" 封装在线定位。二者都把 ROS 2 生命周期、算法对象和输出副作用集中在一个顶层对象中，但算法数据路径不同：SLAM 产生关键帧并更新地图，定位加载既有地图并发布定位状态。

| 系统 | 主要拥有对象 | 输入 | 主要输出/副作用 |
|---|---|---|---|
| `SlamSystem` | LIO、legacy/new backend、G2P5、UI、TF broadcaster | IMU、单/多 LiDAR | 优化位姿、TF、点云/地图 topic、保存地图服务 |
| `LocSystem` | `Localization`、TF broadcaster、telemetry、heartbeat | IMU、单/多 LiDAR、轮速 | 定位结果、TF、诊断/功能安全消息、可选 TUM 轨迹 |

**代码依据：** `src/core/system/slam.h/.cc`、`src/core/system/loc_system.h/.cc`（成员、订阅/发布创建和回调）

# 在线 SLAM 生命周期

1. `Init(yaml)` 读取系统选项，初始化 LIO、所选后端、可选 UI/G2P5，并安装优化回调。
2. `StartSLAM(map_name)` 切换到建图状态；此前到达的数据不会按正常建图流程落图。
3. ROS 回调将 IMU/LiDAR 转成内部消息；LIO 输出新关键帧后交给所选后端，同时驱动地图和发布路径。
4. 当前 SLAM 回调直接驱动 LIO；析构等待新后端 idle 并 `Shutdown()`，最后退出 UI，没有异步前端队列。详见 @ref online_slam_flow "在线建图：ROS 输入到地图服务"。

**代码依据：** `src/core/system/slam.cc`（构造/析构、`Init()`、`StartSLAM()`、`ProcessIMU()`、LiDAR 处理与关键帧派发）

# 在线定位生命周期

在线定位的完整逐帧路径、状态机和发布门控见 @ref online_localization_flow "在线定位端到端流程"；本页只保留系统编排层的共同职责。

1. `LocSystem::Init()` 解析地图路径及 ROS 输出配置，构造 `Localization` 并注册结果/点云回调。
2. ROS 回调只做消息适配、计数和转发；算法队列与定位线程由 `Localization` 管理。
3. 定位回调发布 pose、path、TF、状态与诊断；可选轨迹文件由独立互斥量串行写入。
4. `Finish()` 依次让 sensor、定位和高频输出队列 drain 并 join，再停止地图更新线程、heartbeat 和轨迹文件；析构会再次进入同一幂等路径。

**代码依据：** `src/core/system/loc_system.h`（`Finish()` 和成员），`src/core/system/loc_system.cc`（`Init()`、输入转发、发布与析构），`src/core/localization/localization.cpp`（`Finish()`）

# `AsyncMessageProcess<T>` 的队列契约

这个模板是在线定位数据面的公共线程原语：`Start()` 创建一个 worker；`AddMessage()` 入队并通知；worker 在条件变量上等待，逐个调用已注册处理函数。超过 `max_size_` 时淘汰最旧待处理消息；`Quit()` 会继续处理到队列为空再 join。队列、运行标志和统计量均封装在对象内部。

**代码依据：** `src/core/system/async_message_process.h`（`Start()`、`AddMessage()`、`ProcLoop()`、`Quit()`）

## 线程关系

@dot
digraph online_threads {
  rankdir=LR;
  node [shape=box];
  ros [label="ROS executor callbacks"];
  sensor [label="sensor AsyncMessageProcess"];
  lio [label="LIO / system processing"];
  backend [label="BackendPipeline worker"];
  hba [label="optional HBA thread"];
  loc [label="localization cloud worker"];
  map [label="LidarLoc map-update thread"];
  ros -> sensor -> lio [label="localization"];
  ros -> lio [label="SLAM direct callback"];
  lio -> backend -> hba;
  lio -> loc -> map;
}
@enddot

图表达拥有关系和数据转交，不表示线程会同时处理同一对象。各分支是否存在由入口、在线/离线模式及 YAML 开关共同决定。

# 关键约束

- ROS 回调和算法 worker 的边界必须保留；直接在新增订阅回调里调用耗时算法会改变时序和背压。
- 后端优化回调可能来自 worker/HBA 线程，发布或地图刷新代码必须按跨线程回调处理。
- 退出时先停止生产者、再 drain 消费者、最后释放被回调对象；颠倒顺序会留下悬空回调风险。
- 在线入口的 executor 线程数和回调组设置会影响并发；仅凭系统类不能断言所有 ROS 回调串行。

**推断：** 具体部署中的回调并发度由入口选择的 executor 和 ROS 2 配置决定；本文只陈述源码对象内部的同步，不假定运行时一定单线程。

**代码依据：** `src/app/run_slam_online.cc`、`src/app/run_loc_online.cc`、`src/core/system/*.cc`、`src/core/system/async_message_process.h`
