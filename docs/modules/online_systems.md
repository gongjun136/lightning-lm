@page online_systems 在线系统编排

# 两个系统对象

另有启动时选择的 `InsLocSystem`：单线程 ROS 回调持有 CGI 最新字段、质量门控、有限时长位姿缓冲和多雷达组装；独立心跳线程随对象析构停止。Init 配置失败返回 false，Spin 随 ROS shutdown 退出。它不创建下面的 LIO/地图定位对象，详见 @ref ins_only_operation "纯组合导航生命周期"。

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
4. 当前 SLAM 回调直接驱动 LIO，没有异步前端队列；析构等待新后端 idle 并 `Shutdown()`，或销毁 legacy 回环对象并排空其队列，再退出 UI。回环尾帧仍可调用 G2P5 重绘，因此地图对象保持有效直到回环退出。详见 @ref online_slam_flow "在线建图：ROS 输入到地图服务"。

**代码依据：** `src/core/system/slam.cc`（构造/析构、`Init()`、`StartSLAM()`、`ProcessIMU()`、LiDAR 处理与关键帧派发）

# 在线定位生命周期

在线定位的完整逐帧路径、状态机和发布门控见 @ref online_localization_flow "在线定位端到端流程"；本页只保留系统编排层的共同职责。

1. `LocSystem::Init()` 解析地图路径，校验该包 `map_manifest.yaml` 的几何哈希与变换状态，从清单加载 ROS 地图输出变换，再初始化 `Localization` 并注册结果/点云回调。缺失或未知清单在创建 ROS 节点前返回失败。
2. ROS 回调只做消息适配、计数和转发；算法队列与定位线程由 `Localization` 管理。
3. 定位回调发布 pose、path、TF、状态与诊断；可选轨迹文件由独立互斥量串行写入。
4. `Finish()` 依次让 sensor、定位和高频输出队列 drain 并 join，再停止地图更新线程、heartbeat 和轨迹文件；析构会再次进入同一幂等路径。

**代码依据：** `src/core/system/loc_system.h`（`Finish()` 和成员），`src/core/system/loc_system.cc`（`Init()`、输入转发、发布与析构），`src/core/localization/localization.cpp`（`Finish()`）

# AsyncMessageProcess 的队列契约

定位、legacy 回环和 G2P5 共用 `core/system/async_message_process.h` 中的 @ref lightning::sys::AsyncMessageProcess "AsyncMessageProcess"，只维护这一份实现。`Start()` 创建一个 worker；`AddMessage()` 入队并通知；worker 在条件变量上等待，逐个调用已注册处理函数。超过 `max_size_` 时淘汰最旧待处理消息；`Quit()` 会继续处理到队列为空再 join。队列、运行标志和统计量均封装在对象内部。

@ref lightning::LoopClosing::Init() "LoopClosing::Init()" 显式设置等待容量为 40，另有最多一帧正在处理；离线模式仍由 `AddKF()` 同步调用处理函数。在线回环析构会等待未被容量淘汰的关键帧处理完毕，可能比旧队列退出耗时更长。所有者须在回调依赖对象销毁前显式调用 `Quit()`；公共队列不在析构中隐式停止线程。

`SetThreadInit()` 在启动前设置 worker 初始化回调，worker 进入消息循环前执行一次。
在线定位用它设置传感器、地图匹配及输出线程的 CPU 角色；绑定失败不静默继续。
`run_loc_online` 在日志和 ROS 线程创建前建立默认后台掩码，`LocSystem::Spin()` 再绑定执行线程。
OpenMP 工作线程（含未分到迭代的线程）进入并行区时绑定后台核，主线程保留传感器核。
环境变量、适用硬件、重启与现场核验见 @ref dtc_delivery 。

**代码依据：** `src/core/system/async_message_process.h`（`Start()`、`AddMessage()`、`ProcLoop()`、`Quit()`）

@anchor online_queue_quit_reading
## 退出请求会立即丢掉所有待处理消息吗？

本段只讨论在线定位使用的 `AsyncMessageProcess`。容量淘汰与“正在处理的一帧”先见 @ref lio_sensor_queue_reading "传感器队列导读"；退出时则从 `Quit()` 的这几行读起：

@snippet{lineno} core/system/async_message_process.h async-message-quit-join

先在锁内设置退出标志，再唤醒 worker，最后 `join()` 等它退出。这里没有 `clear()`，也没有取消正在运行的回调。是否排空由 worker 的退出条件决定：

@snippet{lineno} core/system/async_message_process.h async-message-exit-condition

两个条件用的是 `&&`：已经请求退出，而且等待队列为空，才 `break`。若队列还有消息，仍沿 @ref lio_sensor_queue_reading "逐条执行路径" 取出并执行回调。与此同时，新输入被拒绝：

@snippet{lineno} core/system/async_message_process.h async-message-refuse-after-exit

<b>时序例子：</b>假设 A 正在执行、B 和 C 仍在等待，此时所有者调用 Quit。A 不会被中断；A 完成后依次执行 B、C，队列为空才退出并让 join 返回。退出标志已置位后提交的 D 会直接返回，不进入等待队列。此前已因容量限制被淘汰的消息不会恢复。

<b>边界与停止点：</b>这是正常回调能够返回时的时序；没有超时或强制终止，某条回调一直阻塞就会拖住退出。应由 worker 之外的所有者调用 Quit，并先停止生产输入。`ProcessedCount()` 只说明回调完成，不证明算法接受或 ROS 已发布。读到拒收、退出条件和 join 即可，完整源码见 @ref lightning::sys::AsyncMessageProcess "队列模板"。

@anchor online_finish_order_reading
## 为什么先退出传感器队列，再退出地图匹配和高频输出？

`Localization::Finish()` 按数据传递的先后顺序关闭三个队列：

@snippet{lineno} core/localization/localization.cpp localization-finish-order

前一个 Quit 返回时，它的 worker 已完成当前和等待中的回调。传感器尾帧因而仍有机会进入尚未关闭的地图匹配与高频输出队列；随后才分别关闭它们。最后再拿生命周期独占锁，结束地图定位对象和 UI。队列容量、输入时效和各层接受条件仍可能拒绝尾帧，这个顺序不是“保证每帧产生输出”。

外层资源的收尾也有顺序：

@snippet{lineno} core/system/loc_system.cc loc-system-finish-resources

<b>手工推演：</b>假设传感器尾帧生成一个有效定位结果，且后续容量与门控均允许处理，那么算法队列排空时，心跳与轨迹文件尚未由本函数关闭；`loc_->Finish()` 返回后才停止心跳、关闭轨迹流。若先关闭消费者，再让生产者排空，其新消息会命中上一节的退出后拒收分支。

<b>边界与停止点：</b>外层先取消 DTC/Safety 定时器并清除 debug sink，因此尾部处理不代表这些定时诊断仍继续产出。`finished_` 支持完成后的顺序重复调用；这段普通标志不构成多个线程并发调用 Finish 的安全保证，也没有负责停止所有外部生产者。读到三队列 join 与外层资源关闭即可，完整接口见 @ref lightning::loc::Localization::Finish() "Localization::Finish()" 和 @ref lightning::LocSystem::Finish() "LocSystem::Finish()"。

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
