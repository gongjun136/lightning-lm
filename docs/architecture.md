@page architecture 模块地图与入口

# 逻辑模块导航

| 模块 | 主要目录/类型 | 输入 → 输出 | 生命周期与线程 |
|---|---|---|---|
| 程序入口 | `src/app` | flags/ROS/bag → 系统对象与产物 | `main()` 持有顶层对象；决定退出顺序 |
| 公共模型与配置 | `src/common`, @ref lightning::NavState "NavState", `Keyframe` | YAML/传感器状态 → 共享数据结构 | 值对象和共享关键帧；跨前后端传递 |
| LIO 前端 | `src/core/lio`, @ref lightning::LaserMapping "LaserMapping" | IMU + 单/多雷达 → 状态、关键帧、局部图 | 调用者驱动 `RunDetailed()`；内部缓冲受锁保护 |
| 新后端 | `src/core/backend`, @ref lightning::backend::BackendPipeline "BackendPipeline" | 关键帧 → 局部 BA、BTC 约束、PGO/HBA、优化位姿 | 在线 worker + 可选 HBA 线程；析构 `Shutdown()` join |
| 旧后端 | `src/core/loop_closing`, @ref lightning::LoopClosing "LoopClosing" | 关键帧 → 回环与优化位姿 | 兼容 `backend.mode=legacy` |
| 在线 SLAM 编排 | `src/core/system`, @ref lightning::SlamSystem "SlamSystem" | ROS 回调 → LIO/后端/UI/服务 | 系统对象拥有模块；ROS 回调驱动前端；后端停止时等待 idle |
| 定位 | `src/core/localization`, @ref lightning::loc::Localization "Localization" | 地图 + IMU/LiDAR → 定位结果 | `LocSystem` 编排在线资源；全局重定位有 SOLiD/BTC 实现 |
| 地图 | `src/core/maps`, @ref lightning::TiledMap "TiledMap" | 点云/关键帧 → tile、global PCD、PGM/YAML、metadata | 显式导入导出；坐标归一化必须先于最终写出 |
| 2.5D 地图 | `src/core/g2p5`, @ref lightning::g2p5::G2P5 "G2P5" | 位姿/点云 → 2.5D/2D 表示 | 可由系统选项启用 |
| 图优化器 | `src/core/miao` | 顶点/边 → 优化状态 | 作为 `miao.core`/`miao.utils` 链接到主库 |
| ROS/bag 适配 | `src/wrapper`, @ref lightning::RosbagIO "RosbagIO" | ROS 2 bag/message → 内部回调与点云类型 | `Go()` 同步遍历 bag；回调由调用者注册 |
| UI | `src/ui`, @ref lightning::ui::PangolinWindow "PangolinWindow" | 位姿/点云 → 3D 窗口 | 可选；显式 `Init()`/`Quit()` |
| 运行时工具 | `src/utils` | 计时、预算、因果追踪、线程调度 | 多数为进程级辅助状态 |
| 文件/YAML I/O | `src/io` | 文件/YAML ↔ 内部配置/对象 | 同步 I/O |

**代码依据：** `src/CMakeLists.txt`、各目录头文件（编译源清单、公开 API 与成员所有权）

跨模块阅读不要沿物理目录来回跳：在线定位从 @ref online_localization_flow "在线定位端到端流程" 进入，离线建图从 @ref offline_slam_flow "离线 SLAM 端到端流程" 进入；遇到需要理解的类，再回到对应模块页和 API 源码。

# 入口到编排对象

@dot
digraph entry_ownership {
  rankdir=LR;
  node [shape=box];
  offline_slam [label="run_slam_offline\n直接编排"];
  online_slam [label="run_slam_online"];
  online_loc [label="run_loc_online"];
  lio [label="LaserMapping"];
  backend [label="BackendPipeline / LoopClosing"];
  slam [label="SlamSystem"];
  loc [label="LocSystem"];
  offline_slam -> lio;
  offline_slam -> backend;
  online_slam -> slam -> lio;
  slam -> backend;
  online_loc -> loc;
}
@enddot

上图的关键区别是：不能把所有 SLAM 入口都简化成 `main → SlamSystem`。离线入口为控制实验过程、统计和导出，直接编排前后端；在线入口才由系统类封装 ROS 订阅、服务和异步队列。

**代码依据：** `src/app/run_slam_offline.cc`、`src/app/run_slam_online.cc`、`src/core/system/slam.h`（对象构造与成员字段）

# 对象所有权与关闭顺序

- 离线 SLAM：`main` 栈上持有 @ref lightning::LaserMapping "LaserMapping"；后端和 UI 是 `shared_ptr`。bag 播放完成后先 flush 前端，再等待后端空闲，随后导出并退出；后端析构会再次安全调用 `Shutdown()`。
- 新后端：`unique_ptr` 独占 BA/BTC/HBA 算法对象；关键帧用 `shared_ptr` 跨队列与调用者共享。`data_mutex_` 保护优化数据，`queue_mutex_` 保护在线队列，`hba_mutex_` 保护 HBA 请求状态。
- 前端：预处理器、IMU 处理器与 IVox 由 `shared_ptr` 持有；LiDAR/时间/IMU 队列由 `mtx_buffer_` 共同保护。当前帧与地图点云为 PCL shared pointer。
- UI：只有 `Init()` 成功才注入前端；出口显式 `Quit()`，可选 `wait_ui` 仅控制是否等待用户关闭。

**代码依据：** `src/app/run_slam_offline.cc:245-302、421-579`、`src/core/backend/backend_pipeline.h:96-128`、`src/core/lio/laser_mapping.h:335-478`（构造、共享/独占指针、锁和退出流程）

# 并发模型

| 执行上下文 | 触发条件 | 工作 | 同步点 |
|---|---|---|---|
| bag 主线程 | 离线入口 | 顺序触发 IMU/LiDAR 回调和 `drain()` | `rosbag.Go()` 返回代表输入结束 |
| 后端 worker | `BackendPipeline::Init(..., true)` | 从 `queue_` 串行取关键帧 | `queue_cv_`; `WaitUntilIdle()` 等空队列且不处理中 |
| HBA 线程 | `options_.hba.enabled` | 接收合并后的全局优化请求 | `hba_cv_`; idle 条件为无请求且未运行 |
| 在线定位 sensor 队列 | Localization 在线模式 | 串行消费 IMU/LiDAR | 有界队列；停止前 drain/quit；SLAM 回调直接驱动前端 |
| OpenMP 工作线程 | LIO 匹配/预处理 | 点级并行 | 线程数来自配置/环境；共享归约块固定化 |

**推断：** bag 回调本身是否跨线程最终取决于 `RosbagIO::Go()` 的实现；当前 `Go()` 实现按 reader 循环直接调用注册回调，因此本文按单调用线程描述离线入口。

**代码依据：** `src/wrapper/bag_io.cc`、`src/core/backend/backend_pipeline.cc:209-229、444-570`、`src/core/system/async_message_process.h`（回调循环、线程创建、条件变量和退出条件）
