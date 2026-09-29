@page online_slam_flow 在线建图：ROS 输入到地图服务

# 入口与主路径

@ref run_slam_online.cc "run_slam_online.cc" 初始化 ROS，构造栈上 `SlamSystem`，`Init(config)` 后执行 `StartSLAM("new_map")`，再 `Spin`。可选 `--bag` 创建内部回放线程，播放结束等待 `--post_wait_seconds` 后 shutdown；这种模式的回放等待不是所有算法完成的证明。

@dot
digraph online_slam_path { rankdir=LR; node [shape=box]; ros [label="ROS executor"]; lio [label="ProcessIMU / ProcessLidar\nDrainLio"]; kf [label="new keyframe"]; backend [label="Backend worker / legacy"]; out [label="map->odom / G2P5 / UI"]; save [label="SaveMap + final optimization"]; ros -> lio -> kf -> backend -> out; backend -> save; }
@enddot

# 实际线程边界

`SlamSystem::ProcessIMU/ProcessLidar` 直接驱动前端，当前类没有在线定位的 sensor `AsyncMessageProcess`。新后端 `BackendPipeline::Init(..., true)` 创建 worker，HBA 可再有独立线程。优化回调发布 map→odom，并可触发 G2P5 重绘；这一回调可能运行在后端线程。

`DrainLio` 只保存有效且严格递增的 LIO 状态；比较关键帧指针，只有新关键帧才分发给后端、G2P5 和 UI。系统 `with_loop_closing` 与 `backend.mode` 共同决定是否创建后端，单独改 mode 并不保证启动。

# 保存与退出

服务 `lightning/save_map` 使用请求 map_id；`lightning/optimize_backend` 请求全局优化。`SaveMap` 先等待后端（可能运行 HBA），检查关键帧、估计导出坐标，再写地图、轨迹、诊断与数据库。**现有目标目录会被清空重建**，它不是原子替换；应选择新的输出目录。

`--save_map` 在 Spin 退出并 join 回放线程后保存，失败返回 4。析构等待后端 idle、Shutdown，最后退出 UI；并不存在“先 join 异步 LIO 队列”这一步。

**源码依据：** @ref lightning::SlamSystem "SlamSystem" 的 `Init`、`DrainLio`、`SaveMap`、析构；@ref run_slam_online.cc "run_slam_online.cc"。
模块与原理：@ref backend_module "后端、回环与优化"、@ref backend_optimization "体素 BA、回环位姿图与分层优化"；地图契约：@ref map_contracts "地图包、坐标归一化与栅格导出"。
