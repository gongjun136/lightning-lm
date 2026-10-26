@page online_slam_flow 在线建图：ROS 输入到地图服务

# 入口与主路径

@ref run_slam_online.cc "run_slam_online.cc" 初始化 ROS，构造栈上 `SlamSystem`，`Init(config)` 后执行 `StartSLAM("new_map")`，再 `Spin`。可选 `--bag` 创建内部回放线程，播放结束等待 `--post_wait_seconds` 后 shutdown；这种模式的回放等待不是所有算法完成的证明。

@dot
digraph online_slam_path { rankdir=LR; node [shape=box]; ros [label="ROS executor"]; lio [label="ProcessIMU / ProcessLidar\nDrainLio"]; kf [label="new keyframe"]; backend [label="Backend worker / legacy"]; out [label="map->odom / G2P5 / UI"]; save [label="SaveMap + final optimization"]; ros -> lio -> kf -> backend -> out; backend -> save; }
@enddot

# 实际线程边界

`SlamSystem::ProcessIMU/ProcessLidar` 直接驱动前端，当前类没有在线定位的 sensor `AsyncMessageProcess`。新后端 `BackendPipeline::Init(..., true)` 创建 worker，HBA 可再有独立线程。优化回调发布 map→odom，并可触发 G2P5 重绘；这一回调可能运行在后端线程。

`DrainLio` 只保存有效且严格递增的 LIO 状态；比较关键帧指针，只有新关键帧才分发给后端、G2P5 和 UI。系统 `with_loop_closing` 与 `backend.mode` 共同决定是否创建后端，单独改 mode 并不保证启动。

@anchor online_slam_save_reading
# 保存失败时，旧地图一定还在吗？

服务 `lightning/save_map` 使用请求 map_id；`lightning/optimize_backend` 请求全局优化。`SaveMap` 先等待后端（可能运行 HBA），检查关键帧、估计导出坐标，再写地图、轨迹、诊断与数据库。**现有目标目录会被清空重建**，它不是原子替换；应选择新的输出目录。

先读 `SaveMap(path)` 的前置条件：

@snippet{lineno} core/system/slam.cc online-slam-save-preconditions

`WaitUntilIdle(true)` 在获取关键帧之前执行，可能触发最终优化。关键帧为空则返回 false；下一步的地面估计也能在目录修改之前拒绝保存：

@snippet{lineno} core/system/slam.cc online-slam-export-frame-check

但这两道检查通过后，旧目录就会被删除。片段的后半段同时展示了删除之后仍会失败的位置：

@snippet{lineno} core/system/slam.cc online-slam-replace-and-write

<b>状态推演：</b>假设目标目录已有地图。若关键帧为空或地面估计失败，上面两段直接返回，尚未走到删除；若删除重建已完成，而 `ConvertFromFullPCD()` 或全局 PCD 写入失败，则返回 false 时旧地图已被删除，新目录可能不完整。这段函数没有保存旧目录再回滚的分支。

服务把 `SaveMap(path)` 的布尔结果映射成响应值 0/3；命令行 `--save_map` 的失败返回码则是 4。两者是不同接口的返回约定，不能混用。地图保存后还会写轨迹、诊断、数据库与 metadata/manifest；只看到某一个文件也不足以确认整次保存成功。

<b>读到这里：</b>能指出最后一个“尚未修改目录”的失败点和第一个“已经替换目录”的语句即可。等待后端不等于停止一切传感器生产，外部并发调用者仍需遵守系统同步约束。完整代码见 @ref lightning::SlamSystem "SlamSystem 的 SaveMap 两个重载"、@ref slam.cc "slam.cc"；地图内容的坐标与加载校验见 @ref map_contracts "地图包契约"，离线导出对照 @ref offline_slam_export_reading "轨迹重写与清单收尾"。

## 退出

`--save_map` 在 Spin 退出并 join 回放线程后保存，失败返回 4。析构等待后端 idle、Shutdown，最后退出 UI；并不存在“先 join 异步 LIO 队列”这一步。

**源码依据：** @ref lightning::SlamSystem "SlamSystem" 的 `Init`、`DrainLio`、`SaveMap`、析构；@ref run_slam_online.cc "run_slam_online.cc"。
模块与原理：@ref backend_module "后端、回环与优化"、@ref backend_optimization "体素 BA、回环位姿图与分层优化"；地图契约：@ref map_contracts "地图包、坐标归一化与栅格导出"。
