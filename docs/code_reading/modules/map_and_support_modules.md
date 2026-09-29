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

@ref lightning::RosbagIO "RosbagIO" 通过 `AddHandle(topic, callback)` 建立 topic 到回调的映射，`Go()` 用 reader 同步遍历消息并在当前调用线程派发。未注册 topic 不进入算法；回调返回值参与读取流程的继续/失败语义。`YAML_IO` 和 `file_io` 提供同步文件操作，不拥有后台线程。

**代码依据：** `src/wrapper/bag_io.h/.cc`、`src/io/yaml_io.h/.cc`、`src/io/file_io.h/.cc`

# miao 优化器

`src/core/miao` 被编译为 `miao.core` 与 `miao.utils` 两个共享库。调用方通过配置选择算法和线性求解器，向图加入顶点/边，再执行初始化与迭代；SLAM 后端和定位 PGO 都复用这套图优化抽象。顶点、边和鲁棒核主要用 `shared_ptr` 交给图持有。

**代码依据：** `src/core/miao/CMakeLists.txt`、`src/core/miao/core/graph/optimizer.h`、`graph.h`、`src/core/backend/backend_pipeline.cc`（优化器装配和使用）

# UI、公共模型与运行时工具

| 模块 | 作用 | 生命周期/副作用 |
|---|---|---|
| `src/ui` | Pangolin 窗口、车体/点云/轨迹对象 | `Init()` 创建可视资源，`Quit()` 结束；无 UI 模式不得假设窗口存在 |
| `src/common` | `NavState`、`Keyframe`、参数与点类型 | `Keyframe::Ptr` 跨前后端共享；优化位姿可在后端线程更新 |
| `src/utils` | timing、profiling、预算、线程辅助 | 多数采集运行统计；启用与否由配置决定，不应改变算法输入语义 |
| `src/io` | YAML/文件序列化 | 同步磁盘副作用；错误需沿返回值或异常边界处理 |

**代码依据：** `src/ui`、`src/common/keyframe.h`、`src/common/nav_state.h`、`src/utils`、`src/io` 的公开头文件及 `src/CMakeLists.txt` 源清单

# 阅读时的跨模块约束

- 坐标系名称必须和变换方向一起读；不能仅从变量含 `map`/`odom` 猜乘法方向。
- `Keyframe` 同时携带 LIO pose 与 optimized pose；地图导出应使用哪一个由流程阶段决定。
- PCL shared pointer 表示共享所有权，不等于对象内容线程安全；跨线程写点云仍需调用侧同步。
- bag 回放是同步派发，但在线 ROS 输入不是；复用代码时必须保留两种时序差异。
- 文件写出、地图块换入换出、UI 和 ROS 发布都是可观察副作用，文档变更检查不能只看算法类。

**推断：** `src/utils` 中某些诊断工具按设计应只观测而不影响结果；是否完全无时序影响需要性能测试，源码静态阅读只能确认其调用边界。

