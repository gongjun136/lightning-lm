# Lightning-LM 代码阅读首页

Lightning-LM 是一个 ROS 2/ament CMake 工程。生产代码由一个共享核心库、两个内置优化库和十个安装到 ROS 包的程序组成；运行形态分为前端、完整 SLAM、定位、地图转换和辅助工具。本文档按“先运行边界、再数据流、后模块细节”的顺序组织，而不是照目录逐文件复述。

> 文档范围：以当前仓库源码为事实来源。标为 **推断** 的内容是从调用关系或命名推导、但没有被运行时验证的解释。

## 文档导航

### 先建立全局认识

1. @subpage build_and_run "构建、目标与启动"：知道目标如何生成、脚本从哪里启动程序。
2. @subpage architecture "模块地图与入口"：建立模块边界、入口与所有权地图。
3. @subpage script_contracts "Shell 脚本契约"：检查参数、环境、工作目录、调用链和副作用。

### 再沿一条运行流程阅读

- @subpage online_localization_flow "在线定位端到端流程"：topic + map → LIO → 地图匹配/重定位 → PGO → ROS 输出。
- @subpage offline_slam_flow "离线 SLAM 端到端流程"：bag → LIO → 关键帧 → 后端 → 地图产物。

流程页只保留跨模块调用和关键状态转换。后续补充离线定位、在线 SLAM 或离线建图专题时继续放在本层，不把它们混入模块 API 说明。

### 按模块深入

- @subpage laser_mapping_module "LaserMapping 核心模块样板"：前端数据同步、状态更新、关键帧和不变量。
- @subpage backend_module "后端、回环与优化"：关键帧如何进入 BA、BTC、PGO 与 HBA。
- @subpage localization_module "定位与重定位"：定位对象边界、地图加载、结果状态与地图线程。
- @subpage online_systems "在线系统编排"：`SlamSystem`、`LocSystem` 与异步消息处理器的职责。
- @subpage map_and_support_modules "地图与支撑模块"：tile/导航图、2.5D 图、bag/ROS 适配、UI、I/O 和优化器。

### 最后看维护规则

- @subpage documentation_rules "基于代码 diff 的文档更新规则"：判断一次代码修改必须同步哪些文档。

## 从文档跳到 API，再跳到源码

- 文档中的类和函数使用 Doxygen 交叉引用，例如 @ref lightning::LaserMapping "LaserMapping"、@ref lightning::LaserMapping::RunDetailed() "RunDetailed()" 和 @ref lightning::backend::BackendPipeline "BackendPipeline"。
- 点击类或函数进入 API 页面；API 页右侧/底部的定义位置可进入对应源码，文件页提供“转到该文件的源代码”。
- 入口程序不是公共 API，可从 @ref run_slam_offline.cc "run_slam_offline.cc" 的文件页进入带行号源码。

## 运行形态速查

| 形态 | 可执行程序 | 系统编排 | 核心数据路径 |
|---|---|---|---|
| 离线前端 | `run_frontend_offline` | 入口直接持有 @ref lightning::LaserMapping "LaserMapping" | bag → 预处理/融合 → LIO → TUM/PCD |
| 在线前端 | `run_frontend_online` | ROS 2 回调驱动前端 | topic → LIO → topic/UI |
| 离线 SLAM | `run_slam_offline` | 入口直接持有前端和后端 | bag → LIO → keyframe → backend → tiled/global map |
| 在线 SLAM | `run_slam_online` | @ref lightning::SlamSystem "SlamSystem" | topic → 异步前端 → backend/map service |
| 离线定位 | `run_loc_offline` | 离线入口与定位系统 | bag + map → relocalization/tracking → trajectory |
| 在线定位 | `run_loc_online` | @ref lightning::LocSystem "LocSystem" | topic + map → localization → ROS 2 outputs |

在线定位的完整线程和状态路径见 @subpage online_localization_flow "在线定位端到端流程"；离线建图的完整路径见 @subpage offline_slam_flow "离线 SLAM 端到端流程"。

**代码依据：** `src/app/CMakeLists.txt`（`add_executable` 与 `install(TARGETS ...)` 定义生产入口及安装集合）

## 关键约束

- C++17；顶层 CMake 强制 `Release`，并启用 OpenMP、PIC 与调试符号。
- 正常构建不需要 Doxygen；只有 `BUILD_DOCS=ON` 时才查找 Doxygen 和 Graphviz。
- 离线 SLAM 的传感器回调与 `drain()` 同步执行；新后端在离线模式下同步处理关键帧，但 HBA 仍可拥有独立线程。
- `LaserMapping` 的 LiDAR/IMU 队列由同一互斥量保护；队列与时间队列必须保持一一对应。
- 地图导出要求至少产生一个关键帧；程序会在写出前统一转换地图坐标系，并在失败时返回非零状态。

**代码依据：** `CMakeLists.txt`、`src/app/run_slam_offline.cc`、`src/core/lio/laser_mapping.h`、`src/core/backend/backend_pipeline.cc`（构建标志、回调/导出控制流、缓冲区字段和线程创建条件）
