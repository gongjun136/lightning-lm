@page backend_module 后端、回环与优化

# 职责与入口

新后端 @ref lightning::backend::BackendPipeline "BackendPipeline" 接收 `Keyframe::Ptr`，按配置组合局部体素 BA、BTC 回环检测、位姿图优化和分层 BA。`backend.mode` 决定系统使用新后端、旧 @ref lightning::LoopClosing "LoopClosing"，或关闭后端；模式解析失败时默认回退到 `ba_btc_hba`。

**代码依据：** `src/core/backend/backend_pipeline.cc`（`ReadBackendMode()`、`Init()`、`HandleKeyframe()`），`src/core/system/slam.cc`（依据模式构造后端）

# 数据流

@dot
digraph backend_flow {
  rankdir=LR;
  node [shape=box];
  kf [label="Keyframe\nLIO pose + cloud"];
  local [label="local voxel BA\nwindow/stride gate"];
  btc [label="BTC descriptor\n候选与几何验证"];
  pgo [label="miao pose graph\nmotion + loop edges"];
  hba [label="hierarchical BA\ntransactional commit"];
  out [label="optimized pose\nT_map_odom + callback"];
  kf -> local -> btc;
  btc -> pgo [label="accepted &&\noptimization_warranted"];
  pgo -> hba;
  local -> out;
  pgo -> out;
  hba -> out;
}
@enddot

1. `AddKeyframe()` 在离线模式直接调用 `HandleKeyframe()`；在线模式只入 `queue_`，由 worker 串行取出。
2. 关键帧加入 `keyframes_` 后，只有窗口大小和 stride 同时满足才运行局部 BA。
3. BTC 的“检测通过”不等于一定写入图：还要求 `optimization_warranted`。写入的 LiDAR 相对位姿经外参转换成 IMU 相对位姿。
4. 位姿图用相邻 LIO 约束和回环约束；首顶点固定，回环边使用 Cauchy 核并按卡方阈值剔除异常边。
5. HBA 在独立线程运行；若配置要求已有生效回环但条件不满足，代码恢复 HBA 前的全部位姿，拒绝提交。
6. 每轮可能改变优化位姿后更新 `T_map_odom = T_map_imu * inverse(T_odom_imu)`，随后复制并调用通知回调。

**代码依据：** `src/core/backend/backend_pipeline.cc`（`AddKeyframe()`、`HandleKeyframe()`、`OptimizePoseGraph()`、`HbaLoop()`、`UpdateMapToOdomAndNotify()`）

# 关键约束与失败边界

- `AddKeyframe()` 对未初始化 pipeline 或空指针静默返回；调用者不能以返回值判断是否接收，因为接口返回 `void`。
- 在线队列无显式容量上限；吞吐监控应看运行统计和调用侧时延，而不是依赖丢帧语义。
- `WaitUntilIdle(true)` 可能触发一次全局优化，因此它不是纯查询函数。
- BTC 只在约束实际进入图后请求 HBA，避免对不会提交的结果执行昂贵点级优化。
- 优化回调在不持有 `data_mutex_` 时执行，避免回调重新读取后端状态导致自锁。

**推断：** 在线输入持续快于后端时，`queue_` 会增长并增加延迟；源码未给出容量或丢弃策略，但是否在部署中发生需由运行统计验证。

**代码依据：** `src/core/backend/backend_pipeline.cc`（队列插入、HBA 请求门控、回调复制后调用），`src/core/backend/backend_pipeline.h`（队列类型）

# 继续跳转

- @ref lightning::backend::BackendPipeline::AddKeyframe() "AddKeyframe()"
- @ref lightning::backend::BtcLoopDetector "BtcLoopDetector"
- @ref lightning::backend::VoxelBundleAdjuster "VoxelBundleAdjuster"
- @ref lightning::backend::HierarchicalBundleAdjuster "HierarchicalBundleAdjuster"
- @ref lightning::miao::Optimizer "miao::Optimizer"

# 深入阅读

@ref backend_optimization "体素 BA、回环位姿图与分层优化"、@ref pose_graph_theory "位姿图、增量求解与高频平滑"、@ref backend_configuration "后端配置与评测协议"。

# 按需参考：对象管理与并发

修改对象初始化、线程或退出逻辑时核对本节；首次阅读优先掌握上面的主路径与算法对应。

## 生命周期、线程与锁

| 资源 | 创建条件 | 保护/所有权 | 结束条件 |
|---|---|---|---|
| local BA、BTC、HBA 对象 | `Init()` 成功解析 YAML | `unique_ptr` 独占 | 随 pipeline 析构 |
| worker 线程 | `online_mode=true` | `queue_mutex_` + `queue_cv_` | 队列 drain 后设置 `worker_stop_` 并 join |
| HBA 线程 | `hba.enabled=true` | `hba_mutex_` + `hba_cv_` | 当前请求完成后设置 `hba_stop_` 并 join |
| 关键帧/约束/统计 | 首个关键帧到达后 | `data_mutex_` | pipeline 释放 |
| 优化器临界区 | BA/PGO/HBA 执行期间 | `optimization_mutex_` 串行化 | 单次优化返回 |

析构函数调用 `Shutdown()`；`Shutdown()` 先 `WaitUntilIdle(false)`，再唤醒并 join 两类线程。调用者若需要输入结束时的全局优化，必须在关闭前显式调用 `WaitUntilIdle(true)`。

**代码依据：** `src/core/backend/backend_pipeline.h`（成员所有权、锁和线程），`src/core/backend/backend_pipeline.cc`（析构、`WaitUntilIdle()`、`Shutdown()`）
