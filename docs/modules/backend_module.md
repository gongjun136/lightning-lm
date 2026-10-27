@page backend_module 后端、回环与优化

# 职责与入口

新后端 @ref lightning::backend::BackendPipeline "BackendPipeline" 接收 `Keyframe::Ptr`，按配置组合局部体素 BA、BTC 回环检测、位姿图优化和分层 BA。`backend.mode` 决定系统使用新后端、旧 @ref lightning::LoopClosing "LoopClosing"，或关闭后端；模式解析失败时默认回退到 `ba_btc_hba`。

本页验证四个问题：@ref backend_queue_ba_reading "关键帧何时触发局部 BA" → @ref backend_loop_commit_reading "BTC 结果是否真的进入图" → @ref backend_hba_transaction_reading "HBA 为什么会回滚" → @ref backend_map_callback_reading "map→odom 与回调怎样更新"。先对照片段中的条件、容器和锁作用域，目标函数推导放在 @ref backend_optimization "体素 BA、回环位姿图与分层优化"。

<b>模式入口：</b> @ref backend_pipeline.cc "backend_pipeline.cc" 的 `ReadBackendMode()`，以及 @ref slam.cc "slam.cc" 中依据模式构造后端的位置。以下算法主线不适用于 legacy 后端。

legacy 的 @ref lightning::LoopClosing "LoopClosing" 在线使用公共 `sys::AsyncMessageProcess`，等待队列容量固定为 40，超限淘汰最旧待处理关键帧；析构排空保留队列并等待 worker 结束。`SlamSystem` 在释放地图及 UI 前销毁回环对象，以保持尾帧重绘回调有效。离线回环仍同步处理，详见 @ref online_systems "队列契约与退出顺序"。

队列淘汰后，关键帧 ID 不一定连续。`HandleKF()` 保存 ID 到已处理序列下标的映射；`PoseOptimization()` 连接前两个实际处理过的关键帧，候选子图和优化结果通过该映射查找。ID 仍作为图顶点标识，不直接作为数组下标。

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

@anchor backend_queue_ba_reading
# 1. 送入一个关键帧，就会立即执行局部 BA 吗？

先看入口怎样分流。离线调用者直接执行 `HandleKeyframe()`；在线调用者只负责向队尾插入，由 worker 处理。

@snippet{lineno} core/backend/backend_pipeline.cc backend-enqueue-mode

在线 worker 取队头后退出 `queue_mutex_` 的作用域，再调用耗时的 `HandleKeyframe()`。这段代码证明了顺序处理，不包含“新帧替换旧帧”的容量策略：

@snippet{lineno} core/backend/backend_pipeline.cc backend-worker-dispatch

例如 A 正在处理时，B、C 依次入队，等待队列保留 `[B, C]`，随后先处理 B。不要套用在线定位中容量为 1 的定位队列语义。这里也没有显式容量上限；持续积压是否发生，需要运行数据确认。

接下来只找 `local_window` 如何得到内容。以下判断发生在关键帧加入 `keyframes_` 之后：

@snippet{lineno} core/backend/backend_pipeline.cc backend-local-window

读两个 `>=`：总帧数足够形成窗口，而且距上次触发的新增帧数达到 stride，才选择末尾 `window` 帧。实际调用还要通过 `local_ba.enabled`：

@snippet{lineno} core/backend/backend_pipeline.cc backend-local-ba-run

<b>手算：</b> 假设 BA 启用，窗口为 4、stride 为 2，并连续送入不同关键帧。总帧数 1–3 时窗口不足；第 4 帧触发，`last_local_ba_keyframe_count_` 变为 4；第 5 帧距上次只增加 1，不触发；第 6 帧再次触发，使用末尾第 3–6 帧。`local_ba_attempts` 记录调用次数，`local_ba_accepted` 只在求解结果接受时增加，两者不能互换。

<b>读到这里可停：</b> 你已经能确认“何时排队、何时尝试 BA”。`HandleKeyframe()` 还会忽略与最后一帧相同的指针；该去重不是按时间戳全局去重。要理解一次 BA 是否接受步长，再读 @ref backend_optimization "求解与接受"。完整入口：@ref lightning::backend::BackendPipeline::AddKeyframe() "AddKeyframe()"、@ref lightning::backend::BackendPipeline::HandleKeyframe() "HandleKeyframe()"。

@anchor backend_loop_commit_reading
# 2. BTC 检测通过，是否说明回环已经影响地图？

沿 `accepted`、`apply_constraint`、`graph_inlier` 分三步读。第一步在 `HandleKeyframe()`：检测通过还要满足 `optimization_warranted`，才向 `loop_constraints_` 插入约束。

@snippet{lineno} core/backend/backend_pipeline.cc backend-loop-constraint-commit

`loop_results_` 保存检测记录，`loop_constraints_` 保存准备交给位姿图的约束。中间的外参共轭 `T_imu_lidar_ * T_history_lidar_current_lidar * T_imu_lidar_.inverse()` 将 LiDAR 相对变换转换成关键帧使用的 IMU 相对变换。

| 检测结果 | 本次计数变化 | 是否请求后续优化 |
|---|---|---|
| `accepted=false` | `loops_accepted`、`loops_applied` 均不增加 | 不由此候选触发 |
| `accepted=true`，`optimization_warranted=false` | 只增加 `loops_accepted` | 不由此候选触发 |
| 两者均为 true | 两个计数都增加，保存 IMU 约束 | 请求 PGO，并请求 HBA |

第二步看实际触发位置：

@snippet{lineno} core/backend/backend_pipeline.cc backend-loop-optimization-trigger

`OptimizePoseGraph()` 的返回值没有在这里作为 HBA 请求条件；HBA 请求依赖的是 `apply_constraint`。若 HBA 已禁用，这个 `RequestHba(false, ...)` 会直接返回，并不执行 HBA。

第三步进入 @ref lightning::backend::BackendPipeline::OptimizePoseGraph() "OptimizePoseGraph()"。它先构造相邻 LIO 边和回环边，固定首顶点，回环边带 Cauchy 核；第一轮求解后按下面的卡方条件标记内点：

@snippet{lineno} core/backend/backend_pipeline.cc backend-loop-graph-inliers

<b>例子：</b> 配置 `loop_outlier_chi2=25` 时，有限的 `Chi2=30` 会使这条边的 `graph_inlier` 为 false、level 变为 1；若还有其他内点，第二轮只优化 level 0 的边。`loops_applied` 已经在约束插入时增加，所以不能用它证明这条边在图里仍是内点。若没有任何内点，代码跳过第二轮，但仍执行下面的顶点位姿写回；这里没有“整轮 PGO 恢复原位姿”的分支。

<b>边界与停止位置：</b> 这三处足以区分“检测通过、约束入容器、PGO 判为内点”。`OptimizePoseGraph()` 还可能因关键帧不足、找不到回环顶点或初始化失败提前返回；因此保存检测记录也不等于优化已经完成。要追 BTC 的几何确认或修正幅度阈值，再进入 @ref lightning::backend::BtcLoopDetector "BtcLoopDetector"。

@anchor backend_hba_transaction_reading
# 3. HBA 已接受结果，为什么最终仍恢复旧位姿？

先确认 `has_applied_loop` 的来源。HBA 线程在本轮开始时读取累计入图约束计数：

@snippet{lineno} core/backend/backend_pipeline.cc backend-hba-loop-history

<b>这项条件不是最新 PGO 内点数。</b> `runtime_.loops_applied > 0` 说明历史上曾向约束容器加入过回环，并不保证最新 `loops_graph_inliers > 0`。然后看快照、求解和恢复三步在同一个 `optimization_mutex_` 作用域中的顺序：

@snippet{lineno} core/backend/backend_pipeline.cc backend-hba-transaction

读法是：先按请求重跑 PGO，再保存 `poses_before_hba`，然后执行 HBA。只有 `summary.accepted`、`require_applied_loop_for_commit`、`!has_applied_loop` 同时成立时，这个额外提交门控才恢复本轮关键帧快照，并把结果改为拒绝。

<b>状态例子：</b> 假设某个关键帧保存快照时的 x 为 11，HBA 改成 12。若本轮读取的 `has_applied_loop=false`，提交要求开启且 HBA 原本接受，那么恢复到的是快照中的 11。快照保存的位置在可选 PGO 之后，因此这段恢复逻辑不会撤销先前已经完成的 PGO 更新。若关闭 `require_applied_loop_for_commit`，这项回滚条件不成立；是否真正接受仍由 HBA 自身结果决定。数字只说明控制流，不代表实际优化步长。

<b>读到这里可停：</b> 已能确认该门控只撤销本轮 HBA，对快照中的非空关键帧逐一恢复；不应推断它能回滚整个 pipeline 的所有状态。`optimization_mutex_` 串行化 BA/PGO/HBA 写位姿，也不能由此推断所有外部读者都有整批位姿的原子快照。完整实现：@ref lightning::backend::BackendPipeline::HbaLoop() "HbaLoop()"；HBA 内部如何接受各层结果见 @ref lightning::backend::HierarchicalBundleAdjuster "HierarchicalBundleAdjuster"。

@anchor backend_map_callback_reading
# 4. map→odom 怎样更新，回调时到底还持有哪些锁？

读下面完整的短函数，重点看乘法顺序以及内层花括号结束的位置：

@snippet{lineno} core/backend/backend_pipeline.cc backend-map-callback

`GetOptPose()` 是关键帧优化后的 IMU 位姿，`GetLIOPose()` 是同一关键帧的 LIO 位姿。二者组合得到
\f[T_{map,odom}=T_{map,imu}(T_{odom,imu})^{-1}.\f]
例如两者旋转均为单位旋转，x 分别为 12 m 与 10 m，则 `T_map_odom_` 的 x 平移为 2 m；把 odom 中 x=10 的 IMU 点变换后，正好回到 map 中的 x=12。这里始终用后端当前最后一个关键帧配对，不是拿不同帧的两条轨迹相减。

`callback` 在持有 `data_mutex_` 时复制，退出该作用域后调用，因此回调内调用同样只获取 `data_mutex_` 的 `GetMapToOdom()` 不会因这把锁直接自锁。不过上面的 BTC 触发片段和 HBA 提交片段都显示：调用者可能仍持有 `optimization_mutex_`。<b>“释放数据锁”不能扩展成“回调完全无锁、可以任意重入”。</b> 回调中再次进入需要同一优化锁的同步优化路径，或等待本线程负责完成的后端空闲条件，都可能自锁或自等待。

<b>读到这里可停：</b> 可以把回调用于读取当前变换并发出通知，不把回调发生次数当成优化成功次数；`HandleKeyframe()` 末尾即使没有接受新的 BA/回环，也会调用本函数。查看实际订阅方时，再进入 @ref slam.cc "SlamSystem 的 SetOptimizedCallback() 调用处"。完整接口：@ref lightning::backend::BackendPipeline::GetMapToOdom() "GetMapToOdom()"、@ref lightning::backend::BackendPipeline::SetOptimizedCallback() "SetOptimizedCallback()"。

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

`WaitUntilIdle(true)` 会在队列空闲后请求一次全局优化，因此不是纯查询。关闭之前应先结束输入生产；这里的等待不能代替调用者停止持续入队。

<b>代码依据：</b> @ref backend_pipeline.h "backend_pipeline.h"（成员所有权、锁和线程）、@ref backend_pipeline.cc "backend_pipeline.cc"（析构、`WaitUntilIdle()`、`Shutdown()`）。
