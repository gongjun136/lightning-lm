@page pose_graph_theory 位姿图、增量求解与高频平滑

# 两类图优化

建图图优化调整历史关键帧；定位图优化融合地图绝对约束与 LIO/DR 相对运动，并滑动窗口。两者使用 miao，但生命周期和失败处理不同，分别从 @ref backend_module "后端、回环与优化"、@ref localization_module "定位与重定位" 进入。

本页按四个问题阅读：边的残差尺度、入帧与优化失败、重置后的时间边界、平滑的当前输出与下一次因子。每段代码由源码生成，例子用于手算验证，不代表运行测试结果。

@anchor pgo_edge_residual_reading
# 1. 旋转残差究竟是角度，还是四元数的一部分？

先读 @ref lightning::miao::EdgeSE3::ComputeError() "EdgeSE3::ComputeError()" 的最后两次赋值：

@snippet{lineno} core/miao/core/types/edge_se3.h pgo-binary-edge-residual

设测量为 \f$Z_{ij}\f$、估计相对位姿为 \f$T_i^{-1}T_j\f$：
\f[
\Delta=Z_{ij}^{-1}T_i^{-1}T_j,\qquad
e=\begin{bmatrix}\operatorname{translation}(\Delta)\\
\operatorname{vec}(q(\operatorname{rotation}(\Delta)))\end{bmatrix}.
\f]
`from.inverse() * to` 给出估计相对位姿，再左乘测量的逆形成 `delta`。`head<3>()` 取它的平移；`coeffs().head<3>()` 取单位四元数的 xyz 部分。后者对转角 \f$\theta\f$、单位转轴 \f$u\f$ 为 \f$u\sin(\theta/2)\f$，小角度下约为半个旋转向量。

<b>手算例子：</b> 设测量和第一个顶点均为单位位姿，第二个顶点仅绕 z 轴旋转 10°，取正 w 的四元数表示，则平移残差为零，旋转块为 `(0, 0, sin 5°)`，约 `(0, 0, 0.08716)`，不是 `(0, 0, 0.17453)`。

这个结论只对应二元边。定位 PGO 的绝对观测还使用 @ref lightning::miao::EdgeSE3Prior "EdgeSE3Prior"，它有不同的写法：

@snippet{lineno} core/miao/core/types/edge_se3_prior.h pgo-prior-edge-residual

这里平移直接做差，旋转用 `Log(pose.so3().inverse() * measurement.so3())`，只对 SO(3) 取对数。同样令测量为单位位姿、估计绕 z 轴转 10°，旋转残差约为 `(0, 0, -0.17453)`：<b>不仅尺度不同，组合方向也不同</b>。信息矩阵必须与各自的残差定义一起理解，不能把一个边的旋转权重原样解释成另一个边的角度权重。

<b>读到这里：</b> 能指出两种边的最后一次赋值，就已验证公式与代码的对应。小角近似不能扩展成任意转角下严格的二倍关系。顶点更新采用世界平移相加、旋转右乘，继续见 @ref geometry "坐标与扰动约定"；完整代码见 @ref edge_se3.h "二元边源码"、@ref edge_se3_prior.h "先验边源码"。

目标是 \f$\sum_e\rho(e^T\Omega_e e)\f$。线性化后组装稀疏 H、b，求解增量；LM 增加阻尼并比较试探步的代价。固定顶点消除规范自由度；鲁棒核减小异常边影响，但不会让错误闭环自动无害。

@anchor pgo_frame_acceptance_reading
# 2. 一帧地图定位何时能入图，失败是否都清空结果？

沿 @ref lightning::loc::PGO::ProcessLidarLoc() "ProcessLidarLoc()" → `ProcessPGOFrame()` → @ref lightning::loc::PGOImpl::AddPGOFrame() "AddPGOFrame()" 进入。先看入帧检查中的比较符号：

@snippet{lineno} core/localization/pose_graph/pgo_impl.cc pgo-frame-time-check

这里拒绝的是 `< 0` 的时间回退；与上一帧同时间戳不会被<b>这一条</b>判断拒绝。不要把它读成“严格递增”，也不能据此认为同时间戳帧必然能通过其他检查。

接着核对相对与绝对观测：

@snippet{lineno} core/localization/pose_graph/pgo_impl.cc pgo-frame-observation-check

`!interp_lio_success && !interp_dr_success` 意味着两路相对位姿都分配失败才拒绝，任一路成功可以继续。随后 `lidar_loc_set_ && lidar_loc_valid_` 要求地图观测既已设置又有效；有相对运动并不能替代地图约束。这些提前返回尚未进入下面的优化失败处理。

真正调用 `RunOptimization()` 后，如果它返回 false，则执行另一条路径：

@snippet{lineno} core/localization/pose_graph/pgo_impl.cc pgo-optimization-failure

日志可以略读，重点依次看 `result_ = LocalizationResult{}`（默认 `valid_ = false`）、本帧时间戳和 `FAIL`、清图和清输出历史、最后 `return false`。调用者 `ProcessPGOFrame()` 在 `AddPGOFrame()` 返回 false 时也立即返回，跳过本次成功回调。这里没有清除 DR/LIO 队列；它与下一节的完整 `Reset()` 不同。

<b>手算例子：</b> 上一帧时刻为 10 s，本帧为 9.9 s，第一段直接拒绝；若本帧为 10.1 s，LIO 分配成功、DR 失败且地图观测有效，则可以进入优化。假设这次优化失败，结果时刻变成 10.1 s、状态为 `FAIL`，而非保留旧结果冒充本次成功。这三个场景不能统称为同一种失败处理。

<b>阅读边界：</b> `AddPGOFrame()` 的时间比较不等于完整输入合法性校验；外层还有重置水位、静止保持等判断。优化成功后才收集结果、调用 `SlideWindowAdaptively()` 并重建 frame-id 索引。判断“能否入图、失败影响什么”到此即可停止；求解失败的具体来源再读 @ref lightning::loc::PGOImpl::RunOptimization() "RunOptimization()"，完整实现见 @ref pgo_impl.cc "pgo_impl.cc"。

@anchor pgo_reset_watermark_reading
# 3. Reset 清空了队列，为什么还会拒绝随后到达的定位结果？

先在 @ref lightning::loc::PGOImpl::Reset() "PGOImpl::Reset()" 中看清空之前保存了什么：

@snippet{lineno} core/localization/pose_graph/pgo_impl.cc pgo-reset-watermark

`relative_pose_reset_watermark_` 取本次重置时 DR 与 LIO 队尾时间的最大值，空队列不参与。图和帧索引随后被清理，四类历史队列与结果也被清空：

@snippet{lineno} core/localization/pose_graph/pgo_impl.cc pgo-reset-history

再跳到 @ref lightning::loc::PGO::ProcessLidarLoc() "ProcessLidarLoc()" 消费这个值的位置：

@snippet{lineno} core/localization/pose_graph/pgo.cc pgo-reset-watermark-gate

这次比较是 `<=`。它阻止晚到的地图匹配结果重新进入已重置的时间段；队列清空只移除了 PGO 内部已有样本，并不会自动撤销外部队列里尚未到达的结果。

<b>手算例子：</b> Reset 前 DR 队尾为 12.0 s，LIO 队尾为 12.1 s，水位为 12.1 s。此后收到 12.05 s 或恰好 12.1 s 的地图结果，均在此处返回 false；12.2 s 仅通过这一道时间门，仍需要有效地图与重新建立的相对运动支持。

<b>读到这里：</b> 能找到“Reset 写水位 → 接收地图结果读水位”即可确认作用。水位在每次 Reset 中重新从 `-1` 开始计算；若两条相对队列为空，会保留 `-1`，不能把它理解为永不下降的全局时间上限。对外应调用 @ref lightning::loc::PGO::Reset() "PGO::Reset()"，它还清理平滑器和静止保持状态；完整路径见 @ref pgo.cc "pgo.cc" 与 @ref pgo_impl.cc "pgo_impl.cc"。

@anchor pgo_smoothing_reading
# 4. 平滑因子作用于哪一帧，什么时候会直接输出目标？

PGO 地图匹配频率低于 IMU。`PGO::PubResult()` 以已有有效融合结果为基准，在 `dr_extrapolation_enabled_` 启用时调用 `ExtrapolateLocResult()`，再按平滑开关处理并调用高频回调；高频消息的数量不等于独立地图观测数量。

先看 @ref lightning::loc::PoseSmoother::PushPose() "PushPose()" 的入口。每次历史不足三帧时，直接输出输入 `pose`：

@snippet{lineno} core/localization/pose_graph/smoother.h pgo-smoothing-warmup

历史足够后才进入下面的主分支。`pose_queue_.back()` 是上次平滑输出，`dr_queue_` 中的相邻位姿给出相对运动；这两者先合成预测 \f$T_{pred}\f$，再向本次目标 `pose` 插值：

@snippet{lineno} core/localization/pose_graph/smoother.h pgo-smoothing-prediction

\f[
t_{out}=(1-\alpha)t_{pred}+\alpha t_{target},\quad
R_{out}=R_{pred}\operatorname{Exp}\left(\alpha\operatorname{Log}(R_{pred}^{-1}R_{target})\right).
\f]

`motion_effective_` 为 false 时直接采用目标，公式不执行；DR 样本不足两条但该标志为 true 时，使用单位运动，仍会插值。最后才检查输出与目标的平面距离，选择后续因子：

@snippet{lineno} core/localization/pose_graph/smoother.h pgo-smoothing-next-factor

<b>手算例子：</b> 假设已过前三帧、DR 有效、只考虑一维平移，上一输出为 10 m，DR 增量为 1 m，当前目标为 12 m，进入本次调用时因子为 0.2。预测是 11 m，输出为 `11×0.8 + 12×0.2 = 11.2 m`。输出距目标 0.8 m，小于第二门槛 2 m，于是把因子设为 0.01；<b>11.2 m 不会因此重新计算</b>，0.01 留给下一次插值。

<b>边界与停止点：</b> 当前成员值中，平面偏差严格大于 5 m 时直接跳到目标并清空位姿队列，之后再次经历历史不足三帧的路径；偏差大于 2 m 且不超过 5 m 时为下次选 0.2，其余选 0.01。`Reset()` 清历史但不恢复构造时的因子，所以构造参数不是全生命周期常量。理解预测、插值、事后调整的先后顺序后即可停止；只在追查 DR 被拒绝时再读 `PushDRPose()`，完整实现见 @ref smoother.h "smoother.h"。

`PushDRPose()` 的距离门限还依赖时间间隔与两条样本速度模长的最大值：`max(0.3, 1.5*speed*dt+0.1)`。调用者 @ref lightning::loc::PGO::PubResult() "PubResult()" 受 `dr_smoothing_enabled_` 控制，平滑器返回位姿不等于 ROS 已发布；继续核对 @ref output_gate_reading "业务发布门控"。独立的速度影子诊断见 @ref speed_smoothing_shadow "速度平滑影子诊断"。

<b>源码导航：</b> @ref edge_se3.h "edge_se3.h"、@ref vertex_se3.h "vertex_se3.h"、@ref opti_algo_lm.cc "opti_algo_lm.cc"、@ref pgo_impl.cc "pgo_impl.cc"、@ref pgo.cc "pgo.cc"、@ref smoother.h "smoother.h"、@ref pose_extrapolator.cc "pose_extrapolator.cc"。
