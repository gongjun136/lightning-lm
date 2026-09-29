@page pose_graph_theory 位姿图、增量求解与高频平滑

# 两类图优化

建图图优化调整历史关键帧；定位图优化融合地图绝对约束与 LIO/DR 相对运动，并滑动窗口。两者使用 miao，但生命周期和失败处理不同，分别从 @ref backend_module "后端、回环与优化"、@ref localization_module "定位与重定位" 进入。

# 当前边残差不是教科书的 SE(3) 对数

设测量为 \f$Z_{ij}\f$、估计相对位姿为 \f$T_i^{-1}T_j\f$：
\f[
\Delta=Z_{ij}^{-1}T_i^{-1}T_j,\qquad
e=\begin{bmatrix}\operatorname{translation}(\Delta)\\
\operatorname{vec}(q(\operatorname{rotation}(\Delta)))\end{bmatrix}.
\f]
@ref lightning::miao::EdgeSE3 "EdgeSE3" 使用单位四元数的 xyz 部分；小角度下约为半个旋转向量。信息矩阵的旋转权重必须与此尺度一起理解，不能直接套用以 rad 为残差的 `Log` 权重。顶点采用世界平移相加、旋转右乘，见 @ref geometry "坐标系、右扰动与李群雅可比"。

目标是 \f$\sum_e\rho(e^T\Omega_e e)\f$。线性化后组装稀疏 H、b，求解增量；LM 增加阻尼并比较试探步的代价。固定顶点消除规范自由度；鲁棒核减小异常边影响，但不会让错误闭环自动无害。

# 定位 PGO 的实现顺序

@ref lightning::loc::PGOImpl "PGOImpl" 的 `AddPGOFrame` 先检查时间，尝试给该帧分配 LIO/DR 相对位姿，再检查地图观测有效性。优化成功后收集结果并滑窗，重建 frame-id 辅助索引。失败分支清理图和输出队列、将结果标为 FAIL，不能发布上一次结果来冒充本次成功。

`Reset` 清除相对/绝对/输出队列并记录旧时间水位，防止重定位后混入上一段时间的运动约束。图的增量复用只是求解器状态复用，不意味着数据可以乱序输入。

# 高频输出与平滑

PGO 地图匹配频率低于 IMU。`PGO::PubResult` 从最近有效地图结果结合相对运动外推，再检查地图支持和输出门控；高频消息的数量不等于独立地图观测数量。

`PoseSmoother` 在 DR 有效时，以相邻 DR 运动预测 \f$T_{pred}\f$，再向目标位姿插值：
\f[
t_{out}=(1-\alpha)t_{pred}+\alpha t_{target},\quad
R_{out}=R_{pred}\operatorname{Exp}\left(\alpha\operatorname{Log}(R_{pred}^{-1}R_{target})\right).
\f]
前几次直接输出；DR 无效也直接采用目标。正常路径会按平面距离把下一次平滑因子设成 0.01 或 0.2，偏离过大时跳到目标并清队列。因此构造参数不代表全生命周期固定因子。

DR 距离门限还依赖时间间隔和速度：`max(0.3, 1.5*speed*dt+0.1)`。小的平滑因子减少瞬时变化但带来滞后，不能用平滑遮盖失效地图观测。独立的速度影子诊断见 @ref speed_smoothing_shadow "速度平滑影子诊断"，仅记录候选输出。

**源码导航：** @ref edge_se3.h "edge_se3.h"、@ref vertex_se3.h "vertex_se3.h"、@ref opti_algo_lm.cc "opti_algo_lm.cc"、@ref pgo_impl.cc "pgo_impl.cc"、@ref pgo.cc "pgo.cc"、@ref smoother.h "smoother.h"、@ref pose_extrapolator.cc "pose_extrapolator.cc"。
