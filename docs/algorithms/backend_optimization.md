@page backend_optimization 体素 BA、回环位姿图与分层优化

# 先看系统边界

@ref backend_module "后端、回环与优化" 解释调度、锁与提交时机，本页解释目标函数。局部 BA 修正滑窗几何；BTC 提供远距离闭环候选；PGO 调整全局位姿关系；HBA 在更大尺度重新优化点云几何。四者不能用同一个“后端优化”标签代替。

# 平面体素为何可以消去平面参数

将多个关键帧的点变到同一坐标系。设体素内 N 个点的均值和协方差为
\f[
\mu=\frac1N\sum q_j,\quad C=\frac1N\sum q_jq_j^T-\mu\mu^T.
\f]
对单位法向 n，最佳平面经过均值，平均平方点面距离为 \f$n^TCn\f$。最小化该 Rayleigh 商得到最小特征值 \f$\lambda_0(C)\f$，对应特征向量就是法向。于是位姿目标可写成
\f$E=\sum_{v\in\mathcal V}\lambda_0(C_v(T_1,\ldots,T_m))+E_{smooth}\f$。
当前代码对每个有效体素累加一个最小特征值，不额外乘点数 N；不能把“体素均权”误写成“逐点等权”。

`ProcessBuckets` 按当前位姿分桶，要求足够点数与至少指定帧数，使用 \f$\lambda_0/\max(\lambda_1,\epsilon)\f$ 判断平面性；不满足时在深度限制内细分。`BuildFactors` 在一次 `Optimize` 开始构建分桶，迭代重算统计而不每步重新分桶。

# 充分统计与位姿导数

每帧每体素仅保存 \f$N,s=\sum p,M=\sum pp^T\f$。刚体变换后
\f[
s'=Rs+Nt,\quad M'=RMR^T+(Rs)t^T+t(Rs)^T+Ntt^T.
\f]
这些量可直接相加构造整体协方差，避免每轮变换全部点。`PointCluster::Transformed/Covariance` 实现上述公式，分母是 N，而不是无偏估计的 N-1。

简单特征值的微分满足
\f[
d\lambda_0=n^T(dC)n,\quad
d^2\lambda_0=n^T(d^2C)n+
2\sum_{k=1,2}\frac{(n_k^TdCn)^2}{\lambda_0-\lambda_k}.
\f]
`EvaluateFactors` 的 `eigen_coupling` 对应第二项；相近特征值时截断分母。这里不仅是固定法向的点面 Gauss–Newton 近似，还考虑法向随状态改变的耦合。

# 求解与接受

BA 的每帧增量是 **旋转在前、平移在后**。首帧块固定以去掉全局刚体规范自由度。对当前 Hessian H 和梯度 g：
\f[(H+\lambda D)\delta=-g,\quad D_{ii}=\max(|H_{ii}|,10^{-9}).\f]
旋转右乘、平移直接相加。每帧步长先裁剪，再检查相对先验的总修正；只有代价下降且预测下降量为正才接受，接受后减小阻尼，拒绝则增大阻尼。`AddSmoothness` 约束相邻帧修正差；它不是增加 GNSS 真值观测。

关键帧保存 IMU 位姿，BA 用 LiDAR 点：`OptimizeKeyframes` 先乘 `T_imu_lidar`，写回前乘逆外参。忘记这一层会造成杆臂相关的假运动。

# 回环图与 HBA

BTC 提供相对 LiDAR 变换，经过检测和提交门控后转成 IMU 约束。图中邻接 LIO 边维持局部形状，回环边闭合累计漂移。残差参数化见 @ref pose_graph_theory "位姿图、增量求解与高频平滑"；当前 miao 旋转误差不是完整对数。

HBA `BuildNode` 以组内最后一个关键帧作锚，将成员点变到锚 LiDAR 系，优化上层节点后把修正传播回成员。分层降低一次联合优化的变量规模，但不保证全局最优。`BackendPipeline::HbaLoop` 会检查 `require_applied_loop_for_commit`；即使求解成功，也可能恢复原位姿而不提交。

**源码导航：** @ref voxel_bundle_adjustment.cc "voxel_bundle_adjustment.cc"（`PointCluster`、`EvaluateFactors`、`AddSmoothness`、`Optimize`），@ref hierarchical_bundle_adjustment.cc "hierarchical_bundle_adjustment.cc"（`BuildNode`、`ApplyCorrections`），@ref lightning::backend::BackendPipeline "BackendPipeline"。
