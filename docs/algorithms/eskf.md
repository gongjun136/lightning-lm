@page eskf_theory 18 维迭代 ESKF：从预测到观测注入

# 阅读入口与职责

先读 @ref laser_mapping_module "LIO 前端：LaserMapping" 看一次扫描如何到达滤波器，再读 @ref geometry "坐标系、右扰动与李群雅可比" 确认右扰动。@ref lightning::ESKF "ESKF" 维护名义状态 `x_` 和误差协方差 `P_`；它不负责找最近邻，LiDAR 观测由 `LaserMapping::ObsModel` 回调提供，见 @ref lidar_residuals "LiDAR 残差、信息矩阵与配准参数化"。

# 1. 状态不是早期版本的 12 维

\f[
x=(p,R,v,b_g,b_a,g),\qquad
\delta x=(\delta p,\delta\theta,\delta v,\delta b_g,\delta b_a,\delta g)\in\mathbb R^{18}.
\f]
索引依次是 0、3、6、9、12、15，每块三维。`NavState::dim` 和 `full_dim` 都是 18；外参固定。`prediction_gyro_`、停车状态和置信度是元数据，不增加滤波维度。

重力虽然存三维，却不是无约束三自由度：`ApplyGravityDelta` 去掉径向分量，将单次切向修正裁到 0.02 m/s²，再归一化到 9.81 m/s²。它也不是采用独立二维 S² 误差坐标的实现。
\f[
u=g/\|g\|,\quad d=(I-uu^T)\delta g,\quad
g^+=9.81\frac{g+\operatorname{clip}_{0.02}(d)}{\|g+\operatorname{clip}_{0.02}(d)\|}.
\f]
**代码对应：** @ref lightning::NavState "NavState" 的 `boxplus`、`boxminus`、`ApplyGravityDelta`；状态表定义于 @ref nav_state.h "nav_state.h"。

# 2. 名义状态预测

给定去偏角速度 \f$\omega=\omega_m-b_g\f$ 和比力 \f$a=a_m-b_a\f$，代码的一阶传播为
\f[
p_{k+1}=p_k+v_k\Delta t,\quad R_{k+1}=R_k\operatorname{Exp}(\omega\Delta t),\quad
v_{k+1}=v_k+(R_ka+g_k)\Delta t.
\f]
偏置和重力的确定性导数为零。`propagate_velocity_=false` 时只将名义速度导数置零；协方差动力学并未一起删除速度耦合。这是实现中的兼容开关，不能把关闭分支解释成完整一致的常速随机模型。

`ImuProcess` 使用相邻 IMU 输入均值，但 `get_f/df_dx` 仍在区间起点状态求值，因此不是严格的中点状态积分；尤其名义位置没有直接加入半加速度项。去畸变的区间内点位置另外使用二阶平移，见 @ref imu_deskew "IMU 初始化、时间积分与点云去畸变"。

`PredictTo` 只接受有限且严格大于当前 `timestamp_` 的目标时间；积分后才对齐舍入误差。不能把负时间差取绝对值推进状态。

# 3. 误差与协方差预测

右扰动下，连续线性化的非零耦合包括
\f[
\delta\dot p=\delta v,\quad
\delta\dot\theta\simeq-[\omega]_\times\delta\theta-\delta b_g-n_g,
\quad \delta\dot v=-R[a]_\times\delta\theta-R\delta b_a+\delta g-Rn_a.
\f]
偏置随机游走为 \f$\delta\dot b_g=n_{bg},\delta\dot b_a=n_{ba}\f$。噪声排列与 `df_dw()` 的 12 列一致：gyro、acc、gyro bias、acc bias。

实现没有把 \f$-[\omega]_\times\f$ 直接写进 `df_dx`：`Predict` 单独构造旋转块
\f[
\Phi_{\theta\theta}=\operatorname{Exp}(-\omega\Delta t),\quad
\Phi_{\theta b_g}=-J_r(\omega\Delta t)\Delta t.
\f]
向量块用一阶离散，旋转噪声也乘同一个右雅可比。最终
\f[
P^- = \alpha\left(\Phi P\Phi^T+(\Delta t G)Q(\Delta t G)^T\right).
\f]
这里是代码实际的 \f$\Delta t^2Q\f$ 形式，不能直接把连续时间噪声谱密度的 \f$GQ_cG^T\Delta t\f$ 当作同一参数含义；改变 IMU 频率或噪声配置需做离散尺度核验。`predict_cov_inflation_` 默认为 1，按采样次固定膨胀会随采样频率累计。

**代码对应：** @ref lightning::NavState::df_dx "df_dx"、@ref lightning::NavState::df_dw "df_dw"、@ref lightning::ESKF::Predict "Predict"、@ref lightning::ESKF::PredictTo "PredictTo"。

# 4. 迭代更新为什么保留同一个先验

记预测状态为 \f$\hat x\f$，第 i 个线性化点为 \f$x_i\f$。`Update` 开始保存 `start_x` 和 `P_propagated`，每轮从该协方差重新开始，避免把同一帧点云当成多次独立观测重复融合。

计算 \f$d_i=x_i\boxminus\hat x\f$，用旋转块 \f$J_r(d_{i,\theta})\f$ 将先验差和协方差映射到当前切空间。回调返回六维位姿信息 \f$A=H^TWH\f$ 和 \f$b=-H^TWe\f$，误差符号见 @ref lidar_residuals "LiDAR 残差、信息矩阵与配准参数化"。

先对称化 A，再分解 \f$A=U\Lambda U^T\f$。仅保留 \f$\lambda_j>\lambda_{max}\cdot\text{ratio}\f$ 的方向：
\f[
\Pi=U\operatorname{diag}(m_j)U^T,\quad A_e=\Pi A\Pi,\quad b_e=\Pi b.
\f]
这是观测空间退化投影；先验仍会耦合状态，并不等于把退化状态彻底冻结。混合米与弧度的特征值阈值也依赖参数尺度。

令 E 为选择前六维的矩阵，代码在没有显式求逆的情况下，用 LDLT 求解
\f[
M=\sigma^2P_i^{-1}+E^TA_eE,\quad
K_r=M^{-1}E^Tb_e,\quad K_H=M^{-1}E^TA_eE,
\quad \delta_i=K_r+(K_H-I)d_i^{\text{mapped}}.
\f]
`Update(obs, R)` 的 `R` 是这里的标量方差 \f$\sigma^2\f$，不是姿态矩阵。六维观测可经 P 的交叉协方差修正速度和偏置，不能由“只填 6×6 的 HTH”推断“仅更新位姿”。

# 5. 工程门控先于状态注入

- `lidar_update_pose_only_` 保留预测的速度、偏置和重力；`lidar_update_inertial_states_` 单独控制偏置与重力通道。
- LiDAR 惯性修正过大可降级到保留惯性状态；速度步过大可进一步降级为仅位姿。
- 平移、旋转或最终速度步超过限值则恢复 `start_x/P_propagated`。观测无效分支恢复的是 `last_x/P_propagated`；这两个返回路径并不相同。
- 循环索引从 -1 开始，到 `< maximum_iter_`，最多调用 `maximum_iter_+1` 次观测函数。达到迭代上限是退出条件，不是收敛证明。

普通注入为 `x_=x_.boxplus(dx_current)`。Anderson 分支把总增量放回 `start_x` 的局部坐标，解历史残差差分最小二乘：
\f[
F_k=G(u_k)-u_k,\quad \gamma=\arg\min_\gamma\|F_k-\Delta F\gamma\|^2,
\quad u_{k+1}=G(u_k)-\Delta G\gamma.
\f]
列归一化和历史管理见 @ref anderson_acceleration.h "anderson_acceleration.h"。LiDAR 残差统计是**平方残差的上中位数**，并非字段名暗示的均值。加速后残差增大 1% 的分支退出循环，但 `last_x` 已存储加速结果；不要宣称它会重算一个保证下降的非加速解。

# 6. 协方差重置与数值保护

末次注入的旋转右雅可比形成 J。代码先做切空间变换，再执行等价于
\f$P^+=J(P_i-K_HP_i)J^T\f$ 的块矩阵运算；若观测退化，额外膨胀位姿块。
统一保护 `SymmetrizeAndFloorCovariance` 先对称化、修复非有限项和对角下限；若 LDLT 不能确认正定，则特征分解截断特征值，失败时重置为对角下限矩阵。
这些步骤是数值修复，不是几何正确性保证。`last_update_accepted_` 才是调用者判断观测是否提交的接口。

# 7. 轮速是另一条更新路径

`UpdateBodyForwardSpeed` 使用车体前向单位轴 \f$f_b\f$ 在 IMU 系的表达：
\f[
h(x)=(Rf_b)^Tv-v_{offset},\quad r=z-h(x),\quad
H_v=(Rf_b)^T,\quad S=HPH^T+\sigma_v^2.
\f]
这里只对速度求导，并显式清零增益的其他状态行，是条件速度观测，不是对姿态、偏置的完整联合线性化。创新绝对值、\f$r^2/S\f$ 和速度步长分别门控，接受后使用 Joseph 形式
\f$P^+=(I-KH)P(I-KH)^T+K\sigma_v^2K^T\f$。
参考点速度偏移补偿与后轴输出需结合 @ref output_contracts "定位状态、车体参考点与 ROS 输出" 阅读。

# 调试路线

先查 IMU 单位/时间与初始化，再查有效平面数量和退化秩，再查更新拒绝/降级原因，最后看协方差。`LIGHTNING_LM_DEBUG_ESKF_COV=1` 会检查特征值与相关系数，开销不可忽略。
修改预测、残差或注入时，优先运行 `src/test` 中对应 ESKF/IMU 测试，再做固定输入的轨迹回归；公式一致不等于数值实现已验证。
