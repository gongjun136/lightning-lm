@page eskf_theory 18 维迭代 ESKF：从预测到观测注入

# 阅读入口与职责

先读 @ref laser_mapping_module "LIO 前端：LaserMapping" 看一次扫描如何到达滤波器，再读 @ref geometry "坐标系、右扰动与李群雅可比" 确认右扰动。@ref lightning::ESKF "ESKF" 维护名义状态 `x_` 和误差协方差 `P_`；它不负责找最近邻，LiDAR 观测由 `LaserMapping::ObsModel` 回调提供，见 @ref lidar_residuals "LiDAR 残差、信息矩阵与配准参数化"。

# 先抓住公式与代码的对应关系

| 计算步骤 | 核心关系 | 实现入口与阅读重点 |
|---|---|---|
| 状态与增量 | 18 维误差，旋转右乘 | @ref lightning::NavState "NavState" 的 `get_f/boxplus/boxminus`；先核对分块索引 |
| IMU 预测 | 名义状态积分与协方差传播 | @ref lightning::ESKF::Predict "Predict"、`NavState::df_dx/df_dw`；核对时间步、右雅可比与噪声尺度 |
| 构造 LiDAR 观测 | 六维 `HTH/HTr` 正规方程 | @ref lightning::LaserMapping::ObsModel "ObsModel"；对应 @ref lidar_residuals "残差与雅可比" |
| 迭代更新 | 同一个先验，在当前切空间反复线性化 | @ref lightning::ESKF::Update "Update"；对照 `start_x/P_propagated`、退化投影、LDLT 与 `boxplus` |
| 轮速更新 | 前向速度观测与创新门控 | @ref lightning::ESKF::UpdateBodyForwardSpeed "UpdateBodyForwardSpeed"；核对观测参考点与 Joseph 协方差更新 |

首次阅读先核对三个问题：@ref eskf_state_reading "误差增量写进哪些状态" → @ref eskf_prediction_reading "IMU 如何推进名义状态" → @ref eskf_iteration_reading "同一帧为什么需要反复线性化"。每段片段直接来自当前源码；手算例子用于核对语句含义，不是运行测试结果。工程门控要知道何时拒绝或降级；Anderson 加速、协方差数值修复与轮速更新可按需继续读。

@anchor eskf_state_reading
# 1. 18 维增量分别写到哪里，姿态为什么右乘？

\f[
x=(p,R,v,b_g,b_a,g),\qquad
\delta x=(\delta p,\delta\theta,\delta v,\delta b_g,\delta b_a,\delta g)\in\mathbb R^{18}.
\f]
索引依次是 0、3、6、9、12、15，每块三维。`NavState::dim` 和 `full_dim` 都是 18；外参固定。`prediction_gyro_`、停车状态和置信度是元数据，不增加滤波维度。

先看索引定义，再看真正写入状态的 `boxplus()`：

@snippet{lineno} common/nav_state.h eskf-state-layout

@snippet{lineno} common/nav_state.h eskf-state-injection

只追三类语句：`middleRows(index, 3)` 取出对应的三维块；位置、速度与偏置直接相加；旋转通过 `rot_ * SO3::exp(...)` 在右侧注入。`pos_`、`vel_` 表达在世界系，旋转增量表达在当前 IMU 局部轴，单位 rad。`ret = *this` 还说明此函数返回修正后的副本，调用者需要接收返回值；本次观测注入没有推进时间戳。

<b>手算例子：</b>若当前位置为 `(10, 0, 0)` m，仅令位置块增量为 `(1, 0, 0)` m，结果位置就是 `(11, 0, 0)`。若当前姿态为绕世界 z 轴转 90°，再令旋转块为 `(0.01, 0, 0)` rad，增量绕当前 IMU x 轴；用 `R * UnitX = (0, 1, 0)` 就能确认该轴在世界系中指向 y。它不能被解释为绕世界 x 轴左乘修正。

重力虽然存三维，却不是无约束三自由度：`ApplyGravityDelta` 去掉径向分量，将单次切向修正裁到 0.02 m/s²，再归一化到 9.81 m/s²。它也不是采用独立二维 S² 误差坐标的实现。
\f[
u=g/\|g\|,\quad d=(I-uu^T)\delta g,\quad
g^+=9.81\frac{g+\operatorname{clip}_{0.02}(d)}{\|g+\operatorname{clip}_{0.02}(d)\|}.
\f]
例如重力为 `(0, 0, -9.81)`，仅沿 z 方向的增量会被径向投影去掉；不能照着普通向量加法把重力模长改成 8.81。只有需要检查重力约束时再打开 @ref lightning::NavState::ApplyGravityDelta() "ApplyGravityDelta()"。

<b>读到这里：</b>能把六个块、右乘方向和重力特例对上语句即可。误差如何取回见 @ref lightning::NavState::boxminus() "boxminus()"；完整状态定义见 @ref nav_state.h "nav_state.h"。

@anchor eskf_prediction_reading
# 2. 给一段 IMU，位置、姿态与速度怎样推进？

给定去偏角速度 \f$\omega=\omega_m-b_g\f$ 和比力 \f$a=a_m-b_a\f$，代码的一阶传播为
\f[
p_{k+1}=p_k+v_k\Delta t,\quad R_{k+1}=R_k\operatorname{Exp}(\omega\Delta t),\quad
v_{k+1}=v_k+(R_ka+g_k)\Delta t.
\f]
在 @ref lightning::ESKF::Predict() "Predict()" 中先求 `get_f()`，再把得到的导数传给 `oplus(f_, dt)`。核对公式时先读这两个短函数：

@snippet{lineno} common/nav_state.h eskf-state-derivative

@snippet{lineno} common/nav_state.h eskf-state-propagation

`get_f()` 从全零向量开始，只写位置、姿态和速度三个块，所以偏置与重力的确定性导数为零。`res[kPosIdx] = vel_` 使用区间起点速度；随后 `oplus()` 对该导数乘 dt，位置更新里没有另加 `0.5*a*dt*dt`。`rot_ * (acce-ba_) + grav_` 把 IMU 比力旋到世界系再加世界重力，得到速度导数。

<b>手算例子：</b>假设允许传播速度，姿态为单位旋转、偏置为零、陀螺输入为零，世界重力为 `(0, 0, -9.81)` m/s²，加计输入为 `(1, 0, 9.81)` m/s²。此时净加速度为 `(1, 0, 0)` m/s²。若初始位置为零、速度为 `(2, 0, 0)` m/s、dt=0.1 s，本步得到位置 `(0.2, 0, 0)` m、速度 `(2.1, 0, 0)` m/s，而不是把本步位置写成 0.205 m。后者对应另一种位置积分公式。

还要看速度传播开关的位置：

@snippet{lineno} core/lio/eskf.cc eskf-predict-velocity-switch

`propagate_velocity_=false` 时只将名义速度导数置零；`df_dx(acce)` 仍照常调用，协方差动力学并未一起删除速度耦合。这是实现中的兼容开关，不能把关闭分支解释成完整一致的常速随机模型。

`ImuProcess` 使用相邻 IMU 输入均值，但 `get_f/df_dx` 仍在区间起点状态求值，因此不是严格的中点状态积分；尤其名义位置没有直接加入半加速度项。去畸变的区间内点位置另外使用二阶平移，见 @ref imu_deskew "IMU 初始化、时间积分与点云去畸变"。

<b>读到这里：</b>能用 `get_f → oplus` 算出一步状态即可。`PredictTo` 对有限、严格递增目标时间的检查，以及跨扫描时间如何裁剪，见 @ref imu_time_reading "扫描末时间边界"；该检查在包装函数中，直接调用 `Predict(dt, ...)` 的调用者仍须保证 dt 合法。

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

@snippet{lineno} core/lio/eskf.cc eskf-predict-covariance

核对噪声尺度时，只看 `Q` 两边各有一个 `dt * f_w_final`。在固定 G=1、Q=4、膨胀系数为 1 的一维算例中，dt=0.01 的噪声增量为 0.0004，dt=0.02 时为 0.0016。这里仅比较噪声项，不包括 `F_x1_ * P_ * F_x1_.transpose()`；不能据此断言整个协方差也放大四倍。

<b>代码对应：</b>@ref lightning::NavState::df_dx "df_dx"、@ref lightning::NavState::df_dw "df_dw"、@ref lightning::ESKF::Predict "Predict"。

@anchor eskf_iteration_reading
# 4. 同一帧反复线性化，为什么不等于融合多次？

记预测状态为 \f$\hat x\f$，第 i 个线性化点为 \f$x_i\f$。`Update` 开始保存 `start_x` 和 `P_propagated`，每轮从该协方差重新开始，避免把同一帧点云当成多次独立观测重复融合。

先读循环外保存了什么；中间的 `K_r`、`K_H` 和 `dx_current` 是求解用的工作量，本段先略过：

@snippet{lineno} core/lio/eskf.cc eskf-update-save-prior

再看循环内传给观测函数的是当前 `x_`，随后计算相对于固定 `start_x` 的差，并恢复预测协方差：

@snippet{lineno} core/lio/eskf.cc eskf-update-relinearize

@snippet{lineno} core/lio/eskf.cc eskf-update-reuse-prior

这里变化的是线性化点；`start_x` 与 `P_propagated` 没有随迭代覆盖。`P_ = P_propagated` 也不表示丢弃 IMU 预测，它恢复的正是进入本次 LiDAR 更新时的预测结果。随后才做下述切空间映射和信息方程求解。

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

<b>手算同一线性观测的两轮：</b>仅取一维位置，假设先验位置为 0、P=1，观测 z=0.4、H=1、R=1；关闭加速，门控通过，并假设继续进入第二轮。此时两轮都有 `K_H=0.5`：

| 轮次 | 当前位置 | 当前残差 z−x | 相对固定先验的 dx | 本轮增量 `K_r + (K_H−1)dx` |
|---|---|---|---|---|
| 第一轮 | 0 | 0.4 | 0 | 0.2 + (−0.5)×0 = 0.2 |
| 第二轮 | 0.2 | 0.2 | 0.2 | 0.1 + (−0.5)×0.2 = 0 |

第二轮仍回到同一个后验位置 0.2，没有把上一轮的 0.2 再加一遍。真实点云会重新关联和线性化，矩阵不必保持不变；此例只验证“固定先验差项怎样抵消重复修正”。旋转流形、退化投影和修正通道开关不能由这一维算例省略。

<b>停止位置：</b>能区分“保存先验、更新线性化点、恢复先验协方差”后，就已确认迭代的目的。下面只需再认清注入前的拒绝/降级分支；求解器完整实现见 @ref lightning::ESKF::Update() "Update()"，观测矩阵来源见 @ref lidar_residuals "LiDAR 残差导读"。

# 5. 工程门控先于状态注入

- `lidar_update_pose_only_` 保留预测的速度、偏置和重力；`lidar_update_inertial_states_` 单独控制偏置与重力通道。
- LiDAR 惯性修正过大可降级到保留惯性状态；速度步过大可进一步降级为仅位姿。
- 平移、旋转或最终速度步超过限值则恢复 `start_x/P_propagated`。观测无效分支恢复的是 `last_x/P_propagated`；这两个返回路径并不相同。
- 循环索引从 -1 开始，到 `< maximum_iter_`，最多调用 `maximum_iter_+1` 次观测函数。达到迭代上限是退出条件，不是收敛证明。

只有经过前面的增量门控，才到达普通注入分支：

@snippet{lineno} core/lio/eskf.cc eskf-update-plain-injection

这行调用第一节的 `boxplus()`，并用返回值替换当前状态。它仅展示关闭 Anderson 加速时的路径；本轮临时注入也不等同于整次观测已被接受，应核对 @ref lightning::ESKF::LastUpdateAccepted() "LastUpdateAccepted()"。Anderson 分支把总增量放回 `start_x` 的局部坐标，解历史残差差分最小二乘：
\f[
F_k=G(u_k)-u_k,\quad \gamma=\arg\min_\gamma\|F_k-\Delta F\gamma\|^2,
\quad u_{k+1}=G(u_k)-\Delta G\gamma.
\f]
列归一化和历史管理见 @ref anderson_acceleration.h "anderson_acceleration.h"。LiDAR 残差统计是<b>平方残差的上中位数</b>，并非字段名暗示的均值。加速后残差增大 1% 的分支退出循环，但 `last_x` 已存储加速结果；不要宣称它会重算一个保证下降的非加速解。

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
