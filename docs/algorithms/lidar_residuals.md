@page lidar_residuals LiDAR 残差、信息矩阵与配准参数化

先从一个问题进入，再回到公式核对符号。下面的例子都是在明确假设下的手工推演，不代表已运行配准测试。

| 要确认的结论 | 短验证路径 |
|---|---|
| 有很多最近邻，为什么整轮观测仍被拒绝？ | @ref lidar_plane_gate_reading "点面数量的前置门槛" |
| 点面误差怎样变成六维信息，负号在哪里？ | @ref lidar_plane_equation_reading "雅可比、残差与累加" |
| 日志中的 residual mean 到底统计了什么？ | @ref lidar_residual_statistic_reading "平方残差的上中位数" |
| 定位配准的雅可比为何和 LIO 不同？ | @ref lidar_left_update_reading "左乘更新与负梯度" |

# 前端点面残差

@ref lightning::LaserMapping::ObsModel "ObsModel" 将 LiDAR 点变到地图，在 IVox 邻域拟合平面。设 IMU 系点 \f$q=R_{il}p_l+t_{il}\f$，世界平面为 \f$n^Tp+d=0\f$：
\f[
e=n^T(Rq+t)+d,\quad
J=\begin{bmatrix}n^T & (q\times R^Tn)^T\end{bmatrix}.
\f]
旋转雅可比来自 @ref geometry "坐标系、右扰动与李群雅可比" 的右扰动。代码 `A=point_crossmat*C`，C 为 `R.transpose()*normal`；观测残差取 `res=-e`。
\f[
HTH=\sum_j w_j\lambda_{plane}J_j^TJ_j,\qquad
HTr=-\sum_j w_j\lambda_{plane}J_j^Te_j.
\f]
`PointInformationScale` 给出当前信息权重，不能沿用旧笔记把此处写成固定 Cauchy 核；有效平面少于 20 个时，即使有点点匹配也先判本轮观测无效。六维压缩保留正规方程，不保留完整逐点雅可比。

@anchor lidar_plane_gate_reading
## 最近邻很多，为什么还不能构造观测？

在 `ObsModel()` 中，最近邻数量足够只会让点获得后续候选资格。点面项还要经过平面拟合和距离筛选，最终由 `effect_feat_surf_` 计数。构造信息矩阵之前先执行：

@snippet{lineno} core/lio/laser_mapping.cc lio-plane-count-gate

<b>关键语句：</b>`effect_feat_surf_ < 20` 直接设置 `obs.valid_ = false` 并返回，没有用 `effect_feat_icp_` 补足计数。后面的 `enable_icp_part_` 分支还没有执行。

<b>例子与边界：</b>假设只有 19 个有效点面匹配、却有 200 个点点候选，本轮依然无效；20 个点面匹配只能通过这个数量门槛，不能证明六自由度可观，也不能证明 ESKF 最终接受更新。

<b>读到这里：</b>能指出 `return` 位于点点累加之前即可。若需确认调用者如何回退，继续 @ref eskf_iteration_reading "ESKF 迭代与先验"；完整上下文见 @ref lightning::LaserMapping::ObsModel "ObsModel API"、@ref laser_mapping.cc "实现文件"。

@anchor lidar_plane_equation_reading
## 从一个点，怎样读出雅可比和残差的符号？

先只跟踪三个量：`point_this` 是 IMU 系的点 q；`norm_vec` 是世界系的平面法向 n；`Rt` 是当前 IMU 姿态的转置。`corr_pts_[i][3]` 在前面的匹配压缩阶段保存有符号点面距离 e。

@snippet{lineno} core/lio/laser_mapping.cc lio-plane-jacobian

`C = Rt * norm_vec` 把法向转到 IMU 系；`A = point_crossmat * C` 就是 \f$q\times R^Tn\f$。`J << ...` 把平移放在前三维、旋转放在后三维。最后的 `res = -...` 是正规方程右端的负误差，不是把平面的法向反转。

@snippet{lineno} core/lio/laser_mapping.cc lio-plane-information

这里每个点贡献 \f$w\lambda J^TJ\f$ 与 \f$-w\lambda J^Te\f$。`w` 来自测距/角度噪声模型的 `PointInformationScale()`；关闭该模型时返回 1，`plane_icp_weight_` 仍参与加权。这两个局部累加量随后按固定块顺序加到 `obs.HTH_` 和 `obs.HTr_`。

<b>手算一次：</b>假设外参和 R 都为单位变换，q=(2,0,0) m，平面 z=0，当前 t=(0,0,0.1) m。则 n=(0,0,1)，e=0.1 m，`A`=(0,-2,0)，`J`=[0,0,1,0,-2,0]，`res`=-0.1 m。若只放开 z 平移、固定其他量且权重为 1，这一个约束要求 \f$\Delta z=-0.1\f$ m，方向恰好朝平面靠近。

<b>边界与停止点：</b>单个点不足以解六维位姿，上例只验证符号和一维方向。能从 `J`、`res`、两条累加语句写回上面的正规方程即可；完整求解继续 @ref eskf_theory "ESKF"。噪声权重细节在 @ref laser_mapping.cc "PointInformationScale 与 ObsModel 源文件"，不必在第一次阅读中展开全部邻域搜索。

# 点点项及代码差异

可选点点项使用 \f$e=Rp+ t-p_{nearest}\f$ 的三维位置误差、0.5 m 距离门控和 `icp_weight_`。当前旋转块实际写为
\f$J_{rot}=-(R_{wi}R_{il})[p_l]_\times\f$。
这与严格按 IMU 右扰动推导的 \f$-R_{wi}[R_{il}p_l+t_{il}]_\times\f$ 并非一般等价，尤其存在外参旋转或杆臂时。文档保留此实现差异供后续算法核验，不擅自修改公式对应的代码。

计算使用固定 64 点块并行累加，再按块序串行归约，减少线程调度造成的求和顺序变化。`lidar_residual_mean_` 实际是点面平方残差上中位数，供 ESKF 的 AA 判断使用；不是所有残差的平均损失。

@anchor lidar_residual_statistic_reading
## residual mean 是平均误差吗？

先回看上一段的 `observation_residual_squared_[i] = res * res`，再看统计位置：

@snippet{lineno} core/lio/laser_mapping.cc lio-plane-residual-statistic

`begin() + size()/2` 和 `nth_element` 选中排序后零起始索引为 N/2 的值；偶数个数时取两个中间值中的较大者，不求两者平均。数组只由点面项填充，后面的可选点点项也没有进入这个统计。

<b>手算一次：</b>假设通过筛选的 20 个点中，10 个距离的绝对值为 0.01 m，另外 10 个为 0.03 m。平方后，上中位数是 0.0009 m²，算术均值是 0.0005 m²。此时 `lidar_residual_mean_` 和 `lidar_residual_max_` 都为 0.0009 m²；变量名中的 mean 不能据此解释成平均距离或 RMSE。

<b>读到这里：</b>只需连起“平方写入 → N/2 位置 → nth_element → 赋值”四个位置。它用于迭代诊断与 AA 判定，不直接证明定位置信度或输出许可；完整上下文见 @ref lightning::LaserMapping::ObsModel "ObsModel"。

# 定位点面配准使用左扰动

@ref lightning::loc::PointToPlaneRegistration "PointToPlaneRegistration" 独立执行粗细两级迭代。此时世界系变换点为 q，更新是 \f$T^+=\operatorname{Exp}_{SE(3)}(\delta)T\f$：
\f[
e=n^Tq+d,\quad J=[n^T,(q\times n)^T],\quad
(\sum wJ^TJ)\delta=-\sum wJ^Te.
\f]
权重为 Huber IRLS：小残差为 1，大残差为 \f$\delta_H/|e|\f$。求解前要求匹配数、内点率和归一化 Hessian 的最小特征值/条件数通过；再检查相对初值的总修正、收敛以及最终 RMSE。达到迭代上限而未收敛返回失败。

@anchor lidar_left_update_reading
## 定位配准为何使用变换后的点叉乘法向？

这里的 `transformed = pose * source_point` 已经在目标坐标系中。先看 `Linearize()` 中的旋转块和权重，再看 `RunStage()` 中更新乘在哪一侧：

@snippet{lineno} core/localization/point_to_plane_registration.cc loc-plane-left-jacobian

@snippet{lineno} core/localization/point_to_plane_registration.cc loc-plane-left-update

这一实现把正的 e 累加到 `gradient`，到 `solver.solve(-gradient)` 才取负号；前端则在 `res=-e` 时已取负号。两处负号的位置不同，目标都是减小残差。`SE3::exp(increment) * pose` 明确是左乘，所以旋转块使用目标系点 `transformed.cross(normal)`。

<b>手算区别：</b>假设外参为单位变换、R 为单位旋转，源点为 (1,0,0) m，t=(10,0,0) m，n=(0,1,0)。前端右旋转扰动的旋转块 z 分量是 1，定位左扰动使用目标点 (11,0,0)，该分量是 11。左乘旋转同时转动已有平移，右乘姿态注入则保持世界系平移不变；两套增量坐标不能混用。

<b>另验权重：</b>若配置 `huber_delta=0.1` m、已通过对应筛选的误差为 0.2 m，则权重为 0.5。此例只验证 IRLS 分支；其他匹配数、可观性、修正量和收敛条件仍须成立。

<b>读到这里：</b>能对应“叉乘用了哪个坐标系的点、负号出现在哪里、更新在左还是右”即可。完整入口为 @ref lightning::loc::PointToPlaneRegistration "PointToPlaneRegistration API"、@ref point_to_plane_registration.cc "Linearize / RunStage / Refine 源码"；是否接受这次定位还需看 @ref loc_read_acceptance "地图匹配的接受门控"。

# NDT 的分布匹配含义

NDT 在目标体素内估计均值 \f$\mu\f$ 和协方差 \f$\Sigma\f$，以变换点相对该分布的马氏距离
\f$q=(Tp-\mu)^T\Sigma^{-1}(Tp-\mu)\f$ 形成高斯混合式分数。当前 `pclomp` 实现含异常点比例与线搜索等常数，不能用单纯的 \f$\sum q\f$ 代替它的完整目标函数。
定位是否接受还受 `LidarLoc` 置信度、修正量、时间与状态门控影响；“NDT 收敛”不是发布许可。

<b>实现证据：</b> @ref laser_mapping.cc "laser_mapping.cc"（`ObsModel`、`PointInformationScale`）、@ref point_to_plane_registration.cc "point_to_plane_registration.cc"（`Linearize`、`RunStage`、`Refine`）、@ref ndt_omp_impl.hpp "ndt_omp_impl.hpp"（`computeDerivatives`、`computeTransformation`）、@ref lidar_loc.cc "lidar_loc.cc"（`Align` 及定位状态路径）。
