@page lidar_residuals LiDAR 残差、信息矩阵与配准参数化

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

# 点点项及代码差异

可选点点项使用 \f$e=Rp+ t-p_{nearest}\f$ 的三维位置误差、0.5 m 距离门控和 `icp_weight_`。当前旋转块实际写为
\f$J_{rot}=-(R_{wi}R_{il})[p_l]_\times\f$。
这与严格按 IMU 右扰动推导的 \f$-R_{wi}[R_{il}p_l+t_{il}]_\times\f$ 并非一般等价，尤其存在外参旋转或杆臂时。文档保留此实现差异供后续算法核验，不擅自修改公式对应的代码。

计算使用固定 64 点块并行累加，再按块序串行归约，减少线程调度造成的求和顺序变化。`lidar_residual_mean_` 实际是点面平方残差上中位数，供 ESKF 的 AA 判断使用；不是所有残差的平均损失。

# 定位点面配准使用左扰动

@ref lightning::loc::PointToPlaneRegistration "PointToPlaneRegistration" 独立执行粗细两级迭代。此时世界系变换点为 q，更新是 \f$T^+=\operatorname{Exp}_{SE(3)}(\delta)T\f$：
\f[
e=n^Tq+d,\quad J=[n^T,(q\times n)^T],\quad
(\sum wJ^TJ)\delta=-\sum wJ^Te.
\f]
权重为 Huber IRLS：小残差为 1，大残差为 \f$\delta_H/|e|\f$。求解前要求匹配数、内点率和归一化 Hessian 的最小特征值/条件数通过；再检查相对初值的总修正、收敛以及最终 RMSE。达到迭代上限而未收敛返回失败。

# NDT 的分布匹配含义

NDT 在目标体素内估计均值 \f$\mu\f$ 和协方差 \f$\Sigma\f$，以变换点相对该分布的马氏距离
\f$q=(Tp-\mu)^T\Sigma^{-1}(Tp-\mu)\f$ 形成高斯混合式分数。当前 `pclomp` 实现含异常点比例与线搜索等常数，不能用单纯的 \f$\sum q\f$ 代替它的完整目标函数。
定位是否接受还受 `LidarLoc` 置信度、修正量、时间与状态门控影响；“NDT 收敛”不是发布许可。

**实现证据：** @ref laser_mapping.cc "laser_mapping.cc"（`ObsModel`、`PointInformationScale`）、@ref point_to_plane_registration.cc "point_to_plane_registration.cc"（`Linearize`、`RunStage`、`Refine`）、@ref ndt_omp_impl.hpp "ndt_omp_impl.hpp"（`computeDerivatives`、`computeTransformation`）、@ref lidar_loc.cc "lidar_loc.cc"（`Align` 及定位状态路径）。
