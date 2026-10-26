@page geometry 坐标系、右扰动与李群雅可比

# 从哪个问题开始

同一个点在 LiDAR、IMU、里程计和地图坐标系中数值不同。先确定变换方向，再讨论残差和雅可比。本文只保留当前实现使用的数学工具；前端入口见 @ref laser_mapping_module "LIO 前端：LaserMapping"，滤波推导见 @ref eskf_theory "18 维迭代 ESKF：从预测到观测注入"。

# 变换方向与单位

约定 \f$T_{ab}=(R_{ab},t_{ab})\f$ 将 b 系点变到 a 系：
\f[
p_a=R_{ab}p_b+t_{ab},\qquad T_{ac}=T_{ab}T_{bc},\qquad
T_{ba}=(R_{ab}^T,-R_{ab}^Tt_{ab}).
\f]
位置用米，状态与 IMU 时间用秒，角速度用 rad/s；本文旋转向量的模长用 rad。配置角度单位以对应加载函数为准。内部点的 `time` 用毫秒，在去畸变处除以 1000，不能和 ROS header 秒时间混用。

| 量 | 当前实现含义 | 查阅位置 |
|---|---|---|
| `NavState::rot_`, `pos_` | IMU 在世界/局部里程计系的姿态和位置 | @ref lightning::NavState "NavState" |
| `offset_R_lidar_fixed_`, `offset_t_lidar_fixed_` | 主 LiDAR 到 IMU 的固定外参 | @ref lightning::LaserMapping "LaserMapping" |
| `R_lidar_to_primary`, `t_lidar_to_primary` | 副雷达到主雷达 | @ref lightning::MultiLidarSensorConfig "MultiLidarSensorConfig" |
| `T_map_odom` | 将 LIO 的 odom 位姿映射到优化地图系 | @ref lightning::backend::BackendPipeline "BackendPipeline" |

固定外参不在当前 18 维滤波状态中。地图包导出可能再施加重力对齐、原点归一化；导出后的世界系不应与运行时 odom 默认为同一坐标系，见 @ref map_contracts "地图包、坐标归一化与栅格导出"。

# 旋转指数映射

令 \f$[a]_\times b=a\times b\f$，\f$\theta=\|\phi\|\f$：
\f[
\operatorname{Exp}(\phi)=I+\frac{\sin\theta}{\theta}[\phi]_\times+
\frac{1-\cos\theta}{\theta^2}[\phi]_\times^2.
\f]
小角度时用连续极限避免除零。项目同时使用 Sophus 的 `SO3::exp(phi)` 和自定义 `math::exp(v, scale)`；后者构造半角四元数，旋转向量等于 `2*scale*v`，所以积分角速度要传 `0.5*dt`，不能直接照抄 Sophus 的参数。

@anchor geometry_half_angle_reading
## 为什么这里积分角速度要传半个 dt？

先看自定义 `math::exp` 如何构造四元数，注意它的默认 `scale=1`：

@snippet{lineno} core/lightning_math.hpp geometry-half-angle-exp

`cos_sinc_sqrt(scale*scale*norm2)` 返回半角的 cos 与 sinc；`Quatd` 第一个参数是实部，后三个参数是向量部。对有限输入，得到的是旋转向量 `2*scale*vec` 对应的旋转。调用处因此明确使用半个时间间隔：

@snippet{lineno} core/lio/imu_processing.hpp imu-point-rotation

<b>手算：</b>假设区间起点姿态为单位旋转，角速度为 `(0,0,π/2)` rad/s，dt=1 s。传入 `0.5*dt` 得到半角 π/4，旋转为 90°，把 `(1,0,0)` 转为 `(0,1,0)`。若误传 dt，半角变成 π/2，旋转就成为 180°，结果为 `(-1,0,0)`。这是区分接口约定的算例，不是去畸变允许一秒区间的声明；实际点时间仍受 @ref imu_time_reading "扫描时间边界" 限制。

<b>边界与停止点：</b>零向量通过 sinc 的小角展开得到单位旋转；不能把这里的默认参数直接移植到 `SO3::exp(phi)`。读到四元数构造与 `0.5*dt` 即可，完整实现见 @ref lightning::math::exp "math::exp()"、@ref lightning::math::cos_sinc_sqrt "cos_sinc_sqrt()"；点的最终补偿继续 @ref imu_point_reading "点云去畸变参考系"。

# 右扰动与左右雅可比

`NavState::boxplus` 和 `miao::VertexSE3::OplusImpl` 都采用
\f$R^+=R\operatorname{Exp}(\delta\theta)\f$，平移则在世界系直接相加。于是
\f[
R\operatorname{Exp}(\delta\theta)p \simeq Rp-R[p]_\times\delta\theta.
\f]
旋转残差的局部坐标是 IMU/body 切空间；平移增量是世界系。这不是完整左乘 SE(3) 更新。

\f[
J_l(\phi)=I+\frac{1-\cos\theta}{\theta^2}[\phi]_\times+
\frac{\theta-\sin\theta}{\theta^3}[\phi]_\times^2,
\quad J_r(\phi)=J_l(-\phi)=J_l(\phi)^T.
\f]
BCH 的一阶含义是：旋转向量相加与旋转相乘不是同一操作；例如
\f$\operatorname{Exp}(\phi+\delta)\simeq\operatorname{Exp}(\phi)\operatorname{Exp}(J_r(\phi)\delta)\f$。
代码 `A_matrix(phi)` 实现 \f$J_l\f$，其转置用于迭代滤波的切空间转换。

@anchor geometry_jacobian_reading
## A_matrix 的转置为什么换了方向，小角度又实际算了什么？

@snippet{lineno} core/lightning_math.hpp geometry-left-jacobian

先对照两个系数：`(1-cos(norm))/squaredNorm` 乘一次 hat，`(1-sin(norm)/norm)/squaredNorm` 乘两次 hat，与上面的左雅可比一致。hat 矩阵转置变号、hat 的平方转置不变，所以取转置等于把旋转向量取负；ESKF 的使用位置见 @ref eskf_iteration_reading "迭代切空间转换"。

<b>手算：</b>令 `v=(0,0,π/2)`，xy 平面的左雅可比为 `[[2/π,-2/π],[2/π,2/π]]`，转置后的右雅可比为 `[[2/π,2/π],[-2/π,2/π]]`，z 分量均为 1。两者作用在 x 方向的小增量时，y 分量符号相反，不能任意互换。

<b>边界与停止点：</b>代码在模长严格小于 `1e-5` rad 时直接返回单位阵，并未保留 `0.5*hat(v)` 项。例如 `v=(0,0,1e-6)` 时实现返回 I，而解析左雅可比的 xy 非对角一阶项约为 ±`5e-7`；这是该分支的近似。模长恰为阈值时走解析分支。确认一次 hat 的符号和小角分支后即可停止；完整函数见 @ref lightning::math::A_matrix "A_matrix()"，状态右乘与世界平移的注入例子复用 @ref eskf_state_reading "18 维状态与增量"。

# 不同求解器不能混用的约定

| 实现 | 增量顺序 | 更新 |
|---|---|---|
| ESKF 位姿块 / `VertexSE3` | 平移、旋转 | 世界平移相加，旋转右乘 |
| `VoxelBundleAdjuster` | 旋转、平移 | 旋转右乘，世界平移相加 |
| `PointToPlaneRegistration` | 平移、旋转 | `SE3::exp(increment) * pose`，完整左乘 |

因此平面配准中 `transformed.cross(normal)` 与 LIO 中 `point_imu.cross(R.transpose()*normal)` 都有各自成立的扰动前提。调试时先核对参数化，不要只比较矩阵的正负号。

**实现入口：** @ref nav_state.h "nav_state.h"、@ref vertex_se3.h "vertex_se3.h"、@ref voxel_bundle_adjustment.cc "voxel_bundle_adjustment.cc"、@ref point_to_plane_registration.cc "point_to_plane_registration.cc"。
