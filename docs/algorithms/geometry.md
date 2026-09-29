@page geometry 坐标系、右扰动与李群雅可比

# 从哪个问题开始

同一个点在 LiDAR、IMU、里程计和地图坐标系中数值不同。先确定变换方向，再讨论残差和雅可比。本文只保留当前实现使用的数学工具；前端入口见 @ref laser_mapping_module "LIO 前端：LaserMapping"，滤波推导见 @ref eskf_theory "18 维迭代 ESKF：从预测到观测注入"。

# 变换方向与单位

约定 \f$T_{ab}=(R_{ab},t_{ab})\f$ 将 b 系点变到 a 系：
\f[
p_a=R_{ab}p_b+t_{ab},\qquad T_{ac}=T_{ab}T_{bc},\qquad
T_{ba}=(R_{ab}^T,-R_{ab}^Tt_{ab}).
\f]
位置用米，内部时间用秒，角速度用 rad/s；配置字段带 `_deg` 才按度解释。点内 `time` 在去畸变处除以 1000，不能和 ROS header 秒时间混用。

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

**源码：** @ref lightning_math.hpp "lightning_math.hpp"（`exp`、`A_matrix`），@ref imu_processing.hpp "imu_processing.hpp"（`UndistortPcl`）。

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

# 不同求解器不能混用的约定

| 实现 | 增量顺序 | 更新 |
|---|---|---|
| ESKF 位姿块 / `VertexSE3` | 平移、旋转 | 世界平移相加，旋转右乘 |
| `VoxelBundleAdjuster` | 旋转、平移 | 旋转右乘，世界平移相加 |
| `PointToPlaneRegistration` | 平移、旋转 | `SE3::exp(increment) * pose`，完整左乘 |

因此平面配准中 `transformed.cross(normal)` 与 LIO 中 `point_imu.cross(R.transpose()*normal)` 都有各自成立的扰动前提。调试时先核对参数化，不要只比较矩阵的正负号。

**实现入口：** @ref nav_state.h "nav_state.h"、@ref vertex_se3.h "vertex_se3.h"、@ref voxel_bundle_adjustment.cc "voxel_bundle_adjustment.cc"、@ref point_to_plane_registration.cc "point_to_plane_registration.cc"。
