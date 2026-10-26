@page output_contracts 定位状态、车体参考点与 ROS 输出

# 阅读范围与四个问题

本页的 LIO/地图匹配门控适用于 `lidar` 模式。`ins_only` 复用业务消息构造函数，但直接输出已由设备补偿到后轴的导航解，采用独立的严格质量门控；坐标、失效和点云输出契约见 @ref ins_only_operation "CGI-430 纯组合导航模式"。

沿下面四个问题核对输出边界即可：@ref output_callback_reading "结果怎样进入发布函数" → @ref output_gate_reading "为什么有效结果仍会停发" → @ref output_pose_reading "后轴位姿怎样计算" → @ref output_speed_reading "后轴速度怎样计算"。每段先读短片段与例子，需要上下文时再打开完整实现。

| 变量或判断 | 回答的问题 | 不足以证明什么 |
|---|---|---|
| `lidar_loc_valid_` | 这次地图匹配是否有效 | 融合结果已发布 |
| `LocalizationResult::valid_` | 这次融合/DR 结果是否有效 | 最近地图匹配仍有效且新鲜 |
| `PoseOutputsEnabled()` | 当前是否允许业务位姿输出 | 已调用发布器，或订阅端已收到 |

@anchor output_callback_reading
# 1. 算法结果怎样进入发布函数？

在 @ref lightning::LocSystem::Init() "LocSystem::Init()" 中，先找到这个回调连接。这里的 `result` 是实时结果；注册 lambda 时尚未发布任何位姿。

@htmlonly[block]
<details open>
<summary>真实源码：把实时结果交给输出边界</summary>
@endhtmlonly

@snippet{lineno} core/system/loc_system.cc loc-system-init-live-result

@htmlonly[block]
</details>
@endhtmlonly

<b>读法：</b> `SetLocalizationResultCallback` 保存处理函数；之后算法触发该回调，才执行 `PublishLocalizationResult(result)`。`Init()` 紧邻的 `SetGlobalLocalizationResultCallback` 则连接 `CaptureGlobalLocalizationResult()`，用于收集全局校正轨迹。两条回调的目的不同，不能把全局轨迹中的一个样本直接视为一次实时 ROS 发布。

<b>例子与停止位置：</b> 即使回调收到 `timestamp_=12.0` 的结果，也只能确认它进入了发布函数；要判断它能否发出去，继续读下一节的 `return` 条件。若只想确认“谁调用发布函数”，读完这个 lambda 就可以停。完整实现：@ref loc_system.cc "loc_system.cc"；上游调用者：@ref lightning::loc::Localization::Init() "Localization::Init()"。

@anchor output_gate_reading
# 2. 为什么 `valid_=true`，业务位姿仍可能停发？

先按顺序看 @ref lightning::LocSystem::PublishLocalizationResult() "PublishLocalizationResult()" 的三个检查：结果有效、时间戳递增、地图匹配允许输出。任何一处 `return` 都发生在业务消息构造之前。

@htmlonly[block]
<details open>
<summary>真实源码：发布入口的三个检查</summary>
@endhtmlonly

@snippet{lineno} core/system/loc_system.cc output-result-publication-gates

@htmlonly[block]
</details>
@endhtmlonly

`result.valid_` 只通过第一关；`pose_timestamp_gate_.Observe()` 检查有限、正数且严格递增的时间。假设 12.0 已被时间戳门控接受，随后收到另一个 12.0 或 11.9，即使结果有效，也会在 `timestamp_decision.accepted` 处被拒绝。注意这一步早于匹配门控：一个因匹配失效未发布的样本，仍可能已经推进时间戳水位。完整时间规则见 @ref lightning::MonotonicTimestampGate::Observe() "MonotonicTimestampGate::Observe()"。

然后只追一层 `PoseOutputsEnabled()`，直接看最后的布尔表达式：

@htmlonly[block]
<details open>
<summary>真实源码：位姿与地图输出共用的许可条件</summary>
@endhtmlonly

@snippet{lineno} core/system/sany_localization_output.cc output-publication-permission

@htmlonly[block]
</details>
@endhtmlonly

三个条件必须同时成立：曾有有效地图匹配、最近有效匹配未过期、连续失败帧数为零。`MapOutputsEnabled()` 复用同一个判断。这里写的是 `consecutive_lost_frames_ == 0`，所以<b>第一次被门控接受的无效匹配就停止受控位姿和地图输出，无需等到 `lost_frame_threshold`</b>。

匹配状态从哪里来，可以用下面两个短片段接上：`PublishProcessedCloud()` 用 `lidar_loc_valid_` 更新门控；`ObserveLidarMatch()` 的有效分支清零失败计数，无效分支累加。

@htmlonly[block]
<details>
<summary>真实源码：匹配回调如何更新发布许可</summary>
@endhtmlonly

@snippet{lineno} core/system/loc_system.cc output-observe-lidar-match

@snippet{lineno} core/system/sany_localization_output.cc output-match-state-update

@htmlonly[block]
</details>
@endhtmlonly

例如依次收到时间递增的“有效 → 无效 → 有效”匹配，失败计数相应为 `0 → 1 → 0`；若匹配年龄也满足条件，许可相应为“允许 → 停止 → 恢复”。此例假设已有有效匹配，且输入通过 `ObserveLidarMatch()` 开头的时间检查；重复或回退的正时间戳不会更新该状态。

<b>没有新匹配时看什么？</b> 继续看 `MatchStaleLocked()`，同时核对传感器时间差与 steady-clock 到达时间差。后者在传感器时间完全停住时仍会推进。

@htmlonly[block]
<details>
<summary>真实源码：两种时间的过期判断与失效锁存</summary>
@endhtmlonly

@snippet{lineno} core/system/sany_localization_output.cc output-match-freshness

@htmlonly[block]
</details>
@endhtmlonly

设最大年龄为 0.5 s，最近有效匹配时间为 20.0 s，并假设 steady-clock 尚未超时：检查到 20.5 s 时仍满足年龄限制，20.51 s 时超限。若所有传感器都停在 20.0 s，但距有效匹配到达已过 0.501 s，也会失效。`stale_latched_` 保留失效状态；一个更晚的有效匹配先清除锁存，再由上述年龄判断决定是否真正恢复。若已知最新传感器时间为 23.0 s，刚算完的有效匹配却只有 21.4 s，仍会立即重新超时。

读到这里已经可以解释“结果有效，但输出被拒绝”。排查连续失锁时再去看 @ref localization_module "定位模块"；确认订阅端收到消息还需要运行时证据。100 ms 健康定时器调用 @ref lightning::LocSystem::PublishHealthStatus() "PublishHealthStatus()"，与位姿回调分开报告状态，位姿停发不代表诊断也停发。完整实现：@ref sany_localization_output.cc "sany_localization_output.cc"；现有边界用例：@ref sany_localization_output_test.cc "sany_localization_output_test.cc"。

@anchor output_pose_reading
# 3. 主 LiDAR 的位姿如何变成业务地图中的后轴位姿？

先看输出函数中的两个连续步骤。第一个改变车体参考点，第二个改变地图坐标系：

@htmlonly[block]
<details open>
<summary>真实源码：先转后轴，再转业务地图</summary>
@endhtmlonly

@snippet{lineno} core/system/loc_system.cc output-rear-pose-boundary

@htmlonly[block]
</details>
@endhtmlonly

约定 `T_a_b` 把 b 系坐标变到 a 系。输入 `result.pose_` 为 `T_localization_lidar`；`primary_lidar_position_in_body_` 是主 LiDAR 在后轴车体系中的安装位置，单位 m。下面的 `.inverse()` 是阅读重点：`T_localization_rear = T_localization_lidar * inverse(T_rear_lidar)`。

@htmlonly[block]
<details open>
<summary>真实源码：参考点转换中的右乘逆变换</summary>
@endhtmlonly

@snippet{lineno} core/system/sany_localization_output.cc output-rear-pose-transform

@htmlonly[block]
</details>
@endhtmlonly

`MakeRearAxleLidarTransform()` 同时放入安装位置和 `initial_lidar_rotation`；后者经 `MakeLivoxLidarTransform()` 构造成零平移的旋转变换。固定地图对齐则使用左乘：

@htmlonly[block]
<details>
<summary>真实源码：将位姿左乘到业务地图坐标系</summary>
@endhtmlonly

@snippet{lineno} core/system/sany_localization_output.cc output-fixed-map-transform

@htmlonly[block]
</details>
@endhtmlonly

<b>手算例子：</b> 假设安装旋转为单位旋转，LiDAR 位于后轴正前方 2 m，车辆在定位地图中的 yaw 为 90°，LiDAR 地图位置为 `(10, 5, 0)` m。杆臂旋转到地图后为 `(0, 2, 0)`，因此后轴位置为 `(10, 3, 0)`，不能直接从地图 x 减 2。若固定地图变换再沿 x 平移 100 m，业务位置变成 `(110, 3, 0)`。这些数字仅用于说明变换顺序，不是车辆标定值。

读到“右乘改变车体参考点，左乘改变地图坐标系”即可确认此处输出关系。两个函数使用传入位姿生成新值，`result` 以 const 引用进入发布函数；不要把输出对齐后的结果再当作内部估计器位姿反馈。完整实现：@ref lightning::sany_output::MakeMapRearAxlePose() "MakeMapRearAxlePose()"、@ref lightning::sany_output::TransformPoseForOutput() "TransformPoseForOutput()"。

@anchor output_speed_reading
# 4. 后轴速度为何要在旋转之后再减一项？

位姿参考点改变后，速度也要对应后轴。先在发布函数中找 `vehicle_speed`：

@htmlonly[block]
<details open>
<summary>真实源码：旋转、取纵向分量、减杆臂修正，然后构造消息</summary>
@endhtmlonly

@snippet{lineno} core/system/loc_system.cc output-rear-speed-message

@htmlonly[block]
</details>
@endhtmlonly

`result.vel_b_` 在这条输出路径中是 IMU 轴表达的速度。`ImuVelocityToBody()` 只旋转表达坐标系，速度参考点仍是 IMU。设车体系杆臂 r 从后轴指向 IMU，单位 m；角速度为车体系表达、单位 rad/s：
\f[v_{imu}=v_{rear}+\omega\times r,\qquad v_{rear,x}=v_{imu,x}-(\omega\times r)_x.\f]

@htmlonly[block]
<details>
<summary>真实源码：坐标旋转与杆臂速度的计算</summary>
@endhtmlonly

@snippet{lineno} common/imu_body_velocity.h output-imu-velocity-rotation

@snippet{lineno} common/imu_body_velocity.h output-rear-speed-offset

@htmlonly[block]
</details>
@endhtmlonly

看清 `cross(rear_to_imu_body).x()` 与发布处的减号。假设已启用杆臂修正、处于非停车状态、修正量与输出时间一致，旋转后的 IMU 速度为 `(2, 0, 0)` m/s，角速度为 `(0, 0, 1)` rad/s，杆臂为 `(0, 1, 0)` m，则叉乘的 x 分量为 −1 m/s，输出后轴纵向速度为 `2 - (-1) = 3` m/s。若只做旋转，就会遗漏这一差值。

还要核对修正量何时有效：发布代码调用的是结果对象的方法，不是无条件取 `rear_axle_speed_offset_`。

@htmlonly[block]
<details>
<summary>真实源码：停车状态与修正量时间的有效条件</summary>
@endhtmlonly

@snippet{lineno} core/localization/localization_result.h output-speed-offset-validity

@htmlonly[block]
</details>
@endhtmlonly

只有非停车且修正时间与结果时间相差不超过 `1e-6` s，才返回存储的修正量，否则返回零。这里不会自动插值旧修正量。读完这两个条件就可以停止；只有要检查修正量的生产过程，才继续进入 @ref lightning::LaserMapping::GetRearAxleSpeedOffset() "LaserMapping::GetRearAxleSpeedOffset()"。`MakePoseMessage(vehicle_pose)` 直接从同一业务消息派生位姿，因此无需再推导另一套坐标转换。

# 消息与轨迹参考

TUM 记录 `WriteTumPoseLine` 要求时间严格单调。production 与 diagnostic 的记录开关和采样策略不同；调用 `RecordPublishedPoseTum()` 之前还检查业务发布已经发生。空轨迹不能解释为没有内部定位结果，也不能用内部轨迹冒充实际发布的后轴轨迹。

| 实现 | 输入/输出 | 约束 |
|---|---|---|
| `MakeVehiclePoseMessage` | 后轴位姿、速度、时间 → VehiclePose | 位姿与速度参考点一致 |
| `MakePoseMessage` | VehiclePose → PoseStamped | 两种接口使用同一位姿 |
| `MakeCloudMessage` | 点云、起止时间、变换 → ROS 云 | 与对应地图 frame 一致 |
| `LocalizationTelemetryState` | 状态与位姿 → status/fault/path | 容量、采样周期和失效状态独立管理 |

<b>源码：</b> @ref sany_localization_output.h "sany_localization_output.h"、@ref sany_localization_output.cc "sany_localization_output.cc"、@ref imu_body_velocity.h "imu_body_velocity.h"、@ref rear_axle_pose.h "rear_axle_pose.h"、@ref loc_system.cc "loc_system.cc"。
相关原理 @ref eskf_theory "18 维迭代 ESKF：从预测到观测注入"、@ref pose_graph_theory "位姿图、增量求解与高频平滑"；现场排查 @ref localization_troubleshooting "定位失效排查与取证"。
