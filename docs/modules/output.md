@page output_contracts 定位状态、车体参考点与 ROS 输出

# 三层有效性

`LidarLoc` 的 `lidar_loc_valid_` 表示地图匹配状态；`LocalizationResult::valid_` 表示融合结果；`LocalizationPublicationGate` 最后决定 ROS 位姿和地图相关输出是否允许。内部结果有效不等于消息已发送，发送也不等于 UI 已接收。

@ref lightning::LocSystem "LocSystem" 持有发布器、telemetry、heartbeat 和轨迹文件，算法状态由 @ref lightning::loc::Localization "Localization" 管理。关闭顺序见 @ref online_localization_flow "在线定位端到端流程"。

# 位姿和速度不是只改 frame_id

从主 LiDAR 位姿转换到后轴需要旋转和杆臂，`MakeMapRearAxlePose` 负责几何转换。固定地图对齐再左乘 `T_target_localization`，仅用于输出边界；不应反馈到估计器内部再次变换。

速度先从 IMU 局部轴旋转到车体轴，仍位于 IMU 参考点。设车体系杆臂 r 从后轴指向 IMU：
\f[v_{imu}=v_{rear}+\omega\times r,\qquad v_{rear,x}=v_{imu,x}-(\omega\times r)_x.\f]
`RearAxleSpeedOffset` 计算前向偏移；轮速更新和发布必须使用一致的轴、角速度时间与参考点，不能把纯旋转坐标变换当作杆臂补偿。

# 发布门控与时间

`LocalizationPublicationGate` 记录最近有效匹配的传感器时间和 steady-clock 到达时间；检查连续失败帧及匹配年龄，失效后停止受控输出。ROS 时钟/传感器停顿与墙钟推进是不同故障，不能只看一种时间。状态/故障心跳继续承担报告不可用状态的职责。

TUM 记录 `WriteTumPoseLine` 要求时间严格单调。production 与 diagnostic 的记录开关不同，空轨迹不能解释为没有内部定位结果，也不能用内部轨迹冒充实际发布的后轴轨迹。

| 实现 | 输入/输出 | 约束 |
|---|---|---|
| `MakePosResMessage` | 后轴位姿、速度、时间 → PosRes | 位姿与速度参考点一致 |
| `MakePoseMessage/MakeVehiclePoseMessage` | PosRes → 派生消息 | 避免不同接口各算一套位姿 |
| `MakeCloudMessage` | 点云、起止时间、变换 → ROS 云 | 与对应地图 frame 一致 |
| `LocalizationTelemetryState` | 状态与位姿 → status/fault/path | 容量、采样周期和失效状态独立管理 |

**源码：** @ref sany_localization_output.h "sany_localization_output.h"、@ref sany_localization_output.cc "sany_localization_output.cc"、@ref imu_body_velocity.h "imu_body_velocity.h"、@ref rear_axle_pose.h "rear_axle_pose.h"、@ref loc_system.cc "loc_system.cc"。
相关原理 @ref eskf_theory "18 维迭代 ESKF：从预测到观测注入"、@ref pose_graph_theory "位姿图、增量求解与高频平滑"；现场排查 @ref localization_troubleshooting "定位失效排查与取证"。
