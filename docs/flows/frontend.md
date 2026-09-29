@page frontend_flow 纯前端：传感器到 LIO 状态

# 用途与入口

纯前端用于隔离时间同步、初始化、去畸变和滤波问题，不加载定位地图，也不运行完整建图后端。离线入口 @ref run_frontend_offline.cc "run_frontend_offline.cc"，在线入口 @ref run_frontend_online.cc "run_frontend_online.cc"；脚本参数见 @ref script_contracts "Shell 启动与实验脚本契约"。

@dot
digraph frontend_path { rankdir=LR; node [shape=box]; input [label="bag / ROS callbacks"]; sync [label="preprocess / multi-lidar / IMU sync"]; init [label="IMU init + deskew"]; update [label="ESKF + IVox"]; output [label="state / keyframe / trajectory"]; input -> sync -> init -> update -> output; }
@enddot

# 离线顺序

`main` 构造 `LaserMapping`、可选 UI 和输出流。注册的 IMU/LiDAR 回调先输入数据，再调用 `drain`：循环 `RunDetailed` 直到 `kNoData`；已消费但没有输出的帧继续排空，不等于播放结束。只对 `kOutput` 写状态和统计。
bag 播放结束后调用 `FlushMultiLidar` 再 drain，避免尾部组帧滞留。轨迹可以分别输出 IMU、主雷达和后轴坐标，不能无条件直接叠图。`--max_lidar_frames` 限制消费的融合帧，`--playback_rate` 控制节奏，`--wait_ui` 决定结束是否等待窗口。

# 在线顺序与线程

`FrontendNode` 持有 LIO，订阅回调驱动输入和 `Drain`，`rclcpp::spin` 管理节点。当前不是定位的三队列编排；不要把 @ref online_localization_flow "在线定位端到端流程" 的丢帧策略套在这里。点级 OpenMP 仍可在算法调用内并行。

# 调试观察点

先看初始化是否完成、IMU 是否覆盖扫描末、每雷达是否按时组帧，再看有效点数、`LastUpdateAccepted` 与关键帧数。LIO 的全局点云只是前端地图，不具备闭环优化效果。

**深入：** @ref sensor_pipeline "传感器预处理、多雷达组帧与算力预算" → @ref imu_deskew "IMU 初始化、时间积分与点云去畸变" → @ref eskf_theory "18 维迭代 ESKF：从预测到观测注入" → @ref lidar_residuals "LiDAR 残差、信息矩阵与配准参数化"。实现生命周期和缓存约束见 @ref laser_mapping_module "LIO 前端：LaserMapping"。
