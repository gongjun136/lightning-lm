@page sensor_pipeline 传感器预处理、多雷达组帧与算力预算

# 边界和数据流

ROS 消息先由 `PointCloudPreprocess` 转成内部点型，再由 `MultiLidarFrameAssembler` 对齐多雷达扫描，最后 `LaserMapping::SyncPackages` 等待 IMU 时间覆盖。预处理、组帧和 IMU 同步是三个不同阶段；其中任一阶段接收成功都不代表产生 LIO 输出。

| 阶段 | 输入 → 输出 | 状态与失败边界 |
|---|---|---|
| 消息预处理 | Livox/PointCloud2 → 内部点云 | 型号、盲区、抽样和点时间解析 |
| 多雷达组帧 | 带 id 与时间的单雷达云 → `FusedLidarFrame` | 主雷达、周期、容差、乱序窗口和最少雷达数 |
| 外参变换 | 副雷达点 → 主雷达系点 | 保留 `lidar_id`，外参固定 |
| IMU 同步 | 融合云、扫描起止、IMU 队列 → `MeasureGroup` | IMU 必须覆盖需要的积分区间 |
| 观测选择 | 完整云 → 观测子集 | 健康反馈决定预算，输出点云有独立完整性要求 |

# 生命周期与约束

assembler 由 `LaserMapping` 持有，`AddCloud/PopReady` 在调用者同步边界内使用；自身不能由多个线程无锁并发修改。`Flush` 用于离线尾帧，不能在在线流中每帧调用来消除等待。
`LateDropCount/DuplicateDropCount/ToleranceDropCount/InvalidDropCount` 分别说明迟到、重复、容差和无效输入，不能合成一个“算法丢帧率”。统计中 `present/missing_lidar_ids` 是实际组帧证据。

副雷达外参满足 \f$p_{primary}=R_{primary,lidar}p_{lidar}+t_{primary,lidar}\f$。`DownsamplePreservingSource` 保留真实输入点，不平均来源 id。车体自身点过滤使用车辆前左上轴与主雷达原点；不能当作世界轴对齐包围盒。

# 自适应预算

`AdaptiveLidarLoadController` 使用近期处理时长、帧龄和 LIO/地图健康状态选择点数、步长及来源。高负载降级与健康恢复有独立计数窗口；定位失效后的全局搜索不应继续按健康跟踪场景裁剪输入。
`SampleSpatiallyBalanced` 按来源、方位、俯仰和距离分配配额，保留原点、原时间及输入顺序。NDT 的 density 模式按密度分配，避免改写分数分布；空间分区本身不构成可观测性证明。

**源码：** @ref multi_lidar_fusion.h "multi_lidar_fusion.h"、@ref multi_lidar_fusion.cc "multi_lidar_fusion.cc"、@ref pointcloud_preprocess.cc "pointcloud_preprocess.cc"、@ref lightning::LaserMapping::SyncPackages "SyncPackages"。
深入阅读 @ref imu_deskew "IMU 初始化、时间积分与点云去畸变" 与 @ref lidar_residuals "LiDAR 残差、信息矩阵与配准参数化"；运行配置见 @ref compute_budget "点数预算与回归"。
