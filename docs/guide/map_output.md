@page guide_map_output 6. 地图、坐标与输出

将估计器中的状态、定位地图和下游业务坐标连起来。重点是变换方向、扰动约定、地图与数据库的一致性，以及位姿/速度参考点。

**阅读顺序**

1. @subpage geometry "坐标系、扰动与雅可比"：作为各算法的共同约定，按正在阅读的残差核对左右扰动和增量排列。
2. @subpage map_contracts "地图包、坐标归一化与栅格"：区分 IVox、tile 地图与导航栅格，核对导出坐标及数据库引用。
3. @subpage output_contracts "定位状态、参考点与 ROS 输出"：对照后轴位姿、杆臂速度补偿和发布门控。
4. @subpage map_alignment "新旧地图对齐"：完整学习一次性对齐、人工复核及固定输出变换。
5. @subpage map_and_support_modules "地图与运行支撑实现参考"：按需查 I/O、ROS/bag 适配、UI 和优化器等支撑代码。
6. @subpage vehicle_calibration_management "车辆标定档案与验收"：关联车辆、外参版本与地图，核对配置生成和验收记录。

**先核对两条关系**

- 地图归一化：点云、关键帧位姿和描述库使用同一个导出变换；只改点云文件会破坏数据库的坐标关系。
- 业务输出：主 LiDAR/IMU 到后轴的位姿与速度转换需要统一外参和杆臂；固定地图对齐在输出边界施加。

两条关系的公式与实现分别见地图包和输出专题。几何页可随 LIO、配准和图优化阅读反复查阅；无需在进入主线前读完所有雅可比推导。

对照旋转公式时，先核对 @ref geometry_half_angle_reading "自定义旋转指数函数的半角输入" 与 @ref geometry_jacobian_reading "左右雅可比符号"。分析离线结果时，再核对 @ref support_bag_dispatch_reading "回调返回值是否影响派发" 和 @ref support_bag_clock_reading "固定 sleep 与录包时间节拍"，避免把播放速度当成算法的实时性。

<b>带着地图文件核对：</b>

| 想验证的问题 | 源码片段入口 |
|---|---|
| 地图清单检查了哪些文件，哪里拒绝启动？ | @ref map_manifest_validation_reading "manifest 与哈希校验" |
| 新旧地图变换会不会应用两遍？ | @ref map_alignment_once_reading "变换声明与应用状态" |
| 点云、位姿和 PCD 元数据如何使用同一坐标？ | @ref map_frame_normalization_reading "归一化与 VIEWPOINT" |
| tile 离开活动范围，是否已经释放或写盘？ | @ref map_tile_lifecycle_reading "活动地图与持久化" |
| 对齐最高分为什么仍然失败，失败前写过哪些文件？ | @ref alignment_acceptance_reading "独立接受门限"、@ref alignment_persistence_reading "配置与清单的写入顺序" |

完整保存流程见 @ref offline_slam_export_reading "离线地图包导出" 与 @ref online_slam_save_reading "在线地图服务"。

模式差异：@ref ins_only_operation "CGI-430 导航的坐标与参考点契约"。下一章：@ref guide_diagnostics "诊断与排障"。
