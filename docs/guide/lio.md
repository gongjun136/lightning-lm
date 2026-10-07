@page guide_lio 3. 传感器与 LIO 前端

这一章回答“一帧传感器数据怎样变成里程计状态”。重点是输入时间、点云坐标、预测与观测更新之间的关系，以及公式在 `LaserMapping/ImuProcess/ESKF` 中的落点。

**阅读顺序**

1. @subpage laser_mapping_module "LIO 模块与主执行路径"：先认识输入、状态、输出及一次 `RunDetailed()` 的执行顺序。
2. @subpage sensor_pipeline "预处理、多雷达组帧与 IMU 同步"：明确点时间、主雷达、外参和 IMU 覆盖要求。
3. @subpage imu_deskew "IMU 初始化与点云去畸变"：对照逐点坐标变换与 `p_compensate`，区分状态预测和区间内插值。
4. @subpage eskf_theory "ESKF：关键公式与实现"：抓住 18 维状态、预测、迭代更新及轮速观测。
5. @subpage lidar_residuals "LiDAR 残差与配准"：核对残差符号、信息矩阵、左右扰动以及前端与定位的实现差异。
6. @subpage frontend_flow "纯前端运行：离线与在线"：用独立 LIO 入口理解前端行为，再对比完整定位和建图的调用方式。

**公式与代码对照的重点**

| 数学问题 | 代码落点 | 应能解释 |
|---|---|---|
| 扫描内运动补偿 | `ImuProcess::UndistortPcl` | 为什么统一到扫描末时刻，外参在哪一步使用 |
| 名义状态与协方差预测 | `NavState::get_f/df_dx/df_dw`、`ESKF::Predict` | 状态分块、扰动方向、时间步和噪声尺度 |
| 点面观测正规方程 | `LaserMapping::ObsModel` | 点到平面的误差怎样压缩成 `HTH/HTr` |
| 迭代状态更新 | `ESKF::Update` | 为什么每轮保留同一个先验，六维观测如何影响 18 维状态 |

先理解上述关系，Anderson 加速、协方差修复及锁的逐项细节可以第二遍阅读。符号与 API 链接集中在各专题中；修改某个公式时同时核对变量单位、状态索引和调用者的接受条件。

前置约定：@ref geometry "坐标系与扰动"。下一章：@ref guide_localization "地图定位与重定位"。
