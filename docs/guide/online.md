@page guide_online 2. 在线定位主线

以在线定位建立系统心智模型：ROS 输入进入 LIO，去畸变点云与里程计进入地图匹配，定位 PGO 融合相对运动和地图约束，输出边界决定坐标转换和是否发布。

**阅读顺序**

1. @subpage online_localization_flow "在线定位：输入到业务输出"：优先读单帧数据流、跟踪/重定位、PGO 与发布门控。
2. @subpage online_systems "在线编排、队列与时序"：理解传感器、地图匹配、高频输出各自的处理速度与丢弃策略，并与在线建图对比。
3. @subpage ins_only_operation "CGI-430 纯组合导航"：独立学习导航解、坐标迁移、质量门控与点云去畸变。这是另一种定位分支，需按自己的输入与有效性契约阅读。

**按模块继续深入**

| 主路径上的问题 | 对应学习单元 | 首个源码入口 |
|---|---|---|
| 消息何时组成可处理的一帧？ | @ref guide_lio "传感器与 LIO" | `LaserMapping::SyncPackages`、`RunDetailed` |
| LIO 状态怎样由 IMU 预测、LiDAR 修正？ | @ref eskf_theory "ESKF 公式与实现" | `ESKF::Predict/Update`、`LaserMapping::ObsModel` |
| 局部跟踪失效后怎样重新找回位置？ | @ref guide_localization "地图定位与重定位" | `LidarLoc::Align/TryGlobalRelocalization` |
| 匹配结果怎样变成高频业务位姿？ | @ref pose_graph_theory "定位 PGO"、@ref output_contracts "输出契约" | `PGO::PubResult`、`LocSystem::PublishLocalizationResult` |

**阅读重点**

每一步都区分“收到输入”“算法形成结果”“匹配被接受”“业务输出被放行”。高频输出可以来自相对运动外推，消息频率并不代表地图匹配频率。队列与时间约束影响结果新鲜度，应和算法一起理解。

先完成这一章的主路径，再依次进入 @ref guide_lio "LIO 前端"、@ref guide_localization "地图定位" 和 @ref guide_map_output "地图、坐标与输出"。所有权图和关闭顺序留作改线程、重初始化或退出逻辑时的参考。
