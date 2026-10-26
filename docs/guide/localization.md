@page guide_localization 4. 地图定位与重定位

LIO 提供相对运动，地图匹配提供地图坐标下的约束，全局重定位在局部跟踪失效时重新给出候选位姿。定位 PGO 负责这些信息的融合与高频外推。

**阅读顺序**

1. @subpage localization_module "定位模块、跟踪与状态转换"：先区分 `Localization` 的编排和 `LidarLoc` 的地图匹配职责。
2. @subpage relocalization_theory "全局重定位：检索到接受"：对照 SOLiD/BTC、候选相对位姿组合、几何验证与连续帧确认。
3. @subpage pose_graph_theory "定位 PGO 与高频平滑"：理解相对/绝对约束、实际边残差、高频外推和重定位后的重置。
4. @subpage offline_localization_flow "离线定位：bag 与地图到轨迹"：在完整离线入口中对照同一批模块，明确与在线编排的差别。

**重点公式与源码入口**

| 关注点 | 阅读位置 | 代码落点 |
|---|---|---|
| 局部点面配准的残差与左扰动 | @ref lidar_left_update_reading "叉乘坐标与左乘更新" | `PointToPlaneRegistration`；前端 ESKF 使用另一套参数化 |
| 地点相似度与候选航向 | @ref relocalization_theory "SOLiD 描述子" | `SolidDescriptorEngine::Compute/EstimateCandidateFromQueryYaw` |
| 候选相对位姿变成地图绝对位姿 | @ref relocalization_theory "位姿组合与确认" | `SolidRelocalizer/BtcRelocalizer`、`LidarLoc::TryGlobalRelocalization` |
| 图约束的旋转残差 | @ref pgo_edge_residual_reading "半角残差的手算例子" | `EdgeSE3`；区分四元数虚部与旋转对数 |
| 图结果何时可用 | @ref pgo_frame_acceptance_reading "入帧与失败条件"、@ref pgo_reset_watermark_reading "Reset 时间边界" | `PGOImpl`、`PGO` |
| 高频输出为何没有立即贴到优化位姿 | @ref pgo_smoothing_reading "DR 预测与平滑因子" | `PoseSmoother`；确认当前因子与下一次因子的使用顺序 |

先回答“跟踪何时失败、候选何时接受、接受后哪些状态需要重置”。SOLiD-KISS 的检索与粗配准分工见重定位专题；检索分数、配准收敛和最终业务有效性要分别判断。

<b>从代码验证候选到轨迹：</b> @ref relocalization_descriptor_reading "描述子和检索分数" → @ref relocalization_direction_reading "候选变换方向" → @ref relocalization_coarse_acceptance_reading "粗配准接受门槛" → @ref relocalization_database_frame_reading "数据库坐标校验" → @ref loc_read_confirmation "连续确认"。离线对照时再检查 @ref offline_loc_fused_reading "fused 文件的位姿选择" 与 @ref offline_loc_success_reading "退出码和轨迹行的证据边界"。

回到主线：@ref online_localization_flow "在线定位流程"。继续学习 @ref guide_map_output "地图坐标与输出"，再进入 @ref guide_mapping "后端与建图" 理解地图来源。
