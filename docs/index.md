# Lightning-LM 开发与算法手册

这套手册面向熟悉 SLAM、正在学习 Lightning-LM 架构和代码的开发者。先建立系统全貌，以在线定位串起输入、估计、地图匹配和输出，再在每个模块内对照关键公式与实现。建图、离线定位和纯前端有各自的完整阅读入口。

**章节目录**

左侧展开相应章节即可查看专题。每章导览说明阅读顺序、重点问题和源码入口；公式页与对应模块放在一起。

1. @subpage guide_overview "系统全貌与上手"：认识架构、运行模式、工作区、构建与配置。
2. @subpage guide_online "在线定位主线"：先跟踪一帧数据到业务输出，再理解在线编排与 CGI 分支。
3. @subpage guide_lio "传感器与 LIO 前端"：组帧、去畸变、ESKF、残差与纯前端运行。
4. @subpage guide_localization "地图定位与重定位"：局部匹配、全局找回、定位 PGO 与离线定位。
5. @subpage guide_mapping "后端优化与建图"：关键帧、BA、回环、HBA，以及离线/在线建图。
6. @subpage guide_map_output "地图、坐标与输出"：坐标约定、地图包、参考点转换与地图对齐。
7. @subpage guide_diagnostics "诊断与排障"：失效、时间、资源、计算预算与回归证据。
8. @subpage guide_maintenance "文档维护"：新增章节、随代码更新和导航校验。

**推荐学习路线**

第一轮读 @ref architecture "系统架构" → @ref online_localization_flow "在线定位端到端流程"，回答“各模块怎样接力”。
第二轮依次进入 @ref guide_lio "LIO 前端" → @ref guide_localization "地图定位" → @ref guide_map_output "坐标与输出"，回答“关键计算在哪里、结果为什么可信”。
第三轮读 @ref guide_mapping "后端与建图"，理解定位所需地图如何生成；再用下面的运行模式入口对比不同编排方式。

每个模块先看输入输出和主路径，再看公式、参数化及代码对应；对象所有权、锁和关闭顺序作为修改并发或退出逻辑时的参考。当前源码是行为依据，历史性能数字须结合记录版本阅读。

**怎样快速核对代码**

从 @ref online_loc_code_reading "在线定位的问题导读" 开始：先选一个要验证的问题，在正文旁阅读从真实源码提取的短片段，再用给出的例子检查条件和边界。能指出支持结论的语句即可继续主线，需要上下文时再进入完整 API/源码。页面顶部显示生成时的源码提交、未提交修改状态和时间；这些信息描述手册的来源，不代表当前已安装程序的版本。

| 正在核对的内容 | 可直接进入的验证路径 |
|---|---|
| 配置是否真正生效 | @ref config_multilidar_enable_reading "多雷达启用分支"、@ref config_extrinsic_matrix_reading "矩阵读取与校验边界" |
| 几何与回放约定 | @ref geometry_half_angle_reading "旋转增量的半角"、@ref geometry_jacobian_reading "左右雅可比"、@ref support_bag_dispatch_reading "bag 回调返回值"、@ref support_bag_clock_reading "回放节拍" |
| 输入时间与多雷达融合 | @ref sensor_point_time_reading "点时间"、@ref sensor_frame_ready_reading "组帧等待"、@ref sensor_imu_sync_reading "IMU 覆盖" |
| LIO 内部计算 | @ref imu_point_reading "去畸变"、@ref eskf_iteration_reading "迭代先验"、@ref lidar_plane_equation_reading "点面正规方程" |
| 找回地图位置 | @ref relocalization_direction_reading "候选变换方向"、@ref relocalization_coarse_acceptance_reading "几何接受"、@ref loc_read_confirmation "连续确认" |
| 融合与业务输出 | @ref pgo_smoothing_reading "PGO 平滑"、@ref output_gate_reading "发布门控"、@ref output_pose_reading "参考点转换" |
| 地图来源与一致性 | @ref offline_slam_code_reading "离线建图到地图包"、@ref map_manifest_validation_reading "地图清单校验" |
| 离线结果是否足以作结论 | @ref offline_loc_fused_reading "轨迹来源"、@ref offline_loc_success_reading "退出成功的含义" |
| 独立运行分支 | @ref frontend_online_delivery_reading "纯前端发布队列"、@ref ins_input_acceptance_reading "纯 INS 输入接受"、@ref ins_output_validity_reading "导航恢复与发布" |
| 退出与多文件保存 | @ref online_finish_order_reading "队列关闭顺序"、@ref alignment_persistence_reading "对齐配置与清单的部分失败" |
| 进程健康与安全状态 | @ref dtc_heartbeat_evidence_reading "心跳是否代表工作推进"、@ref dtc_failure_recovery_reading "故障保持与恢复"、@ref dtc_safety_decision_reading "Safety 的真实含义" |

**运行模式入口**

| 运行模式与阅读入口 | 程序 | 编排者 | 主要边界 |
|---|---|---|---|
| @ref online_localization_flow "在线定位（优先）" | run_loc_online | LocSystem + Localization | 传感器、地图匹配、高频输出队列 |
| @ref ins_only_operation "CGI-430 纯组合导航" | run_loc_online（ins_only） | InsLocSystem | 后轴导航解、严格门控、导航驱动的点云去畸变 |
| @ref offline_localization_flow "离线定位" | run_loc_offline | main 持有 LIO/LidarLoc/PGO | bag + 地图，同步编排与轨迹核验 |
| @ref frontend_flow "纯前端（离线/在线）" | run_frontend_offline / online | main / FrontendNode | LIO 里程计与局部图 |
| @ref offline_slam_flow "离线建图" | run_slam_offline | main 直接持有前后端 | 同步回放、尾帧处理、优化与地图导出 |
| @ref online_slam_flow "在线建图" | run_slam_online | SlamSystem | ROS 回调驱动 LIO，后端 worker/HBA 与地图服务 |

**代码依据：** @ref run_frontend_offline.cc "run_frontend_offline.cc"、@ref run_frontend_online.cc "run_frontend_online.cc"、@ref run_slam_offline.cc "run_slam_offline.cc"、@ref run_slam_online.cc "run_slam_online.cc"、@ref run_loc_offline.cc "run_loc_offline.cc"、@ref run_loc_online.cc "run_loc_online.cc"。

**从公式进入源码**

例如 @ref laser_mapping_module "LIO 主路径" → @ref eskf_theory "ESKF 公式与代码对应" → @ref lightning::ESKF "ESKF API" → 定义位置与带行号源码。左侧章节用于学习，类/文件索引和搜索用于定位符号。
