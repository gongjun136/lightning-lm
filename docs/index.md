# Lightning-LM 开发与算法手册

这套手册面向框架维护者。先看系统怎样运行，再沿流程定位模块，最后按需进入公式与 API；当前源码是行为依据，历史性能数字不能代替当前版本验证。

## 第一层：建立全局认识

- @subpage architecture "架构、入口与对象所有权"
- @subpage build_and_run "构建与启动"
- @subpage workspace_layout "工作区与实际源码位置"
- @subpage script_contracts "脚本参数、环境与副作用"
- @subpage configuration "配置入口与运行证据"

## 第二层：沿输入到输出阅读

- @subpage frontend_flow "纯前端：IMU/LiDAR → LIO"
- @subpage offline_slam_flow "离线建图：bag → 地图包"
- @subpage online_slam_flow "在线建图：ROS → 后端与地图服务"
- @subpage offline_localization_flow "离线定位：bag + map → 轨迹"
- @subpage online_localization_flow "在线定位：队列 → 匹配 → 发布门控"

## 第三层：按模块查实现

- @subpage sensor_pipeline "传感器、多雷达与预算"
- @subpage laser_mapping_module "LIO、状态更新与关键帧"
- @subpage backend_module "后端、回环与线程"
- @subpage localization_module "定位与重定位生命周期"
- @subpage online_systems "在线系统与队列契约"
- @subpage map_contracts "地图包、坐标与栅格"
- @subpage output_contracts "业务输出、参考点与有效性"
- @subpage map_and_support_modules "I/O、ROS 适配、UI 与运行时支撑"

## 第四层：进入算法原理

- @subpage geometry "坐标系、右扰动与李群雅可比"
- @subpage eskf_theory "18 维迭代 ESKF 与轮速更新"
- @subpage imu_deskew "IMU 初始化、积分与去畸变"
- @subpage lidar_residuals "点面/点点残差、NDT 与配准"
- @subpage backend_optimization "体素 BA、回环图与 HBA"
- @subpage relocalization_theory "SOLiD/BTC 检索与几何确认"
- @subpage pose_graph_theory "定位 PGO、外推与平滑"

## 日常运行和维护

- @subpage backend_configuration "后端配置与评测协议"
- @subpage localization_troubleshooting "定位失效排查"
- @subpage resource_diagnostics "资源与端到端延迟"
- @subpage causal_diagnostics "因果链与时间关联"
- @subpage compute_budget "点数预算与回归"
- @subpage timing_contracts "时间一致性检查"
- @subpage map_alignment "新旧地图对齐"
- @subpage speed_smoothing_shadow "速度影子诊断"
- @subpage incident_20260925 "当前失效保护相关的故障证据"
- @subpage documentation_rules "随代码维护文档"
- @subpage documentation_structure "目录设计与清理原则"

## 从说明到源码

以 ESKF 为例：@ref laser_mapping_module "LIO 前端：LaserMapping" → @ref eskf_theory "18 维迭代 ESKF：从预测到观测注入" → @ref lightning::ESKF "ESKF" → API 的定义位置 → 带行号源码。反向链接放在关键类的文档注释中。网页左侧导航用于主题阅读，Classes/Files 用于定位符号。

## 运行形态速查

| 程序 | 编排者 | 主要边界 |
|---|---|---|
| run_frontend_offline / online | main / FrontendNode | 仅前端，无地图定位和完整后端 |
| run_slam_offline | main 直接持有前后端 | 同步回放、尾部 flush、地图导出 |
| run_slam_online | SlamSystem | 回调驱动 LIO，后端 worker/HBA |
| run_loc_offline | main 持有 LIO/LidarLoc/PGO | 离线直接编排 |
| run_loc_online | LocSystem + Localization | sensor、定位、高频输出队列 |

**代码依据：** @ref run_frontend_offline.cc "run_frontend_offline.cc"、@ref run_frontend_online.cc "run_frontend_online.cc"、@ref run_slam_offline.cc "run_slam_offline.cc"、@ref run_slam_online.cc "run_slam_online.cc"、@ref run_loc_offline.cc "run_loc_offline.cc"、@ref run_loc_online.cc "run_loc_online.cc"。
