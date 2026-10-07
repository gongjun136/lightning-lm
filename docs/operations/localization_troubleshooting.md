@page localization_troubleshooting 定位失效排查与取证

# 先确定失效发生在哪一层

| 现象 | 先查的证据 | 对应实现 |
|---|---|---|
| 无 LIO 输出 | 传感器 topic、组帧 drops、初始化和 IMU gap | @ref sensor_pipeline "传感器预处理、多雷达组帧与算力预算" |
| LIO 正常但地图匹配失败 | 原始修正量、静态地图一致性、分数与地图版本 | @ref localization_module "定位与重定位" |
| 全局恢复失败 | 候选数、细化拒绝原因、时间确认计数 | @ref relocalization_theory "全局重定位：检索、几何验证与时间确认" |
| 地图结果有效但融合失败 | PGO 相对位姿赋值、求解、reset 日志 | @ref pose_graph_theory "位姿图、增量求解与高频平滑" |
| 内部结果有效但不发布 | 连续失败、最近有效匹配年龄、发布门控 | @ref output_contracts "定位状态、车体参考点与 ROS 输出" |
| UI 位置异常 | 实际消息、参考点、固定地图变换和接收时间 | @ref map_alignment "新旧地图对齐" |

高分或 matcher 的收敛标记不能独立证明全帧匹配正确。排查时使用平滑之前的配准修正，避免小比例平滑掩盖大跳变。

# 采集与复现

使用 `scripts/run_sany_lidar_loc.sh` 和 @ref script_contracts "Shell 启动与实验脚本契约" 约定的环境，保存最终配置、源码/二进制身份、地图 metadata、传感器 bag 和完整日志。production 的低频日志不能重建每条实际发布消息；有需要时显式诊断采集并量化额外开销。
回放只输入传感器话题，避免同时回放旧 PosRes/path。压缩点云须先经过匹配的解压节点；不要把类型不匹配误判为算法不处理。

从 `analyze_localization_health_log.py` 检查失效时间线，再用 `analyze_localization_run.py` 和实际发布轨迹交叉核验。只有日志时标明缺少点云、CAN 或真值的证据限制。

# 当前需要保留的故障依据

@ref incident_20260925 "2026-09 定位异常证据与修复边界" 记录了静止保持、匹配器重建、迭代耗尽与失效后发布问题及其证据边界。配套 JSON 保留在同目录供追溯；它解释当前保护逻辑，不是每次定位异常的通用根因。
