@page guide_diagnostics 7. 诊断与排障

完成算法主线后，按故障所在层查证据。这一章是实践与回归入口，历史现场记录须结合其代码版本和验证范围阅读。

需要从源码确认健康状态时，先核对 @ref dtc_heartbeat_evidence_reading "心跳序列与工作序列"、@ref dtc_failure_recovery_reading "FAIL 保持与 UNKNOWN"、@ref dtc_safety_decision_reading "Safety 优先级和覆盖开关"。这三条路径区分进程仍在运行、诊断证据通过和位姿实际发布。

**按问题阅读**

1. @subpage localization_troubleshooting "定位失效：按流水线排查"：先定位输入、LIO、地图匹配、融合还是最终发布层。
2. @subpage timing_contracts "时间一致性与回归"：核对传感器、轮速、消费及输出的时间边界。
3. @subpage resource_diagnostics "资源、处理时间与排队延迟"：区分算法耗时和等待时间，检查主 Orin 运行证据。
4. @subpage causal_diagnostics "因果链与时间关联"：沿输入、队列、状态和输出寻找关联，使用可复现对照。
5. @subpage compute_budget "点数预算与计算回归"：理解观测点选择、负载降级与恢复条件。
6. @subpage speed_smoothing_shadow "速度平滑影子诊断"：区分记录的候选输出和实际业务输出。
7. @subpage dtc_delivery "DTC 与安全状态接口交付"：核对源 DTC、兼容安全状态、标定与验收边界。
8. @subpage incident_20260925 "定位失效案例与修复证据"：按记录日期学习当前失效保护的来源和回归边界。
9. @subpage incident_20261008_can_schema_startup "CAN 消息版本冲突与启动静默退出"：查看接口冲突、输入预检与现场验证范围。
10. @subpage incident_20261020_ptp "雷达诊断接口与 PTP 同步故障"：区分接口不兼容与真实同步中断，查看现场证据、优先级试验及回退结果。
11. @subpage incident_20261020_bounded_sync_policy "有界同步告警与雷达停车开关"：核对有界放行条件、独立屏蔽开关和现场验证边界。
12. @subpage lidar_ins_accuracy "建图与纯雷达定位精度评估"：旁路录制 CGI 参考，以独立数据冻结地图配准，比较完整三维位姿和采样覆盖。
13. @subpage lio_scheduling_20261021 "动态联调的 LIO 调度方案与复测"：区分计算耗时与调度等待，核对单线程候选的适用条件及现场验证边界。
14. @subpage incident_capture_guide "进程内故障采集与快速排查"：了解触发、资源边界和事件包使用。

先记录配置、提交、地图版本、输入时间和实际发布结果，再用 @ref online_localization_flow "在线定位流程" 定位边界。公式对应模块仍从 @ref guide_lio "LIO" 和 @ref guide_localization "地图定位" 回查，故障专题保存运行证据与方法。

**Safety 上报与排查留痕**：`diagnostics.safety_fault_reporting_enabled: false` 将 Safety.runstate 的所有非零故障上报置为 0（当前为 2、4、-9）；`heartbeat.active_fault` 独立保留，其他消费者行为需分别验证。默认 reporting 为 true。开关仅在 LocSystem 初始化时读取，修改后须重启定位进程；实际生效配置保存在本次 `run/config.yaml`。

`SAFETY_DIAGNOSTIC` 独立于上报开关及 `radar_time_sync_shadow_enabled`，可对照 `raw_runstate`、`selected_runstate`、`published_runstate` 和 `suppressed`。协方差无效或位姿不可观测时 sigma 记为 NaN，并提供有效性标志；旧 `SAFETY_TIME_SYNC_SHADOW` 保留兼容。新日志需要相应版本，尚未部署到 Orin。

SANY 不确定度调试门限将横向 sigma 从 0.3 放宽至 0.5 m、yaw sigma 从 2° 放宽至 5°，主轴 sigma 保持 0.5 m。运动一致性、不可观测/未知及输出超时保护保留；门限调整不证明旧告警误报或真实定位精度达标。

下一章：@ref guide_maintenance "文档维护"。
