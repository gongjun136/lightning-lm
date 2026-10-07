@page guide_diagnostics 7. 诊断与排障

完成算法主线后，按故障所在层查证据。这一章是实践与回归入口，历史现场记录须结合其代码版本和验证范围阅读。

**按问题阅读**

1. @subpage localization_troubleshooting "定位失效：按流水线排查"：先定位输入、LIO、地图匹配、融合还是最终发布层。
2. @subpage timing_contracts "时间一致性与回归"：核对传感器、轮速、消费及输出的时间边界。
3. @subpage resource_diagnostics "资源、处理时间与排队延迟"：区分算法耗时和等待时间，检查主 Orin 运行证据。
4. @subpage causal_diagnostics "因果链与时间关联"：沿输入、队列、状态和输出寻找关联，使用可复现对照。
5. @subpage compute_budget "点数预算与计算回归"：理解观测点选择、负载降级与恢复条件。
6. @subpage speed_smoothing_shadow "速度平滑影子诊断"：区分记录的候选输出和实际业务输出。
7. @subpage incident_20260925 "定位失效案例与修复证据"：按记录日期学习当前失效保护的来源和回归边界。

先记录配置、提交、地图版本、输入时间和实际发布结果，再用 @ref online_localization_flow "在线定位流程" 定位边界。公式对应模块仍从 @ref guide_lio "LIO" 和 @ref guide_localization "地图定位" 回查，故障专题保存运行证据与方法。

下一章：@ref guide_maintenance "文档维护"。
