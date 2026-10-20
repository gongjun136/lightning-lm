# 有界同步告警与雷达停车开关

本变更承接 [PTP 故障取证](20261020_ptp_sync.md)。PTP TX 时间戳超时的根因尚未修复；不能用停车策略变化证明时钟故障消失。

## 实际停车策略

用户批准采用既有 ±1 ms 边界。`diagnostics.radar_bounded_sync_warning_enabled: true` 时，仅新鲜驱动证据满足以下全部条件的原生 `0x04060002` 告警可不再单独产生 Safety=2：

- PTP 类型为 1，至少连续两次同步时间推进，最近推进不超过 2.5 s。
- 设备报告偏差在闭区间 [-1000000, 1000000] ns；同步时间不倒退，设备时间与主机单调时钟的增量差不超过 1 ms。
- HMS、点云、工作状态均新鲜且雷达可用，无失锁、冲突或其他当前同步错误/致命告警。夹带无法识别的原始 HMS 码时拒绝有界放行。
- 定位一致性、不确定性检查通过，无其他阻断故障。

原始 DTC 保留 FAIL。仅停车判定消费独立证据；不会直接采用忽略全部同步故障的 shadow 值。证据失效立即恢复严格判定，重新放行需重新累计两次推进。接收延迟也计入证据龄期；主机调度抖动可能保守地拒绝放行。

这不是定位精度认证，也不放行持续 26 秒失同步。实际偏差超过 1 ms 仍停车。

## 独立屏蔽开关

```yaml
diagnostics:
  safety_fault_reporting_enabled: true
  radar_bounded_sync_warning_enabled: true
  radar_fault_mask_enabled: false
```

`radar_fault_mask_enabled` 默认 false。显式设为 true 并重启定位后，只屏蔽雷达分支的 `runstate=2`，包括普通雷达故障、未知或过期雷达证据；保留定位故障/未知的 4、严重定位故障或全雷达失效的 -9。它不会清除任何 DTC，不等于故障恢复。

该开关优先于有界策略；为 true 时即使偏差超限，雷达分支的 2 也会被屏蔽。全局 `safety_fault_reporting_enabled: false` 则连 4/-9 都屏蔽，语义不同。本次部署保持全局上报 true、独立屏蔽 false。

启动日志记录开关；`SAFETY_TIME_SYNC_SHADOW` 审计同时记录 strict、shadow、bounded、实际值、角色掩码和屏蔽状态。即使关闭 shadow 审计开关，只要独立屏蔽开启，审计仍记录。

## 消息与部署

`LidarDiagnosticState` 新增 `bounded_time_sync_warning`、`time_sync_offset_ns`、`time_sync_evidence_age_us`。驱动语义版本为 `livox-dtc-v3-bounded-sync`；旧版本不能授权放行。消息结构已改变，须按 workspace.repos 在主从 Orin 重建消息、驱动/聚合器和定位消费者，再协调切换，不能混用旧二进制。

回归覆盖：边界 ±1 ms、超限、时钟跳变、停止推进、过期、原生错误、组合故障、旧接口、定位未知及严重故障；独立屏蔽开关只影响 2。源代码测试通过不等于线上验证，线上验收需实际收到四路新版本诊断和定位/Safety，并记录自然告警分支的覆盖情况。
