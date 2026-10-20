@page incident_20261020_bounded_sync_policy 有界同步告警与雷达停车开关

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

## 现场验证结果

2026-10-20 已部署实现提交 `49b2d36`，依赖 `common_msgs=79c2a68`、`livox_ros_driver2=430d803`，完整哈希见 `workspace.repos`。主机重建消息、驱动/聚合器、定位及全部安装程序，从机重建消息和后雷达驱动；使用新的 `build_bounded_20261020` / `install_bounded_20261020`，未改动的依赖沿用原安装。安装文件与测试候选逐字节一致，已核查定位进程加载的实际库路径。

- 接口校验 8 项、清单校验 6 项通过。主机定位/心跳两项 CTest、驱动三项 CTest 通过；从机驱动三项 CTest 及隔离 ROS 域 166 的聚合器失效/恢复集成测试通过。
- 主从雷达及主机定位已重启。主机运行目录 `/home/nvidia/Sany/runs/sany_lidar_20261020_121811`，启动元数据提交为 `49b2d36`；有界开关 true、雷达全屏蔽 false、全局上报 true。
- 300.009 s 实际接收窗口：四路新版本诊断各约 3000 条，定位 41731 条、最大接收间隔 102.608 ms；Safety=0 收到 13906 条、Safety=2 收到 1057 条。接口 errors 为空，**不代表全程无故障**。
- 12:23:11，前雷达 eth3 的 ptp4l 报 TX timestamp timeout，MASTER→FAULTY；12:23:33 回到 MASTER。12:23:12.014765–12:23:37.012802，发布审计实际/严格/有界为 `2/2/2`，影子为 0，前雷达掩码 1；约 25 s 后恢复为 0。定位一致性、不确定性均 PASS。
- 接收窗口末条前雷达诊断仍为同步 FAIL、usability=1，报告偏差 -11392 ns，但 `bounded_time_sync_warning=false`、证据龄期为 UINT64_MAX（无合格证据）。说明**小的残留偏差值不能代替持续同步证据**；真实失锁未被新策略放行。接收窗口在恢复前结束，后续恢复时间由审计记录确认。
- 自启动至取证共 337 条审计样本，实际发布值均符合所选策略；包含启动 4/-9、健康 0 和故障 2。审计是限频记录，不能当作 Safety 消息计数。
- 补采恢复后的 15.003 s：实际接收 Safety=0 共 752 条，定位 1872 条、最大接收间隔 65.559 ms，四路诊断各 150 条，errors 为空；确认接收端也已恢复。原始数据为 `runtime_bounded_recovery.json`。

本轮自然事件验证了真实故障仍停车及其恢复。四路 `bounded_time_sync_warning` 次数均为 0，未遇到满足放行条件的轻微告警；该放行分支和独立全屏蔽开关由离线测试覆盖，不能宣称已有自然事件放行验证，也不能宣称 PTP 根因已修复。

结构化统计见 [验证结果](20261020_bounded_sync_validation.json)。原始证据位于现场 `/home/nvidia/Sany/validation_safety_20261020/` 及知识库 `outputs/sany_safety_20261020/`：`runtime_bounded_300s.json`、`bounded_online.log`、`bounded_audit_summary.json`、`ptp_bounded_fault.txt`、主从构建测试日志。旧完整安装包均保留在现场 `pre_bounded_install/`；回退必须协调停止相关节点并同时恢复主从匹配的消息/驱动及主机定位，不能单独替换消息库。

@anchor device_heartbeat_dtc_evidence_20261020
## 同日后续设备心跳改造

16:10 后的四路自身心跳 DTC 改造及部署记录见 @ref dtc_device_heartbeat_20261020 "四路自身心跳 DTC 部署验证"；
[机器可读验证记录](20261020_device_heartbeat_dtc_validation.json) 包含精确版本、安装哈希和实际接收计数。
该次主机沿用现场已开启的雷达屏蔽开关，且定位存在间歇失锁，不能与上文较早的现场窗口混为同一次验收。
