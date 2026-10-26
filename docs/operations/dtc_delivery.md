@page dtc_delivery DTC 与安全状态接口交付

接口依据：
https://sanygroup.feishu.cn/wiki/XoT4wYusciDuhKkuw5Fcfo1Bn7T
及 https://rs2flu7c17.work.sany.com.cn/sheets/shtk5ejY5B0puSUC3dbBtmCE0ec?sheet=1dfnNH 。

下面三段先核对当前 `lidar` 定位进程的运行证据；接口参数与部署记录保留在后文。手算例子假定同一进程启动周期、没有序列号回绕，不是现场验证结果。`ins_only` 的有效性和恢复计数见 @ref ins_only_operation "纯组合导航模式"。

@section dtc_heartbeat_evidence_reading 1. 心跳在跳、状态为 Running，就能证明定位正常吗？

先比较心跳发送时对两个序列号的操作：

@snippet{lineno} utils/functional_safety_heartbeat.cc dtc-heartbeat-publish

<b>只读 `Advance` 和 `load` 的区别</b>：每次发送推进 `heartbeat_seq`，`work_seq` 只是读取已有值。后者由另一个调用推进：

@snippet{lineno} utils/functional_safety_heartbeat.cc dtc-heartbeat-record-work

在 `LocSystem::PublishLocalizationResult()` 中，它位于业务 `VehiclePose` 的发布调用之后；进入此前还需通过 @ref output_gate_reading "位姿发布门控"：

@snippet{lineno} core/system/loc_system.cc dtc-business-work-progress

再看定时健康发布中的节点状态设置：

@snippet{lineno} core/system/loc_system.cc dtc-heartbeat-state-from-work

<b>状态例子</b>：已有一次业务发布，`work_seq=1`。之后位姿门控关闭，但心跳线程继续运行，后续心跳可能是 `(heartbeat_seq, work_seq)=(20,1),(21,1),(22,1)`；健康定时器仍按非零工作序列设置 `Running`。这组数值只能证明新的心跳产生过和此前完成过一次发布，不能证明新的定位结果仍在推进。

| 观察到的证据 | 在本实现中能确认什么 | 仍需继续检查什么 |
|---|---|---|
| 同一 `boot_id` 的 `heartbeat_seq` 推进 | 心跳发送路径继续执行 | 算法进度、输入/输出年龄 |
| `work_seq` 推进 | 业务 `VehiclePose` 的 `publish()` 调用已返回 | 接收端是否收到、位姿是否准确、DTC 是否通过 |
| `state=Running` | 健康回调看到非零工作序列 | 最近一次工作是否过期，不能代替原始诊断 |
| `active_fault` 为空 | 当前进程没有保留活动故障码 | 是否尚未完成检测、证据是否仍新鲜 |
| 新鲜的 DTC PASS | 对应检查在其模型与门限下通过 | 位姿发布门控与独立实测精度 |

<b>边界与停止点</b>：启动时 UNKNOWN 可以对应空故障集合；重启后的新 `boot_id` 不能直接续接旧序列。ROS `publish()` 返回也不是订阅端确认。读到两个计数器及 `WorkSequence()` 分支即可停止；要查“为何不出位姿”，转到输出门控，不能只看心跳频率。完整源码：@ref lightning::functional_safety::HeartbeatPublisher "HeartbeatPublisher"、@ref lightning::LocSystem::PublishLocalizationResult() "PublishLocalizationResult()"、@ref lightning::LocSystem::PublishHealthStatus() "PublishHealthStatus()"。

@section dtc_failure_recovery_reading 2. 检测变成 UNKNOWN，为什么故障码和停车状态还没恢复？

`Monitor::ObservePose()` 分别记住平移、航向和不确定度三个检查的明确故障：

@snippet{lineno} core/localization/localization_dtc_monitor.cc dtc-confirmed-failure-latch

<b>注意没有 UNKNOWN 分支</b>：FAIL 写入 `true`，对应子项的 PASS 才写回 `false`；缺证据不能清除此前故障。`Evaluate()` 将三个标志的或值作为 `confirmed_failure`。转交心跳时，又将当前检查结果分成两个集合：

@snippet{lineno} utils/localization_dtc_reporting.h dtc-observations-to-heartbeat

心跳活动集合的更新规则如下：

@snippet{lineno} utils/functional_safety_heartbeat.cc dtc-heartbeat-fault-update

先移除 `passed`，再加入 `failed`；未出现的码保持原状态，同批矛盾时 FAIL 优先。`optional_can_unavailable` 是一个有条件的例外：只有必需航向、不确定度通过且没有待恢复故障时，才允许清除一致性码，条件见下一节。

<b>状态例子</b>：先出现明确平移 FAIL，一致性码进入活动集合。后续位姿评估因 CAN 缺失将平移判为 UNKNOWN，即使航向和不确定度为 PASS，平移的 `failed_checks_` 仍为 `true`，不能靠断开 CAN 恢复。恢复完整有效 CAN 窗口并得到平移 PASS 后，该标志才清除；若其余检查也通过，心跳才收到用于清除一致性码的通过证据。

<b>边界与停止点</b>：活动集合是本进程内的保持状态，不表示每个码此刻都有新 FAIL，也不跨重启持久化。原始三态、新鲜度和保持标志共同参与 Safety；不能把“没有新报错”当作恢复。读到两个集合如何加入/移除即可停止，外部诊断中心的消抖不在这几段里。完整源码：@ref lightning::loc::dtc::Monitor::ObservePose() "ObservePose()"、@ref lightning::loc::dtc::Monitor::Evaluate() "Evaluate()"、@ref localization_dtc_reporting.h "定位 DTC 转交"、@ref lightning::functional_safety::HeartbeatPublisher::UpdateFaults() "UpdateFaults()"。

@section dtc_safety_decision_reading 3. Safety 怎样选择 -9、4、2、0，哪些开关会改变含义？

先看门控丢失如何进入诊断，再看输出年龄怎样成为严重条件。这里的年龄使用单调时钟；计时点是监测器接受位姿观测，不是下游收到 ROS 消息的时刻。

@snippet{lineno} core/system/loc_system.cc dtc-publication-gate-loss

@snippet{lineno} core/localization/localization_dtc_monitor.cc dtc-output-severity

曾发布过位姿后，发布门控关闭会令 `lost=true`；还没见过位姿时则使用启动宽限。最终判定按下面的提前返回顺序执行：

@snippet{lineno} core/localization/localization_dtc_monitor.cc dtc-safety-priority

<b>从上往下读</b>：严重定位条件或全部启用雷达不可用优先给 `-9`；待恢复故障、一致性 FAIL、不确定度非 PASS 给 `4`；其余一致性未通过若没有明确允许条件也给 `4`；最后才检查雷达分支并决定 `2/0`。`mask_radar_fault` 仅跳过最后的雷达分支，不能越过前面的 `-9/4`。

CAN 可选不是把 UNKNOWN 改成 PASS，它要求以下所有条件同时成立：

@snippet{lineno} core/localization/localization_dtc_monitor.cc dtc-optional-can-evidence

必需航向也变成 UNKNOWN 时，只能再检查有限的旧通过证据：

@snippet{lineno} core/localization/localization_dtc_monitor.cc dtc-unknown-evidence-age

`grace` 限制的是最后一次航向 PASS、输出观测、IMU 到达和 IMU 时间戳差的总年龄；不是每次 UNKNOWN 后重新获得一段宽限。所有条件以 `&&` 连接。

<b>手算例子</b>：假设无严重条件和待恢复故障、不确定度 PASS、启用雷达全部通过。仅 CAN 平移 UNKNOWN、航向 PASS 时，`optional_can_unavailable=true`，可以走到 `0`；若此前平移 FAIL 未恢复，则先得到 `4`。再假设 `grace=0.15 s`，航向 UNKNOWN 且四项年龄均为 `0.12 s`，可有限容忍；最后一次航向 PASS 年龄达到 `0.16 s` 时，即使其他项仍是 `0.12 s`，也会给 `4`。若曾有输出、之后超过配置的 `output_timeout_sec`，则严重条件优先给 `-9`。

最后核对实际发布的选择，避免把内部判定值当成对外消息值：

@snippet{lineno} core/system/loc_system.cc dtc-safety-output-switches

<b>边界与停止点</b>：`safety_fault_reporting_enabled_=false` 会直接发布 `0`；雷达屏蔽和有界同步选项也会改变所选雷达判定。因而单独看到 `runstate=0`，不能证明位姿在发布或全部诊断健康。Safety 是独立输出；`4` 也不等于位姿发布门控必然关闭。读到实际 `msg.runstate` 和 `publish()` 即完成本节，具体配置及历史验证继续看 @ref dtc_interfaces "接口与职责"、@ref dtc_unknown_and_cpu_20261020 "证据年龄与部署边界"。完整源码：@ref lightning::loc::dtc::SafetyValue() "SafetyValue()"、@ref lightning::loc::dtc::Monitor::Evaluate() "Evaluate()"、@ref lightning::LocSystem::PublishSafety() "PublishSafety()"。

@section dtc_interfaces 接口与职责

2026-10-21 动态联调后的 CAN 可选输入修订：CAN 未接入、冻结、过期或历史窗口缺口，
本身不再触发 `runstate=4`。平移轮速一致性仍记录 UNKNOWN，航向独立使用 IMU 检查；
必须有当前航向 Pass（含预热窗口）及有效不确定度，且不存在待恢复的明确故障。
短时 IMU 证据缺失仍受下文 150 ms 总证据年龄约束，输出超时、门禁丢失、全雷达失效的 -9 条件不变。
无 CAN 时缺少独立平移验证，不能把 UNKNOWN 解读成平移精度已通过。
SANY 模板及代码默认诊断版本为 `localization-dtc-v3-optional-can`；更新后核对启动版本及配置快照。

平移、航向、不确定度分别保存明确故障。CAN 断开不能解除已确认平移故障，
需恢复完整有效 CAN 窗口并重新通过；航向故障可在无 CAN 时由独立航向 Pass 恢复。
对外保留原 DTC 码及消息定义；心跳清除一致性码时要求没有任何待恢复子项且必需检查已通过。
日志增加 `translation_condition`、`yaw_condition`、`optional_can_unavailable`、`yaw_pass_age_sec`。
CAN 单独不可用不计入 `unknown_tolerated`，该字段只表示必需 IMU 证据的有限容忍。

启动入口使用 `check_localization_inputs --optional-wheel-speed`，缺 CAN 不阻止启动，
仍检查实际 LiDAR/IMU 收包；单独手动使用旧 `--wheel-speed` 仍可做严格 CAN 检查。
融合保留既有 250 ms 传感器时间年龄门限、未来时间容差和创新门限，不因本次修订放宽。
CAN 恢复后仅有效且符合时间门限的样本参与融合，不能以新接收时间替代旧 Header 时间。
LIO 共享工作线程等待及两核调度对照方案见 [动态联调 LIO 调度说明](lio_scheduling_20261021.md)。

本轮本地 WSL Ubuntu 22.04 / ROS Humble 验证：`run_loc_online`、启动探针及三个回归目标构建成功；
`localization_dtc_test`、`localization_pgo_test`、`heartbeat_dtc_test` 3/3 通过。
真实 ROS 输入探针验证四项：可选 CAN 缺失通过、严格 CAN 缺失失败、必需 IMU 缺失失败、必需 LiDAR 缺失失败。
启动脚本回归 8/8、依赖接口回归 9/9 通过；真实三雷达进程 smoke 通过，Safety steady 平均约 50.1 Hz。
挂载盘启动脚本测试使用 `LIGHTNING_TEST_LAUNCH_TIMEOUT_SEC=120`（默认测试预算仍为 30 秒）；
初次 30 秒预算发生 TimeoutExpired，未作为通过记录，不修改实车超时参数。
这些验证包含构造的无 CAN、5 m/s 输入、故障保持/恢复及停走回归，不是 ARM 实车负载验收或整包原始数据回放。

对外 DTC 统一通过 diagnostic_monitor_interfaces/msg/NodeHeartbeat.active_fault 上报。
使用 workspace.repos 固定的 common_msgs 既有定义，不增加公共字段。
定位 ID 6 的 /diagnostics/heartbeat/lightning_slam 负责
0x30811813（一致性异常）、0x40450198（不确定度超限）；
四台雷达各自在自己的心跳 active_fault 上报自身九类当前设备 DTC：
前 /diagnostics/heartbeat/livox_front（ID 7）、后 livox_rear（ID 12）、
左 livox_left（ID 8）、右 livox_right（ID 9），均使用 /diagnostics/heartbeat/ 前缀。
每台驱动的检测快照直接更新自己的心跳，设备码不集中上报。
不再创建未登记的聚合节点 ID 13 或 /diagnostics/heartbeat/livox_dtc。
0x10F00001/0x10F00002 没有已登记的节点归属，不自行挂到定位或其他设备心跳；
多台/全雷达不可用仍在定位 Safety 内部组合判定，不依赖这些外部汇总码。
心跳默认 10 Hz、best_effort、volatile、KeepLast(1)，每次携带排序去重后的完整活动集合。
FAIL 加入，明确且新鲜的 PASS 才移除，UNKNOWN、缺项或过期保留此前活动状态。
零码和历史码不进入 active_fault；雷达四角色的 base+12/base+13 仅保留在内部观测。
源端不重复实现诊断中心的 1 秒确认和 2 秒恢复消抖。
HeartbeatPublisher::UpdateFaults 原子更新一批判定，发布线程复制快照，发送不推进 work_seq。

旧 /diagnostics/source_dtc/localization 和 /diagnostics/source_dtc/livox 不再发布。
定位按本次启用的物理角色消费 /sensors/livox/{front,rear,left,right}/diagnostic_state，
每台为 livox_ros_driver2/msg/LidarDiagnosticState，10 Hz、reliable、volatile、KeepLast(10)。
每路保留十一项 DTC 观测、设备可用性和独立同步证据；
physical_condition：0 UNKNOWN，1 PASS，2 FAIL；保留检测序号、证据年龄和历史观测。
已移除独立聚合进程和 /internal/diagnostics/livox_dtc；
旧 diagnostics.internal_radar_topic、radar_source_id 与聚合器配置不再使用。
升级须同步主从驱动和定位，并停止旧聚合进程；使用新的安装目录以避开残留旧可执行文件。
现有 /localization/fault_status 留作兼容遥测，统一诊断节点只消费心跳 DTC，
不得把内部观测或兼容状态再次接入故障确认/动作链路。

活动集合是进程内状态，不跨重启持久化。启动 UNKNOWN 不生成虚构故障码；
空集合不能证明检测已完成或系统健康。boot_id 改变后的历史故障保持及就绪判定仍由
诊断中心负责，本次不修改其实现，也不假设 state=3/4 会阻止恢复。
保留旧 FAIL 的兼容方案不能表达“当前证据未知”，中心仍可能据此累计确认时间；
这是现有公共接口的表达限制。Safety 继续使用原始三态及新鲜度，不使用保留后的活动集合。

内部观测的 comm_header.seq 表示发送方本轮健康评估；evaluation_seq 表示该项检测所依据的
新证据/检测进度；evaluation_age_us 是最近实际评估的稳态时钟年龄。
从未完成检测时为 UINT64_MAX。boot_id 与同节点心跳共用启动标识。
接收方应同时检查来源、版本、完整码集合、启动标识、序号及证据年龄。
UNKNOWN、缺项、超时、重启不能当作 PASS，也不能清除已确认故障。
定位源 ID 默认 6，四台雷达 ID 为前 7、后 12、左 8、右 9；
诊断中心保持上述已有节点 ID 与各自设备 DTC 的对应关系，不新增汇总节点。

定位另发布 /safety/throstle，复用既有 geosun_msgs/msg/Safety，字段 int8 runstate，50 Hz。
Safety.msg 定义保持不变，不新增兼容消息类型；header.stamp 为发布时的 ROS 时间，
timesync、errorcode、strategies、reserved 保持默认 0，header.frame_id 默认空。
Topic 拼写 throstle 保留；接收端订阅 Safety 并读取 runstate。
`config/sany/localization.yaml` 的 `diagnostics.safety_fault_reporting_enabled`
控制该兼容接口的故障上报，默认及缺省值为 `true`。设为 `false` 并重启定位节点后，
@ref lightning::LocSystem::PublishSafety() 始终以 50 Hz 发布 `runstate=0`，
包括启动 UNKNOWN、定位失锁、超时及雷达故障；时间戳仍正常更新。
该开关只覆盖 Safety 输出，不改变源 DTC、心跳或 `/localization/fault_status`。
开关为 `true` 时使用以下判定规则。
优先级 -9 > 4 > 2 > 0：定位失锁/有效输出超时或本次启用的雷达全部不可用给 -9；
定位 DTC FAIL 或尚未恢复的明确故障给 4；其他雷达故障给 2。
定位检测 UNKNOWN 不满足 @ref dtc_unknown_and_cpu_20261020 的有限证据条件时给 4；
满足条件的一致性 UNKNOWN 继续检查雷达分支，全部分支允许时给 0，原始观测仍为 UNKNOWN。
仅雷达检测 UNKNOWN/过期给 2。
已恢复通信历史（四个角色码的 base+12）与同步历史（base+13）仍保留在源 DTC 中，
但新鲜的历史 FAIL 不单独参与当前 Safety 故障判定；历史 UNKNOWN、过期、其他当前 DTC
以及本次启用设备的实际可用性仍须校验。因此历史事件本身不会让正常设备永久保持 runstate=2。
已确认的全雷达不可用在重启 UNKNOWN 或通信中断时维持 -9，直至至少一路启用雷达恢复新鲜可用证据。
该兼容接口执行本地实时停车判定；内部原始观测及诊断中心的等级与恢复策略仍各自保留。
兼容接口同时读取启用雷达 diagnostic_state 的设备可用性，因此单台未映射 native fatal
或实际不可用也给 2；设备可用性未知时不能给 0，不新增未定义的系统 DTC。

@section dtc_lidar 雷达实现与部署

标准启动脚本通过 `/opt/conf/vehicle_config.yaml` 解析本次有效配置。
沿用现有约定：`LIDAR: MASTER AND SLAVE` 为前、左、右三雷达，`LIDAR: MASTER` 为前、左、右、后四雷达。
Safety 从同一运行快照的 `vehicle_configuration` 和 `multi_lidar.topics` 确定必检角色，
不在运行中重新读取车辆文件，也不根据话题在线状态或自适应点云抽样减少必检雷达。
三雷达配置不订阅后雷达诊断，不因后雷达缺失给 2；三路全部不可用仍给 -9。
配置矛盾、缺失所需话题或未经解析的 SANY 模板会拒绝启动。
不含 SANY 车辆快照的旧配置保留四路检查。角色 ID 和同步故障掩码仍按前、后、左、右排列，不能按启用数量重新编号。
启动日志逐项输出 `LOCALIZATION_RADAR_REQUIRED role=... required=...`，用于核对实际检查范围。

运行检查优先指定同次定位使用的有效配置：
`python3 scripts/check_runtime_interfaces.py --workspace ~/Sany/lightning_lm_ws --config /path/to/run/config.yaml`。
未传 `--config` 时，脚本用共享模板和当前 `/opt/conf/vehicle_config.yaml` 解析所需雷达；
若启动后改过车辆文件，必须传运行快照，避免检查范围与正在运行的进程不同。
检查结果记录 `required_radar_roles` 和车辆配置来源，仍要求每个启用雷达收到真实诊断和心跳数据。

2026-10-21 静止包 `sany_lidar_20261021_111553`（旧二进制源码 `d187c2e`）仍固定检查四路，
实际为三雷达，缺少后雷达会阻止 Safety 恢复 0；前、右雷达另有同步告警，不能将全部 2 都归为拓扑问题。
用户提供的 PID 100549 核验与包内资源元数据一致：关键线程绑核正确，稳定阶段 894 个资源样本中未发现线程掩码越界；
但内核 `isolated_cpus=[]`、启动参数缺失隔核项、默认 IRQ 掩码 `fff`，这台域控的隔核尚未生效。
约 15 分钟保留日志中，启动阶段约 5.16 s 为 4；随后仅在 11:28:00.527773 记录一次 4，
到 11:28:00.536057 恢复 2，日志切换间隔约 8.28 ms。该次一致性为 UNKNOWN、非确认 FAIL，
CAN/IMU 接收年龄约 107 ms，位姿输出年龄约 0.37 ms；同秒执行器/传感器线程累计运行队列等待约 214/216 ms。
附近一次雷达处理墙钟耗时 188 ms、线程 CPU 耗时约 15 ms，支持排查调度及等待，但这不能替代对停车分级的审查。
旧版本代码中的直接触发链为：接收年龄超过 `max_input_gap_sec=0.1` → 一致性 UNKNOWN →
未获有限容忍的 UNKNOWN 映射为 4；本次 `consistency_unknown_grace_sec=0`。
旧版本有限容忍还要求 CAN/IMU 接收年龄均不超过 100 ms，因此仅开启该参数也不能覆盖本次 107 ms 事件。
同事件附近高频 TUM 的相邻源时间戳从 `1792553280.310244` 跳到 `1792553280.416021`，
间隔约 105.8 ms，也会触发旧代码清除此前一致性 Pass 的分支；源时间戳间隔并不等于 ROS 实际发布间隔。
这条日志证明的是短时失去一致性校验依据，并未检测到明确定位不一致；应分别定义校验有效性与持续失去校验依据的停车边界。
静止包不足以单独标定自动行驶的容忍上限，仍需运行车速边界及实际负载验证；明确 FAIL 与原有 -9 条件保持。
用户补充自动行驶上限为 3–5 m/s，且控制端收到 runstate=2 同样停车；不能将 4 降为 2 作为继续行驶的缓冲。
5 m/s 下 107 ms 对应约 0.54 m 行程，这不是已测得的定位误差，也不能直接用行程代替误差预算。
资源采样同时出现名为 `ros2` 的进程，具体命令不明，不能认定其为唯一原因。
因果追踪未发现超过 100 ms 的系统时钟跳变。保留日志未见 -9，但有 13 条异步日志丢失（含 1 条 critical）
和 32 条因果记录丢失，不能作为无损的控制端接收验收。
后续先按本机硬件核验并应用隔核、重启，确认 `isolated=10-11` 和默认 IRQ 掩码 `3ff`，
再在相同生产负载和采样设置下复测。本次保留定位 -9 条件，将 SANY 的一致性 UNKNOWN
总证据年龄上限设为待实车验证的 150 ms；不能据此宣称行驶安全或此前所有停车均已消除。

每台驱动只绑定一个配置设备，读取 Livox SDK2 状态 JSON 和周期内部信息查询。
hms_code 必须是完整 8 个 uint32；高 16 位为 HMS ID，低 8 位为原生等级。
同类别多原因按 OR 合并，完整全零快照才能解除 HMS 原因。
设备 native fatal、work_state 故障、点包断流、连续无有效点云算设备不可用。
冻结点云源时间戳不刷新有效性；非零有限 XYZ 点数门限可配置。
低温加热 0x0201 不生成 DTC；0x0401、0x0405 是历史提示，不能单独算当前设备失效。
`livox-dtc-v2-sync-source` 增加同步源验证：仅在类型为 PTP、最后同步时刻非零且持续前进、
设备时间与同步时刻相差不超过 2.5 s、绝对 PTP 偏差不超过 1 ms，且验证证据仍新鲜时，
不将未使用 GPS/PPS 源的 Warning（0x0404、0x0407–0x0409）计入当前授时故障。
这些值分别采用驱动 `dtc.state_timeout_ms` 和 `dtc.max_sync_offset_ns`；
缺字段、查询超时、同步时刻冻结/回退时不应用该排除。PTP 故障、低精度、PTP/gPTP 冲突、
原生 Error/Fatal 保持上报。原始 HMS 历史信息不清除，也不把真实断流当作历史故障。
未知 native fatal 仍计入不可用，不杜撰新的系统 DTC。
多台与全部启用雷达故障可同时成立，按启用的物理角色去重统计。

master_lidar.launch.py 启动前、左、右三台驱动，master_4lidar 间接复用它；
slave_lidar 只运行后雷达。四雷达部署须确保四个 diagnostic_state Topic 跨主机可达；
三雷达部署仅要求前、左、右三路。
驱动 dtc.* 参数
可设置角色、超时、有效点数和配置版本。默认点包/点云 500 ms、HMS 2500 ms、
内部状态 500 ms、启动等待 8 s，设备状态每 1 s 查询。
遗漏设备在启动宽限后按不可用处理，因此三雷达实车必须显式补齐后雷达输入才能正常。
主、从驱动与定位须使用 `livox-dtc-v3-bounded-sync` 语义并重新构建、重启后生效；
混合配置版本的设备状态会被拒绝，不得只更新主机就宣称后雷达完成修复。

2026-10-15 现场只读查询确认前、左、右三台 `work_state=1`、`lidar_diag_status=0`，
PTP 偏差约 -5 微秒，最后同步时刻连续前进；三台均保留 0x04050002，左雷达另有
0x04090002（GPS 信号异常）。原实现未区分同步源，并将历史事件参与 Safety，导致误报。
ROS 快照中四个角色的可用性均为 USABLE，后雷达点云时间戳也在前进。
现场接线描述与 ROS 图不一致时，应核对发布者、配置和持续前进的源时间戳，不能直接跳过后雷达检查。
随后经 `sany-03-slave` 只读查询确认后雷达配置指向 `192.168.4.143`：
设备同样为 Sampling、`lidar_diag_status=0`、PTP，偏差约 -8 微秒，最后同步时刻连续前进，
HMS 为 0x04050002 与 0x04090002，属于上述同类误报。12 s 内后角色 121 条诊断均为 USABLE。
从机现场驱动源码为 `fbf39f93`，需随主机一起升级；另有一个处于 T（暂停）状态的旧驱动进程，
仍占用同一套 UDP 端口。当前只有另一实例持续发布；部署时应结束旧实例后启动唯一驱动，
不要恢复暂停实例，否则可能重新引入端口竞争及 boot_id 防重放冲突。
该次诊断没有重启雷达或修改现场授时配置。
修复在本地完成驱动和定位程序构建；`lidar_dtc_monitor_test`、`sany_ptp_guard_test`、
`localization_dtc_test` 与实际聚合节点集成测试均通过，真实 Safety 节点 smoke 约 49.9 Hz。
回归覆盖健康 PTP 配合 GPS Warning、证据冻结/回退/过期、原生 Error/Fatal、历史告警、
当前授时异常和设备不可用；上述结果是本地验证，不表示已部署到车辆。

@section dtc_localization 定位检测与协方差

一致性检查使用原始有符号 CAN 轮速和经 IMU 到车体旋转后的原始角速度，
积分到与最终后轴位姿相同的时间窗口，检查纵向、横向及航向误差。
同时检查平滑前位姿和最终发布位姿；相邻增量检查覆盖窗口预热时的突跳。
输入间隙、无效值、重复/回退时间戳、缺少完整窗口均为 UNKNOWN。
CAN 缺失只使平移子项 UNKNOWN，航向继续独立检查；UNKNOWN 不自动等于停车，详见上方 CAN 可选输入修订。
重定位/图重置推进 pose_epoch，跨 epoch 不比较，重新预热检测。
rpm 比例、轮速符号、时间对齐、陀螺偏差和安装方向须用实车独立参考标定；
平面运动近似需在项目道路坡度范围内验证。

不确定度取实际 PGO 约束在优化点的 J'ΩJ，使用鲁棒 rho' 权重，
不加入 LM 阻尼；Schur 消去其他状态，得到最新位姿完整 6×6 边缘协方差。
窗口达到五帧时，将全部历史状态边缘化到最新状态作为完整先验，
移除已消费的观测，避免先验与历史观测重复计数。
先验残差旋转符号与 VertexSE3 更新符号的差异通过变换信息矩阵处理。
高频 DR 传播使用运动 Jacobian 和配置过程噪声，平滑偏移加入保守误差项。
静止保持继续用地图观测估计后验，保持发布位置固定，并计入保持位置到融合均值的偏移。
后轴杆臂及固定地图变换传播保留平移/旋转交叉项。
最终检测 XY 椭圆长轴、车辆横向及航向标准差。
数值无效或缺失协方差给 UNKNOWN，可观测性丢失给 FAIL。
模型方差倍率与标准差下限防止相关观测和模型误差导致过度自信。

当前默认 `covariance_variance_scale=1.0`，版本为
`localization-dtc-v2-model-scale`。原台架默认值 4.0 会将标准差放大两倍，
而 PGO 已使用 0.3 m / 1° 的观测噪声，并传播 DR 过程噪声、后轴杆臂和输出偏移。
这与 0.3 m 横向门限不匹配，连无误差的正常 PGO 轨迹也会持续触发 FAIL。
移除未经独立参考标定的额外放大；0.5 m 长轴、0.3 m 横向、2° 航向门限，
0.1 m / 0.5° 下限，以及失锁、超时、不可观和 UNKNOWN 的处理均保持不变。
额外方差倍率仍可配置，但必须根据独立误差覆盖率标定，不能仅用轨迹看起来平稳来确定。

升级后须重新由 `config/sany/localization.yaml` 生成运行配置并重启定位。
既有 `runs/*/config.yaml` 是不可变运行快照；更新代码不会把其中显式写入的 4.0 改成 1.0。
排查 `runstate=4` 时同时查看 `LOCALIZATION_DTC` 日志及定位心跳的 active_fault：
`consistency=1 uncertainty=2 severe=0` 表示不确定度超限；条件 0 表示证据不足。
启动日志 `LOCALIZATION_DTC_CONFIG` 显示实际生效的版本、倍率和三个门限；
周期日志 `can_stamp_lag_sec` / `imu_stamp_lag_sec` 显示最后位姿时间减相应最新输入时间，
缺输入时为 NaN。持续的大偏差应检查时钟同步，不能通过增加运动容差消除。
`Safety.errorcode=0` 是兼容消息的默认值，不能据此判断没有 DTC。

@section dtc_calibration 标定、测试与验收

config/sany/localization.yaml 的 diagnostics 和 localization_pgo 为台架起始值，
需在联调后冻结 config_version。不得将建议阈值当作已验收的安全响应指标。

1. 采集静止、低速、倒车、弯道、颠簸、各坡度及地图边缘场景，
   同步独立 RTK/全站仪参考、输入时间戳、LOCALIZATION_DTC 日志、定位及启用雷达的心跳与设备观测。
2. 正常数据统计纵/横/航向残差分布，按速度、转弯、坡度分层标定基础容差和比例容差。
   时间不同步或参考自身无效的片段单独处理，不能用于放宽容差。
3. 用独立参考误差标定协方差方差倍率、下限及过程噪声，检查二维椭圆和航向覆盖率；
   窗口相关性及静止保持偏差必须纳入。推荐以 95% 覆盖率为起点，项目另行确定验收目标。
4. 注入单设备 HMS 多原因及清除、未知 fatal、断网、空点云、冻结时间戳、
   部分/全部启用雷达失效、驱动退出、内部消息重放、进程重启及乱序。
   分别验证三/四雷达布局下各路十一项完整性、各自心跳九类当前故障的归属、UNKNOWN 不解除故障及全失效 Safety。
5. 注入位姿/航向跳变、轮速/IMU 停止或冻结、协方差超限/无效、重定位及定位输出停止。
   验证 0/2/4/-9 优先级、故障发生与恢复时间、50 Hz 输出及接收端消息兼容。
6. 用稳态时钟测量故障注入到源 DTC、兼容状态和诊断中心安全动作的各段延迟。
   默认上界预算：设备 HMS 检出约 1 s 查询 + 100 ms 评估 + 网络；
   断流约 500 ms + 100 ms 评估 + 网络；定位有效输出超时 300 ms + 20 ms 状态周期；
   运动窗口 500 ms，突跳可由相邻增量提前发现。
   对外 DTC 改由心跳发送后，还需计入最多一个心跳周期（默认 100 ms）；
   检测变化只更新缓存，不额外插发心跳。Safety 仍使用独立的 20 ms 发布周期。
   这些是配置预算，最终最大响应时间必须在目标硬件负载及网络环境中实测验收。

自动回归：雷达 lidar_dtc_monitor_test、heartbeat_dtc_test、driver_dtc_reporting_test；
定位 localization_dtc_test、radar_dtc_state_test、heartbeat_dtc_test、localization_pgo_test。
驱动回归通过生产发布函数向真实四路 ROS 心跳注入九类 HMS，覆盖角色隔离、历史过滤、
UNKNOWN/过期保持、逐码恢复和有界告警仍保留原 DTC。
定位心跳回归使用真实运动一致性和协方差检测器产生两项 DTC，验证真实 ROS 接收及独立恢复。
PGO 回归直接将实际优化器后验送入 DTC，覆盖五帧窗口边缘化后的健康轨迹；
定位回归另外验证 0.31 m 横向、2.1° 航向超限仍为 FAIL。
radar_dtc_state_test 验证三/四雷达配置、无效配置拒绝、未启用雷达不参与判定，以及设备快照的角色/码集合/版本校验、重放和过期拒绝、
重启 UNKNOWN 保持严重状态、历史过滤、真实可用性以及有界 PTP 边界与证据年龄。
src/test/dtc_online_smoke.py 检查实际 run_loc_online 的心跳、50 Hz Safety、按配置直接订阅（`--sany-lidars 3` 或 `4`）
以及旧聚合/旧源 DTC 话题停用。历史 2026-10-19 的聚合器验证记录仅适用于当时实现。

历史 SANY bag 的 `SpeThrCAN4` 可能只有 `Header+x+y`，新版多了 `TopicCommHeader`。
`OnlineBagPlayer` 原样发布 CDR 字节，由在线输入已有的新旧布局解码器校验，
避免把旧 header 误读成通信头，进而破坏传感器时间戳和一致性检测。
回放健康状态时，若 bag 没有设备 HMS/诊断消息，不能把缺失雷达证据认定为正常；
隔离测试可明确注入健康雷达夹具，但该结果只验证定位和状态组合，不证明真实雷达健康。

2026-10-15 修复回归保留了以下证据，不能将其表述为实车安全验收：

- sany-03-master 现场日志中，一致性 PASS、不确定度 FAIL 的 589 个周期样本，
  横向标准差中位数约 0.498 m；使用相同协方差和下限、仅移除额外 4 倍方差放大后，
  对应中位数约 0.263 m，589 个样本均不再超过三个不确定度门限。
- 宝马展复刻场地 `sany_3lidar_20261007_182426` 的 111.88 s 录包，
  CAN / 主 IMU 相对录包时刻的中位偏差分别约 -0.34 / -0.64 ms。
  同一地图和配置仅改变倍率进行对照；统计首个位姿后 2 s 至最后位姿的完整区间，
  旧倍率的 5842 个 Safety 样本全部为 4；新倍率得到 5662 个 0、48 个 4、52 个 -9。
  新倍率的不确定度未出现 FAIL，但仍观测到输入证据短缺、一次输出超时和约 2.72° 的
  航向一致性残差，故不能声称整段回放全程正常，也未为消除这些状态而放宽门限。
- 筑工搅拌站 `sany_4lidar_can_on_20260826_174223` 的 33.52 s 静止录包，
  同一统计区间得到 1238 个 0、2 个 4、0 个 -9；不确定度和运动残差没有 FAIL，
  两个 4 对应一致性 UNKNOWN。上述回放均使用明确标记的健康雷达夹具。
- 筑工搅拌站 `sany_4lidar_can_on_20260826_174022` 的 95.46 s 含运动诊断录包，
  同一统计区间得到 2807 个 0、1332 个 4、519 个 -9。
  周期日志仍有横向标准差超过 0.3 m（最高约 0.362 m）、航向一致性超限和输出超时。
  该诊断录包不能当作“所有输入与定位都健康”的正例；当前修复不把这些实际触发的保护改成 PASS。
  本次数据没有独立精度参考，因此不足以证明移动工况的误报率或误差覆盖率满足实车验收要求。
- 筑工搅拌站 2026-08-19 两段 CAN 与雷达/IMU 相差约 37 s 的录包已按用户要求移除，
  原包及本次临时校时副本均移入回收站，不用于正常工况验收。

本次重新执行的七项定位相关 CTest、雷达 DTC 与 PTP 启动保护两项 CTest 均通过。
真实节点 smoke 测试约 50.1 Hz，保留启动 UNKNOWN=4、无定位输出超时=-9。
smoke 夹具通过 `export_map_manifest` 生成合法地图清单；需同时构建该程序和 `run_loc_online`。

@section dtc_validation 开发验证与部署准备

2026-10-19 在 ROS 2 Humble / Ubuntu 22.04 完成既有 dev 消息定义、雷达驱动/聚合器、
在线/离线定位程序及相关测试构建。九项 CTest 通过：两个包各自的 heartbeat_dtc_test，
雷达 DTC、PTP 启动保护、定位 DTC、后轴位姿、三一定位输出、PGO、定位输入并发。
另有三项实际节点 ROS 验证通过，均使用 Cyclone DDS 和隔离 domain：
test/dtc_ros_integration.py 通过完整快照、两台/四台失效、恢复、重放、HMS 过期和断流，
并验证 UNKNOWN、过期 PASS、设备重启不清除已活动码，历史码不进入心跳、旧话题无发布者；
src/test/dtc_online_smoke.py 启动真实 run_loc_online，验证心跳 DTC 和旧定位 DTC 话题停发，
且 Safety 只订阅新内部雷达观测话题。默认 Safety 故障上报平均 49.7 Hz，
启动 UNKNOWN 为 4，超时后为 -9；关闭故障上报时平均 50.0 Hz，runstate 始终为 0。
两次 smoke 均校验消息类型、时间戳、心跳序号和无业务输入时 work_seq 不递增。
文档源链接/导航检查通过（59 页）；全量 HTML 校验仍有原有 sany_configuration.md 的
resolve 链接告警，不属于本次 DTC 代码回归。上述均为本地验证，未部署车辆。
频率测试使用与 C++ WallTimer 相同的 ROS 稳态时钟；header 为 ROS 发布时刻。
本机 WSL 的校正时钟与原始稳态时钟曾出现约 10% 的速率差，header 频率仅作为诊断值。
ROS 稳态时钟实现参考 https://github.com/ros2/rcutils/blob/humble/src/time_unix.c 。
后者使用无传感器的临时空地图，故正常运行的 0/2/4/-9 组合由定位 DTC 单元测试覆盖。

在 lightning_lm_ws 根目录完成构建后，可复现节点验证：
@verbatim
source /opt/ros/humble/setup.bash
source install/setup.bash
RMW_IMPLEMENTATION=rmw_cyclonedds_cpp ROS_DOMAIN_ID=173 ctest --test-dir build/livox_ros_driver2_core --output-on-failure -R '^(driver_dtc_reporting_test|heartbeat_dtc_test|lidar_dtc_monitor_test)$'
RMW_IMPLEMENTATION=rmw_cyclonedds_cpp ROS_DOMAIN_ID=172 python3 src/lightning-lm/src/test/dtc_online_smoke.py build/lightning_lm/bin/run_loc_online
@endverbatim

部署时须同步 common_msgs 的三个更新消息包、驱动及定位程序，并重新构建消费者。
使用一致版本的安装目录，检查 lightning_lm/lib 内历史复制的 ROSIDL 库是否遮蔽新消息包。
升级 NodeHeartbeat 后须将该目录内旧 diagnostic_monitor_interfaces 库副本备份移出，
用 ldd 确认程序实际加载 install/diagnostic_monitor_interfaces/lib 中的新库；
仅重新生成头文件不能修复被旧类型支持库遮蔽的运行环境。
目标硬件的 Livox 状态查询、真实传感器数据、跨主机 DDS 及故障响应时间仍按上述流程联调验收。

@section dtc_device_heartbeat_20261020 四路自身心跳 DTC 部署验证（2026-10-20）

本次移除未登记的聚合节点 ID 13，四路设备及定位统一沿用既有 NodeHeartbeat 类型，
各自的 active_fault 只包含自身当前 DTC。多台/全雷达不可用由定位直接消费四路
diagnostic_state 后组合 Safety，不另造聚合心跳或外部系统 DTC。
完整的版本、安装哈希、接收计数及验证边界见
@ref device_heartbeat_dtc_evidence_20261020 "设备心跳验证记录"。

本地从固定版本的新消息源码重新生成两项变更相关消息包并构建消费者。
驱动四项、定位七项 CTest、两次真实 run_loc_online smoke、九项接口检查及六项依赖门禁回归通过。
四路 ROS 回归逐一覆盖 36 个当前设备故障码的归属、UNKNOWN 保持、逐码恢复及历史过滤；
定位心跳用真实一致性/协方差检测器验证两项 DTC。两次 smoke 的 Safety 均为约 50.0 Hz。
公司 dev 消息包与已验证公共定义合并后另外完成两包全新编译，保留同事接口及公司 CI。

sany-03-master 完成 ARM 驱动四项、定位七项 CTest 和两次真实节点 smoke；
sany-03-slave 完成后雷达 ARM 四项 CTest。目标机使用已安装的 Fast DDS，
回归在 localhost-only 的隔离域执行。主从安装包切换后已核对实际运行的二进制 SHA256，
并确认定位从活动 install 加载本次 liblightning_lm.libs.so。
直接用 CMake 安装新前缀时还须保留 colcon 的依赖索引和库路径 hook，
合并 package.dsv 中对应 hook 行；本次补齐该元数据后，地图工具和现场入口均已成功启动。
不能仅因新前缀有 .so 文件就认为工作区 setup 会加入它的库目录。

60.02 s 生产网络实收中，四路设备诊断及五路心跳各约 600 条，Safety 为 3005 条。
十一项接口均只有一个发布者，四组设备诊断与自身心跳 boot_id 一致，旧聚合心跳无发布者。
定位心跳实际收到 0x40450198；四台雷达本次九类当前观测均为 PASS，心跳为空符合判定结果。
真实车辆未注入硬件故障，非空雷达心跳由上述隔离 ROS 回归验证。

完整运行检查并未全通过：位姿只收到 195 条，最大间隔约 32.48 s，结束时不新鲜，
Safety 为 -9/4/0 分别 2942/28/35 条。切换前的 10 s 基线同样没有位姿且 Safety 全为 -9，
旧定位日志已有 pose_enabled=0、stale=1；新日志显示匹配失败及重复重定位。
本次未改 NDT/PGO 源码或放宽定位门限，不能把 DTC 接口正确表述为定位连续性或整车健康验收通过。
主机原有 radar_fault_mask_enabled=true 已原样保留，也不能用 Safety=0 证明雷达健康。

@section dtc_unknown_and_cpu_20261020 Unknown、明确失败与关键线程隔离

Safety 的判定顺序为：

| 条件 | runstate | 本次边界 |
|---|---|---|
| 定位发布门禁丢失、输出超时，或全部雷达不可用 | -9 | 保留原条件、阈值和优先级 |
| 一致性或不确定度明确 Fail，或已确认故障尚无对应 Pass 恢复证据 | 4 | Unknown 不得清除明确故障 |
| 不确定度 Unknown，或必需航向 Unknown 不满足有限容忍条件 | 4 | 启动、重定位预热、无效协方差仍停车 |
| 仅可选 CAN 平移检查 Unknown，必需航向和不确定度 Pass，无待恢复故障 | 继续检查雷达分支 | 不为 CAN 缺失设置停车计时器，观测保留 Unknown |
| 仅一致性 Unknown，且满足下述证据条件 | 继续检查雷达分支 | 原始观测仍为 Unknown，不伪造 Pass |
| 雷达故障或雷达证据未全部通过，且未启用雷达屏蔽 | 2 | 保留原策略 |
| 其余情况 | 0 | 不代表独立验证的整车安全结论 |

`diagnostics.consistency_unknown_grace_sec` 的代码默认值为 0，必需航向 Unknown 不使用旧证据。
2026-10-21 SANY 模板采用 0.15 s 作为部署调试候选：接收年龄超过 100 ms 仍将一致性判为 UNKNOWN，
但停车判断独立使用下述有限证据条件。150 ms 高于本次约 107 ms 的接收年龄，并在原 300 ms 输出急停之前截止；
是否可以容忍仍取决于最后一次 Pass 的年龄及全部条件，不能只比较接收年龄。
这是初始调试值，不是仅凭静止数据完成的行驶验收。非零值必须小于 `output_timeout_sec`。
CAN 可选修订后的计时基准是本定位 epoch 中**最后一次明确航向一致性 Pass 的单调时间**，
不是首次 Unknown 或定时器每次评估时间；Unknown 自身不能续期。
输出年龄、IMU 到达年龄和 IMU 与位姿时间戳差均不得超过该总上限；CAN 不参与这一必需证据计时；
不确定度当前 Pass、没有待恢复明确故障且没有 severe 时，才可使用旧的通过证据。
epoch 切换和非法 IMU 输入会取消旧航向证据；位姿源时间增量与单调到达时间增量相差超过 100 ms 也取消旧证据。
单纯源时间间隔超过 100 ms 不再独立等同于时钟跳变，但旧 Pass 仍会在总上限到期后失效。
明确 Fail 必须由该检查后续 Pass 清除。
300 ms 输出超时、定位门禁丢失和全雷达不可用均不受容忍参数影响。

运行日志新增 `confirmed_failure`、`unknown_tolerated`、
`consistency_pass_age_sec`、`output_age_sec`、CAN/IMU 到达年龄、
`output_stale` 和 `localization_lost`，区分证据缺失、明确失败与急停来源。
无证据的年龄为 infinity，不应解释为零。

2026-10-21 本地 WSL Ubuntu 22.04 / ROS Humble 构建及验证：
`localization_dtc_test`、`radar_dtc_state_test`、`thread_scheduling_test` 全部通过；
覆盖 107 ms 缺口、150 ms 总上限、延迟位姿不续期、持续缺失、时钟异常、明确故障保持，
以及三/四雷达的缺失、恢复和全失效。107 ms 用例及 5 m/s 输入为构造的边界测试，不是本次数据包的完整回放。
真实 `run_loc_online` 三/四雷达 smoke 均通过，ROS steady 平均 Safety 频率分别为 49.5/50.1 Hz，
并验证实际订阅集合、启动 UNKNOWN=4 和无有效输出后的 -9。
两个 smoke 的 header 时间间隔中位数折算均为 55.6 Hz；这是本地 WSL ROS 时钟观测，不能据此验收域控授时。
车辆配置解析测试 9 项、工作区接口测试 9 项通过，依赖源码与 `workspace.repos` 一致。
此验证不包含新域控 ARM64 部署、内核隔核或实际行驶验收。

2026-10-20 的 `sany_lidar_20261020_144934` 包标记提交为 `2b27e3a`。
与 `49daa6d` 对比，定位 DTC 算法未变化，设备诊断接入有变化。
`scripts/analyze_safety_timing.py <run-directory>` 可复现以下**日志证据**：
30 次 wall/steady 偏移突增约 3.584 s；33 段 -9 中 31 段开始于最近一次
跳变后的约 2.3–4.2 wall seconds；4 的 1917 条审计日志中 1898 条为
一致性 Unknown、不确定度 Pass。日志按变化加周期记录，计数不是故障时间占比。
wall 时间段含跳变，不可直接据其分位数标定容忍时间，也不能仅凭相关性断言 Chrony 是根因。
首段 -9 前有 stale lidar correction 约 3.1 s 和发布门禁关闭记录；
另外两段分别在 15:00:35 与 15:04:46，门禁匹配年龄约 0.50/0.60 s，
不能用同一时钟跳变解释所有急停。当前仍保持容忍参数为 0，等待稳定授时下的负载验证。

Orin 12 核部署候选使用 `config/sany/orin_cpu_roles.bash`：
CPU 10 为 `loc_sensor`，CPU 11 为 `loc_executor`/`loc_output`，
CPU 0–9 为地图匹配、OpenMP 工作线程、日志和其他后台线程。
进程先以 0–11 启动，初始化前将默认继承掩码缩为 0–9；
仅上述具名关键角色可使用启动时验证过的专用掩码。
LIO 的 OpenMP 每个工作线程（包括没有分到迭代的线程）在 parallel 区域中转至后台核，
主线程保留传感器核。`OMP_WAIT_POLICY=PASSIVE` 避免空闲团队持续自旋。
绑定失败或读回不一致会报错，不能静默视为成功。

`scripts/configure_orin_localization_isolation.py` 默认只预览，
`--apply` 需要 root，备份 extlinux 配置并仅修改 primary APPEND：
`isolcpus=domain,managed_irq,10-11 irqaffinity=0-9`；脚本不重启。
现场内核未启用 CONFIG_NO_HZ_FULL，因此不承诺无调度滴答或无中断。
回滚为恢复脚本输出的备份并重启；上线前还需保存旧启动入口及旧安装前缀。
重启后运行 `scripts/verify_localization_affinity.py <pid> --require-isolation`，
核对 /proc 中每条线程的实际 CPU 掩码及 /sys 中的隔离核。
该检查不能替代 IRQ 增量、实际运行 CPU、输入/发布延迟和实车负载测试。

当前候选验证记录：独立 C++ 核心回归、ARM CTest 两项、九项接口检查、
两项 extlinux 修改测试通过。ARM 真实节点在 localhost-only、ROS domain 172
验证了启动 Unknown、超时急停与约 50 Hz Safety 发布。
测试须显式检查 `ldd`：加载旧 install 库会造成结构 ABI/符号不匹配，
不能把该错误归因于判定算法，也不能把未核实库路径的测试视为候选验证。
现场最小候选只包含 `run_loc_online`、两个启动检查工具及同版内部库，
避免把旧离线可执行文件与新内部 ABI 混装；它不是完整离线工具发行包。

切换前 20 s 实收基线（编译期间、现场原 production 配置）收到 918 条 Safety，
全部为 0；位姿最大接收间隔约 0.074 s，设备身份/发布者/实收检查通过。
第一次探针没有加载现场 Fast DDS XML，未发现发布者，因此该次结果不是车辆无输出证据；
加载与运行进程相同的配置后才得到以上基线。
这些结果不能证明原问题已消失，也不能作为 CPU 隔离后的效果。

2026-10-20 现场候选启动补充：本地正式 CMake/CTest 两项均通过。
CAN 恢复过程中短暂出现旧 SpeThrCAN4（无 comm_header）与新布局不匹配；
随后现场切换到 wheel_loader/ROS2_CAN 桥接，实收 source_id=3、有效时间戳，
五路输入检查通过（CAN-IMU 时间差约 -13 ms）。不要将“进程已启动”视作接口兼容证据。
候选运行目录 safety_cpu_validation_can_20261020，启动 PID 167457。
逐线程亲和性读回通过，包括 LIO OpenMP 与 SOLID 工作线程；内核隔离列表仍为空，尚待重启。
启动后地图初始化尚未通过，日志 pose_enabled=0、output_age_sec=inf，当前 -9
属于没有有效定位输出，不能据此评价运行中调度改善；等待现场确认初始位姿流程。
容忍参数仍为 0，急停触发条件未放宽。
车辆换位置后，原候选已自行恢复定位（pose_enabled=1、Safety=0）。
随后应现场要求正常重启定位，运行目录 safety_cpu_new_position_20261020_1833，PID 180402。
五路输入检查通过，CAN-IMU 时间差约 -14 ms；启动 18:33:45，18:33:51.269 Safety 转为 0。
首段 25 s 实收 1252 条 Safety（142 条 4、1110 条 0），无 -9；
日志中的 4 出现在初始化完成前。收到 3303 条位姿，最大接收间隔约 55 ms。
重启后 24 条线程亲和性核验通过，内核 isolated 列表仍为空。
初始化后另采 40 s：2003 条 Safety 全部为 0；5900 条位姿，最大实收间隔约 92 ms。
Safety 与位姿各 1 个发布者。该结果仅代表当前短时运行，不能替代隔核重启及实车负载验证。


整机重启后验证（2026-10-20）：boot ID 105b42b2-61f5-4a99-803a-e78253fbb34a，
内核 isolated=10-11，默认 IRQ mask=3ff。运行 safety_cpu_postboot_20261020，PID 7058。
24 条线程通过 --require-isolation；30 s 采样关键角色的实际 CPU 与分配一致。
隔离核上仅观察到 arch_timer/IPI 增量，未观察到其他设备中断；不承诺无中断。
40 s 稳定实收：Safety 2003 条全部为 0，位姿 5171 条，最大接收间隔约 71 ms。
启动五路输入通过；后雷达诊断在该窗口没有实收，雷达屏蔽仍为用户已有 true。
日志仍有雷达诊断告警，不能将 Safety=0 解释成全雷达健康。
启动初期 PTP 等待 Chrony，随后门禁通过且三路服务 active，没有绕过门禁。
另发现 BottomLevel_X 两个订阅端为 lightning_interfaces/msg/VehiclePose，
定位发布端为 lightning/msg/VehiclePose；本轮不修改控制侧，控制端实收仍需单独验收。
证据在 localization_safety_stage_20261020/postboot-*.json（本地同名副本）。
内核隔核及线程绑定已验证，动态负载和非零 Unknown 容忍仍未验收。
