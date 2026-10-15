@page dtc_delivery DTC 与安全状态接口交付

接口依据：
https://sanygroup.feishu.cn/wiki/XoT4wYusciDuhKkuw5Fcfo1Bn7T
及 https://rs2flu7c17.work.sany.com.cn/sheets/shtk5ejY5B0puSUC3dbBtmCE0ec?sheet=1dfnNH 。

@section dtc_interfaces 接口与职责

雷达发布 /diagnostics/source_dtc/livox，包含 44 个单雷达 DTC 和
0x10F00001（至少两台不可用）、0x10F00002（四台不可用），共 46 条。
定位发布 /diagnostics/source_dtc/localization，包含
0x30811813（一致性异常）、0x40450198（不确定度超限）。
两者类型均为 diagnostic_monitor_interfaces/msg/DtcObservationArray，
完整快照按码排序，10 Hz；已观察到的状态改变尽快追加发布。
QoS 为 reliable、volatile、KeepLast(10)。
physical_condition：0 UNKNOWN，1 PASS，2 FAIL。
这里只报告内部物理判断；故障分级、确认/恢复消抖、历史、锁存和后处理由诊断节点负责。
现有 /localization/fault_status 留作兼容，统一诊断节点应只消费新接口，避免重复处理。

comm_header.seq 表示发送方本轮健康评估；evaluation_seq 表示该项检测所依据的
新证据/检测进度；evaluation_age_us 是最近实际评估的稳态时钟年龄。
从未完成检测时为 UINT64_MAX。boot_id 与同节点心跳共用启动标识。
接收方应同时检查来源、版本、完整码集合、启动标识、序号及证据年龄。
UNKNOWN、缺项、超时、重启不能当作 PASS，也不能清除已确认故障。
定位源 ID 默认 6，四台雷达 ID 为前 7、后 12、左 8、右 9；
新增聚合节点 ID 默认 13，部署前须与统一诊断节点登记，保持唯一。

定位另发布 /safety/throstle，复用既有 geosun_msgs/msg/Safety，字段 int8 runstate，50 Hz。
Safety.msg 定义保持不变，不新增兼容消息类型；header.stamp 为发布时的 ROS 时间，
timesync、errorcode、strategies、reserved 保持默认 0，header.frame_id 默认空。
Topic 拼写 throstle 保留；接收端订阅 Safety 并读取 runstate。
优先级 -9 > 4 > 2 > 0：定位失锁/有效输出超时或四雷达不可用给 -9；
定位 DTC FAIL 给 4；其他雷达故障给 2；全部检测有效且正常才给 0。
定位检测 UNKNOWN 给 4；仅雷达检测 UNKNOWN/过期给 2。
已恢复通信历史（四个角色码的 base+12）与同步历史（base+13）仍保留在源 DTC 中，
但新鲜的历史 FAIL 不单独参与当前 Safety 故障判定；历史 UNKNOWN、过期、其他当前 DTC
以及四台设备实际可用性仍须校验。因此历史事件本身不会让正常设备永久保持 runstate=2。
已接收的四雷达 FAIL 在重启 UNKNOWN 或通信中断时维持 -9，直至新鲜 PASS。
该兼容接口表达实时原始状态，不承担诊断中心的等级与恢复消抖策略。
兼容接口同时读取四台 diagnostic_state 的设备可用性，因此单台未映射 native fatal
或实际不可用也给 2；设备可用性未知时不能给 0，不新增未定义的系统 DTC。

@section dtc_lidar 雷达实现与部署

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
两台以上与四台全故障可同时成立，按四个物理角色去重统计。

master_lidar.launch.py 启动一次聚合节点，master_4lidar 间接复用它；
slave_lidar 只运行后雷达。主从部署须确保四个 diagnostic_state Topic 跨主机可达。
配置为 livox_ros_driver2_core/config/dtc_aggregator.yaml，驱动 dtc.* 参数
可设置角色、超时、有效点数和配置版本。默认点包/点云 500 ms、HMS 2500 ms、
内部状态 500 ms、启动等待 8 s，设备状态每 1 s 查询。
遗漏设备在启动宽限后按不可用处理，因此三雷达实车必须显式补齐后雷达输入才能正常。
主、从驱动与聚合节点须一起升级到 `livox-dtc-v2-sync-source` 并重新构建、重启后生效；
混合配置版本的设备状态会被拒绝，不得只更新主机就宣称后雷达完成修复。

2026-10-15 现场只读查询确认前、左、右三台 `work_state=1`、`lidar_diag_status=0`，
PTP 偏差约 -5 微秒，最后同步时刻连续前进；三台均保留 0x04050002，左雷达另有
0x04090002（GPS 信号异常）。原实现未区分同步源，并将历史事件参与 Safety，导致误报。
ROS 快照中四个角色的可用性均为 USABLE，后雷达点云时间戳也在前进。
现场接线描述与 ROS 图不一致时，应核对发布者、配置和持续前进的源时间戳，不能直接跳过后雷达检查。
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
排查 `runstate=4` 时同时查看 `LOCALIZATION_DTC` 日志及定位 DTC Topic：
`consistency=1 uncertainty=2 severe=0` 表示不确定度超限；条件 0 表示证据不足。
启动日志 `LOCALIZATION_DTC_CONFIG` 显示实际生效的版本、倍率和三个门限；
周期日志 `can_stamp_lag_sec` / `imu_stamp_lag_sec` 显示最后位姿时间减相应最新输入时间，
缺输入时为 NaN。持续的大偏差应检查时钟同步，不能通过增加运动容差消除。
`Safety.errorcode=0` 是兼容消息的默认值，不能据此判断没有 DTC。

@section dtc_calibration 标定、测试与验收

config/sany/localization.yaml 的 diagnostics 和 localization_pgo 为台架起始值，
需在联调后冻结 config_version。不得将建议阈值当作已验收的安全响应指标。

1. 采集静止、低速、倒车、弯道、颠簸、各坡度及地图边缘场景，
   同步独立 RTK/全站仪参考、输入时间戳、LOCALIZATION_DTC 日志及两个 DTC Topic。
2. 正常数据统计纵/横/航向残差分布，按速度、转弯、坡度分层标定基础容差和比例容差。
   时间不同步或参考自身无效的片段单独处理，不能用于放宽容差。
3. 用独立参考误差标定协方差方差倍率、下限及过程噪声，检查二维椭圆和航向覆盖率；
   窗口相关性及静止保持偏差必须纳入。推荐以 95% 覆盖率为起点，项目另行确定验收目标。
4. 注入单设备 HMS 多原因及清除、未知 fatal、断网、空点云、冻结时间戳、
   两台/四台失效、驱动退出、内部消息重放、进程重启及乱序。
   验证 44+2 条完整性、同设备不重复计数、UNKNOWN 不解除故障。
5. 注入位姿/航向跳变、轮速/IMU 停止或冻结、协方差超限/无效、重定位及定位输出停止。
   验证 0/2/4/-9 优先级、故障发生与恢复时间、50 Hz 输出及接收端消息兼容。
6. 用稳态时钟测量故障注入到源 DTC、兼容状态和诊断中心安全动作的各段延迟。
   默认上界预算：设备 HMS 检出约 1 s 查询 + 100 ms 评估 + 网络；
   断流约 500 ms + 100 ms 评估 + 网络；定位有效输出超时 300 ms + 20 ms 状态周期；
   运动窗口 500 ms，突跳可由相邻增量提前发现。
   这些是配置预算，最终最大响应时间必须在目标硬件负载及网络环境中实测验收。

自动回归：lidar_dtc_monitor_test、localization_dtc_test、localization_pgo_test。
PGO 回归直接将实际优化器的后验送入 DTC，覆盖五帧窗口边缘化后的健康轨迹；
定位回归另外验证 0.31 m 横向、2.1° 航向超限仍为 FAIL。
test/dtc_ros_integration.py 验证实际聚合节点的 ROS 输入、恢复、重放和断流行为。

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

ROS 2 Humble / Ubuntu 22.04 下已完成消息包、驱动、定位程序及相关测试构建。
八项相关 CTest 通过：雷达 DTC、定位 DTC、轮速消息、雷达跟踪状态、后轴位姿、
三一定位输出、PGO 和定位输入并发回归。
两项实际节点 ROS 验证使用 Cyclone DDS：
test/dtc_ros_integration.py 通过完整快照、两台/四台失效、恢复、重放、HMS 过期和断流；
src/test/dtc_online_smoke.py 启动真实 run_loc_online，验证定位 DTC 元数据和完整性，
Safety.runstate 约 50 Hz，启动 UNKNOWN 为 4，超时后为 -9；消息类型和时间戳一并校验。
频率测试使用与 C++ WallTimer 相同的 ROS 稳态时钟；header 为 ROS 发布时刻。
本机 WSL 的校正时钟与原始稳态时钟曾出现约 10% 的速率差，header 频率仅作为诊断值。
ROS 稳态时钟实现参考 https://github.com/ros2/rcutils/blob/humble/src/time_unix.c 。
后者使用无传感器的临时空地图，故正常运行的 0/2/4/-9 组合由定位 DTC 单元测试覆盖。

在 lightning_lm_ws 根目录完成构建后，可复现节点验证：
@verbatim
source /opt/ros/humble/setup.bash
source install/setup.bash
RMW_IMPLEMENTATION=rmw_cyclonedds_cpp ROS_DOMAIN_ID=173 python3 src/livox_ros_driver2/test/dtc_ros_integration.py
RMW_IMPLEMENTATION=rmw_cyclonedds_cpp ROS_DOMAIN_ID=172 python3 src/lightning-lm/src/test/dtc_online_smoke.py build/lightning_lm/bin/run_loc_online
@endverbatim

部署时须同步 common_msgs 的三个更新消息包、驱动及定位程序，并重新构建消费者。
使用一致版本的安装目录，检查 lightning_lm/lib 内历史复制的 ROSIDL 库是否遮蔽新消息包。
目标硬件的 Livox 状态查询、真实传感器数据、跨主机 DDS 及故障响应时间仍按上述流程联调验收。
