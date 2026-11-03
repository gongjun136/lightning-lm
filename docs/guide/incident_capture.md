@page incident_capture_guide 进程内故障采集与快速排查

SANY 配置默认开启 `diagnostics.incident_capture`，独立于全程录包开关 `SANY_RECORD_BAG`。采集器复用现有回调中的点云和主 IMU 引用，不新增 DDS 订阅或进程；独立线程约每 100 ms 检测状态，写线程使用 nice 10。采集器变更需要重新编译并重启定位后生效；现场验收应核对实际运行版本和采集结果。

## 开启、关闭与路径

在实际使用的定位 YAML 中设置：

```yaml
diagnostics:
  incident_capture:
    enabled: true  # 改为 false 关闭
```

SANY 启动脚本通过 `--incident_dir` 指向本次运行的 `$run_dir/results/incidents`，并在 `run_metadata.txt` 记录该路径；检查本次运行的配置快照和启动日志，确认设置与路径。手工启动时可显式传入绝对路径 `--incident_dir=/绝对路径/incidents`；未提供该参数时，优先回退到 `--output_published_tum` 父目录下的 `incidents`，其次使用 `--output_tum` 的父目录；如果都未提供，日志会提示禁用采集，算法继续运行。事件各有目录；全程 bag 与事件采集分别开关。使用此采集器时可设置 `SANY_RECORD_BAG=0`、`SANY_ENABLE_POSE_VEL_WATCHDOG=0`，避免同时启用全程录包和独立看门狗；原始数据默认缓存且不增加重复 DDS 订阅。

标准启动器默认把本次匹配失败快照放在 `$run_dir/match_failures`，把事件放在 `$run_dir/results/incidents`。两类目录均在默认位置时，需要交付取证数据可在确认保存完成后复制整个 `$run_dir`，一并携带这两类目录、日志和运行元数据；仅复制单个事件目录不包含整个运行及匹配快照。复制后元数据中仍可能保留原机器上的绝对路径，应对照复制后的运行目录检查文件。旧版可能使用 `$out_root/match_failures` 与 `$out_root/incidents`，诊断保留对旧路径的兼容，路径变更不会迁移或删除历史文件；旧运行或手工指定运行目录外路径时，应按元数据和实际位置另行携带所需目录。

本文的 `$run_dir` 表示本次运行目录，不一定是当前 shell 中已设置的变量。可从启动或正常结束时输出的运行目录路径查找，也可查看 `--diagnose --json` 输出的 `resolved_run`；实际操作时用该次运行的绝对路径替换，注意多个实例时核对选中的进程和时间。

## 现场单命令分析

现场部署入口 `~/Sany/run_loc.sh` 的 `lidar` 分支须将参数原样传给仓库内的 `scripts/run_sany_lidar_loc.sh`；确认入口已按此更新后执行：

```bash
bash ~/Sany/run_loc.sh --diagnose
bash ~/Sany/run_loc.sh --diagnose /路径/事件目录
bash ~/Sany/run_loc.sh --diagnose --all --color always
bash ~/Sany/run_loc.sh --diagnose --json > diagnosis.json
```

开发联调也可从仓库根目录直接执行 `bash scripts/run_sany_lidar_loc.sh --diagnose [path]`，并附加相同选项。仓库入口在启动定位、录包或创建运行目录之前转入 `scripts/diagnose_localization.sh`；独立 shell 和 Python 脚本仍可直接使用。支持 `--pid`、`--root`、`--log-mb`、`--bag-mb`、`--all` 和 `--color auto|always|never`。默认文本显示最近 12 段及按来源分类的全部汇总；同一真实故障可能同时出现在日志段和事件时间线段，不能将这些统计相加当作事件数。仅有入口日志的新运行也优先于旧事件检查；报告区分 runtime_loader 缺库、map_preflight 拒绝和 input_preflight 超时，若运行目录创建之前未留下日志则不能推断原因。诊断开头优先给出结论和简短问题链，同时保留最近明确航向一致性 FAIL 的段数、时间、残差及恢复情况，避免较早失锁掩盖后续状态 4；`--all` 展开所有已扫描判据，JSON 包含全部已分析段，日志读取仍受预算限制，不能保证覆盖全量日志。`auto` 遵守 TTY、`NO_COLOR` 和 `TERM=dumb`；JSON 始终无 ANSI。红色表示最新记录故障，黄色表示历史或缺失证据，绿色仅表示最新记录为 0，青色为信息。

`path` 可以是事件目录、匹配失败快照目录、`summary.json`、`timeline.csv`、运行目录、归档根目录、裸 MCAP 或 db3 文件，也可以是 rosbag 目录；不支持压缩归档。显式指定运行目录时检查该次运行日志，仍受 `--log-mb` 预算限制；显式指定事件目录时按该事件时间窗检查。无进程自动查找仍选择最新包。例如：

```bash
bash ~/Sany/run_loc.sh --diagnose /路径/事件目录
```

自动检查当前定位时，会调用雷达诊断的 `--json --no-refresh` 汇总当前上游证据；显式指定历史定位包时不拼入当前雷达状态。雷达独立检查可使用 `bash ~/Sany/run_lidar.sh --diagnose`，其退出码 0 只表示有证据，并不表示健康。

不带路径时，脚本通过 `/proc/<pid>/exe` 识别 `run_loc_online`，不依赖可能显示为 `loc_executor` 的线程名。多个实例时选择启动最晚的实例，并提示可用 `--pid` 指定；通过 `run_metadata` 和启动时间关联本次运行，避免 PID 复用或误选其他运行的事件。本次运行尚无事件时分析本次运行日志。无定位进程时，在 `LIGHTNING_LM_OUT_ROOT`、`SANY_WS/runs` 等默认路径寻找最新包；前者可来自环境或 `~/.config/sany/workspace.bash` 中的导出值。`--root` 仅覆盖无进程时的搜索根目录。

脚本默认用标准库读取摘要、时间线和运行日志尾部 32 MiB，另读取日志头部 1 MiB 以获取配置；`--log-mb` 可调整尾部范围，报告会标明覆盖范围。新版本时间线包含细分 Safety 判据；旧包需要同次运行日志补充判据证据；日志无法恢复未记录的物理根因。指定事件只关联其时间窗内的日志；自动检查正在运行的实例时，还会显示当前运行日志尾部的历史故障并标明来源时间。

裸 bag 缺少摘要和时间线时，脚本通过 `rosbag2_py` 离线筛选 Safety，不创建 ROS 节点、不解点云；默认文件预算 128 MiB，可用 `--bag-mb` 调整，读取还受 15 秒和 30000 次状态变化上限约束。只有上报码而缺少原始细分判据时，报告会说明证据不足。脚本不订阅 DDS、不启动定位、不修改阈值或输入数据；文本或 `--json` 报告打印到标准输出，可自行重定向保存。选中完整匹配快照时默认另在快照包中写入俯视图和分析文件；`--no-match-image` 只跳过画图，仍读取和显示可用的快照摘要。

脚本及 ROS 环境准备辅助进程以 nice 10 运行，从自身允许的 CPU 集合中排除系统隔离核；若没有可用的非隔离核则拒绝分析。退出码 `0` 表示读到了可分析的 Safety 记录，不能理解为没有故障；无证据或选择失败返回 `2`。检查报告中的事件时间、原始及发布故障码、具体字段及阈值或状态、证据位置和缺失证据。判据触发原因不能直接当作已证实的物理根因，证据不足时应继续调查。

记录按来源、设备和判据时段保留，区分启动未输出与运行超时，并保留首次异常、较坏样例及附近重定位、CAN、IMU 日志。Safety 日志 sequence 不连续时报告会提示；未知的新故障码会提示规则缺失，不按健康处理。状态优先依据最新时间戳；旧记录或超过 5 秒的记录不能证明当前健康。摘要和时间线范围会与故障起始时刻交叉检查；未找到覆盖时提示可能仍在导出、被抑制或已轮转，范围重叠也不保证逐条完整。事件摘要最多 256 条，采集统计区分累计计数和窗口丢弃增量。

关联日志受同源 Safety 恢复为 0 和观测缺口的边界约束，避免相邻故障时段串联；保留最近 1000 条普通日志及最早 128 条关键边界日志，去重并在超限时明确提示；保留首次局部匹配失败和禁用输出证据，避免长期重定位日志刷掉故障起点。独立 runtime 保护统计覆盖扫描范围或指定事件范围内的过期校正丢弃和自适应降载，不能作为 Safety 故障次数。匹配失败日志记录地图点数、XY 范围内点的比例、重力余弦及其阈值（`map_points`、`map_inside_xy`、`map_gravity_cos`、`map_min_inside_xy`、`map_min_overlap`、`map_min_gravity_cos`）；解析优先使用日志中的实际阈值，旧包配置值仅作线索。`map_checked=1` 或 overlap 不足不能证明车辆驶出地图。航向检查比较位姿航向增量与 IMU 机体系 z 角速度积分，CAN 参与平移检查；日志缺少积分累计转角或窗口分支时，不能重建确切动态阈值，也不能判定物理根因。首次输出或 epoch 重建后尚无 yaw PASS 时，应检查窗口建立与历史 IMU 覆盖，不能仅凭 UNKNOWN 区分这两种情况。启动阶段无有效输出时，优先核查初始化；摘要自动关联首次 -9 起点附近的 sensor queue lag 与过期 LIO 校正证据，区分首次 -9、门控关闭与解除、原始码恢复 0，并包含状态 2 的同步方位与恢复；lag 是最新入队与最新处理传感器时间之差，不能直接归因于 CPU 或 DDS。并存雷达条件仅作为对应历史 -9 记录的旁证，不虚报状态 2，也不能据此证明当前仍异常或认定故障因果。事件覆盖故障起点不等于覆盖整段故障。上述规则的现场验证以实际复验结果为准。

诊断按日志的 `policy_version` 区分 v8 与历史 v7、v6、v5。v8 和 v7 的 `raw=0` 同时运动一致性 UNKNOWN 时，文本单独显示“运动证据提示”，JSON 的 `motion_evidence` 保留 `status` 和 `unknown_advisory`；该提示不计为故障段，也不把 UNKNOWN 描述为 PASS。v8 或 v7 的 `raw=4` 若只有运动 UNKNOWN 而缺少其他支持判据，应报告判据缺失，不能声称 UNKNOWN 触发了 4；还需核查初始化、此前 FAIL 是否恢复、协方差、输出及匹配新鲜度与门控等证据。历史 v5/v6 日志继续按各自策略解释，不以 v8/v7 的规则倒推旧记录。

已初始化后若 `output_stale=1`、`last_published_pose=0` 且 `output_age_sec=inf`，`--diagnose` 报告尚未完成首次有效位姿发布。对照事件时间线中的 `prepared_age` 和 `publish_begin_age`：前者是距结果进入发布校验的时间，不证明结果有效；后者为 `-1` 表示采集器尚未观测到业务发布开始，并不能单独确定原因。匹配成功、雷达组帧正常或 TUM 为空都不能单独证明发布状态；TUM 还受配置、缓冲刷新和采集完整性影响。航向 UNKNOWN 本身不单独触发 4。

一种可能的代码时序是静止 DR 请求先到，随后地图结果在下一条 DR 前设置静止标志；旧版可能因标志已置位而跳过静止输出缓存初始化，导致返回无效结果。修复后，在静止且缓存无效时，下一条被接受的 DR 从有效融合位姿补全缓存，结果仍须通过有效性、时间戳递增及匹配质量门控。`event=LOC_OUTPUT_RESULT_INVALID` 表示结果无效或时间戳非正而未进入业务发布，每 200 次节流记录一次；字段 `valid`、`timestamp`、`is_parking`、`pose_epoch`、`covariance_valid` 描述该结果。`event=PGO_STATIC_HOLD_SNAPSHOT_RECOVERED` 表示静止标志已置位但缓存缺失时补全了缓存，字段 `dr_stamp`、`fused_stamp`、`pose_epoch`、`covariance_valid` 描述补全时刻和来源；时间戳单位为秒。缺少事件可能是运行版本、日志范围或节流所致；补全事件本身不证明业务发布、下游接收或 Safety 已恢复，仍需核对 `publish_begin_age`、`publish_end_age` 与 Safety 记录。

报告对照运行库的 `LOCALIZATION_DIAGNOSTIC_BUILD` 标记与脚本合同；缺失或不匹配时为未验证。标记包含合同 ID、source revision、dirty 状态和构建时间；`run_metadata` 的 `git_commit` 不能证明运行二进制版本。

## IMU、雷达与 CAN 输入时序

定位进程对主 IMU、每路已配置雷达以及启用轮速观测时的 CAN 分别统计相邻消息的源时间间隔和订阅回调的单调到达间隔。源时间取消息头时间戳并转为 double 统计，日志保留 17 位有效数字，不保证保留原始硬件时间戳的纳秒精度；到达间隔使用回调开始时采样的单调时钟，不受墙钟校时影响。CAN 消息无法解码时仍记录回调到达，源时间标为无效，不将其误认为回调静默。

统计复用已有回调和健康定时器：每条消息在现有输入统计锁内做 O(1) 更新，每路使用固定 16 条异常记录缓存；健康定时器最多每秒排出一次记录，正常运行约每 10 秒输出一次累计摘要。正常退出时只排空已有时序证据，不将停止订阅后的收尾等待生成新的 `callback_silence`。缓存满时记录丢弃数，不保证每个异常都留存。该功能不新增运行时监控进程、线程、DDS 订阅或 ROS 消息字段，不改变 Safety 2、4、-9 的判定策略；原有故障采集仍独立保留。

日志标记为 `INPUT_TIMING`，正文为“输入时序观测”，当前 `schema=1`。各行按 `stream=imu|lidar|can` 和 `topic` 区分输入；`kind=start|interval|silence|summary` 分别表示首个样本、异常相邻间隔、回调静默及累计摘要。常用字段如下，时间单位均为秒：

| 字段 | 含义 |
|---|---|
| `sample_seq` | 该路从进程启动以来已观测的回调样本序号；摘要为 0 表示尚未观测到该路样本。正常样本不逐帧输出，不能以该序号的间隔判断日志丢弃。 |
| `event_seq` | 该路首样本、异常间隔和静默等事件的单调编号；事件缓存满时也递增，`summary` 仅记录当前编号而不递增。 |
| `source_stamp_sec`、`previous_source_stamp_sec`、`source_gap_sec` | 在 `interval` 中表示当前和上一条源时间戳及其间隔；`silence` 中的 `source_stamp_sec` 是最近收到的源时间戳。上一条或当前源时间无效时不能计算间隔；回退时的间隔为负，需结合原始时间戳解释。 |
| `arrival_gap_sec`、`arrival_steady_sec`、`emitted_steady_sec` | 回调间隔、记录的单调时钟时刻和批量输出时刻。静默记录的间隔是采样时距上次回调的时间；单调时钟只用于本进程内比较，不能当作墙钟或跨进程直接比较。 |
| `threshold_sec` | 时序取证阈值：IMU 和 CAN 为 0.1 秒，每路雷达为 0.2 秒；间隔严格大于该值才计为缺口。这不是 Safety 阈值，也不等于已证实丢包。 |
| `source_gaps`、`arrival_gaps`、`duplicates`、`rollbacks`、`invalid` | 该路进程生命周期内源间隔超阈、到达间隔超阈、源时间重复、回退和无效的累计次数；对应最大间隔见 `max_source_gap_sec`、`max_arrival_gap_sec`。 |
| `events_lost` | 该路原生异常记录缓存满导致的累计丢弃数。日志本身也可能因异步日志丢弃或读取预算而不完整，不能以该值为 0 证明证据完整。 |

现场可用 `bash ~/Sany/run_loc.sh --diagnose` 查看“上游输入时序”，JSON 输出对应 `input_timing`。报告对照相邻事件记录和累计摘要的 `event_seq` 识别扫描范围内的未覆盖事件，`event_sequence_gaps` 表示观察到的事件序号缺口；`events_lost` 表示原生缓存丢弃，两者可能重叠，不能相加当作丢弃总数。报告关联相应时间窗口的 Safety 段，但不将时序记录增加为 Safety 故障次数；回调墙钟时间由日志墙钟与单调时钟间隔估算，批量排出、调度和日志写入可能带来偏差，不应当作精密授时。具体分类及含义如下：

| 分类 | 观察到的条件和解释 |
|---|---|
| `callback_silence` | 该路曾有样本，之后回调超过阈值未到达；只说明观察点暂未收到新样本。下一条消息尚未到达时，没有新的源时间间隔可比较。 |
| `delivery_delay` | 同一对消息的源间隔有效且未超过取证阈值，而回调间隔超限，表明有延迟交付线索，不能保证逐帧无缺样，也不能确定延迟环节。JSON 的 `source_within_threshold` 表示源间隔是否在阈值内；`relative_delivery_change_sec` 为同一对消息的回调间隔减源间隔，仅在源间隔为正且有效时有值，表示相对节奏变化，不是绝对单向网络延迟。 |
| `source_gap` | 同一对消息的源间隔超过取证阈值，说明观察到超过该阈值的源时间缺口，并不能覆盖所有物理丢样。源端停发、传输或接收过程丢失及源时钟跳变都可能产生类似现象，不能直接认定硬件停发。 |
| `source_rollback`、`duplicate_source_stamp`、`invalid_source_stamp` | 源时间回退、重复或无效；需结合原始消息、解码状态和时钟状态进一步检查。 |

这些间隔不是绝对端到端时延，也不能单独区分硬件故障和网络丢弃。旧版运行库、缺少对应日志或受日志读取预算限制时，报告会缺少时序证据；没有证据不能证明输入正常或发生断流。新功能需重新编译并重启后采集，再核对运行库标记及实际现场结果。

## 定位日志与事件码

定位程序自身输出的自然语言日志正文使用中文，诊断使用稳定的 `event=` 大写事件码识别事件。机器字段、单位和枚举值，以及 ROS、DDS 和第三方库原始错误 detail 保留原文。`SAFETY_DIAGNOSTIC`、`INPUT_TIMING`、`LOCALIZATION_DIAGNOSTIC_BUILD`、`LOCALIZATION_DTC_CONFIG`、`GLOBAL_RELOCALIZATION[backend]`、`SOLID_KISS_REJECT`、`SOLID_KISS_ACCEPT`、`COMPUTE_BENCH_*` 和 `LIO_BENCH_FRAME` 等结构化标签保留；字段名、数值、时间戳和日志级别保持原有契约。长中文日志超过原有 8192 字节上限时，在完整 UTF-8 字符边界截断并标记 `async_truncated`，因此不能假定正文完整。

`/localization/debug_message` 的 `lightning/msg/DebugMessage.message` 人类可读正文也使用中文，覆盖节点就绪、首次定位、全局重定位开始/失败/接受及 LIO、ESKF、PGO、平滑器异常与恢复等事件。原有数值字段名、单位、枚举值及原始 `reason`/`detail` 保留，便于诊断；`header` 时间戳、`node_name`、ROS 话题和消息结构，以及内部限频键和状态键不变。上述日志的 `event=` 规则不新增消息字段：`DebugMessage` 仍为 `header`、`node_name` 和 `message`，正文中已有的 `event=` 文本仍可能出现。

常见事件码包括 `LOC_INITIAL_POSE_UNSET`、`LOC_MATCH_FAILED`、`LOC_INIT_MATCH_FAILED`、`LOC_OUTPUT_DISABLED`、`LOC_OUTPUT_STOPPED`、`LOC_OUTPUT_RECOVERED`、`LOC_SENSOR_QUEUE_LAG`、`LOC_SENSOR_QUEUE_RECOVERED`、`LIO_STALE_CORRECTION_DROPPED`、`LIO_ADAPTIVE_LOAD_CHANGED`、`LOC_RELOCALIZATION_ACCEPTED`、`LOC_TRACKING_LOST`、`LOC_CAN_SPEED_IGNORED`、`IMU_LONG_INTERVAL`、`SOLID_CLOUD_LOAD_FAILED`、`SOLID_TARGET_BUILD_FAILED` 和 `SOLID_QUERY_PREPARE_FAILED`。队列延迟使用 `lag_sec`，自适应负载级别变化使用 `from_step`（变化前级别）和 `to_step`（变化后级别）；`age`、`hard_deadline`、`processing_ms`、`latency_ms`、`tracking_healthy` 等原有字段仍按原名称读取。

例如，以下仅为格式示意，并非实际采集结果，也不能据此判断物理根因：

```text
event=LOC_MATCH_FAILED 激光定位匹配失败 score: 2.28, map_checked=1, map_overlap=0.72
```

历史英文日志没有事件码时，`scripts/diagnose_localization.sh` 调用的诊断逻辑继续按旧规则识别；日志存在显式大写事件码时，即使该码暂不被识别，也不再根据正文猜测事件。新增或修改源码日志时，应同步维护诊断解析与回归用例，并按下方开发同步检查审查合同。新增匹配取证字段不改变 Safety 消息、事件触发条件或时间线 CSV；匹配快照 schema 2 和诊断 JSON 可追加向后兼容的字段。现场应以实际运行版本和证据为准。

## 匹配失败点云证据

排查匹配拒绝与 Safety `4` 的关系时，先按时间对照 `SAFETY_DIAGNOSTIC`、`LOC_MATCH_FAILED`（跟踪）和 `LOC_INIT_MATCH_FAILED`（初始化）。失败日志中的 `rejection_reason` 和 `registration_reason` 记录实际判定项；诊断逐项解释 NDT、原始里程计差异、地图边界、重合率、重力一致性或前置证据缺失，并显示已记录的数值和阈值。跟踪日志还记录原始里程计差异的阈值 `raw_odom_limit`、传感器时间戳和连续拒绝帧数及阈值。两类日志的 `snapshot_recorder_present` 仅表示采集器对象是否存在，不能据此判断 YAML 开关或写线程是否成功启动；`snapshot_capture_admitted` 表示本帧是否成功入队。单帧拒绝本身不等于 Safety `4` 或确认失锁；判断是否相关仍需核对同一时段的 Safety 判据、匹配新鲜度和状态变化。日志中的判定原因也不等于已查明的物理根因。

基础配准日志及新匹配快照的 `registration` 记录 NDT 是否实际执行和收敛、迭代次数与上限、收敛步长阈值、迭代预算是否耗尽、变换与分数的数值有效性，以及 `confidence` 与 `min_confidence`。`ndt_budget_exhausted=1` 表示本次实际执行、尚未收敛且迭代次数达到上限；最后一次已收敛时不能仅凭达到上限认定预算耗尽。未执行 NDT 时快照的 `ndt_metrics_available` 为 `false`，新增的 NDT 收敛、数值有效性、迭代上限、epsilon 和耗尽状态字段为 `null`；保留的通用分数和迭代次数字段及日志中的默认布尔值不能作为 NDT 实测结果。旧日志和旧 schema 2 快照缺少字段时不能将缺项当作 `false`。诊断只依据实际记录解释，缺少必要字段时报告证据不足。

普通 NDT 基础配准接受要求 `hasConverged()` 为真、confidence 和变换均为有限值，且 confidence 不低于当前跟踪或初始化阈值。零增量或有限步长小于 `transformation_epsilon` 可以置收敛；仅看到迭代次数达到 20 次上限，不能判断是否收敛，最后一次也可能满足收敛条件。基础配准通过后，候选仍可能因相对与绝对里程计差异或地图一致性检查被拒绝；单帧拒绝也不等于已确认定位失锁。

例如，旧 13:48 包的五帧失败记录中 `iterations=20`，score 为 2.848–3.039，高于当时 tracking 阈值 0.4，且 `map_checked=0`；日志没有 `registration_reason`。这些记录表明分数阈值不是已知的拒绝原因，也未记录到地图一致性检查；不能仅凭 20 次迭代断言具体失败原因或物理根因。

普通 tracking 拒绝、启动 seed 和全局重定位候选均可采集接受前的扫描，以及当时可用的候选和局部地图，保存到本次运行的匹配快照目录。SANY 配置默认开启 `diagnostics.match_failure_capture`；启动器导出 `SANY_MATCH_FAILURE_DIR=$run_dir/match_failures` 并将原运行路径写入 `run_metadata.txt`。应用也支持 `--match_failure_dir`，未显式指定时使用输出 TUM 旁的目录。此采集器独立于 Safety 触发，不能保证每次 `4` 或 `-9` 都有对应的匹配快照。

schema 2 记录 phase、stage、rejection_reason、检索指标、候选以及原始、精修和实际检查位姿；只有实际可用的位姿有值，无候选为 null。tracking 包含 `summary.json`、配置、扫描及局部地图，可按三维最近邻原判据点复算分数。候选在地图范围外、局部地图为空时仍保存地图边界、扫描和实际位姿；单传感器坐标扫描图不伪叠加地图，也不能据此推定车辆在地图外。`run_loc.sh --diagnose [包路径]` 展示证据；支持 `--json` 和 `--no-match-image`。例如 `bash ~/Sany/run_loc.sh --diagnose /绝对路径/匹配失败包`；默认在包内 `birdseye/comparison.png`、`birdseye/index.html` 和 `birdseye/analysis.json` 输出图与分析。历史保存图的时间可能不同于最新 Safety，须按时间分别解释。红色未命中点不直接证明动态物体。

一次连续失败最多入队三包：尚无入队记录时可尝试记录失败；首包入队后，仅相对上一次入队记录出现阶段或原因变化、从无候选变为有候选，或候选显著改善才再次尝试。改善以本轮已入队记录的最佳有效重合率或检索分数为基准，分别要求增加至少 0.05 或 0.10；如果之前没有对应的有效值，首次取得有效值也可触发尝试。实际入队仍受空槽等条件限制。普通 tracking 的初始拒绝也计入本轮，连续成功 2 秒且相邻成功匹配通知间隔不超过 0.5 秒后才重置；停车复用旧位姿或少点时没有有效匹配，不算恢复，也不刷新匹配新鲜度。不会每 30 秒周期补包。旧 `cooldown_seconds` 键保留兼容，不再触发周期补包。其余轻量最新观测和状态继续更新。

日志中的入队成功只表示工作交给后台写入线程，不证明写盘成功。`MATCH_FAILURE_CAPTURE saved=` 日志证明当时完成保存；实际排查还应核对对应目录中状态为 complete 的 `summary.json` 及其列出的文件。已保存日志也不能保证该目录现在仍存在；当前未找到目录或上传包中缺少摘要，不能据此断言当时从未保存。诊断 JSON 对可关联的故障项提供 `matching_rejection`，顶层 `match_snapshots` 区分已读取的保存日志与当前文件可用性。复制整次运行目录可一并携带默认位置的两类采集数据；如果只取得事件子目录，需另查对应快照，未找到时应报告证据缺失，不能用事件包或轨迹代替。

采集不增加 DDS 订阅或独立进程；两个预分配槽最多保存 50000 个展示点及 50000 个原 global 检查点，后台线程避开隔离核 10/11、nice 10，限速 10 MiB/s。单包上限 64 MiB，匹配快照归档目录的默认预算为 1 GiB，磁盘预留 2 GiB，保护期至少 7 天；保护期内空间不足时拒绝新包并报告丢弃原因。新默认位置下该预算分别作用于每次运行的匹配快照目录，不是所有运行合计的全局上限。功能需要部署后的新采集，不能恢复已轮转的 2026-10-30 16:29 原始点云；测试和现场验收以实际结果为准。

## 航向一致性检查与证据（v8）

策略及配置版本为 `localization-dtc-v8-yaw-native-imu`。航向检查是定位之外的短窗口参考：用现有 IMU 到机体系的旋转转换主 IMU 角速度，再积分 z 分量，与位姿的航向增量比较；该旁路不扣除 ESKF 估计的陀螺零偏。定位发布路径本身使用三维 ESKF、逐输入时间步传播并扣除零偏，还受雷达、地图、PGO 和平滑修正影响；旁路 z 轴积分不是航向真值，与发布位姿之间出现差异并不自动说明定位错误。CAN 参与平移检查，不是航向参考。

v8 的航向参考按 IMU 原生样本的源时间戳分段，窗口两端落在样本之间时插值截取，逐段对角速度作梯形积分，不跳过中间样本；绝对转角按分段线性角速度曲线的绝对值面积计算，过零时分别计算两侧面积。平移参考仍沿用约 20 ms 中点采样积分。航向积分端点或原生源时间覆盖不足时，该次航向检查为 UNKNOWN；平移仍按 CAN 与 IMU 中点证据独立判定，中点证据不足时也可能为 UNKNOWN，明确平移 FAIL 在无更高优先级状态时仍可报 4。输入源间隔、最近有效回调的到达年龄和时间戳不推进等检查仍保留。SANY 的 `motion_window_sec=0.5`、`max_input_gap_sec=0.1` 是可配置值；运动检查使用后者约束相邻源样本间隔和最近有效回调的到达年龄，严格超过时证据不足；重复或回退的源时间戳不能刷新有效性。积分末端超过最新 IMU 样本时，允许在该间隔内保持末值，保持时长见下表。这个门限不应与上文输入时序取证使用的 0.1 秒阈值混为一谈。v7 的航向参考使用约 20 ms 中点采样，遇到高频角速度振荡时可能产生混叠，历史日志应按当时的版本解释，不能当作 v8 重算结果。

对相邻位姿及约 0.5 秒窗口分别检查发布位姿和未平滑位姿；未平滑位姿不可用时该分支回退为发布位姿，并非独立观测。残差为 `abs(wrap(位姿航向增量 - IMU 积分转角))`，其中 `wrap` 将角差归到 -180° 至 180°，用于处理航向跨越 ±180° 的情况。任一已评估分支的残差严格大于实际阈值时为 FAIL；相邻检查可在窗口建立前发现 FAIL。相邻检查没有 FAIL（包括 UNKNOWN）时由窗口结果决定，因此窗口 PASS 可给出航向 PASS；窗口尚未建立时为 UNKNOWN。之前的航向 FAIL 在后续整次航向检查明确 PASS 时才解除。阈值为 `yaw_tolerance_deg + angle_tolerance_ratio × 绝对转角`，角度项须使用相同单位。SANY 现场 YAML 的 `yaw_tolerance_deg` 为 3.0 度，`angle_tolerance_ratio` 为 0.1；通用 `Config` 默认仍为 2.0 度和 0.1。实际运行值以配置快照和启动日志为准，检查旧记录时使用其当时的配置。

运行日志的 `SAFETY_DIAGNOSTIC` 沿用约 1 Hz 及状态变化时的记录方式，追加下列航向证据，不保证记录每次位姿检查。角度字段单位为弧度，时间和时长字段单位为秒。

| 字段 | 含义 |
|---|---|
| `yaw_evidence_available`、`yaw_comparison` | 前者为 1 表示有可解释的航向比较证据；后者为 `adjacent` 或 `window`，表示按相邻／窗口规则选中的比较。证据可用不等于 PASS；前者为 0 时，其余默认值不能当作实测结果。 |
| `yaw_from_stamp`、`yaw_to_stamp` | 所选比较的两端位姿源时间戳，也是积分的起止时间。 |
| `yaw_gyro_delta_rad`、`yaw_gyro_turn_rad`、`yaw_limit_rad` | 分别为 IMU 有符号积分转角、累计绝对转角和包含比例项后的实际航向阈值。 |
| `yaw_published_delta_rad`、`yaw_raw_delta_rad` | 发布和未平滑位姿在所选比较中的航向增量。 |
| `yaw_published_error_rad`、`yaw_raw_error_rad` | 发布和未平滑位姿分别相对 IMU 参考的包角绝对残差；与 `yaw_limit_rad` 比较可识别超限分支。 |
| `yaw_raw_independent` | 为 0 时，两端中至少一端的未平滑位姿回退为发布位姿，不能将该分支视为独立证据；为 1 时两端均有未平滑位姿。 |
| `yaw_imu_segments` | 积分分段数，包含裁剪的原生样本区间及末端保持段（若有），不等于 IMU 消息数。 |
| `yaw_max_source_gap_sec`、`yaw_held_sec` | 前者为所用相邻原生样本的最大完整源时间间隔，即使区间被窗口端点截取也记录完整间隔；后者为使用最新角速度保持的时长。没有遇到相邻样本区间时前者可为 0，需结合分段数和保持时长解读。 |

这些字段对应被选中的相邻或窗口比较；旧字段 `error_yaw_rad` 仍为本次已检查的分支和窗口中的最大残差，可能不同于所选比较的残差，不能混用两者与动态阈值。航向 FAIL 且证据可用时，`run_loc.sh --diagnose` 解释实际记录的数值和超限分支；JSON 也保留对应字段。事件时间线未记录这些新增字段时，应同时检查同次运行日志；历史日志缺失字段不代表数值为零，也不能据此重建精确阈值或断定物理根因。排查时还须对照同一时段的输入时序、匹配、输出和协方差证据。

本次修改需重新编译并重启定位后采集，现场仍须检查运行库标记、实际配置、航向证据及触发和恢复情况。验收应覆盖高频振荡、源时间缺口和冻结、窗口边界及末端保持、两条位姿分支超限和 FAIL 后 PASS 恢复，并核对各自判据及实际上报码；旧运行不能用来证明新版本的行为或验收结果。

## 运动证据与雷达输入组帧策略（v8，兼容 v7）

v8 延续 v7 的运动 UNKNOWN 规则：运动一致性证据不足时，DTC 保持 UNKNOWN；UNKNOWN 本身不单独触发 Safety `runstate=4`。例如 IMU 150–208 ms 未更新可能使运动一致性证据不足，但须结合其他判据判断 Safety，不能由此推定定位已经通过运动一致性检查。此前的明确 FAIL 仍需相应检查明确 PASS 才能解除，UNKNOWN 不会解除；在无更高优先级状态时未恢复仍为 4。运动一致性明确 FAIL、协方差 FAIL/UNKNOWN、初始化未完成、输出或匹配过期及门控关闭仍按各自规则处理为 4 或 -9。输入源间隔和回调到达间隔继续记录，供排查使用。

`consistency_unknown_grace_sec` 配置和 `unknown_tolerated` 日志字段继续保留历史观测兼容，但 v8 和 v7 不依赖有限宽容期判断运动 UNKNOWN 是否单独停车。必须按同次运行的版本、配置和原始判据解释历史日志；较早的 v6 运动 UNKNOWN 宽容判据以及下节 v5 的设备 DTC 判据仍仅适用于各自版本。

v8 沿用 v7/v6 的原始雷达输入组帧策略，以雷达输入持续无法满足组帧条件作为 Safety `runstate=2` 的判据。定位 C++ 进程内现有的 20 ms Safety wall timer 使用 `steady_clock` 独立检查，不新增 Python 运行检测或服务。默认 `diagnostics.radar_frame_warning_sec=0.2` 秒为预警、`diagnostics.radar_frame_timeout_sec=0.3` 秒为故障；`diagnostics.radar_frame_recovery_frames=3` 要求连续三个新鲜合格输入组帧恢复。启动宽限复用 `diagnostics.startup_grace_sec`，默认 8 秒，从 `SetInitPose` 开始。上述阈值尚不能据此声称现场安全。

检查位于 `LocSystem` 原始雷达订阅入口，在预处理、限流和 LIO 排队之前，复用组帧器的时间戳、来源及最低来源数规则，不合并点云、不执行算法计算。仅接受非空原始消息，要求主雷达；完成新鲜输入组帧条件检查时更新进度，原始主雷达接收索引用于防止旧输入缓存续命。算法排队、处理或地图匹配滞后不参与该计时和恢复，不能作为报 2 的原因。`cloud_publish_min_lidars=3` 仅控制对外点云发布，不是报 2 的门槛。200 ms 为两个 100 ms 时间桶周期，现场同桶容差为 20 ms、乱序窗口为 100 ms；20 ms 约束雷达之间的相对组帧时间跨度，不保证绝对授时或与 IMU 对齐。诊断应结合 frame age、seen、started、warning、fault 和 recovery count 判断，并检查实际运行配置与版本。

优先级为：确认定位失锁先给 `-9`；雷达输入组帧超时及恢复确认期间给 `2`，优先于其后果性输出过期及无输入导致的启动过期；其余启动、质量或输出降级按各自判据给 `4` 或 `-9`；其他情况给 `0`。雷达 DTC 和独立心跳继续保留，严格 DTC 评估仅作对照；设备 DTC 不直接触发 2 或 -9。兼容接口 `errorcode` 默认仍为 0，业务状态看 `runstate`。部分雷达故障仍能满足实际组帧最低来源规则时不报 2；无法满足时可以报 2，这不等于全部雷达硬件损坏。历史 v5/v6 记录仍按原策略解析。

实施后须回归断流、旧缓存、时间戳不前进、来源不足、启动宽限、恢复及状态优先级，并在目标设备验证触发与恢复延迟、日志和实际发布值。本文不表示新策略已部署或通过实车验收。

启动边界：定位 Monitor 在 `Init` 时开始计时，帧看门狗在 `SetInitPose` 后启动。Monitor 已启动超时、但尚未见到合格雷达输入组帧且帧看门狗宽限未结束时保持 4；帧看门狗启动后满 8 秒仍无首帧则报 2。若已见过合格雷达输入组帧但初始定位始终未建立，保留原有启动超时 -9。

## 历史 v5：Safety 码与确认失锁

策略和配置版本为 `localization-dtc-v5-confirmed-loss`。输出超过 0.5 秒未更新时为 4；有效匹配超过 0.5 秒未更新且尚未确认失锁时也为 4，评估新鲜度同步采用 0.5 秒。启动后从未成功定位的前 8 秒为 4，严格超过 8 秒为 -9；曾成功定位后，一旦确认失锁或进入重定位立即为 -9，不享有启动宽限。单帧匹配拒绝不等于确认失锁：原有连续五帧、状态 FAIL 或已定位后进入 INITIALIZING 等确认条件仍适用。原算法匹配接受阈值不变，不可信位姿的输出门控仍关闭。

部分雷达错误时，只有新鲜有效匹配和持续输出、持续输出至少 `motion_window=0.5` 秒，且当前运动一致性的必需证据明确 PASS（合法的可选 CAN 除外）、不确定度 PASS 等条件均满足，才保留 raw 2 和 strict 证据而将 selected 设为 0；真实质量失败或 UNKNOWN 证据不足仍为 4，`all_required_lidar_failed` 仍为 -9。不要把 selected 0 当作雷达故障已消失。

## 现场三步

1. 找到本次 `$run_dir/results/incidents` 下对应时间的事件目录，先打开本地 `summary.html`。查看触发原因、首个观测异常阶段、原始和上报码、图表与丢失计数；旧运行按元数据和实际路径查找。
2. 查看 `summary.json` 和 `timeline.csv`，按时间关联输入回调与时间戳、队列、匹配和发布；对照运行配置。检查运行元数据中的 `match_failure_dir` 及对应时段的匹配失败包，结合 `saved=` 日志核对实际保存情况。新默认路径可通过复制整个 `$run_dir` 一并移交；诊断阶段是观测顺序，不等于硬件根因。
3. 需要复现时再检查或回放事件中的分段 MCAP，例如 `ros2 bag info /路径/segment_0`；事件包不含地图，复现需要对应地图。采集只包含实际启用雷达点云、主 IMU、开启融合时的 CAN、`pose_vel` 及 Safety、heartbeat、状态、流水线和参与雷达诊断等低带宽证据；不录重复 raw 位姿和不需要的 TF。是否足够回放取决于实际记录及丢弃情况。

## 触发与报告解读

原始 Safety 码变化到 2、4 或 -9 时均触发事件采集，不限于某种具体的 4/-9 原因；安全上报关闭仍有效。启动无输出、发布调用完成超时、输入回调断流或主 IMU/雷达 header 时间戳不推进的状态新出现时也触发。持续故障不会每 100 ms 重复创建事件。活跃且未关闭的事件可合并并延长，最长 60 秒；导出收尾或默认 10 秒冷却中的新触发会被抑制并计数，不保证每次触发都有独立目录。它与上述匹配失败快照的触发和保存条件分别判断。

报告每个 stream 的 `first_relative_sec`、`last_relative_sec` 和 `written` 描述相对触发时间的首末记录及写入数；`callback_age_at_trigger`、`stamp_age_at_trigger` 描述触发时回调和时间戳推进的年龄，`-1` 表示未观测。采集器源码的并发修复将已观测时间戳晚于采样时刻的微小负年龄钳位为 0，避免误判为未观测及误触发启动无输出；修复需重新编译并重启后生效，旧版本事件可能受此问题影响。`omitted_messages` 是导出容量限制省略的消息；`collector_totals` 分别记录 busy、budget、oversize 丢弃及 suppressed、writer_errors。`capture_loss_seen` 表示本进程生命周期中曾发生采集丢弃，不限于当前事件；有累计丢弃时状态保守标记 partial。complete 只表示事件导出结束，不保证每个话题的触发前 30 秒完整。

发布完成只证明本进程 `publish` 调用返回，不能证明下游 DDS 消费者收到。正常 `Finish` 会结束写入；如果后窗口尚未结束，会缩短窗口并标记 partial。`kill -9` 可能丢失未触发的内存历史和当前未关闭片段；已关闭的独立片段通常可读，不能承诺强杀后完整导出。

## 默认资源边界

目标窗口为触发前 30 秒、后 10 秒，容量不足时可能更短。默认 2 GiB 内存预算估算被保留的消息和记录，并非进程 RSS 上限；序列化单条消息、库缓冲和页缓存等另占资源。单条默认上限 64 MiB；单事件序列化硬上限 2 GiB。分段目标为 32 MiB 或约 2 秒关闭，单条消息可以超过分段目标。

事件归档目录的默认保留目标 `disk_mb` 为 8 GiB，最多保留 10 个事件，每个分段关闭时检查实际磁盘大小；最后分段、单条大消息及元数据可带来额外空间占用。新默认位置下限额分别作用于每次运行的事件目录，不是所有运行合计的全局上限。新段前预留 `max(segment_mb, 本条大小) + 1 MiB` 和 `min_free_mb`（默认 2 GiB），其他进程仍可能同时写盘，因此不是绝对空闲空间保证。默认平均写速限制 20 MiB/s，不压缩；nice 10 降低 CPU 调度优先级，不保证 IO 调度优先级。调整窗口或预算后需结合现场数据速率、内存和磁盘实际验证。

## 重新编译后的现场验收

现场验收应在重新编译并重启后执行：核对运行库编译标记与脚本合同一致，再分别检查有进程、无进程回退和指定事件包。核对 2、4、-9 的实际字段，以及关闭上报时原始故障仍可见；检查覆盖范围、缺口、导出和抑制提示，并确认分析进程避开 CPU 10、11。匹配拒绝场景还应核对跟踪及初始化日志中的实际判据、`--diagnose` 的解释、快照入队与后台保存结果，并分别验证有快照和缺失快照的报告。记录实际结果和版本，不能以旧运行结果代替新版本验收。

采集器导出收尾或冷却期抑制新触发的行为尚未修复。并发年龄修复需在编译并重启后生效；旧格式与缺少构建标记的报告应按未验证处理。

## 开发同步检查

修改诊断相关源码后，在仓库根目录执行：

```bash
python3 scripts/check_diagnostic_contract.py
python3 scripts/test_diagnose_localization.py
```

源码摘要覆盖 Monitor、雷达和拓扑、Safety 发布及初始化/丢失门控、采集器、diagnostics 配置段和 `workspace.repos`。源码变化时先更新并审查解析逻辑和回归用例，再执行 `python3 scripts/check_diagnostic_contract.py --update`；不能只更新摘要。CMake 配置和 CTest 也会检查合同。该合同不是所有算法行为或物理根因的形式化保证。
