@page causal_diagnostics 因果链与时间关联

# 主 Orin 定位诊断：资源、处理长尾与状态变化

本页说明如何联合读取资源概况、处理耗时与业务事件，分析定位长尾、速度变化和静止保持。三类证据的采集入口与输出如下：

| 证据 | 采集实现 | 输出文件 | 控制开关 |
|---|---|---|---|
| 资源概况 | `scripts/monitor_process_resources.py` 独立采样进程 | `results/process_resources.jsonl` | `LIGHTNING_LM_RESOURCE_PROFILE` |
| 计算与队列耗时 | C++ 中的逐帧计时与窗口汇总；独立有界队列后台写入 | `results/compute_profile.log` | `LIGHTNING_LM_COMPUTE_PROFILE` |
| 业务事件追踪 | C++ 调用方填写 `CausalEvent`，共享 `Sink` 异步写入 | `results/causal_trace.csv` | `LIGHTNING_LM_CAUSAL_TRACE` |

资源文件提供进程、线程及整机 CPU 背景，CSV 提供业务事件与部分耗时记录。详细计算采样直接写入独立文件，主 stderr 保留低频端到端健康摘要；该摘要依赖合格 tracking 样本，不是停流心跳。各采集开关可以独立覆盖；`LIGHTNING_LM_RESOURCE_INTERVAL_SEC` 只控制资源采样周期。资源字段的完整单位与统计口径见 @ref resource_diagnostics "资源与端到端延迟"，启动及文件收集流程见 @ref run_sany_lidar_loc.sh "运行脚本源码"。

## 部署范围

这里的诊断开关控制证据采集；CAN融合权重、速度步长门槛、静止检测条件、线程预算与保护阈值由各自配置决定。分析已有run时，需核对其代码版本与配置快照。当前输出坐标约定由所选地图包的 `map_manifest.yaml` 加载，入口为 @ref lightning::LocSystem::Init() "LocSystem::Init()"；配置归属见 @ref configuration "配置入口"。

采集代码更新后需重新编译C++并部署脚本。`scripts/run_sany_lidar_loc.sh` 在 `LIGHTNING_LM_RUN_MODE=diagnostic` 时默认开启详细计算、资源和因果采集，`production` 时默认关闭。显式环境变量覆盖优先：

```bash
export LIGHTNING_LM_RUN_MODE=diagnostic
export LIGHTNING_LM_CAUSAL_TRACE=1          # 业务事件追踪
unset LIGHTNING_LM_COMPUTE_PROFILE          # 清除显式覆盖，使用 diagnostic 默认 1
export LIGHTNING_LM_RESOURCE_PROFILE=1      # 进程/线程资源采样
export LIGHTNING_LM_RESOURCE_INTERVAL_SEC=1 # 资源采样默认1秒，范围0.5–60秒
```

运行脚本自动将异步追踪限制在本次run的 `results/causal_trace.csv`，清除继承的旧路径；不复用旧文件。直接调用二进制的离线回放可显式设置 `LIGHTNING_LM_CAUSAL_TRACE_PATH`，必须是不存在的新文件。没有设置路径时不创建后台线程或文件。

详细计算采样不节流，使用独立队列与本次 run 的新文件 `results/compute_profile.log`；每条记录带成功入队序号 `profile_seq`，不再从 stderr 复制。序号连续本身不能证明无丢失，须结合 `received`、`dropped`、footer 和完整性状态检查。`profile_seq`、主日志 `async_seq` 与因果 CSV 的 `sequence` 各自独立，不能横向比较；应按可用的 `frame_id`、时间戳和 TID 关联。文件写入在后台，但现有 `FormatTimingSample` 字符串格式化仍在调用线程执行，详细采样有开销。正常结束后 `scripts/compute_profile_report.py` 检查 header、footer、序号、计数、丢弃、I/O 错误及截断，并写 `results/compute_profile_status.json`。metadata 记录存储方式、文件和完整性；`--diagnose` 显示证据状态。旧二进制缺少独立文件时报告缺证，不自动从 stderr 提取。采集损失告警不改变算法原退出码，也不能直接判定算法故障。手动复核用 `python3 scripts/compute_profile_report.py --run /path/to/run`，会更新状态文件并输出 metadata 键，不修改原始 profile。正常退出排空已接收记录；SIGKILL 可能缺 footer，断电持久性不受保证，完整采集也不能保证足以定位所有根因。

直接运行安装了 `OnlineLogging` 的 `run_loc_online`（包括 `--bag` 回放）可设置 `LIGHTNING_LM_COMPUTE_PROFILE_PATH`，缺省为当前目录的 `compute_profile.log`；文件已存在时启动失败。运行脚本清除继承路径并使用本次 run 的新文件。计算采样与因果追踪独立；关闭因果追踪会失去对应时序关联证据。

## 资源概况：process_resources.jsonl

`scripts/monitor_process_resources.py` 在本次算法进程组中找到实际 `run_loc_online` PID，读取 Linux `/proc`；`sample_interval()` 生成区间统计，当前metadata的 `schema_version=2`。默认每1秒采样，CPU增量使用两次采样的实际单调时间间隔作为分母；PSS默认每30秒读取一次。

| JSON字段 | 可回答的问题与口径 |
|---|---|
| `cpu_core_equivalents`、`threads[].cpu_core_equivalents` | 定位进程及各线程在采样区间消耗多少CPU；3.2表示平均消耗约3.2个逻辑核的计算时间。 |
| `memory`、`io_delta` | 定位进程的内存占用与读写计数增量；内存单位为kB，I/O计数需按字段区分字节与调用次数。 |
| `threads[].sched_wait_ns_delta`、上下文切换增量 | 线程等待CPU调度及切换情况；调度等待与业务队列、互斥锁等待的口径不同。 |
| `threads[].wchan` | 读取瞬间线程位于哪个内核等待通道；与线程CPU、切换次数、schedstat联合检查，不是区间内完整等待栈。 |
| `host_cpus` | 整机各逻辑核的busy/iowait背景，与定位进程占用分别读取。 |
| `system_process_cpu.top` | 每个采样区间CPU最多的20个已配对进程：PID、start_ticks、进程名、进程组、CPU核等效量。用于检查其他进程负载与定位长尾的时间关联。 |

`system_process_cpu.top` 和 `threads[].wchan` 都是资源JSON中的字段。采样器本身的耗时、缺失计数及首个区间的未配对值也需检查；详细定义见 @ref resource_diagnostics "资源与端到端延迟"。

## 事件追踪：causal_trace.csv

### 记录结构

文件第一行为CSV列名，后续每行是一条事件，正常退出时附加以 `#` 开头的计数行。一个传感器处理周期可能产生多个事件。

@ref lightning::profiling::CausalEvent "CausalEvent" 定义可复用的业务字段；具体 `kind`、`source`、`reason` 以及实际填写哪些字段，由业务调用点决定。事件类型目前采用字符串。@ref lightning::profiling::RecordCausalEvent() "RecordCausalEvent()" 经共享 `Sink` 补充记录时间与线程ID，再交给 @ref lightning::BoundedAsyncQueue "有界异步队列" 分配成功入队序号并写出。

| 公共列 | 含义 |
|---|---|
| `sequence` | 从0开始的成功入队序号；跨线程业务事件发生顺序可能与入队顺序不同。丢弃发生在分配序号之前，因此序号连续也不能证明没有丢记录。 |
| `wall_ns` | 调用线程记录时的系统时钟时间，单位ns，用于关联外部日志；系统校时可能使其跳变。 |
| `steady_ns` | 调用线程记录时的单调时钟时间，单位ns，用于比较同一运行中的记录时间间隔。 |
| `tid` | 产生事件的Linux线程ID，用于结合资源采样中的线程记录。 |

这两个本机时刻是在调用 `RecordCausalEvent()` 时采样，并非所有操作的开始时刻；例如 `sensor_process` 在消费回调完成后记录。

业务列按用途阅读；不是每种事件都会填写全部列：

| 列 | 含义与单位 |
|---|---|
| `kind`、`source`、`reason` | 事件类型、产生位置、接受/拒绝或状态变化原因；具体取值由调用点决定。 |
| `frame_id` | 调用点指定的关联帧标识，默认0；不是 `sequence`。 |
| `stamp`、`input_stamp`、`prior_stamp` | 关联状态或传感器、输入、此前状态的时间戳，单位s；具体含义随事件变化。 |
| `before_v`、`after_v`、`measured_v`、`innovation`、`stddev` | 更新前后、量测速度、速度创新及观测标准差，单位m/s；速度的具体坐标和含义由调用点决定。 |
| `nis`、`torque` | 归一化新息平方（无量纲）及关联电机扭矩；CAN调用点的扭矩单位为N·m。 |
| `duration_ms`、`queue_wait_ms`、`thread_cpu_ms` | 对应阶段墙钟经过时间、入队前采样到消费入口的等待、调用线程CPU时间，单位ms。 |
| `position_delta`、`yaw_delta` | 前后位置变化幅度（m）和航向变化量（rad）。 |
| `accepted`、`hold` | 接受状态或静止保持状态：1为是，0为否，−1为未记录或不适用。 |

字段默认值及各调用点的实际赋值见 @ref lightning::profiling::CausalEvent "字段定义" 和下表；CSV列顺序与数值格式以 `src/utils/causal_trace.cc` 中的表头和 `Sink::WriteBatch()` 为准。字符串直接写入CSV，不进行引号转义，调用方应保证其不含逗号或换行。

### 各事件的业务历元

下表中的时间戳单位均为秒；“未填”表示默认 `NaN`。这些业务历元不能与本机 `wall_ns` 或 `steady_ns` 直接相减。

| `kind` | `stamp` | `input_stamp` | `prior_stamp` |
|---|---|---|---|
| `can_input` | CAN样本Header历元 | 未填 | 未填 |
| `can_update` | 目标滤波器状态历元 | 所选CAN样本历元 | 未填 |
| `lidar_update` | 更新后的激光帧末状态历元 | 未填 | 更新前预测状态历元 |
| `hf_rebuild` | 重建后的高频状态历元 | 作为重建起点的LIO状态历元 | 重建前高频状态历元 |
| `static_decision` | 决策使用的IMU样本历元 | 所选CAN样本历元（有样本时） | 未填 |
| `static_transition` | 保持切换后的高频状态历元 | 未填 | 切换前状态历元 |
| `sensor_process` | IMU样本或点云Header历元 | 未填 | 未填 |
| `lock_wait` | 调用方提供的关联传感器或状态历元 | 未填 | 未填 |

CSV没有导出保证唯一的CAN消息标识。`input_stamp` 有助于关联输入和重放，但仅按时间戳去重不能保证得到独立CAN消息数；追踪也可能丢事件。需要精确输入数时应查原始输入证据。

### 事件类型与用途

| `kind` | 可回答的问题 |
|---|---|
| `can_input` | 实际收到的CAN换算速度、扭矩、Header历元及本机接收追踪时刻；即使关闭速度融合仍记录。 |
| `can_update` | 实际送入ESKF的CAN样本时刻、滤波器历元、更新前后纵向速度、创新、NIS、观测标准差、接受/拒绝原因。source区分LIO预测中的 `lio_predict`、实时高频IMU预测中的 `imu_predict` 和高频重建重放中的 `imu_replay`。 |
| `lidar_update` | 同一LIO观测更新前后纵向速度、位置变化范数、yaw变化（rad）、是否接受及frame_id。 |
| `hf_rebuild` | 从LIO状态重建高频状态的前后速度/历元、位置与yaw变化；包含内部预测和CAN重放的总变化，不可与内部事件重复相加。 |
| `static_decision` | 静止进入/退出所选的CAN样本时刻和速度；退出reason是最后一个防抖样本触发的CAN/LIO/IMU条件组合，不代表三个样本都由同一条件触发。 |
| `static_transition` | 高频滤波器实际保持切换前后速度、状态时刻。进入时归零，退出不直接恢复历史速度。 |
| `sensor_process` | 每个LiDAR消息及慢IMU消息的Header、frame_id、入队等待、消费回调墙钟/调用线程CPU。IMU详细记录条件为排队≥5ms或处理≥20ms，不能用这些详细记录计算全量IMU分位数；原全量窗口摘要仍保留。 |
| `lock_wait` | IMU/LiDAR入口的生命周期锁与处理锁等待≥0.2ms的记录；未涵盖所有第三方库或OpenMP内部锁。 |

事件调用点可从 @ref lightning::loc::Localization::ProcessSensorInput() "传感器处理"、@ref lightning::loc::Localization::ProcessWheelSpeed() "CAN输入" 和 @ref lightning::LaserMapping::ApplyWheelSpeedObservation() "CAN融合" 继续核对。

不同事件用 `stamp/input_stamp/prior_stamp` 区分输出状态、所选输入和更新前状态历元，单位为s；它们与本机记录时刻属于不同时间口径。浮点字段默认以NaN表示未测量，未适用布尔为−1，不是零；`sensor_process.queue_wait_ms` 未取得入队采样时为−1，入队采样依赖 `LIGHTNING_LM_COMPUTE_PROFILE` 开关。两个高频重建状态可能处于不同历元，先比较历元才能解释速度差。跨事件比较速度还要核对 `kind/source`、观测轴、参考点及运行配置：CAN更新在启用杆臂模式时会考虑后轴速度偏移，激光更新、高频重建和静止切换采用车体前向投影，不能直接假设各列对应同一参考点。重放滤波器分支中的CAN接受次数不等于独立CAN消息数。

CAN速度更新只改变速度状态；reason区分stale、invalid_measurement、nonfinite_prediction、innovation_gate、invalid_variance、nis_gate、velocity_step_gate、accepted。正式配置 `config/sany/localization.yaml` 中 `system.wheel_speed_dr_max_velocity_step_mps` 当前为0.35m/s，用于拒绝异常更新，不是平滑器，也不是下游发布相邻速度步长的绝对上限；具体run需核对其配置快照。

## 联合阅读资源、耗时与事件

1. 在 `causal_trace.csv` 中定位异常附近的事件，先按 `kind/source` 区分CAN融合、激光更新、高频重建和静止保持，再检查业务时间戳及更新前后状态。
2. 延迟问题先看 `sensor_process.queue_wait_ms`、`duration_ms` 与 `thread_cpu_ms`，并结合 `compute_profile.log` 的逐帧/窗口记录细分阶段。墙钟耗时与调用线程CPU耗时的差值还可能包含调度、锁、I/O及并行工作线程完成前的等待，不能直接解释为某一种等待。
3. 使用事件的 `wall_ns` 与资源sample的 `wall_time_ns` 关联时间窗口，在同一进程和run中结合 `tid` 查找对应线程的CPU、调度等待和 `wchan`，并结合 `start_ticks` 留意线程ID复用，再看整机各核与Top20进程负载。两种系统时间都可能受校时影响，遇到跳变需各自结合单调时间检查；传感器Header时间不用于计算本机排队耗时。

例如（假设数据）：某条慢帧事件显示排队20 ms、回调耗时30 ms、调用线程CPU 5 ms，且同一资源窗口中该TID有调度等待。这支持继续调查等待，但不能把25 ms差值全部归因于调度；还需结合锁、并行任务及其他证据。

资源采样描述区间负载与读取瞬间状态，事件追踪描述应用已记录的操作。1秒资源采样不足以完整还原毫秒级调度经过，CPU繁忙与慢帧同时出现也不能单独证明抢占原因；分析时还需检查下面的采样、丢弃与I/O边界。

## 开销与完整性

- C++追踪使用待接收容量为8192条的双缓冲队列、try-lock入队，后台处理中的批次另占一组缓冲。独立线程每轮等待20ms后取批次写出并flush，实际间隔还包含回调执行和线程调度时间；入队不主动唤醒，停止时提前唤醒。估计线程不等待文件I/O或缓冲空间。满缓冲/争用会丢诊断记录，不丢传感器消息。字段数值格式化在写线程执行，队列机制见 @ref lightning::BoundedAsyncQueue "BoundedAsyncQueue"。
- 正常结束末尾输出 `# received=... written=... dropped=... io_errors=...`，运行元数据会摘录它。`received` 为入队尝试次数（含随后丢弃），`dropped` 为入队失败次数，`written` 为事件 `fprintf()` 成功的行数，`io_errors` 为事件写入或批次刷新错误数。追踪关闭或I/O失败后直接跳过的事件不计入 `received`；`written` 也不保证磁盘持久化。分析必须核对footer、事件行数、丢弃和I/O错误。SIGKILL、断电或磁盘错误可能缺末尾，不能将缺失事件当成未发生。
- 进程采样只读 `/proc/*/stat`、定位线程wchan及已有资源字段；不读取全机命令行、其他进程环境、内核栈。新/退出进程没有配对CPU差分会记录未配对/消失数；短命进程可能遗漏。没有权限的记录不等于零CPU。
- 采样器自身耗时仍保留。Top20是截断榜单，不代表所有进程CPU之和。进程名可能被截断或重复，结合PID及start_ticks判断。
- 资源采样器与事件追踪器不修改调度优先级、CPU亲和性或内核参数。`kernel_sched_schedstats=0`时sched等待零值无解释力；wchan也仅是瞬时线索。若仍不能确定原因，再在域控获准进行短时调度跟踪，不能承诺仅靠应用日志区分所有抢占、锁和OpenMP等待。

## CAN A/B 回放约定

同一输入/地图/参数/二进制，顺序回放两组，仅切换 `system.enable_wheel_speed_dr_observation`。保持 `enable_wheel_speed_observation=true`，让两组静止检测仍用同样的CAN输入。不能用 `SANY_ENABLE_CAN_OBSERVATION=0` 替代此开关，否则同时改变静止检测。

对照中的CAN换算速度来自电机转速，不是独立车速真值。时间分析先画未对齐原始曲线，再计算动态段偏移；正偏移定义为estimate(t+lag)最接近CAN(t)，即估计曲线相对CAN滞后。不能把相关性偏移直接当成CAN网络传输时延，也不能由此断言应修改消息时间戳。

本机回放使用localhost隔离ROS域，不能将回放定位输出发送到车辆控制域；本机处理时间不作为Orin实时性验收。

## 测试

`python3 scripts/test_causal_trace.py -v`：Linux独立编译并发写入测试，核对footer/丢弃计数、不覆盖已有文件以及关闭路径。

`python3 scripts/test_resource_diagnostics.py -v`：CPU单位、PID/TID复用、系统进程排序、正式配置及真实/proc采样。
