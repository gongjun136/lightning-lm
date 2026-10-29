@page incident_capture_guide 进程内故障采集与快速排查

SANY 配置默认开启 `diagnostics.incident_capture`，独立于全程录包开关 `SANY_RECORD_BAG`。采集器复用现有回调中的点云和主 IMU 引用，不新增 DDS 订阅或进程；独立线程约每 100 ms 检测状态，写线程使用 nice 10。采集器并发修复需要重新编译并重启定位后生效；新编译版本的域控现场验收待执行。

## 开启、关闭与路径

在实际使用的定位 YAML 中设置：

```yaml
diagnostics:
  incident_capture:
    enabled: true  # 改为 false 关闭
```

SANY 启动脚本通过 `--incident_dir` 指向 `$out_root/incidents`；检查本次运行的配置快照和启动日志，确认设置与路径。手工启动时可显式传入绝对路径 `--incident_dir=/绝对路径/incidents`；未提供该参数但提供 `--output_tum` 时，回退到 TUM 父目录下的 `incidents`；若两者均未提供，日志会提示禁用采集，算法继续运行。事件各有目录；全程 bag 与事件采集分别开关。使用此采集器时可设置 `SANY_RECORD_BAG=0`、`SANY_ENABLE_POSE_VEL_WATCHDOG=0`，避免同时启用全程录包和独立看门狗；原始数据默认缓存且不增加重复 DDS 订阅。

## 现场单命令分析

现场使用仓库外已有的部署入口 `~/Sany/run_loc.sh`；其 `lidar` 分支将参数原样传给仓库内的 `scripts/run_sany_lidar_loc.sh`。执行：

```bash
bash ~/Sany/run_loc.sh --diagnose
bash ~/Sany/run_loc.sh --diagnose /路径/事件目录
bash ~/Sany/run_loc.sh --diagnose --all --color always
bash ~/Sany/run_loc.sh --diagnose --json > diagnosis.json
```

开发联调也可从仓库根目录直接执行 `bash scripts/run_sany_lidar_loc.sh --diagnose [path]`，并附加相同选项。仓库入口在启动定位、录包或创建运行目录之前转入 `scripts/diagnose_localization.sh`；独立 shell 和 Python 脚本仍可直接使用。支持 `--pid`、`--root`、`--log-mb`、`--bag-mb`、`--all` 和 `--color auto|always|never`。默认文本显示最近 12 段及按来源分类的全部汇总；同一真实故障可能同时出现在日志段和事件时间线段，不能将这些统计相加当作事件数。诊断开头优先给出结论和简短问题链；`--all` 展开所有已扫描判据，JSON 包含全部已分析段，日志读取仍受预算限制，不能保证覆盖全量日志。`auto` 遵守 TTY、`NO_COLOR` 和 `TERM=dumb`；JSON 始终无 ANSI。红色表示最新记录故障，黄色表示历史或缺失证据，绿色仅表示最新记录为 0，青色为信息。

`path` 可以是事件目录、`summary.json`、`timeline.csv`、运行目录、归档根目录、裸 MCAP 或 db3 文件，也可以是 rosbag 目录；不支持压缩归档。显式指定运行目录时检查该次运行日志，仍受 `--log-mb` 预算限制；显式指定事件目录时按该事件时间窗检查。无进程自动查找仍选择最新包。例如：

```bash
bash ~/Sany/run_loc.sh --diagnose /路径/事件目录
```

不带路径时，脚本通过 `/proc/<pid>/exe` 识别 `run_loc_online`，不依赖可能显示为 `loc_executor` 的线程名。多个实例时选择启动最晚的实例，并提示可用 `--pid` 指定；通过 `run_metadata` 和启动时间关联本次运行，避免 PID 复用或误选其他运行的事件。本次运行尚无事件时分析本次运行日志。无定位进程时，在 `LIGHTNING_LM_OUT_ROOT`、`SANY_WS/runs` 等默认路径寻找最新包；前者可来自环境或 `~/.config/sany/workspace.bash` 中的导出值。`--root` 仅覆盖无进程时的搜索根目录。

脚本默认用标准库读取摘要、时间线和运行日志尾部 32 MiB，另读取日志头部 1 MiB 以获取配置；`--log-mb` 可调整尾部范围，报告会标明覆盖范围。新版本时间线包含细分 Safety 判据；旧包需要同次运行日志补充判据证据；日志无法恢复未记录的物理根因。指定事件只关联其时间窗内的日志；自动检查正在运行的实例时，还会显示当前运行日志尾部的历史故障并标明来源时间。

裸 bag 缺少摘要和时间线时，脚本通过 `rosbag2_py` 离线筛选 Safety，不创建 ROS 节点、不解点云；默认文件预算 128 MiB，可用 `--bag-mb` 调整，读取还受 15 秒和 30000 次状态变化上限约束。只有上报码而缺少原始细分判据时，报告会说明证据不足。脚本不订阅 DDS、不启动定位、不修改阈值或输入数据；默认不写报告，`--json` 将 JSON 打印到标准输出，可自行重定向保存。

脚本及 ROS 环境准备辅助进程以 nice 10 运行，从自身允许的 CPU 集合中排除系统隔离核；若没有可用的非隔离核则拒绝分析。退出码 `0` 表示读到了可分析的 Safety 记录，不能理解为没有故障；无证据或选择失败返回 `2`。检查报告中的事件时间、原始及发布故障码、具体字段及阈值或状态、证据位置和缺失证据。判据触发原因不能直接当作已证实的物理根因，证据不足时应继续调查。

记录按来源、设备和判据时段保留，区分启动未输出与运行超时，并保留首次异常、较坏样例及附近重定位、CAN、IMU 日志。Safety 日志 sequence 不连续时报告会提示；未知的新故障码会提示规则缺失，不按健康处理。状态优先依据最新时间戳；旧记录或超过 5 秒的记录不能证明当前健康。摘要和时间线范围会与故障起始时刻交叉检查；未找到覆盖时提示可能仍在导出、被抑制或已轮转，范围重叠也不保证逐条完整。事件摘要最多 256 条，采集统计区分累计计数和窗口丢弃增量。

关联日志受同源 Safety 恢复为 0 和观测缺口的边界约束，避免相邻故障时段串联；保留最近 1000 条普通日志及最早 128 条关键边界日志，去重并在超限时明确提示；保留首次局部匹配失败和禁用输出证据，避免长期重定位日志刷掉故障起点。独立 runtime 保护统计覆盖扫描范围或指定事件范围内的过期校正丢弃和自适应降载，不能作为 Safety 故障次数。匹配失败日志记录地图点数、XY 范围内点的比例、重力余弦及其阈值（`map_points`、`map_inside_xy`、`map_gravity_cos`、`map_min_inside_xy`、`map_min_overlap`、`map_min_gravity_cos`）；解析优先使用日志中的实际阈值，旧包配置值仅作线索。`map_checked=1` 或 overlap 不足不能证明车辆驶出地图。首次输出或 epoch 重建后尚无 yaw PASS 时，应检查窗口建立与历史 IMU 覆盖；状态 4 的 UNKNOWN 无法可靠区分这两种情况，应说明证据不足。启动阶段无有效输出且状态 4 的判据为 UNKNOWN 时，优先核查初始化；文本须区分切换为 -9 与恢复为 0。并存雷达条件仅作为对应历史 -9 记录的旁证，不虚报状态 2，也不能据此证明当前仍异常或认定故障因果。事件覆盖故障起点不等于覆盖整段故障。上述规则的现场验证以实际复验结果为准。

报告对照运行库的 `LOCALIZATION_DIAGNOSTIC_BUILD` 标记与脚本合同；缺失或不匹配时为未验证。标记包含合同 ID、source revision、dirty 状态和构建时间；`run_metadata` 的 `git_commit` 不能证明运行二进制版本。

## 现场三步

1. 找到 `$out_root/incidents` 下对应时间的事件目录，先打开本地 `summary.html`。查看触发原因、首个观测异常阶段、原始和上报码、图表与丢失计数。
2. 查看 `summary.json` 和 `timeline.csv`，按时间关联输入回调与时间戳、队列、匹配和发布；对照运行配置。诊断阶段是观测顺序，不等于硬件根因。
3. 需要复现时再检查或回放事件中的分段 MCAP，例如 `ros2 bag info /路径/segment_0`；事件包不含地图，复现需要对应地图。采集只包含实际启用雷达点云、主 IMU、开启融合时的 CAN、`pose_vel` 及 Safety、heartbeat、状态、流水线和参与雷达诊断等低带宽证据；不录重复 raw 位姿和不需要的 TF。是否足够回放取决于实际记录及丢弃情况。

## 触发与报告解读

原始 Safety 码变化到 2、4 或 -9 时触发，安全上报关闭仍有效；启动无输出、发布调用完成超时、输入回调断流或主 IMU/雷达 header 时间戳不推进的状态新出现时也触发。持续故障不会每 100 ms 重复创建事件。活跃且未关闭的事件可合并并延长，最长 60 秒；导出收尾或冷却中的新触发会被抑制并计数，不保证每次触发都有独立目录。

报告每个 stream 的 `first_relative_sec`、`last_relative_sec` 和 `written` 描述相对触发时间的首末记录及写入数；`callback_age_at_trigger`、`stamp_age_at_trigger` 描述触发时回调和时间戳推进的年龄，`-1` 表示未观测。采集器源码的并发修复将已观测时间戳晚于采样时刻的微小负年龄钳位为 0，避免误判为未观测及误触发启动无输出；修复需重新编译并重启后生效，旧版本事件可能受此问题影响。`omitted_messages` 是导出容量限制省略的消息；`collector_totals` 分别记录 busy、budget、oversize 丢弃及 suppressed、writer_errors。`capture_loss_seen` 表示本进程生命周期中曾发生采集丢弃，不限于当前事件；有累计丢弃时状态保守标记 partial。complete 只表示事件导出结束，不保证每个话题的触发前 30 秒完整。

发布完成只证明本进程 `publish` 调用返回，不能证明下游 DDS 消费者收到。正常 `Finish` 会结束写入；如果后窗口尚未结束，会缩短窗口并标记 partial。`kill -9` 可能丢失未触发的内存历史和当前未关闭片段；已关闭的独立片段通常可读，不能承诺强杀后完整导出。

## 默认资源边界

目标窗口为触发前 30 秒、后 10 秒，容量不足时可能更短。默认 2 GiB 内存预算估算被保留的消息和记录，并非进程 RSS 上限；序列化单条消息、库缓冲和页缓存等另占资源。单条默认上限 64 MiB；单事件序列化硬上限 2 GiB。分段目标为 32 MiB 或约 2 秒关闭，单条消息可以超过分段目标。

档案保留目标 `disk_mb` 为 8 GiB，最多保留 10 个事件，每个分段关闭时检查实际磁盘大小；最后分段、单条大消息及元数据可带来额外空间占用。新段前预留 `max(segment_mb, 本条大小) + 1 MiB` 和 `min_free_mb`（默认 2 GiB），其他进程仍可能同时写盘，因此不是绝对空闲空间保证。默认平均写速限制 20 MiB/s，不压缩；nice 10 降低 CPU 调度优先级，不保证 IO 调度优先级。调整窗口或预算后需结合现场数据速率、内存和磁盘实际验证。

## 重新编译后的现场验收

sany-03-master 的新编译版本现场验收待用户重新编译并重启后执行：核对运行库编译标记与脚本合同一致，再分别检查有进程、无进程回退和指定事件包。核对 2、4、-9 的实际字段，以及关闭上报时原始故障仍可见；检查覆盖范围、缺口、导出和抑制提示，并确认分析进程避开 CPU 10、11。记录实际结果和版本，不能以旧运行结果代替新版本验收。

采集器导出收尾或冷却期抑制新触发的行为尚未修复。并发年龄修复需重新编译后生效；旧格式与缺少构建标记的报告应按未验证处理。

## 开发同步检查

修改诊断相关源码后，在仓库根目录执行：

```bash
python3 scripts/check_diagnostic_contract.py
python3 scripts/test_diagnose_localization.py
```

源码摘要覆盖 Monitor、雷达和拓扑、Safety 发布及初始化/丢失门控、采集器、diagnostics 配置段和 `workspace.repos`。源码变化时先更新并审查解析逻辑和回归用例，再执行 `python3 scripts/check_diagnostic_contract.py --update`；不能只更新摘要。CMake 配置和 CTest 也会检查合同。该合同不是所有算法行为或物理根因的形式化保证。
