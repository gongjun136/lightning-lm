@page incident_capture_guide 进程内故障采集与快速排查

SANY 配置默认开启 `diagnostics.incident_capture`，独立于全程录包开关 `SANY_RECORD_BAG`。采集器复用现有回调中的点云和主 IMU 引用，不新增 DDS 订阅或进程；独立线程约每 100 ms 检测状态，写线程使用 nice 10。原采集器版本 `dc951ed` 已由用户部署到 Orin 并观察到运行；本次诊断脚本已在现场验证。采集器并发修复需要重新编译并重启定位后生效，本次未重启定位。

## 开启、关闭与路径

在实际使用的定位 YAML 中设置：

```yaml
diagnostics:
  incident_capture:
    enabled: true  # 改为 false 关闭
```

SANY 启动脚本通过 `--incident_dir` 指向 `$out_root/incidents`；检查本次运行的配置快照和启动日志，确认设置与路径。手工启动时可显式传入绝对路径 `--incident_dir=/绝对路径/incidents`；未提供该参数但提供 `--output_tum` 时，回退到 TUM 父目录下的 `incidents`；若两者均未提供，日志会提示禁用采集，算法继续运行。事件各有目录；全程 bag 与事件采集分别开关。使用此采集器时可设置 `SANY_RECORD_BAG=0`、`SANY_ENABLE_POSE_VEL_WATCHDOG=0`，避免同时启用全程录包和独立看门狗；原始数据默认缓存且不增加重复 DDS 订阅。

## 现场单命令分析

域控部署的快捷入口为 `~/Sany/diagnose_loc.sh`；其他机器不自动带此链接，可使用仓库中的 `scripts/diagnose_localization.sh`。例如：

```bash
bash ~/Sany/diagnose_loc.sh                     # 自动选择当前实例或最新事件
bash ~/Sany/diagnose_loc.sh /路径/事件目录       # 指定事件或包
bash ~/Sany/diagnose_loc.sh --json > diagnosis.json  # 保存 JSON 输出
```

完整用法为：

```text
diagnose_localization.sh [path] [--pid PID] [--root PATH] [--json] [--log-mb 32] [--bag-mb 128]
```

`path` 可以是事件目录、`summary.json`、`timeline.csv`、运行目录、归档根目录、裸 MCAP 或 db3 文件，也可以是 rosbag 目录；不支持压缩归档。显式路径优先，不选择当前进程。例如：

```bash
bash ~/Sany/lightning_lm_ws/src/lightning-lm/scripts/diagnose_localization.sh /路径/事件目录
```

不带路径时，脚本通过 `/proc/<pid>/exe` 识别 `run_loc_online`，不依赖可能显示为 `loc_executor` 的线程名。多个实例时选择启动最晚的实例，并提示可用 `--pid` 指定；通过 `run_metadata` 和启动时间关联本次运行，避免 PID 复用或误选其他运行的事件。本次运行尚无事件时分析本次运行日志。无定位进程时，在 `LIGHTNING_LM_OUT_ROOT`、`SANY_WS/runs` 等默认路径寻找最新包；前者可来自环境或 `~/.config/sany/workspace.bash` 中的导出值。`--root` 仅覆盖无进程时的搜索根目录。

脚本默认用标准库读取摘要、时间线和运行日志尾部 32 MiB，另读取日志头部 1 MiB 以获取配置；`--log-mb` 可调整尾部范围，报告会标明覆盖范围。新版本时间线包含细分 Safety 判据；旧包需要同次运行日志补充判据证据；日志无法恢复未记录的物理根因。指定事件只关联其时间窗内的日志；自动检查正在运行的实例时，还会显示当前运行日志尾部的历史故障并标明来源时间。

裸 bag 缺少摘要和时间线时，脚本通过 `rosbag2_py` 离线筛选 Safety，不创建 ROS 节点、不解点云；默认文件预算 128 MiB，可用 `--bag-mb` 调整，读取还受 15 秒和 30000 次状态变化上限约束。只有上报码而缺少原始细分判据时，报告会说明证据不足。脚本不订阅 DDS、不启动定位、不修改阈值或输入数据；默认不写报告，`--json` 将 JSON 打印到标准输出，可自行重定向保存。

脚本及 ROS 环境准备辅助进程以 nice 10 运行，从自身允许的 CPU 集合中排除系统隔离核；若没有可用的非隔离核则拒绝分析。退出码 `0` 表示读到了可分析的 Safety 记录，不能理解为没有故障；无证据或选择失败返回 `2`。检查报告中的事件时间、原始及发布故障码、具体字段及阈值或状态、证据位置和缺失证据。判据触发原因不能直接当作已证实的物理根因，证据不足时应继续调查。

域控已验证 live、指定事件和裸 MCAP 分析；无进程回退通过向诊断程序注入空进程列表模拟，未停止真实定位。普通 SSH 环境下的一次默认运行耗时 2.8027 秒、CPU user 2.2159 秒及 system 0.4363 秒、峰值 RSS 119768 KiB；这是单次现场观测，不是性能承诺。测试前三个关键线程 `loc_sensor`、`loc_executor`、`loc_output` 分别绑定 CPU 10、11、11，测试后保持不变。

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

## 验证范围

本次验证中，域控 Python 单测 24/24 通过；本机 `run_loc_online` 与库编译成功，CTest 7/7 通过，覆盖并发时间回归和 CSV 判据列及值检查；依赖与接口检查 11 项通过，63 页文档链接检查无错误。对本地 10 月 27 日旧运行包解析出 3940 条诊断记录，并识别出 2、4、-9 状态。域控定位进程 PID 1374884 在验证期间持续运行，未重启；本次 C++ 采集器修复仅在源码中，需后续重新编译并重启定位后生效。

以下是此前的采集器与 DDS 集成验证记录，并非本次重新运行的测试：

独立 C++ 采集测试已覆盖上报为 0 时的短时 raw 4、前后窗口 MCAP、冻结时间戳、启动与发布开始后未完成的超时，以及容量、合并、归档保护、锁竞争和错误隔离。慢盘与小预算压力测试验证预算限制及回调不等待写盘；这是本机测试，不代表现场性能。SIGKILL 测试中，一个已关闭分段的 10 帧原始点云读回完整，报告仍处于 capturing 状态。

主程序 `run_loc_online` 和库编译成功，CTest 7/7、启动脚本测试 11/11、接口测试 11/11、依赖检查及 63 页文档链接检查通过。ROS 2 Humble、Cyclone DDS、`ROS_LOCALHOST_ONLY=1` 的 production 测试中：三雷达、Safety 上报及 shadow 关闭时，301 条 Safety 约 50.0 Hz；一个事件的 479 条 MCAP 消息包含三路原始点云和主 IMU，逐条读回校验模拟 payload，raw 4/-9 和启动超时保留，输入订阅者数量未增加。四雷达默认 Safety 上报场景 301 条 Safety 约 50.1 Hz。此结果仅覆盖 Cyclone DDS，不代表其他 RMW 已验收，也不是实车精度或负载验收。 本机 WSL 的独立 ROS 域对照中，显式使用 `rmw_fastrtps_cpp`、`ROS_LOCALHOST_ONLY=1` 且关闭采集器时，探针 20 秒内仍未收到 Safety，而进程日志显示持续发布诊断；Fast DDS 因此尚未通过联调，也不能据此归因采集器或认定 DDS 根因。HTML 的离线 JSON 结构和 JavaScript 语法检查通过，尚未完成浏览器视觉渲染验收。

测试依赖固定版本如下；`check_workspace_repos` 与接口检查通过，本次没有消息接口或依赖提交变动。

| 依赖 | 提交 |
|---|---|
| common_msgs | `a6e3045097e5433f35370694d858609ef6a17a11` |
| cgi430_driver | `6bdedcaf317aca7d0f74119005871dc8a4989271` |
| livox_ros_driver2 | `1718003576eed30e1801e3e90e2f2a7e16d1ecf7` |
| Livox-SDK2 | `92a02baa44d87aa53562af4db60cfed684835d89` |
