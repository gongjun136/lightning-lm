@page incident_capture_guide 进程内故障采集与快速排查

SANY 配置默认开启 `diagnostics.incident_capture`，独立于全程录包开关 `SANY_RECORD_BAG`。采集器复用现有回调中的点云和主 IMU 引用，不新增 DDS 订阅或进程；独立线程约每 100 ms 检测状态，写线程使用 nice 10。尚未部署到 Orin，现场性能尚未实测。

## 开启、关闭与路径

在实际使用的定位 YAML 中设置：

```yaml
diagnostics:
  incident_capture:
    enabled: true  # 改为 false 关闭
```

SANY 启动脚本通过 `--incident_dir` 指向 `$out_root/incidents`；检查本次运行的配置快照和启动日志，确认设置与路径。手工启动时可显式传入绝对路径 `--incident_dir=/绝对路径/incidents`；未提供该参数但提供 `--output_tum` 时，回退到 TUM 父目录下的 `incidents`；若两者均未提供，日志会提示禁用采集，算法继续运行。事件各有目录；全程 bag 与事件采集分别开关。使用此采集器时可设置 `SANY_RECORD_BAG=0`、`SANY_ENABLE_POSE_VEL_WATCHDOG=0`，避免同时启用全程录包和独立看门狗；原始数据默认缓存且不增加重复 DDS 订阅。

## 现场三步

1. 找到 `$out_root/incidents` 下对应时间的事件目录，先打开本地 `summary.html`。查看触发原因、首个观测异常阶段、原始和上报码、图表与丢失计数。
2. 查看 `summary.json` 和 `timeline.csv`，按时间关联输入回调与时间戳、队列、匹配和发布；对照运行配置。诊断阶段是观测顺序，不等于硬件根因。
3. 需要复现时再检查或回放事件中的分段 MCAP，例如 `ros2 bag info /路径/segment_0`；事件包不含地图，复现需要对应地图。采集只包含实际启用雷达点云、主 IMU、开启融合时的 CAN、`pose_vel` 及 Safety、heartbeat、状态、流水线和参与雷达诊断等低带宽证据；不录重复 raw 位姿和不需要的 TF。是否足够回放取决于实际记录及丢弃情况。

## 触发与报告解读

原始 Safety 码变化到 2、4 或 -9 时触发，安全上报关闭仍有效；启动无输出、发布调用完成超时、输入回调断流或主 IMU/雷达 header 时间戳不推进的状态新出现时也触发。持续故障不会每 100 ms 重复创建事件。活跃且未关闭的事件可合并并延长，最长 60 秒；导出收尾或冷却中的新触发会被抑制并计数，不保证每次触发都有独立目录。

报告每个 stream 的 `first_relative_sec`、`last_relative_sec` 和 `written` 描述相对触发时间的首末记录及写入数；`callback_age_at_trigger`、`stamp_age_at_trigger` 描述触发时回调和时间戳推进的年龄，`-1` 表示未观测。`omitted_messages` 是导出容量限制省略的消息；`collector_totals` 分别记录 busy、budget、oversize 丢弃及 suppressed、writer_errors。`capture_loss_seen` 表示本进程生命周期中曾发生采集丢弃，不限于当前事件；有累计丢弃时状态保守标记 partial。complete 只表示事件导出结束，不保证每个话题的触发前 30 秒完整。

发布完成只证明本进程 `publish` 调用返回，不能证明下游 DDS 消费者收到。正常 `Finish` 会结束写入；如果后窗口尚未结束，会缩短窗口并标记 partial。`kill -9` 可能丢失未触发的内存历史和当前未关闭片段；已关闭的独立片段通常可读，不能承诺强杀后完整导出。

## 默认资源边界

目标窗口为触发前 30 秒、后 10 秒，容量不足时可能更短。默认 2 GiB 内存预算估算被保留的消息和记录，并非进程 RSS 上限；序列化单条消息、库缓冲和页缓存等另占资源。单条默认上限 64 MiB；单事件序列化硬上限 2 GiB。分段目标为 32 MiB 或约 2 秒关闭，单条消息可以超过分段目标。

档案保留目标 `disk_mb` 为 8 GiB，最多保留 10 个事件，每个分段关闭时检查实际磁盘大小；最后分段、单条大消息及元数据可带来额外空间占用。新段前预留 `max(segment_mb, 本条大小) + 1 MiB` 和 `min_free_mb`（默认 2 GiB），其他进程仍可能同时写盘，因此不是绝对空闲空间保证。默认平均写速限制 20 MiB/s，不压缩；nice 10 降低 CPU 调度优先级，不保证 IO 调度优先级。调整窗口或预算后需结合现场数据速率、内存和磁盘实际验证。

## 验证范围

独立 C++ 采集测试已覆盖上报为 0 时的短时 raw 4、前后窗口 MCAP、冻结时间戳、启动与发布开始后未完成的超时，以及容量、合并、归档保护、锁竞争和错误隔离。慢盘与小预算压力测试验证预算限制及回调不等待写盘；这是本机测试，不代表现场性能。SIGKILL 测试中，一个已关闭分段的 10 帧原始点云读回完整，报告仍处于 capturing 状态。

主程序 `run_loc_online` 和库编译成功，CTest 7/7、启动脚本测试 11/11、接口测试 11/11、依赖检查及 63 页文档链接检查通过。ROS 2 Humble、Cyclone DDS、`ROS_LOCALHOST_ONLY=1` 的 production 测试中：三雷达、Safety 上报及 shadow 关闭时，301 条 Safety 约 50.0 Hz；一个事件的 479 条 MCAP 消息包含三路原始点云和主 IMU，逐条读回校验模拟 payload，raw 4/-9 和启动超时保留，输入订阅者数量未增加。四雷达默认 Safety 上报场景 301 条 Safety 约 50.1 Hz。此结果仅覆盖 Cyclone DDS，不代表其他 RMW 已验收，也不是实车精度或负载验收。 本机 WSL 的独立 ROS 域对照中，显式使用 `rmw_fastrtps_cpp`、`ROS_LOCALHOST_ONLY=1` 且关闭采集器时，探针 20 秒内仍未收到 Safety，而进程日志显示持续发布诊断；Fast DDS 因此尚未通过联调，也不能据此归因采集器或认定 DDS 根因。HTML 的离线 JSON 结构和 JavaScript 语法检查通过，尚未完成浏览器视觉渲染验收。

测试依赖固定版本如下；`check_workspace_repos` 与接口检查通过，本次没有消息接口或依赖提交变动。

| 依赖 | 提交 |
|---|---|
| common_msgs | `a6e3045097e5433f35370694d858609ef6a17a11` |
| cgi430_driver | `6bdedcaf317aca7d0f74119005871dc8a4989271` |
| livox_ros_driver2 | `1718003576eed30e1801e3e90e2f2a7e16d1ecf7` |
| Livox-SDK2 | `92a02baa44d87aa53562af4db60cfed684835d89` |
