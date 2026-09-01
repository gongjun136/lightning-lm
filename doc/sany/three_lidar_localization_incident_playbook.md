# SANY 三/四雷达域控定位丢失排查与取证方案

> 当前现场拓扑为主 Orin 三雷达，唯一操作入口为 `bash scripts/run.sh`。本文中的四雷达、双 Orin、Zstd 和旧故障目录只用于历史取证；不得直接复制历史启动命令。三雷达默认重定位后端为 `solid_kiss`，纯 `solid` 为回退。

## 结论

三次故障不是 `/PosRes` 发布器自身退出，而是同一条故障链在不同阶段的表现：输入处理产生大幅时间倒退，LIO 随后发散，连续 5 帧激光定位无效后发布门控主动停止 `/PosRes`。算法仍在运行并继续尝试全局重定位。

第三次故障还暴露了恢复路径中的独立缺陷：全局重定位在 21:51:38 成功后，PGO 重置保留了旧的 DR/LidarOdom 相对位姿队列，旧纪元数据被带入新图，随后连续出现 `Cholesky failure, solve failed`。本分支已清除这些队列，并拒绝倒退的 LidarOdom 输入。

现有三次运行没有录制传感器包，因此能确认软件内的故障链，不能仅凭日志唯一确定最初的外部诱因是 DDS/网络、驱动、系统时钟还是算力抢占。新增的压缩 MCAP 和 10 Hz 管线诊断正是用于补齐这段证据。

## 2026-08-18 四雷达故障结论

五个 `sany_4lidar_solid_*` 目录补齐了压缩点云和管线诊断，可以把问题分成两个同时存在、但边界不同的故障：

1. **算法消费速率低于输入速率。** 五组原始运行的 LIO 单帧平均耗时为 48.6–82.0 ms、P95 为 96.0–193.3 ms。10 Hz 雷达仅 LIO 就占用 49%–82% 的单核墙钟预算，再叠加 200 Hz IMU、全局定位、DDS 解压和调度开销后，长期平均利用率可超过 100%。因此“某一帧小于 100 ms”不等于系统实时；队列会先缓慢积累，几分钟后才满并丢帧。五组旧运行的最大 pending 为 12,454–16,534，最大传感器滞后为 107–235 s。队列满是过载的后果，不是最初原因。
2. **内核更新后算法端 Best-Effort 点云订阅更容易丢失大而分片的 DDS 消息。** 压缩器/录包使用 Reliable，算法原来使用 `SensorDataQoS`（Best Effort）。以 bag 中实际录到的帧为基准，旧算法回调覆盖率在主雷达约为 72%–93%，后雷达 143 最差仅约 28%。同时，部分 bag 在录制前就已经缺少后雷达帧并出现 14–23 s 间隙，这一部分属于第二块 Orin、Livox 驱动、网卡或两 Orin 之间链路，算法无法补回。

逐帧 Zstd 解压比较了新数据 99,313 帧和旧数据 24,477 帧。四台雷达的有效点数 P5/P50/P95 均约为 `19968/19968/20064`，无空帧、无解压失败，布局一致。因此本次问题不是内核更新后点云点数突然增大。

正式四雷达配置采用三层处理：LiDAR QoS 改为 Reliable/KeepLast(16)；LIO 从 10 Hz 降为 5 Hz而 IMU 传播和 `/PosRes` 保持原频率；传感器时延超过 0.5 s 时临时拒绝新雷达帧，恢复到 0.1 s 后继续，并始终保留 IMU。传感器队列上限为 1,000，激光定位队列只保留最新待处理帧。

本地 1.0× 完整回放五个最新 bag 后，`/PosRes` 为 195.7–198.0 Hz，传感器队列峰值 17–28、每组结束均为 0，地图输出开启后没有再次关闭。最差乱序样本 `200009` 在单调时间前沿修复后，最大实时 backlog 为 0.15 s、过载保护丢帧为 0。重定位 `data1..data5 × 10` 冷启动矩阵为 50/50，通过全部安全门限；累计处理 P95 为 4.60 s，最低地图重叠率为 0.9879，与历史同起点结果的最大位置/航向差为 0.082 m/0.53°。

### 最终丢定位应看哪里

不要只找一条 `ERROR`。完整因果链分布在三个位置：

1. `F:\tmp\runs\sany_4lidar_solid_*\logs\posres_watchdog.csv` 是最终外部症状时间线。`LOST` 表示 `/PosRes` 静默，`RECOVERED` 表示恢复；例如 `200009` 最后一次不可恢复的 `LOST` 为 `2026-08-18 20:05:12 +08:00`。
2. `F:\tmp\runs\sany_4lidar_solid_*\bag\localization_incident\localization_incident_0.mcap` 中的 `/localization/pipeline_diagnostics` 是根因主日志。先看 `current_sensor_lag_sec`、`sensor_queue_pending/dropped`；再看四路 `lidar_callback_counts`/静默时间；最后看 `consecutive_lost_frames` 达到 5 且 `map_outputs_enabled=false`。这条诊断能证明“积压/断流 → 陈旧或不完整点云 → 连续定位失败 → 输出门控关闭”。
3. `F:\tmp\runs\sany_4lidar_solid_*\logs\run_loc_online.stderr.log` 是算法细节。搜索 `传感器顺序队列 exceeds largest size` 可定位队列开始实际丢帧的时刻，搜索 `GLOBAL_RELOCALIZATION` 可查看丢定位后的候选接受/拒绝原因。以 `193432` 为例，顺序队列从 `19:37:50` 起反复满 10,000；它晚于积压开始，只是故障已经恶化的证据。

因此，**最终触发地图定位输出关闭**的是诊断中的 `consecutive_lost_frames >= 5`；**最早且最有因果价值的信号**是 MCAP 诊断中的 lag/pending 持续增长，而不是几分钟后的 queue-full 日志。

## 三次日志的共同证据

| 运行 | 初次定位 | 明显时间倒退 | `/PosRes` 对应的定位丢失 | 后续恢复 |
|---|---:|---:|---:|---|
| 2026-08-07 09:41 | 09:42:10 | 09:53:13，最大约 12.54 s | 09:55:38 | 未恢复；LIO 数值持续发散，重定位均拒绝 |
| 2026-08-10 20:13 | 20:13:49 | 20:17:40，最大约 34.27 s | 20:17:42 | 未恢复；点云有效面特征不足，LIO 数值发散 |
| 2026-08-10 21:03 | 21:03:18 | 21:51:04，最大约 33.55 s | 21:51:11 | 21:51:38 重定位成功；随后 PGO Cholesky 持续失败 |

在旧测试分支中，ROS 输入回调与后端处理共用一把大锁。后端计算变慢时，输入回调也会被阻塞，单线程 executor 更容易积压。最新分支已经把生命周期、输入和后端处理锁拆开；本次又补上了时间倒退拒绝和重定位后的 PGO 历史清理。

## `/PosRes` 为什么停止，以及是否会立即重定位

`PublishProcessedCloud()` 每帧观察 `lidar_loc_valid_`。连续失败达到 YAML 中的 `relocalization.lost_frame_threshold`（当前为 5）后，地图坐标输出门控关闭，因此 `/PosRes`、`/slamPoseRaw_topic` 和地图点云停止；进程、LIO 和激光定位线程仍继续运行。

`LidarLoc` 在进入 LOST 后会继续逐帧尝试全局重定位。当前三雷达链路由 SOLiD 召回地点，KISS-Matcher 风格模块执行鲁棒无初值粗配准，候选再通过地图预检、点到面精配准、最终地图一致性和连续两帧确认。新的有效匹配到达后门控自动重新打开。外层脚本只记录和告警，不重启算法，也不触发车辆安全停机。

## 新增在线诊断

话题 `/localization/pipeline_diagnostics` 以 10 Hz 发布：

- 每路 LiDAR 的 topic、累计消息数、最后传感器时间和墙钟静默时长；
- 主 IMU 的相同指标；
- 传感器顺序队列和激光定位队列的 pending、dropped、processed；
- 最新入队/处理时间、当前/最大处理落后量、超过 1 s 的时间倒退次数及最坏值；
- 全局重定位尝试/成功次数、候选、得分、耗时和最近原因；
- 连续失败帧数、地图输出门控、最后定位结果和最后 `/PosRes` 时间。

现场判断顺序：

| 现象 | 首要检查 | 更可能的问题 |
|---|---|---|
| 某一路 `lidar_silence_sec` 持续大于 0.5 s，或 IMU 大于 0.1 s | 对应计数是否停止、`ip -s link` | 雷达/驱动/网络/DDS 上游断流 |
| 输入计数正常增长，`sensor_queue_pending`、`current_sensor_lag_sec` 持续增长 | `tegrastats`、CPU、内存、温度 | LIO 后端变慢或调度/锁竞争 |
| 传感器队列正常，`localization_queue_pending` 增长 | 重定位耗时、NDT/BTC 日志 | 激光定位或全局检索跟不上输入 |
| 队列正常，`map_outputs_enabled=false`，重定位尝试持续增长 | 最近 `last_relocalization_reason` | 场景退化、地图变化、点云质量或初值问题 |
| `relocalization_accept_count` 增长，但 `/PosRes` 不恢复 | PGO/solver 日志和最后两个时间戳 | 融合恢复路径或输出链异常 |

阈值是排查起点，不是新的算法判定条件；最终应结合车辆速度、LiDAR 实际频率和 MCAP 回放调整。

## 域控部署与运行

在 `feature/gj_change_2025.11.19` 完成 Release 构建并部署后，先启动主 Orin 上的前 184、左 108、右 133 三路 Livox；开启 CAN 观测时同时启动 `/SpeThrCAN4_topic`。新终端只执行统一入口：

```bash
set +e
set +u
cd /home/nvidia/project/gj_ws/lightning-lm
bash scripts/run.sh
```

当前 `run.sh` 固定 `SANY_LIDAR_LAYOUT=3`、`SANY_RECORD_BAG=0`、`SANY_ENABLE_POSE_VEL_WATCHDOG=0`，并清除外部 `LIGHTNING_LM_CONFIG`。因此：

- 操作员不得直接调用 `run_sany_online_diagnostics.sh`；它只是执行器；
- 当前运行不要求 Zstd 压缩器或 `/zstd` Topic；
- 当前地图目录必须同时包含 `solid_relocalization/database.yaml` 和 `btc_relocalization/submap_*.pcd`；
- 启动日志必须显示 `coarse_registration_backend=kiss_matcher`、`validation_top_k=8`、`retrieval_pool_size=64`；
- CAN 隔离只允许在同一入口下执行 `SANY_ENABLE_CAN_OBSERVATION=0 bash scripts/run.sh`；
- 若要启用 MCAP 或 `/localization/pose_vel` watchdog，必须先修改统一入口、形成新的代码提交并重新验收，不能绕过入口临时启动。

每次运行的事实源是 `runs/sany_3lidar_<timestamp>/run_metadata.txt`、冻结的 `config.yaml` 和 `logs/run_loc_online.stderr.log`。正常重定位日志应从 `GLOBAL_RELOCALIZATION[solid_kiss] ... count=1/2` 进入 `accepted ... confirmations=2`。

## 故障包回放

压缩包必须先解压回配置中的原始 PointCloud2，再喂给定位：

```bash
# 终端 1：解压节点
source /opt/ros/humble/setup.bash
source /home/nvidia/project/gj_ws/sdk/livox_sdk/install/setup.bash
ros2 run livox_ros_driver2 pointcloud_zstd_decompressor

# 终端 2：启动待验证版本定位（不要用嵌入式 --bag 模式）
source /home/nvidia/project/gj_ws/lightning-lm/install/setup.bash
ros2 run lightning run_loc_online \
  --config=/absolute/path/config.yaml \
  --map=/absolute/path/to/map \
  --output_tum=/tmp/replay_global.tum \
  --output_high_frequency_tum=/tmp/replay_high_frequency.tum

# 终端 3：只播放传感器输入，避免回放旧输出与新输出冲突
ros2 bag play /path/to/localization_incident \
  --topics \
  /livox/lidar_192_168_3_184/zstd \
  /livox/lidar_192_168_1_108/zstd \
  /livox/lidar_192_168_2_133/zstd \
  /livox/lidar_192_168_4_143/zstd \
  /livox/imu_192_168_3_184
```

验证至少包括：全部配置雷达的解压消息计数与 bag 一致；时间戳乱序应被统计但不能使已处理前沿倒退；两级队列无持续增长，传感器队列最终 pending/lag 为 0；失配时重定位尝试增加；接受重定位后 `map_outputs_enabled` 和 `/PosRes` 自动恢复；日志中不出现连续 Cholesky failure。

## 本地验收结果（2026-08-11）

- Ubuntu 22.04 / ROS 2 Humble Release 构建成功；
- CTest 全量 14/14 通过；
- 使用 `F:/datasets/SANY/3lidar_online_loc/data1` 的 124.4 s 三雷达 MCAP，以 1× 速度、同批三雷达地图和 BTC 数据库完成在线回放；
- BTC 全局重定位接受 1 次；PGO reset 丢弃 3 个早于 reset 相对位姿水位线的积压定位帧，随后自动恢复；
- 全局轨迹 860 帧，时间严格单调，最大间隔 0.100397 s，最大相邻位移 0.453932 m；
- 高频轨迹 16119 帧，时间严格单调，最大间隔 0.080250 s，最大相邻位移 0.335778 m；
- 两条轨迹均到达最后一帧 LIO；Cholesky failure、LidarOdom 时间倒退、段错误均为 0；
- DR 有 2 次微小逆序（0.106 ms、0.0088 ms），属于跨流时间抖动，远小于现场的 12–34 s；
- 运行中订阅 `/localization/pipeline_diagnostics` 成功，三路 LiDAR 和 IMU 指标齐全；采样时队列 dropped 为 0、严重时间倒退为 0。

三次现场目录缺少原始或压缩传感器 bag，因此本地回放验证的是代码恢复路径和正常 SANY 三雷达数据兼容性，不等价于复现三次现场的最初诱因。下一次现场故障必须带回本方案录制的 MCAP，才能完成外因闭环。

## 每次故障需要带回的最小材料

完整复制单次 `sany_loc_diag_*` 目录，包括：

- `bag/localization_incident/`；
- `logs/run_loc_online.stderr.log`、`tegrastats.log`、`posres_watchdog.csv`；
- `snapshots/`；
- `run_metadata.txt`、`git_status.txt`、`config.yaml`；
- `results/trajectory_global.tum` 和 `trajectory_high_frequency.tum`（正常 Ctrl-C 后生成）。

地图应使用不可变版本目录，并同时记录地图版本或校验值。没有对应地图时，即使传感器包完整也无法做等价重定位回放。
