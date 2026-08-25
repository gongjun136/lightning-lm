# Lightning-LM 脚本说明

仓库沿用 `scripts/` 目录。根目录脚本是稳定运行入口和公共辅助工具；`reproduction/` 保存已交付实验的复现、审计和评价脚本。

## 稳定运行入口

| 脚本 | 用途 |
|---|---|
| `run_frontend_offline.sh` | 单个 SQLite3 ROS2 bag 的统一离线前端入口。单雷达/多雷达由 YAML 控制，并输出轨迹、地图、帧统计、处理耗时、资源监控和验收元数据。 |
| `run_frontend_offline_batch.sh` | 对同一个 bag 依次运行多个 YAML；默认运行 SANY C0—C6，可用 `--repeats` 重复实验。每个子任务仍调用 `run_frontend_offline.sh`。 |
| `run_frontend_online.sh` | 启动在线前端节点，供实时传感器或 ROS2 bag 在线回放使用。 |
| `run_slam_offline.sh` | 运行完整离线 SLAM。 |
| `run_slam_online.sh` | 运行完整在线 SLAM。 |
| `run_loc_offline.sh` | 使用已有地图进行离线定位。 |
| `run_loc_online.sh` | 使用已有地图进行在线定位。 |
| `save_default_map.sh` | 调用保存地图服务。 |
| `install_dep.sh` | 安装 Ubuntu 22.04 下的基础依赖。 |

### run_frontend_offline.sh

run_frontend_offline.sh 输出文件说明：

* `results/trajectory_lidar114.tum`：主 LiDAR 位姿轨迹，脚本主要拿它做完整性检查、时间戳单调性检查、轨迹间隔检查。

* `results/trajectory_imu.tum`：IMU 坐标系轨迹。

* `results/trajectory_rear_axle.tum`：后轴坐标系轨迹，通常用于车辆轨迹评估或和车体参考点对齐。

* `results/map_lio.pcd`：LIO 建出的点云地图。

* `results/frame_stats.csv`：每个融合帧的统计，包括时间范围、是否 partial、融合点数、各雷达是否到齐、每个雷达点数。

* `run_metadata.txt`：这次运行的核心摘要，包括 bag/config/binary 哈希、输出路径、验收结果、最后轨迹时间、资源参数等。

  `run_metadata.txt` 是这次离线运行的“总账本”：它把输入数据、配置、二进制版本、输出文件、资源限制、验收结果都记录下来，方便复现实验和定位失败原因。

  **重点字段**

  | 字段                                             | 作用                                                         |
  | ------------------------------------------------ | ------------------------------------------------------------ |
  | `method=lightning_lm`                            | 标记本次运行的方法/算法名。                                  |
  | `sequence=SANY_4lidar`                           | 数据序列标签，是你命令里传的 `--sequence`，用于区分是哪组数据。 |
  | `repeat=1`                                       | 重复实验编号，是你命令里传的 `--repeat`，用于多次重复运行时区分结果。 |
  | `bag=...`                                        | 输入 rosbag 路径。                                           |
  | `config=...`                                     | 本次使用的 YAML 配置文件路径。                               |
  | `primary_lidar_topic=/livox/lidar_192_168_1_114` | 脚本认为用于验收的主 LiDAR topic。                           |
  | `sensor_duration_s=120.222785744`                | 这包数据的传感器时间跨度，约 120.22 秒。                     |
  | `expected_last_lidar_end_s=1782898362.7087536`   | 主 LiDAR 在 bag 里的最后一帧时间，也是脚本期望轨迹至少接近到达的终点。 |
  | `completion_tolerance_s=0.25`                    | 允许轨迹尾部比主 LiDAR 最后一帧早多少秒，默认 0.25 秒。      |
  | `completion=incomplete`                          | 完整性验收结果。这次表示输出轨迹没有跑到主 LiDAR 末尾附近。  |
  | `algorithm_rc=0`                                 | 算法进程返回码。`0` 表示程序本身正常退出。                   |
  | `watchdog_status=completed`                      | watchdog 结果。`completed` 表示没有超时、没有被强杀。        |
  | `trajectory_lines=1083`                          | 输出的主 LiDAR 轨迹行数，也就是有效位姿数量。                |
  | `last_stamp=1782898360.9002664`                  | 输出轨迹最后一个时间戳。失败主要就是因为它距离 `expected_last_lidar_end_s` 还差约 1.81 秒。 |
  | `invalid_count=0`                                | 轨迹格式/数值非法的行数。这里为 0，说明轨迹格式没问题。      |
  | `nonmonotonic_count=0`                           | 时间戳非递增的行数。这里为 0，说明轨迹时间顺序正常。         |
  | `maximum_output_gap_s=0.100343...`               | 输出轨迹中相邻两帧最大时间间隔。                             |
  | `excessive_output_gap_count=0`                   | 超过允许间隔的 gap 数量。这里为 0，说明中间没有异常大断帧。  |

* `bag_contract.json`：脚本从 rosbag 和配置里推导出的“应跑到哪里”的合同，比如主雷达 topic、第一帧/最后一帧时间、消息数量。

* `logs/algorithm.stderr.log`：算法日志，虽然叫 stderr，但 glog 的 INFO 也在这里。

* `results/processing_timing.csv` 和 `results/processing_timing_summary.json`：处理耗时明细。CSV 便于横向比较各算法阶段；JSON 同时记录整次离线运行的墙钟耗时、输出轨迹吞吐率，以及完整处理全包时的实时因子和处理倍速。

  阶段统计目前包括离线初始化、bag 回放、尾部冲刷、地图导出，以及预处理、点云去畸变、LiDAR 观测匹配、增量建图和 IVox 插点。字段包括 `average_ms`、`median_ms`、`p95_ms`、`retained_samples`。底层计时器每个阶段最多保留最近 2000 个样本；`sample_window_limit_reached=true` 表示已经达到该上限，统计值可能只代表尾部窗口。不同阶段可能嵌套或重叠，例如 bag 回放包含逐帧算法阶段，因此 `retained_window_estimated_total_ms` 不能相加后当作整次墙钟耗时。

  只有完整处理到 bag 主 LiDAR 末帧、进程和 watchdog 正常结束且未等待 UI 时，JSON 才计算 `realtime_factor`（墙钟耗时/传感器时长）和 `processing_speed_x`（传感器时长/墙钟耗时）。限帧、尾部不完整或 `--wait-ui true` 的运行会把不可靠的性能字段置为 `null`，并通过 `performance_ratio_suppressed_reasons` 说明原因。`playback_rate>0` 表示墙钟指标包含主动限速，比较算法最大吞吐能力时应使用 `--playback-rate 0`。

* `resource_samples.csv` 和 `resource_summary.json`：CPU/RSS 资源监控采样和汇总。

   `resource_summary.json` 相关字段：

  | 字段                         | 作用                      | 例子         |
  | ---------------------------- | ------------------------- | ------------ |
  | `samples`                    | 资源采样次数              | `677` 次     |
  | `duration_s`                 | 资源监控持续时间          | `139.06s`    |
  | `mean_cpu_cores`             | 平均用了多少个 CPU 核     | `1.24` 核    |
  | `peak_cpu_cores`             | 峰值用了多少个 CPU 核     | `2.62` 核    |
  | `p95_cpu_cores`              | 95 分位 CPU 使用核数      | `1.95` 核    |
  | `mean_cpu_pct_of_allocation` | 相对分配 CPU 的平均占用率 | `15.52%`     |
  | `peak_cpu_pct_of_allocation` | 相对分配 CPU 的峰值占用率 | `32.81%`     |
  | `mean_rss_mb`                | 平均内存占用 RSS          | `997.42 MB`  |
  | `peak_rss_mb`                | 峰值内存占用 RSS          | `2990.04 MB` |

  `resource_samples.csv` 相关字段：

  | 字段                    | 作用                             |
  | ----------------------- | -------------------------------- |
  | `elapsed_s`             | 从运行开始到该采样点经过了多少秒 |
  | `processes`             | 当前被监控到的进程数量           |
  | `cpu_cores`             | 该时刻等效用了多少个 CPU 核      |
  | `cpu_pct_of_allocation` | 相对分配 CPU 的占用百分比        |
  | `rss_mb`                | 该时刻 RSS 内存占用，单位 MB     |

* `watchdog_status.json`：watchdog 是否超时。

  以一下内容为例：

  ```
  {
    "margin_s": 300.0,
    "playback_rate": 0.0,
    "sensor_duration_s": 120.222785744,
    "status": "completed",
    "timeout_s": 421.0,
    "updated_wall_ns": 1783925136865566082
  }
  ```

  | 字段                | 含义                                          | 这次的值说明                                                 |
  | ------------------- | --------------------------------------------- | ------------------------------------------------------------ |
  | `margin_s`          | watchdog 额外宽限时间，单位秒                 | `300.0`，在数据时长之外额外允许跑 300 秒                     |
  | `playback_rate`     | 回放倍率，用来估算 watchdog 超时时间          | `0.0`，表示离线尽快处理，不按传感器时间节奏慢速回放          |
  | `sensor_duration_s` | bag 里传感器数据的时间跨度，单位秒            | `120.222785744`，这包数据约 120.22 秒                        |
  | `status`            | watchdog 对进程的最终判断                     | `completed`，说明算法在超时前正常结束，没有被 watchdog 杀掉  |
  | `timeout_s`         | watchdog 给本次运行设置的最大墙钟时间，单位秒 | `421.0`，约等于 `120.22 + 300` 后取整                        |
  | `updated_wall_ns`   | watchdog 最后更新状态的墙钟时间戳，单位纳秒   | `1783925136865566082`，主要用于机器记录/排查日志时间，不影响算法结果 |







## 公共辅助工具

| 脚本 | 用途 |
|---|---|
| `inspect_rosbag2_sqlite.py` | 读取 SQLite3 ROS2 bag 和 YAML，核对 Topic 合同、主 LiDAR 末帧、数据时长及输入哈希。由离线前端入口调用。 |
| `monitor_process_tree.py` | 采样算法进程树的 CPU 和 RSS，生成 `resource_samples.csv` 与 `resource_summary.json`。 |
| `extract_frontend_timing.py` | 从离线前端日志提取阶段耗时，生成 `processing_timing.csv` 与 `processing_timing_summary.json`。 |
| `analyze_frontend_topic_bag.py` | 审计在线前端录制的五个公开 Topic，包括消息类型、字段、频率、时间戳和 `lidar_id`。 |
| `plot_posres_statistics.py` | 可视化 `/PosRes` 接收间隔、端到端时延和每秒接收数量。输入由 `pos_res_recorder` 生成的两个 CSV。 |

## 公共复现工具

| 脚本 | 用途 |
|---|---|
| `reproduction/common/build_reproduction_manifest.py` | 递归记录复现目录中每个文件的相对路径、大小和 SHA-256。 |

## SANY 四雷达复现工具

目录：`reproduction/multi_lidar/sany_4livox/`

| 脚本 | 用途 |
|---|---|
| `run_sany_formal_matrix.py` | 运行历史正式 7×3 随机化矩阵，生成实验指纹、attempt 和完成清单。日常复现优先使用根目录的批量入口。 |
| `validate_sany_matrix_posthoc.py` | 对正式矩阵的轨迹、组帧统计、地图字段、来源 ID 和证据文件做独立严格审计。 |
| `evaluate_sany_formal_matrix.py` | 将矩阵轨迹与 114 Voxel-SLAM 代理轨迹比较；代理轨迹不是独立真值。 |
| `analyze_sany_map_cross_source.py` | 对单个六字段 PCD 做跨 LiDAR 平面残差、覆盖率和来源 ID 检查。 |
| `evaluate_sany_map_matrix.py` | 批量调用地图分析器并汇总 21 个正式矩阵单元。 |
| `plot_sany_formal_results.py` | 根据轨迹和地图评价 JSON 生成正式实验图。 |
| `run_sany_online_contract.sh` | 以 1× 回放运行在线前端，录制并验收五个公开 Topic。 |
| `create_sany_fault_bag.py` | 从 SQLite3 ROS2 bag 生成单个确定性丢包/中断故障包及故障合同。 |
| `create_sany_formal_fault_bags.sh` | 一次生成正式定义的四个故障输入包。 |
| `run_sany_formal_fault_matrix.sh` | 顺序执行四个故障场景，并调用在线合同脚本验收。 |
| `analyze_sany_fault_recovery.py` | 对故障场景的初始化、降级、恢复和公开 Topic 连续性进行审计。 |
| `summarize_sany_final_evidence.py` | 汇总原正式实验的冻结证据目录；它绑定历史 attempt/final/retry 布局，不用于新批次目录。 |

`create_sany_fault_bag.py` 依赖 ROS 2 Python 模块，直接调用前需先执行 `source /opt/ros/humble/setup.bash`；两个故障批处理 shell 会在其运行流程中加载相同环境。

## 常用命令

单雷达使用统一单次入口：

```bash
bash scripts/run_frontend_offline.sh \
  --bag /path/to/rosbag2_directory \
  --config config/reproduction/single_lidar/m3dgr/lightning_m3dgr_mid360_benchmark.yaml \
  --output-dir /path/to/new_run_directory
```

SANY 七配置批量运行：

```bash
bash scripts/run_frontend_offline_batch.sh \
  --bag /path/to/sany_rosbag2_directory \
  --output-root /path/to/new_batch_directory \
  --repeats 3 -- --wait-ui false
```

所有输出目录都应预先不存在或不含同名运行产物。正式实验前执行 `git status --short`，确保工作树干净。

### 记录并可视化在线 `/PosRes`

`run_sany_online_diagnostics.sh` 会在运行期间实时写入三份定位轨迹：

* `results/trajectory_global.tum`：PGO 全局校正轨迹。

* `results/trajectory_high_frequency.tum`：定位器内部的高频原始位姿。

* `results/trajectory_published_rear_axle.tum`：经过后轴外参和固定地图变换后，实际发布到
  `/PosRes` 的位姿。该文件逐条刷新，可供下游联调时实时读取。

诊断入口 `run_sany_online_diagnostics.sh` 默认使用 `LIGHTNING_LM_RUN_MODE=diagnostic` 和
`LIGHTNING_LM_COMPUTE_PROFILE=1`，将 LIO、NDT/重定位、PGO、IMU/DR、ROS
发布和诊断 I/O 等热点的墙钟/线程 CPU/进程 CPU 计时写入算法 stderr，并在退出后提取到
`results/compute_profile.log`。高频链路按一秒窗口输出 count、mean、P50/P95/P99/max；LiDAR
帧链路逐帧输出。`run_metadata.txt` 会记录开关值和计时记录条数。

现场统一从 `scripts/run.sh` 启动；地图、雷达布局、CAN、bag、诊断计时和减负模式都固定在这一个
入口中，新增运行选项时直接更新该文件，不再增加模式包装脚本。当前现场基线关闭 bag，保持
`LIGHTNING_LM_RUN_MODE=diagnostic`、`LIGHTNING_LM_COMPUTE_PROFILE=1` 和
`LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD=0`。需要生产减负时也只修改该入口，将运行模式设为
`production`，并按现场验证结果关闭计时、开启非必要开销裁剪。

`LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD` 只裁剪诊断输出；用于下次启动恢复的
`recover_pose` 始终更新，不受该开关影响。

SANY 正式定位 YAML 通过统一的 `compute_budget` 限制 LIO、NDT 和 SOLiD ICP 工作池。域控调参时可
用 `LIGHTNING_LM_LIO_THREADS`、`LIGHTNING_LM_NDT_THREADS`、
`LIGHTNING_LM_SOLID_ICP_WORKERS` 临时覆盖，并用 `LIGHTNING_LM_CPU_AFFINITY=0-9` 将算法限制到
已为定位预留的核心。实际核心列表必须按域控的 IRQ、驱动和其他进程分配填写；脚本会在启动前
校验 SOLiD worker 核心集合是定位进程有效核心集合的子集，不满足时拒绝启动。

构建并加载工作空间后启动接收节点；按 `Ctrl-C` 停止时会保留最后一个不足一秒的统计窗口：

```bash
ros2 run lightning pos_res_recorder --ros-args \
  -p topic:=/PosRes \
  -p output_dir:=./posres_record
```

输出包括 `posres_trajectory.tum`（可由 evo 读取）、逐消息的 `posres_timing.csv`，以及一秒窗口的
`posres_rate.csv`。生成统计图：

```bash
python3 scripts/plot_posres_statistics.py \
  --timing ./posres_record/posres_timing.csv \
  --rate ./posres_record/posres_rate.csv \
  --output ./posres_record/posres_statistics.png
```

`latency_ms` 是接收节点 ROS 时钟减消息 `header.stamp`；两台机器运行发布端和接收端时，应先做
NTP/PTP 时钟同步。`interarrival_receive_ms` 使用单调时钟，不受系统时间校准影响。
