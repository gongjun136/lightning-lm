# Lightning-LM 脚本说明

仓库沿用 `scripts/` 目录。根目录脚本是稳定运行入口和公共辅助工具；`reproduction/` 保存已交付实验的复现、审计和评价脚本。

## 稳定运行入口

| 脚本 | 用途 |
|---|---|
| `run_frontend_offline.sh` | 单个 SQLite3 ROS2 bag 的统一离线前端入口。单雷达/多雷达由 YAML 控制，并输出轨迹、地图、帧统计、资源监控和验收元数据。 |
| `run_frontend_offline_batch.sh` | 对同一个 bag 依次运行多个 YAML；默认运行 SANY C0—C6，可用 `--repeats` 重复实验。每个子任务仍调用 `run_frontend_offline.sh`。 |
| `run_frontend_online.sh` | 启动在线前端节点，供实时传感器或 ROS2 bag 在线回放使用。 |
| `run_slam_offline.sh` | 运行完整离线 SLAM。 |
| `run_slam_online.sh` | 运行完整在线 SLAM。 |
| `run_loc_offline.sh` | 使用已有地图进行离线定位。 |
| `run_loc_online.sh` | 使用已有地图进行在线定位。 |
| `save_default_map.sh` | 调用保存地图服务。 |
| `install_dep.sh` | 安装 Ubuntu 22.04 下的基础依赖。 |

## 公共辅助工具

| 脚本 | 用途 |
|---|---|
| `inspect_rosbag2_sqlite.py` | 读取 SQLite3 ROS2 bag 和 YAML，核对 Topic 合同、主 LiDAR 末帧、数据时长及输入哈希。由离线前端入口调用。 |
| `monitor_process_tree.py` | 采样算法进程树的 CPU 和 RSS，生成 `resource_samples.csv` 与 `resource_summary.json`。 |
| `analyze_frontend_topic_bag.py` | 审计在线前端录制的五个公开 Topic，包括消息类型、字段、频率、时间戳和 `lidar_id`。 |

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
