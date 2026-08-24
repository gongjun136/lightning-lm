# SANY 四雷达 SOLiD 重定位

## 数据与坐标约定

- 地图输入：`F:\datasets\SANY\4lidar_lm\mapping\data1`。
- 冷启动验证：`F:\datasets\SANY\4lidar_lm\relocalization\data1` 至 `data5`。
- 四个 Mid-360 点云先按主配置中的外参变换到前雷达 `lidar_0` 坐标系，再生成一个融合 SOLiD 描述子。不要为每个雷达独立检索后再投票，否则雷达遮挡和不同视场会产生互相冲突的候选。
- 单个 Mid-360 的标称视场为水平 360°、垂直 -7° 至 52°。外参变换后的四雷达融合云覆盖范围更宽，因此描述子的垂直角范围配置为 -90° 至 90°，并对分箱索引做边界保护。

## 重定位链路

1. 从优化后的地图子图生成 SOLiD range/angle 描述子数据库。
2. 查询端维护 1、3、5、10 帧滚动子图，允许车辆静止或运动后完成初始化。
3. SOLiD 提供地点相似度和航向初值；平移由对应优化子图上的局部 ICP 恢复。
4. 对重复场景并行评估 0°、90°、-90°、180° 四个确定性航向种子，每个地点只保留 ICP fitness 最优解。
5. 候选继续通过现有地图边界、重叠率、重力方向、NDT 和连续两帧一致性检查。描述子/ICP 结果不能绕过这些安全门限。

建图与定位使用独立的正式 YAML，两者都保存在 `config/` 下：

- 建图：`config/reproduction/multi_lidar/sany_4livox/sany_4lidar_mapping.yaml`；
- SOLiD 定位：`config/reproduction/multi_lidar/sany_4livox/sany_4lidar_localization_solid.yaml`。

`runs/` 只保存地图、轨迹、日志、分析结果和临时变体，不是域控正式配置的来源。需要从建图配置重新冻结定位配置时，执行：

```bash
python3 scripts/reproduction/multi_lidar/sany_4livox/prepare_relocalization_config.py \
  --backend solid \
  --compute-backend cpu \
  --top-k 20 \
  --retrieval-pool-size 200 \
  --output config/reproduction/multi_lidar/sany_4livox/sany_4lidar_localization_solid.yaml
```

生成数据库：

```bash
source /opt/ros/humble/setup.bash
source install/setup.bash
./install/lightning/lib/lightning/build_solid_database \
  --config=config/reproduction/multi_lidar/sany_4livox/sany_4lidar_localization_solid.yaml \
  --map_path=runs/sany_4lidar_mapping_data1_20260811_full/data/new_map
```

运行 5×10 冷启动矩阵并分析：

```bash
bash scripts/reproduction/multi_lidar/sany_4livox/run_phase_a_relocalization_matrix.sh \
  /mnt/f/datasets/SANY/4lidar_lm/relocalization \
  config/reproduction/multi_lidar/sany_4livox/sany_4lidar_localization_solid.yaml \
  runs/sany_4lidar_mapping_data1_20260811_full/data/new_map \
  runs/sany_4lidar_relocalization_latency_matrix_20260812/worker12

python3 scripts/reproduction/multi_lidar/sany_4livox/analyze_phase_a_relocalization_matrix.py \
  --runs-root runs/sany_4lidar_relocalization_latency_matrix_20260812/worker12 \
  --reference-runs-root runs/sany_4lidar_relocalization_latency_matrix_20260812/batch8 \
  --output-json runs/sany_4lidar_relocalization_latency_matrix_20260812/analysis/worker12.json \
  --output-csv runs/sany_4lidar_relocalization_latency_matrix_20260812/analysis/worker12.csv
```

## 2026-08-12 耗时优化结论

新口径将 SOLiD 检索、多航向 ICP 和接受前处理纳入累计耗时。早期 `P95≈0.2002 s` 只是从播放起点到接受帧的传感器时间延迟，且旧检索计时在 ICP 前结束；它不是端到端 CPU 计算耗时，不能支持“完整重定位亚秒”的结论。

本轮对 `data1` 至 `data5` 各选取 `0、3、5、7、10、15、20、25、30、35 s` 十个启动偏移，每个配置执行 50 次冷启动。历史时延验收配置保留四航向种子和 8 个候选，并将 ICP 工作线程从 8 提高到 12：

```yaml
relocalization:
  solid:
    icp_batch_size: 8
    icp_workers: 12
    icp_yaw_hypothesis_offsets_deg: [0.0, 90.0, -90.0, 180.0]
```

该 12-worker 数值是独占算力下的历史实验条件。在线部署现由顶层 `compute_budget` 统一控制
LIO、NDT 和 SOLiD，SANY 正式定位配置默认 `solid_icp_workers: 4`，避免重定位与持续跟踪争抢全部
域控核心；需要复现实验时可用 `LIGHTNING_LM_SOLID_ICP_WORKERS=12` 显式覆盖。

| 配置 | 成功率 | 累计处理 P95 | 累计检索/ICP P95 | 结论 |
| --- | ---: | ---: | ---: | --- |
| batch=1 | 48/50 | — | — | data2 的 7 s、10 s 启动失败，召回不可接受 |
| batch=2 | 50/50 | 12.614 s | 10.845 s | 尾延过大 |
| batch=4 | 50/50 | — | — | data4 局部出现约 10.6 s 尾延 |
| batch=8, workers=8 | 50/50 | 5.544 s | 4.732 s | 成功率基线 |
| batch=8, workers=12 | 50/50 | 4.683 s | 3.870 s | 最终 Orin CPU 默认 |
| batch=8, workers=16 | 50/50 | 3.766 s | 2.988 s | 更快，但 P95 峰值 CPU 约 17.65 核，不作默认 |

workers=12 相比 workers=8 将处理 P95 降低约 `15.5%`，检索/ICP P95 降低约 `18.2%`；位姿结果与 workers=8 基线一致，最低地图重叠率为 `0.9868`。资源采样的 P95 峰值约为 `13.26` 个 CPU 核、`752.2 MiB` RSS。

已否决的方向包括：渐进式 batch 4→8 使困难样本恶化到约 13.3 s；单航向在 data2/7 需 16 次尝试、22.3 s；将点云下采样从 0.20 m 改为 0.30 m 在两个困难点反而慢约 7.8%–8.9%。因此本轮的主要优化是“保留召回候选与四航向安全覆盖，只提高有界并行度”，而不是通过减少候选或关闭几何门限换取耗时。

## Orin 构建边界

- 目标基线：Ubuntu 22.04、ROS 2 Humble、Jetson Linux/L4T R36.4.3、CUDA 12.6。
- 当前验证路径为标准 C++/Eigen/PCL CPU 实现，不包含 x86 专用指令，也不要求 CUDA 参与编译。`compute_backend: auto` 当前选择 CPU。
- `compute_backend: gpu` 是保留值；未编译 CUDA 后端时初始化会明确失败，不会静默回落到 CPU。
- 建议先在 Orin 上执行 Release 构建和单元测试，再用同一地图包做至少一次完整 data1-data5 回放。由于本机不是 aarch64，不能把本机通过等同于已经完成 Orin 实机构建。

```bash
source /opt/ros/humble/setup.bash
colcon build --packages-select lightning --cmake-args -DCMAKE_BUILD_TYPE=Release
source install/setup.bash
./bin/solid_descriptor_test
```
