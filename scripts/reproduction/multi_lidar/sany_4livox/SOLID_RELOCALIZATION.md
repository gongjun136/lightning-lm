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

生产配置默认仍为 `relocalization.backend: btc`，便于回退。使用 SOLiD 时应通过脚本冻结独立配置：

```bash
python3 scripts/reproduction/multi_lidar/sany_4livox/prepare_relocalization_config.py \
  --backend solid \
  --compute-backend cpu \
  --top-k 20 \
  --retrieval-pool-size 200 \
  --output runs/sany_4lidar_20260803_validation/configs/sany_4lidar_solid_cpu_multiyaw.yaml
```

生成数据库：

```bash
source /opt/ros/humble/setup.bash
source install/setup.bash
./install/lightning/lib/lightning/build_solid_database \
  --config=runs/sany_4lidar_20260803_validation/configs/sany_4lidar_solid_cpu_multiyaw.yaml \
  --map_path=runs/sany_4lidar_20260803_validation/mapping_data1_phase_a/data/new_map
```

运行 5×10 冷启动矩阵并分析：

```bash
bash scripts/reproduction/multi_lidar/sany_4livox/run_phase_a_relocalization_matrix.sh \
  /mnt/f/datasets/SANY/4lidar_lm/relocalization \
  runs/sany_4lidar_20260803_validation/configs/sany_4lidar_solid_cpu_multiyaw.yaml \
  runs/sany_4lidar_20260803_validation/mapping_data1_phase_a/data/new_map \
  runs/sany_4lidar_20260803_validation/solid_relocalization_matrix_cpu_multiyaw_final

python3 scripts/reproduction/multi_lidar/sany_4livox/analyze_phase_a_relocalization_matrix.py \
  --runs-root runs/sany_4lidar_20260803_validation/solid_relocalization_matrix_cpu_multiyaw_final \
  --reference-runs-root runs/sany_4lidar_20260803_validation/relocalization_matrix_strict \
  --output-json runs/sany_4lidar_20260803_validation/solid_relocalization_matrix_cpu_multiyaw_final/acceptance.json \
  --output-csv runs/sany_4lidar_20260803_validation/solid_relocalization_matrix_cpu_multiyaw_final/acceptance_trials.csv
```

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
