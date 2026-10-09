@page map_alignment 新旧地图对齐

# SANY 新旧地图一次性对齐

## 目标与坐标约定

定位器始终加载最新四雷达地图及其 SOLiD 数据库，内部定位、NDT、PGO 和航迹推算均保持在新地图坐标系中。离线工具只计算并固化一个六自由度变换：

```text
T_old_body = T_old_new * T_new_body
```

其中 `T_old_new` 在 YAML 中使用 `target_from_localization` 约定。在线发布前左乘该变换，使外部系统继续收到老地图坐标系中的结果。对已有且不包含 `output.fixed_map_transform` 的配置，功能默认关闭，行为与原来一致。

## 构建

```bash
source /opt/ros/humble/setup.bash
MAKEFLAGS="-j4" colcon build --packages-up-to lightning_lm --cmake-args -DCMAKE_BUILD_TYPE=Release
source install/setup.bash
```

## 一次性对齐

在 WSL Ubuntu 22.04 的仓库根目录执行：

```bash
bash scripts/align_maps_offline.sh \
  --old_map=runs/sany_3lidar_mapping_blind5_ground_z0_btc_20260725/data/new_map \
  --new_map=runs/sany_4lidar_mapping_data1_20260811_full/data/new_map \
  --input_bag=/mnt/f/datasets/SANY/sites/zhugong_mixing_plant/raw/rig_4lidar/mapping/data1 \
  --config=config/reproduction/multi_lidar/sany_4livox/sany_4lidar_localization_solid.yaml \
  --output_dir=runs/sany_map_alignment_3lidar_old_4lidar_new_20260813
```

工具从新地图 BTC 关键帧序列的起始段自动检测静止窗口，默认静止门限为相对首帧平移不超过 0.15 m、旋转不超过 2 度，并在检测到连续运动后停止取样。输入 bag 的 metadata 用于校验数据来源及时间范围。

候选生成优先使用老地图 SOLiD；若老地图尚无 SOLiD 数据库，工具会利用其 BTC 子图生成一次。BTC 和全局 PCA 几何候选作为交叉验证与后备。候选经过多尺度 ICP、全局双向重叠、局部分区覆盖、多个静止子图一致性及歧义检查。

只有全部自动门限通过时，工具才会备份并更新定位 YAML 的 `output.fixed_map_transform`。失败时不改配置，并返回非零退出码。可用 `--update_config=false` 只生成判断材料而不写 YAML。

## 人工复核材料

无论自动判断成功或失败，输出目录都会保留：

- `alignment_report.yaml`：变换、静止窗口、完整指标、门限和 PASS/FAIL。
- `alignment_candidates.csv`：所有收敛候选及其质量，便于检查重复场景歧义。
- `alignment_birdseye.png`：老地图灰色、新地图红色、静止数据蓝色的俯视叠加。
- `alignment_overlay_rgb.pcd`：可在 CloudCompare 或 RViz 中查看的彩色叠加点云。
- `old_map_review.pcd`、`new_map_aligned_review.pcd`、`static_window_aligned.pcd`：分开的降采样复核点云。
- `MANUAL_REVIEW.md`：人工复核提示。
- `localization_config.before.yaml`、`localization_config.after.yaml`：配置修改前后快照；仅成功写入时产生。

## 在线输出范围

### 2026-10-09 宝马展复刻场地地图

正式三雷达配置 `config/reproduction/multi_lidar/sany_3livox/sany_3lidar_localization_solid.yaml`
现已启用从 `bauma_3lidar_20261009_113953` 到旧资产地图
`sany_3lidar_baoma_sim_20260901` 的固定变换。使用该配置时，应把
`SANY_MAP_PATH` 指向新地图包；配置中的变换和地图包必须成对切换。
旧地图目录继续保留，便于回退。PGM 和三维地图文件仍处于新地图内部坐标系。

变换平移 XYZ 为 `[2.34416675567627, -3.11785960197449, -0.0604841262102127]` m，
四元数 XYZW 为
`[0.000469976801548861, -0.00108713562520467, 0.0238052228746426, 0.999715913958474]`。
全局双向 0.5 m / 1 m 重叠率分别为 80.0% / 91.3%；起始静止段的 0.5 m
平均重叠率为 95.4%。8 个静止子图的变换最大差异为 0.07076 m 和 0.10476 度。
全部自动门限通过，用户已复核叠加图。以上是点云配准指标，不是测量真值精度。
完整算法报告见 [2026-10-09 对齐报告](map_alignment/bauma_20261009/alignment_report.yaml)。

新地图来自 `sany-03-master:~/data_gj/rosbag2_2026_10_09-11_39_53`，
数据与新旧地图已按场地归档在
`F:/datasets/SANY/sites/bauma_2026_replica/`。
新 PGM 使用 0.05 m 分辨率、X `[-10.8, 28.7]` m、Y `[-2.8, 22.6]` m，
障碍物高度为归一化地图中的 Z >= 0.3 m。

### 发布结果

启用固定变换后，只改变下列对外发布结果：

- `/slamPoseRaw_topic`
- `/PosRes`
- `/tf`（保持 `map -> base_link`）
- `/LidarDataInL`
- `/localization/path`

这些消息的 `frame_id` 保持现有的 `map`。`/LidarDataInv`、定位器内部状态、内部新地图位姿导出以及离线 `run_loc_offline` 均不应用该变换。
