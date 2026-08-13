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
colcon build --packages-select lightning --cmake-args -DCMAKE_BUILD_TYPE=Release
source install/setup.bash
```

## 一次性对齐

在 WSL Ubuntu 22.04 的仓库根目录执行：

```bash
bash scripts/align_maps_offline.sh \
  --old_map=runs/sany_3lidar_mapping_blind5_ground_z0_btc_20260725/data/new_map \
  --new_map=runs/sany_4lidar_mapping_data1_20260811_full/data/new_map \
  --input_bag=/mnt/f/datasets/SANY/4lidar_lm/mapping/data1 \
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

启用固定变换后，只改变下列对外发布结果：

- `/slamPoseRaw_topic`
- `/PosRes`
- `/tf`（保持 `map -> base_link`）
- `/LidarDataInL`
- `/localization/path`

这些消息的 `frame_id` 保持现有的 `map`。`/LidarDataInv`、定位器内部状态、内部新地图位姿导出以及离线 `run_loc_offline` 均不应用该变换。

## 本次部署结果（2026-08-13）

自动对齐通过，写入配置的 `T_old_new` 为：

```yaml
translation_xyz: [-2.41949081420898, 18.3118782043457, -0.130619361996651]
quaternion_xyzw: [0.00079188585717829, -0.0124272375022884, -0.305301087578141, 0.952174449672147]
```

本次静止窗口为 16.0001 s，共验证 8 个子图；全局双向 1 m 重叠率为 0.9324，静止子图平均 1 m 重叠率为 0.9896，静止窗口间最大变换离散为 0.0444 m / 0.2584 度，未发现近质量的异解。
