@page map_alignment 新旧地图对齐

# SANY 新旧地图一次性对齐

## 目标与坐标约定

定位器始终加载最新四雷达地图及其 SOLiD 数据库，内部定位、NDT、PGO 和航迹推算均保持在新地图坐标系中。离线工具只计算并固化一个六自由度变换：

```text
T_old_body = T_old_new * T_new_body
```

其中 `T_old_new` 在 YAML 中使用 `target_from_localization` 约定。在线发布前左乘该变换，使外部系统继续收到老地图坐标系中的结果。在线启动以所选地图包的 `map_manifest.yaml` 为准；历史地图缺少清单时须先补齐，不能从共用配置猜测变换。

## 构建

```bash
source /opt/ros/humble/setup.bash
MAKEFLAGS="-j4" colcon build --packages-up-to lightning_lm --cmake-args -DCMAKE_BUILD_TYPE=Release
source install/setup.bash
```

## 一次性对齐

在 WSL Ubuntu 22.04 的仓库根目录执行：

```bash
python3 scripts/resolve_sany_config.py \
  --base-config config/sany/localization.yaml \
  --vehicle-config /absolute/vehicle_config.yaml \
  --output /absolute/new_alignment_run/config.yaml
bash scripts/align_maps_offline.sh \
  --old_map=runs/sany_3lidar_mapping_blind5_ground_z0_btc_20260725/data/new_map \
  --new_map=runs/sany_4lidar_mapping_data1_20260811_full/data/new_map \
  --input_bag=/mnt/f/datasets/SANY/sites/zhugong_mixing_plant/raw/rig_4lidar/mapping/data1 \
  --config=/absolute/new_alignment_run/config.yaml \
  --output_dir=runs/sany_map_alignment_3lidar_old_4lidar_new_20260813
```

工具从新地图 BTC 关键帧序列的起始段自动检测静止窗口，默认静止门限为相对首帧平移不超过 0.15 m、旋转不超过 2 度，并在检测到连续运动后停止取样。输入 bag 的 metadata 用于校验数据来源及时间范围。

候选生成优先使用老地图 SOLiD；若老地图尚无 SOLiD 数据库，工具会利用其 BTC 子图生成一次。BTC 和全局 PCA 几何候选作为交叉验证与后备。候选经过多尺度 ICP、全局双向重叠、局部分区覆盖、多个静止子图一致性及歧义检查。

只有全部自动门限通过时，工具才会备份并更新传入运行 YAML 的 `output.fixed_map_transform`，并更新新地图的 `map_manifest.yaml`，记录源/目标地图版本、参考点云 SHA-256 与完整矩阵。保存清单失败返回非零退出码。可用 `--update_config=false` 只生成判断材料，不更新运行 YAML 或地图清单。审核后将配置与对应地图成套发布；不要直接修改共用模板来启用某个场地的变换。

## 人工复核材料

进入自动门限判断后的正常路径，会尝试在输出目录保存以下复核材料；输入、候选生成等更早阶段抛出异常时，不能保证这些文件齐全，应先看 `fatal_error.txt`：

- `alignment_report.yaml`：变换、静止窗口、完整指标、门限和 PASS/FAIL。
- `alignment_candidates.csv`：所有收敛候选及其质量，便于检查重复场景歧义。
- `alignment_birdseye.png`：老地图灰色、新地图红色、静止数据蓝色的俯视叠加。
- `alignment_overlay_rgb.pcd`：可在 CloudCompare 或 RViz 中查看的彩色叠加点云。
- `old_map_review.pcd`、`new_map_aligned_review.pcd`、`static_window_aligned.pcd`：分开的降采样复核点云。
- `MANUAL_REVIEW.md`：人工复核提示。
- `localization_config_before.yaml`、`localization_config_after.yaml`：配置更新前后的快照。前快照在替换前写入；失败时可能只有前快照，不能以它存在证明配置已更新。

@anchor alignment_acceptance_reading
## 最高分候选为什么仍然会被拒绝？

`quality` 用于给已细化候选排序，最高分被选为 `best`；接受条件要另外查看。先读这组独立门限：

@snippet{lineno} app/align_maps_offline.cc alignment-metric-gates

每个 `gate(false, ...)` 都把原因追加到 `failed_gates`。全局重叠、局部分区、静止子图残差及变换离散度分别检查，不会用高总分抵消其中一项失败。候选还要经过歧义检查：

@snippet{lineno} app/align_maps_offline.cc alignment-distinct-ambiguity

按质量降序找到第一个相对最佳解平移超过 1 m 或旋转超过 5° 的不同解，若它的质量距最高分不超过 0.05，就记为歧义。最后再核对支持证据和接受语句：

@snippet{lineno} app/align_maps_offline.cc alignment-support-and-acceptance

描述子支持或更强的几何后备条件只满足其中一个门；它们不能绕过前面已经失败的指标。只有 `failed_gates.empty()` 才通过自动验收。

<b>手算例子：</b>假设其他门限全通过，最佳质量为 0.90，另一个相差 2 m 的解质量为 0.87，则 `0.87 >= 0.90-0.05`，仍因歧义拒绝。另一个独立例子是全局重叠都很高，但 `regional_p10_overlap_10=0.20`，小于 0.25，局部覆盖门仍失败。这些是假设输入的分支推演，不是对某张地图的实测。

<b>边界与停止点：</b>自动 PASS 表示这些代码门限通过，不代表获得了独立测量真值。候选全部细化失败会在进入这里之前抛出异常。读到失败原因集合与 success 赋值即可；变换方向和只应用一次的规则复用 @ref map_alignment_once_reading "地图输出变换导读"。完整入口见 @ref align_maps_offline.cc "align_maps_offline.cc"。

@anchor alignment_persistence_reading
## 关闭配置更新会写什么？配置与清单能一起回滚吗？

配置更新和地图清单写入分别受相同开关控制，先看配置：

@snippet{lineno} app/align_maps_offline.cc alignment-config-write-gate

再看后续清单。前一步若把 success 改为 false，本步骤就会跳过：

@snippet{lineno} app/align_maps_offline.cc alignment-manifest-write-gate

`update_config=false` 会跳过这两块，`config_updated` 与 `map_manifest_updated` 都保持 false；仍会尝试写复核点云、候选表、报告等到输出目录，因此它不是全程只读。接受门限通过时，两个 updated 都为 false 也可以返回 0。

`PatchConfig()` 自身先备份，再写临时文件、解析 YAML、替换原配置，最后保存后快照：

@snippet{lineno} app/align_maps_offline.cc alignment-config-backup

@snippet{lineno} app/align_maps_offline.cc alignment-config-replace

<b>故障推演：</b>假设全部自动门限通过、开关为 true、PatchConfig 完成，但随后 SaveAlignmentManifest 失败，那么配置已经替换，报告记录 `config_updated=true`、`map_manifest_updated=false`、`success=false`，正常错误路径返回 1。这段没有恢复旧配置的调用。甚至后快照复制失败也发生在 rename 之后，PatchConfig 返回 false 不能一概理解为原配置从未变化。

<b>边界与停止点：</b>临时文件替换只涉及单份 YAML，不是配置、清单和复核产物的整体事务；部分复核 writer/stream 也未检查写入结果，退出 0 不能单独证明全部材料完整。更早或抛出的异常走外层捕获并返回 2。读到两个写入条件、rename 的位置和失败返回即可；在线实际使用哪份变换继续 @ref map_manifest_validation_reading "地图启动校验"。

## 在线输出范围

### 2026-10-09 宝马展复刻场地地图

2026-10-09 的历史三雷达配置启用了从 `bauma_3lidar_20261009_113953` 到旧资产地图
`sany_3lidar_baoma_sim_20260901` 的固定变换。2026-10-11 共用配置
`config/sany/localization.yaml` 保留该数值，按用户要求默认关闭。当前若恢复该场地变换，应先为对应新地图包生成并校验审核后的对齐清单，再把
`SANY_MAP_PATH` 指向该地图包；在线启动从地图清单加载变换，仅启用运行 YAML 的历史矩阵并不足够，见 @ref map_alignment_once_reading "变换的实际来源"。
旧地图目录继续保留，便于回退。PGM 和三维地图文件仍处于新地图内部坐标系。

变换平移 XYZ 为 `[2.34416675567627, -3.11785960197449, -0.0604841262102127]` m，
四元数 XYZW 为
`[0.000469976801548861, -0.00108713562520467, 0.0238052228746426, 0.999715913958474]`。
全局双向 0.5 m / 1 m 重叠率分别为 80.0% / 91.3%；起始静止段的 0.5 m
平均重叠率为 95.4%。8 个静止子图的变换最大差异为 0.07076 m 和 0.10476 度。
全部自动门限通过，用户已复核叠加图。以上是点云配准指标，不是测量真值精度。
完整算法报告见 <a href="alignment_report.yaml">2026-10-09 对齐报告</a>。

新地图来自 `sany-03-master:~/data_gj/rosbag2_2026_10_09-11_39_53`，
数据与新旧地图已按场地归档在
`F:/datasets/SANY/sites/bauma_2026_replica/`。
新 PGM 使用 0.05 m 分辨率、X `[-10.8, 28.7]` m、Y `[-2.8, 22.6]` m，
障碍物高度为归一化地图中的 Z >= 0.3 m。

### 发布结果

启用固定变换后，只改变下列对外发布结果：

- `/slamPoseRaw_topic`
- `/localization/pose_vel`
- `/tf`（保持 `map -> base_link`）
- `/LidarDataInL`
- `/localization/path`

这些地图消息的 `frame_id` 保持现有的 `map`。`/LidarDataInv`、定位器内部状态、内部新地图位姿导出以及离线 `run_loc_offline` 的原生轨迹均不应用该变换。实际 publisher 与参考点转换见 @ref output_pose_reading "业务输出的源码导读"；历史 `/PosRes` 名称不能代替当前创建的 `/localization/pose_vel` publisher。
