@page map_contracts 地图包、坐标归一化与栅格导出

# 地图的三种职责

LIO 的 IVox 是最近邻查询用局部结构；`TiledMap` 是持久地图的分块加载和定位支撑；G2P5/导航栅格是下游二维表达。它们的分辨率、更新时机和生命周期不同。

| 模块 | 主要操作 | 阅读源码 |
|---|---|---|
| IVox | 插点、体素邻域检索、容量管理 | @ref ivox3d.h "ivox3d.h" |
| TiledMap / Chunk | 转换全图、读取索引、位姿附近加载/卸载 | @ref tiled_map.cc "tiled_map.cc"、@ref tiled_map_chunk.cc "tiled_map_chunk.cc" |
| map_frame | 起点地面估计、统一导出变换、metadata 校验 | @ref map_frame.h "map_frame.h"、@ref map_frame.cc "map_frame.cc" |
| map_manifest | 地图版本、文件 SHA-256、业务坐标声明与启动校验 | @ref map_manifest.h "map_manifest.h"、@ref map_manifest.cc "map_manifest.cc" |
| 导航地图导出 | 多帧观测、障碍高度带、射线/栅格输出 | @ref navigation_map_export.cc "navigation_map_export.cc" |
| G2P5 | 关键帧投影、局部/全局栅格、回环后重绘 | @ref g2p5.cc "g2p5.cc"、@ref g2p5_map.cc "g2p5_map.cc" |

# 持久化契约

分块地图使用 `index.txt` 和 chunk PCD；完整导出可含 `global.pcd`、关键帧/帧轨迹、`map_frame.yaml`、导航 `map.pgm/map.yaml`，以及配置启用的 BTC/SOLiD 数据库。不是所有运行都保证产生所有可选文件。

在线/离线 SLAM 和 Voxel-SLAM 转换导出同时写 `map_manifest.yaml`。它记录 `map_id`、`site_id`、实际 `coordinate_frame_id`、根目录几何文件 SHA-256、地面归一化，以及 `legacy_alignment` 与可直接使用的 `output` 配置。`map_frame.yaml` 的 SLAM→导出地图归一化与新→老业务变换分开记录。

| `legacy_alignment.status` | `needs_transform_to_old` | 在线定位行为 |
|---|---|---|
| `unknown` | null | 拒绝启动；不得把未知状态当作无需变换 |
| `not_required` | false | 明确使用该地图原生坐标，输出单位变换 |
| `aligned`，`applied_to_map_files=false` | true | 只在业务输出边界左乘 `T_old_new` |
| `aligned`，`applied_to_map_files=true` | false | 几何已经转换，保留历史矩阵但禁用再次变换 |

新图默认 `unknown`；可在建图配置的 `map_export.alignment_status` 中明确声明 `not_required`，并用 `map_export.site_id/map_id` 指定身份。不能为新图声明 `aligned` 来沿用另一张地图的矩阵。`align_maps_offline` 全部门限通过且 `update_config=true` 时写入与源/参考几何绑定的变换；`update_config=false` 保持地图清单及配置不变。

`LocSystem::Init()` 必须加载所选地图的清单，校验全部登记的几何 SHA-256、矩阵/四元数一致性和应用状态，然后覆盖共用 YAML 的输出变换。清单缺失、状态未知、几何变动、字段矛盾均返回失败；纯 INS 和离线原生轨迹输出保留各自流程。仅有 `global.pcd` 的可视化导出也可有清单，但 `map_format=global_pcd_only`，不能用来启动地图定位。

### 为历史地图补齐或检查清单

`export_map_manifest` 是安装到 `lib/lightning_lm` 的 C++ 工具；参数由 @ref export_map_manifest.cc "main()" 解析。写入采用临时文件后 rename；不会改变点云或数据库。

```bash
# 明确使用原生坐标。
ros2 run lightning_lm export_map_manifest --map=/absolute/map \
  --site_id=new_vehicle_test_ground --status=not_required
# 从已通过的报告导入，检查报告所指的源图和参考图全局点云哈希。
ros2 run lightning_lm export_map_manifest --map=/absolute/new_map \
  --site_id=zhugong_mixing_plant --old_map=/absolute/old_map \
  --alignment_report=/absolute/accepted/alignment_report.yaml
# 只读检查；非零退出说明不能启动定位。
ros2 run lightning_lm export_map_manifest --map=/absolute/map --check_only
```

物理转换过的历史点云额外指定 `--files_already_transformed=true`。`--apply_to_config=/absolute/run/config.yaml` 只在校验成功后更新该次运行快照的地图输出字段和清单哈希；保留车辆标定、参考点和其他算法参数。现场入口在传感器检查之前调用它，错误记录在 `logs/map_manifest.log`。旧 `SANY_ENABLE_FIXED_MAP_TRANSFORM` 已弃用，不覆盖所选地图的清单；现场脚本若仍设置它，会记录忽略提示。

若启用归一化，所有相关数据必须使用同一个变换：
\f$p_{export}=T_{export,slam}p_{slam},\quad T_{export,i}=T_{export,slam}T_{slam,i}.\f$
metadata 的 transform id/矩阵引用用于核对数据库与地图是否来自同一次导出。仅移动点云但保留旧关键帧数据库会制造表面“匹配成功”的错误坐标。

全局地图、分块地图和对齐预览 PCD 的 XYZ 已经位于其目标坐标系中，因此通过 `pcd_io::SaveWorldCloudBinaryCompressed` 写出单位 `VIEWPOINT 0 0 0 1 0 0 0`。CloudCompare 导入 PCD 时会应用 `VIEWPOINT`；若再写入地面归一化或地图对齐位姿，就会造成重复变换。此写出接口保留 XYZ、强度、时间、雷达编号等点属性，以及内存中的传感器元数据。传感器局部坐标系的 BTC/SOLiD 子地图仍由原有写出流程和 manifest 位姿定义其坐标。

`EstimateStartGroundFrame` 利用起点附近点、高度预期、内点和倾角等约束估计地面；它有明确失败分支，不是任意场景都能可靠归零。先检查地面估计再创建/清理输出目录。

# 二维栅格的几何含义

`georeference_map` 已将几何物理转换到审核后的目标 ENU，因此在 `tiled/` 写入 `not_required` 清单，避免再应用历史业务地图变换；地理转换依据继续保存在上级 `georeference.yaml`。几何 SHA-256 使用 OpenSSL Crypto，构建依赖为 `libssl-dev`。

平面坐标通常按 \f$i=\lfloor(x-x_0)/r\rfloor,j=\lfloor(y-y_0)/r\rfloor\f$ 落格，r 是米/像素。PGM 的图像行方向与世界 y 方向必须结合 map.yaml 的原点/分辨率读取。高度筛选、多帧观测门槛和障碍膨胀是独立参数；不能将空白像素直接等同于已观测自由空间。

导航导出与 G2P5 运行时地图是两条实现路径，修改任一路径时分别核对坐标、未知格、射线起点和回环重绘。外参影响多雷达射线起点，不能都用主雷达原点替代。

# 维护与工具入口

`convert_voxel_slam_map` 转换外部地图；`build_solid_database` 生成检索数据；`align_maps_offline` 求地图间固定变换。对齐输出施加于 ROS 边界还是改写地图文件，必须按 @ref map_alignment "新旧地图对齐" 区分。
`SlamSystem::SaveMap` 会重建已有目标目录，不提供整体事务回滚。定位的活动 tile 由地图线程更新，匹配器替换与使用由 `LidarLoc` 协调；调用者不要并发直接操作底层 chunk。
