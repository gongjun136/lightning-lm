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

下面四个问题按“启动接受什么地图 → 输出变换如何选择 → 导出几何如何变换 → 运行时如何选 tile”阅读。例子用于手工核对分支，不代表已经运行的地图测试。

@anchor map_manifest_validation_reading
# 1. 地图能打开，为什么定位仍拒绝启动？

`LocSystem::Init()` 通过 `LoadOutputConfig()` 读取所选地图的清单，再使用清单中的输出变换。先看地图格式和状态这两个门：

@snippet{lineno} core/maps/map_manifest.cc map-manifest-startup-gates

<b>先读 `ReadAndVerify`，再读两个状态判断</b>：有 `global.pcd` 不等于有可定位的分块地图；`map_format` 必须为 `tiled`。`unknown` 明确拒绝，只有 `aligned`、`not_required` 进入后续校验。异常最终转为调用者可处理的失败结果：

@snippet{lineno} core/maps/map_manifest.cc map-manifest-startup-result

再核对 `ReadAndVerify()` 怎样防止“清单还在，但点云已换过”。下面两遍检查分别覆盖已登记文件和目录中新出现的几何文件：

@snippet{lineno} core/maps/map_manifest.cc map-manifest-file-verification

第一遍的 `FileSha256(...) != ...` 检查实际文件内容；文件无法读取也会抛出异常。第二遍的 `!files[...]` 拒绝未登记的新文件。这里的“几何文件”由以下条件限定：

@snippet{lineno} core/maps/map_manifest.cc map-manifest-geometry-scope

<b>状态例子</b>：清单登记了 `index.txt` 和 `7.pcd`，且状态为 `not_required`。改动 `7.pcd` 会在第一遍失败；不改清单而新增 `8.pcd` 会在第二遍失败；将状态改为 `unknown`，即使哈希都相同也会失败。

<b>边界与停止点</b>：检查范围是根目录的 `global.pcd`、`index.txt`、`map_frame.yaml` 和纯数字文件名的 `.pcd`；`7_dyn.pcd`、轨迹和数据库子目录不属于这组哈希检查。清单缺失、schema/身份不合法、矩阵与四元数不一致、坐标声明矛盾还有各自失败分支。只核对“为何此图不能启动”，读到对应报错和 `return false` 即可；哈希通过不能证明地图与现场相符。纯 INS 和离线原生轨迹输出保留各自流程。

完整入口：@ref lightning::LocSystem::Init() "LocSystem::Init()" → @ref lightning::map_manifest::LoadOutputConfig() "LoadOutputConfig()"；校验辅助函数在 @ref map_manifest.cc "map_manifest.cc"。

@anchor map_alignment_once_reading
# 2. 新图到老图的变换为什么有时启用、有时禁用？

关键不是历史矩阵是否存在，而是它是否已应用到地图文件。启动时只用以下表达式决定输出端是否还要应用：

@snippet{lineno} core/maps/map_manifest.cc map-manifest-alignment-once

`enabled` 需要同时满足 `status == "aligned"` 和 `!applied`；`needs_transform_to_old` 必须与这个计算结果相同。实际输出配置由 `OutputNode()` 生成：

@snippet{lineno} core/maps/map_manifest.cc map-manifest-effective-output

<b>逐句核对</b>：`enabled ? transform : SE3()` 决定写入原变换还是单位变换，随后还写入 `fixed["enabled"]`。因此禁用不是保留一个等待外部猜测是否应用的非单位输出矩阵。源/目标名称固定为 `localization_map`→`map`。

| `legacy_alignment.status` | `needs_transform_to_old` | 在线定位行为 |
|---|---|---|
| `unknown` | null | 拒绝启动；不得把未知状态当作无需变换 |
| `not_required` | false | 明确使用该地图原生坐标，输出单位变换 |
| `aligned`，`applied_to_map_files=false` | true | 只在业务输出边界左乘 `T_old_new` |
| `aligned`，`applied_to_map_files=true` | false | 几何已经转换，保留历史矩阵但禁用再次变换 |

<b>手算例子</b>：假设 `T_old_new` 只有 x 方向 `+100 m` 平移。地图文件仍在新图坐标中，定位得到 `x=3 m`，则 `applied=false`，输出变为 `103 m`。若点云等几何已物理转换到老图坐标，定位本身得到 `103 m`，则 `applied=true`，输出仍为 `103 m`。若把已转换几何整套错误声明为新图坐标，并在输出再次平移，就会得到 `203 m`；只改一个造成字段矛盾的标志，会被启动校验拒绝。

<b>边界与停止点</b>：`SaveAlignmentManifest(..., already_transformed, ...)` 只记录已完成的转换，不会改写点云；`--files_already_transformed=true` 同样是声明，不能代替转换操作。启动还要求 `coordinate_frame_id`、输出位姿和应用状态一致，但不能凭这些字段判断调用者的物理坐标声明是否真实。读完 `enabled` 与单位变换分支即可确认本节；后轴位姿如何左乘该变换见 @ref output_pose_reading "输出坐标核对"。

新图默认 `unknown`；可在建图配置的 `map_export.alignment_status` 中明确声明 `not_required`，并用 `map_export.site_id/map_id` 指定身份。不能为新图声明 `aligned` 来沿用另一张地图的矩阵。`align_maps_offline` 全部门限通过且 `update_config=true` 时写入与源/参考几何绑定的变换；`update_config=false` 保持地图清单及配置不变。

完整源码：@ref lightning::map_manifest::LoadOutputConfig() "LoadOutputConfig()"、@ref lightning::map_manifest::SaveAlignmentManifest() "SaveAlignmentManifest()"、@ref map_manifest.cc "map_manifest.cc"。启动读取的是当前地图包的哈希；`target_global_pcd_sha256` 在此仅检查长度，不能将启动成功理解为重新读取并核验了老地图目录。

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

@anchor map_frame_normalization_reading
# 3. 地面归零会怎样改变点、位姿和 PCD 的 VIEWPOINT？

`EstimateStartGroundFrame()` 先从起点附近估计单位法向量 `ground_normal` 和平面偏移 `ground_offset`，满足 `n·p+d=0`。估计通过后，变换由以下语句构造：

@snippet{lineno} core/maps/map_frame.cc map-ground-plane-transform

<b>读旋转与平移两行</b>：`FromTwoVectors(n, UnitZ())` 将地面法向量转到 z 轴；平移是 `(0,0,d)`。因此地面点的新 z 为 `n·p+d=0`。这不是把车辆起点强制放到零高度，倾斜地面还会引起旋转。

点云和位姿分别这样使用同一变换：

@snippet{lineno} core/maps/map_frame.cc map-ground-transform-cloud

@snippet{lineno} core/maps/map_frame.h map-ground-transform-pose

点云循环只赋值 XYZ；之后另外变换内存里的 `sensor_origin_` 和 `sensor_orientation_`。位姿使用左乘，让传感器/关键帧与点云落在同一导出坐标系。所有相关数据必须使用同一个变换：

\f$p_{export}=T_{export,slam}p_{slam},\quad T_{export,i}=T_{export,slam}T_{slam,i}.\f$

<b>手算例子</b>：假设地面水平且位于 SLAM 坐标 `z=-2 m`，则 `n=(0,0,1)`、`d=2 m`，旋转为单位旋转。地面点 `(1,2,-2)` 变成 `(1,2,0)`；位于 `z=0` 的传感器位姿变成 `z=2 m`。`metadata.normalized` 的含义是启用该变换，不是记录某个点云已处理；对同一点云调用两次，地面点会从 `z=-2` 变为 `0`，再变为 `2`。

全局地图、分块地图和对齐预览 PCD 的 XYZ 已在目标坐标系中，写出时使用下面的 helper：

@snippet{lineno} common/pcd_io.h map-world-pcd-viewpoint

`toPCLPointCloud2` 序列化点属性；writer 的最后两个参数显式指定零原点、单位四元数，形成 `VIEWPOINT 0 0 0 1 0 0 0`。它不会再变换 XYZ，也不会修改 `const` 输入点云的传感器元数据。CloudCompare 会在导入时应用 `VIEWPOINT`，所以不能再把刚才的 `+2 m` 放进世界点云的文件头，否则已经位于 `z=0` 的地面又会被移到 `z=2 m`。BTC/SOLiD 的传感器局部子地图仍使用各自坐标约定和常规 PCL writer。

<b>边界与停止点</b>：地面估计有高度预期、内点、倾角和残差门限，可以失败；上述例子以前面估计已成功为条件。metadata 的 transform id/矩阵引用用于核对数据库与地图是否使用相同的导出变换，不能证明它们来自同次导出，也不能自动修复只变换点云、未变换关键帧的混合坐标。读到“点与位姿使用相同左乘、世界 PCD 写单位 VIEWPOINT”即可停止；导出目录何时重建继续读 @ref offline_slam_flow "离线 SLAM 流程"。

完整入口：@ref lightning::map_frame::EstimateStartGroundFrame() "EstimateStartGroundFrame()"、@ref lightning::map_frame::TransformCloudInPlace() "TransformCloudInPlace()"、@ref lightning::map_frame::TransformPose() "TransformPose()"、@ref pcd_io.h "pcd_io.h"。

@anchor map_tile_lifecycle_reading
# 4. tile 离开附近范围后，会释放内存或写回磁盘吗？

先区分“点云已读入”和“被选入当前活动范围”。`LoadMapIndex()` 解析每个静态索引项后就执行：

@snippet{lineno} core/maps/tiled_map.cc map-tile-index-preload

因此静态点云不是等到车辆靠近才首次尝试读取。随后 `LoadOnPose()` 按位姿更新活动集合；若仍在上次的同一网格，会直接返回：

@snippet{lineno} core/maps/tiled_map.cc map-tile-same-grid

发生网格变化后，附近项按以下条件加入 `loaded_chunks_`：

@snippet{lineno} core/maps/tiled_map.cc map-tile-load-radius

<b>先算 `d`，再看 `n`</b>：`d` 是 tile 与当前网格的二维索引差，`n=abs(dx)+abs(dy)` 是曼哈顿距离，单位是格。只有 `n` 不大于 `load_map_size_` 才加入活动集合，点云指针为空时才再次读文件；已分配但没有点的 cloud 不满足这个重读条件。离开范围走另一门限：

@snippet{lineno} core/maps/tiled_map.cc map-tile-active-removal

<b>手算例子</b>：假设 `load_map_size_=2`、`unload_map_size_=3`，偏移 (1,1) 的 tile 距离为 2，可以加入；偏移 (2,1) 的距离为 3，不会新加入，但如果此前已活动就继续保留；偏移 (2,2) 的距离为 4，会被移出活动集合。这是两条不同门限形成的保留区间。

这里执行的是 `loaded_chunks_.erase`；静态 `Unload()` 调用仍被注释，静态点云留在缓存。动态图层还受策略控制：`SHORT` 清空对应点云，`LONG` 在此保留，`PERSISTENT` 进入以下分支：

@snippet{lineno} core/maps/tiled_map.cc map-tile-persistent-unload

`save_dyn_when_unload_` 控制写 `<id>_dyn.pcd`，`delete_when_unload_` 单独控制 `Unload()`。例如当前动态 chunk 已加载且非空、两个开关分别为 `true/false`，移出活动范围会写盘，但仍保留内存；只有保存开关不能证明内存已经释放。

<b>边界与停止点</b>：同网格的 early return 不重新检查这些条件；`PERSISTENT` 的空点云会在保存和释放之前 `continue`。本段 writer 的返回值未被检查，看到调用只能确认尝试写盘。`SaveToBin(false)` 则写静态索引/点云后返回；`SaveToBin(true)` 才写动态图层，它们都不负责更新 manifest。读到“活动集合、内存对象、磁盘文件分别由哪些语句改变”即可停止，不能将它们统称为一次卸载完成。

完整源码：@ref lightning::TiledMap::LoadMapIndex() "LoadMapIndex()"、@ref lightning::TiledMap::LoadOnPose() "LoadOnPose()"、@ref lightning::TiledMap::SaveToBin() "SaveToBin()"、@ref tiled_map_chunk.cc "MapChunk 的读取/释放实现"。定位的退出保存条件见 @ref online_localization_flow "在线定位流程"。

# 二维栅格的几何含义

`georeference_map` 已将几何物理转换到审核后的目标 ENU，因此在 `tiled/` 写入 `not_required` 清单，避免再应用历史业务地图变换；地理转换依据继续保存在上级 `georeference.yaml`。几何 SHA-256 使用 OpenSSL Crypto，构建依赖为 `libssl-dev`。

平面坐标通常按 \f$i=\lfloor(x-x_0)/r\rfloor,j=\lfloor(y-y_0)/r\rfloor\f$ 落格，r 是米/像素。PGM 的图像行方向与世界 y 方向必须结合 map.yaml 的原点/分辨率读取。高度筛选、多帧观测门槛和障碍膨胀是独立参数；不能将空白像素直接等同于已观测自由空间。

导航导出与 G2P5 运行时地图是两条实现路径，修改任一路径时分别核对坐标、未知格、射线起点和回环重绘。外参影响多雷达射线起点，不能都用主雷达原点替代。

# 维护与工具入口

`convert_voxel_slam_map` 转换外部地图；`build_solid_database` 生成检索数据；`align_maps_offline` 求地图间固定变换。对齐输出施加于 ROS 边界还是改写地图文件，必须按 @ref map_alignment "新旧地图对齐" 区分。
`SlamSystem::SaveMap` 会重建已有目标目录，不提供整体事务回滚。定位的活动 tile 由地图线程更新，匹配器替换与使用由 `LidarLoc` 协调；调用者不要并发直接操作底层 chunk。
