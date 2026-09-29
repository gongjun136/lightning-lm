@page map_contracts 地图包、坐标归一化与栅格导出

# 地图的三种职责

LIO 的 IVox 是最近邻查询用局部结构；`TiledMap` 是持久地图的分块加载和定位支撑；G2P5/导航栅格是下游二维表达。它们的分辨率、更新时机和生命周期不同。

| 模块 | 主要操作 | 阅读源码 |
|---|---|---|
| IVox | 插点、体素邻域检索、容量管理 | @ref ivox3d.h "ivox3d.h" |
| TiledMap / Chunk | 转换全图、读取索引、位姿附近加载/卸载 | @ref tiled_map.cc "tiled_map.cc"、@ref tiled_map_chunk.cc "tiled_map_chunk.cc" |
| map_frame | 起点地面估计、统一导出变换、metadata 校验 | @ref map_frame.h "map_frame.h"、@ref map_frame.cc "map_frame.cc" |
| 导航地图导出 | 多帧观测、障碍高度带、射线/栅格输出 | @ref navigation_map_export.cc "navigation_map_export.cc" |
| G2P5 | 关键帧投影、局部/全局栅格、回环后重绘 | @ref g2p5.cc "g2p5.cc"、@ref g2p5_map.cc "g2p5_map.cc" |

# 持久化契约

分块地图使用 `index.txt` 和 chunk PCD；完整导出可含 `global.pcd`、关键帧/帧轨迹、`map_frame.yaml`、导航 `map.pgm/map.yaml`，以及配置启用的 BTC/SOLiD 数据库。不是所有运行都保证产生所有可选文件。

若启用归一化，所有相关数据必须使用同一个变换：
\f$p_{export}=T_{export,slam}p_{slam},\quad T_{export,i}=T_{export,slam}T_{slam,i}.\f$
metadata 的 transform id/矩阵引用用于核对数据库与地图是否来自同一次导出。仅移动点云但保留旧关键帧数据库会制造表面“匹配成功”的错误坐标。

`EstimateStartGroundFrame` 利用起点附近点、高度预期、内点和倾角等约束估计地面；它有明确失败分支，不是任意场景都能可靠归零。先检查地面估计再创建/清理输出目录。

# 二维栅格的几何含义

平面坐标通常按 \f$i=\lfloor(x-x_0)/r\rfloor,j=\lfloor(y-y_0)/r\rfloor\f$ 落格，r 是米/像素。PGM 的图像行方向与世界 y 方向必须结合 map.yaml 的原点/分辨率读取。高度筛选、多帧观测门槛和障碍膨胀是独立参数；不能将空白像素直接等同于已观测自由空间。

导航导出与 G2P5 运行时地图是两条实现路径，修改任一路径时分别核对坐标、未知格、射线起点和回环重绘。外参影响多雷达射线起点，不能都用主雷达原点替代。

# 维护与工具入口

`convert_voxel_slam_map` 转换外部地图；`build_solid_database` 生成检索数据；`align_maps_offline` 求地图间固定变换。对齐输出施加于 ROS 边界还是改写地图文件，必须按 @ref map_alignment "新旧地图对齐" 区分。
`SlamSystem::SaveMap` 会重建已有目标目录，不提供整体事务回滚。定位的活动 tile 由地图线程更新，匹配器替换与使用由 `LidarLoc` 协调；调用者不要并发直接操作底层 chunk。
