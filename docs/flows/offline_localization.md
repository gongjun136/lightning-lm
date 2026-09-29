@page offline_localization_flow 离线定位：bag 与地图到可核验轨迹

# 对象组成

@ref run_loc_offline.cc "run_loc_offline.cc" 的 `main` 直接创建 `LaserMapping`、`loc::LidarLoc` 和 `loc::PGO`。它不经过在线 `LocSystem` 的完整 ROS/诊断编排，也不直接实例化在线 `Localization` 的三队列流程。

输入包括 `--input_bag`、`--config`、`--map_path`。地图须包含可用分块索引及配置要求的重定位数据库；`--use_config_initial_pose` 控制是否采用配置初值，`--force_relocalization_frame` 用于显式恢复实验。

# 一帧的执行顺序

1. IMU 回调送入 LIO，并把 DR 状态转交 `LidarLoc::ProcessDR` 和 `PGO::ProcessDR`，再 drain。
2. 单/多 LiDAR 回调预处理、组帧并 drain；`RunDetailed` 返回 `kOutput` 才进入定位输出路径。
3. LIO 状态交给 `LidarLoc::ProcessLO` 和 `PGO::ProcessLidarOdom`；点云交给 `ProcessCloud`。
4. 获取地图匹配结果，传到 `PGO::ProcessLidarLoc`，分别记录原始地图定位和融合结果，并反馈定位健康状态给点数预算控制。
5. 播放结束 `FlushMultiLidar`、drain、`LidarLoc::Finish`，再关闭输出与 UI。

离线主调用顺序确定，不代表所有模块都无后台线程：地图更新仍有自己的生命周期。`--playback_rate=0` 的吞吐结果不能直接证明在线 10 Hz 无丢帧。

# 输出与证据边界

| 选项 | 内容/用途 |
|---|---|
| `--output_tum` | 融合定位轨迹 |
| `--output_lidar_loc_tum` | 原始 LiDAR 地图匹配轨迹 |
| `--output_csv` | 每次定位的状态和匹配统计 |
| `--output_frame_stats_csv` | 多雷达融合帧统计 |
| `--output_bag` | 直接写出 SANY 五类输出消息 |
| `--publish_topics` | 是否创建发布路径，不能用来推断线下文件内容 |

比较轨迹前核对 IMU/后轴参考点、地图变换和时间基准。输入 bag 已有旧定位输出时，只回放传感器输入，避免把旧消息当作新算法结果。

**深入：** @ref localization_module "定位与重定位"、@ref relocalization_theory "全局重定位：检索、几何验证与时间确认"、@ref pose_graph_theory "位姿图、增量求解与高频平滑"、@ref output_contracts "定位状态、车体参考点与 ROS 输出"。在线线程及发布门控另见 @ref online_localization_flow "在线定位端到端流程"。
