@page offline_localization_flow 离线定位：bag 与地图到可核验轨迹

# 对象组成

@ref run_loc_offline.cc "run_loc_offline.cc" 的 `main` 直接创建 `LaserMapping`、`loc::LidarLoc` 和 `loc::PGO`。它不经过在线 `LocSystem` 的完整 ROS/诊断编排，也不直接实例化在线 `Localization` 的三队列流程。

输入包括 `--input_bag`、`--config`、`--map_path`。地图须包含可用分块索引及配置要求的重定位数据库；`--use_config_initial_pose` 控制是否采用配置初值，`--force_relocalization_frame` 用于显式恢复实验。

先验证 @ref offline_loc_fused_reading "fused 文件中的位姿从哪里来" 和 @ref offline_loc_success_reading "退出成功究竟证明了什么"。两段直接展示入口中的选择、写入和计数条件；例子是源码条件的手工推演。

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
| `--output_tum` | 优先采用足够新的 PGO 回调结果，否则回退原始地图结果；见下节 |
| `--output_lidar_loc_tum` | 原始 LiDAR 地图匹配轨迹 |
| `--output_csv` | 每次定位的状态和匹配统计 |
| `--output_frame_stats_csv` | 多雷达融合帧统计 |
| `--output_bag` | 直接写位姿、点云及健康/路径消息，各自仍受对应条件控制 |
| `--publish_topics` | 是否创建发布路径，不能用来推断线下文件内容 |

比较轨迹前核对 IMU/后轴参考点、地图变换和时间基准。输入 bag 已有旧定位输出时，只回放传感器输入，避免把旧消息当作新算法结果。

@anchor offline_loc_fused_reading
# fused 轨迹每一行都来自本次 PGO 吗？

先找 `latest_final_result` 的写入者。入口将普通全局回调和高频回调都接到同一个缓存：

@snippet{lineno} app/run_loc_offline.cc offline-loc-capture-pgo

然后看每帧地图匹配后的选择逻辑。重定位接受时，前面的分支会调用 `pgo.Reset()` 并将 `latest_final_result_set` 清为 false；下段调用后才读取缓存：

@snippet{lineno} app/run_loc_offline.cc offline-loc-fused-selection

`final_result` 先取 `loc_result`，只有缓存已设置且时间戳不早于当前地图时间减 1 微秒，才替换成缓存结果。这里没有检查 `ProcessLidarLoc()` 的返回值，也没有要求缓存时间与当前地图时间严格相等。

<b>手算例子：</b>当前地图结果时刻为 10.1 s，PGO 缓存仍为 10.0 s，则保留地图结果；不存在缓存也一样。若调用后缓存达到 10.1 s，才选择它。这个条件只有下界：若缓存为 10.2 s，也会被选中；所以还必须继续看写入时使用的时间戳。

@snippet{lineno} app/run_loc_offline.cc offline-loc-trajectory-output

融合 TUM 使用的是 `final_result.pose_`，而时间参数传入 `loc_result.timestamp_`。因此不能只凭文件名或一行时间戳，断言该行来自时间完全相同的一次成功图优化。两条轨迹偶尔完全重合，也可能只是走了原始地图结果回退分支。

<b>读到这里：</b>能连起“两个回调写缓存 → 默认地图结果 → 时间条件替换 → TUM 参数”即可。真正的 PGO 失败路径见 @ref pgo_frame_acceptance_reading "入图与失败条件"；完整入口见 @ref run_loc_offline.cc "run_loc_offline.cc"。

@anchor offline_loc_success_reading
# 有 TUM 输出或退出码 0，能证明整段定位成功吗？

先核对上段反复调用的判定函数：

@snippet{lineno} app/run_loc_offline.cc offline-loc-usable-result

这里是两个有效标志的 `||`，没有要求 `status_ == GOOD`，也没有检查是否通过 ROS 发布门控。`loc_valid_frames` 数的是满足该条件的<b>原始地图结果</b>，不是 PGO 成功次数，更不是已经发布位姿的消息数。

实际写一行 TUM 还经过下面的函数：

@snippet{lineno} app/run_loc_offline.cc offline-loc-tum-write

它检查流已打开、时间大于 0 且比上一条写入时间新。例如上一条为 10.0 s，本次仍是 10.0 s，即使有效标志为 true 也不会再写一行。此 helper 没有显式的 `isfinite()` 检查；只能从这里确认上述门槛，不能宣称它完成了所有数值合法性检查。

最后的退出码条件更宽：

@snippet{lineno} app/run_loc_offline.cc offline-loc-exit-status

<b>计数例子：</b>假设其他初始化/文件步骤均未提前失败，产生 10 次定位结果，其中仅 1 次满足 `valid_ || lidar_loc_valid_`，最终仍返回 0。即便没有配置 TUM 输出路径，计数也独立进行。退出码 0 证明本程序满足这个结束条件，不能代替持续跟踪成功率、误差评测或业务消息验收。

<b>阅读边界：</b>TUM 写入、ROS 发布和直接写 bag 是不同路径；`--publish_topics=false` 控制 ROS publisher 的创建，配置 `--output_bag` 仍可创建 writer。具体位姿消息还经过 `OfflineLocalizationPublisher::PublishPose()` 的有效性和 publication gate 检查，不能用 TUM 行数推算消息数。完整实现见 @ref run_loc_offline.cc "输出辅助函数与 OfflineLocalizationPublisher"；参考点差异见 @ref output_pose_reading "位姿输出参考点"。

<b>深入：</b> @ref localization_module "定位与重定位"、@ref relocalization_theory "全局重定位：检索、几何验证与时间确认"、@ref pose_graph_theory "位姿图、增量求解与高频平滑"、@ref output_contracts "定位状态、车体参考点与 ROS 输出"。在线线程及发布门控另见 @ref online_localization_flow "在线定位端到端流程"。
