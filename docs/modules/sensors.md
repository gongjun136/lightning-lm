@page sensor_pipeline 传感器预处理、多雷达组帧与算力预算

# 边界和阅读路径

ROS 消息先由 `PointCloudPreprocess` 转成内部点型，再由 `MultiLidarFrameAssembler` 组帧，最后 `LaserMapping::SyncPackages` 等待 IMU 时间覆盖。接收单雷达消息、排出融合帧、取出同步量测是不同事件，任一步成功都不能直接当作 LIO 已输出。

| 想确认的问题 | 短验证路径 |
|---|---|
| 点的 time 是绝对秒还是相对毫秒？ | @ref sensor_point_time_reading "预处理的时间换算" |
| 后一帧已收齐，为什么还取不出来？ | @ref sensor_frame_ready_reading "时间桶、容差和顺序排出" |
| 融合后坐标与点时间分别变成什么？ | @ref sensor_fusion_coordinates_reading "外参与时间重基准" |
| 有一条晚于扫描末的 IMU，就会把它取走吗？ | @ref sensor_imu_sync_reading "覆盖检查与消费边界" |
| 降载是否会把完整输出云一起裁掉？ | @ref sensor_budget_reading "观测子集、健康条件和点数上限" |

下面的数值例子用于手工核对语句，均带假设，不代表现场配置或已执行的回归结果。

@anchor sensor_point_time_reading
# 1. 点的 time 以哪个时刻为零，单位是什么？

先限定输入类型：这里核对 Livox 的 `PointCloud2` 和 `CustomMsg` 两条路径。两者最终都把点时间写成<b>相对本扫描开始的毫秒</b>，但读入的字段和尺度不同。

`LivoxPointCloud2Handler()` 先取消息头的绝对秒，再用阈值判断缩放后的逐点 `timestamp` 应当如何解释：

@snippet{lineno} core/lio/pointcloud_preprocess.cc sensor-livox-time-reference

@snippet{lineno} core/lio/pointcloud_preprocess.cc sensor-livox-pointcloud2-time

按三步读：`livox_point_time_scale_` 先把输入字段换成秒；超过 `1e6` 的值减去 `head_time`，其余按相对秒解释；最后乘 `1e3` 写成毫秒。这是实现中的数值判断，不是自动识别任意驱动时间协议。负时间的拒绝阈值 `-1e-3` 此时也是毫秒。

<b>手算例子：</b>消息头为 1700000000.000 s，点时间为 1700000000.025 s、尺度为 1，得到约 25 ms；相对字段为 0.025 s 时也得到 25 ms。若驱动给的是绝对纳秒，应先配置尺度 `1e-9`。只改变消息头，不能替代点字段的单位换算。

`CustomMsg` 不走上述尺度与绝对时间判断，直接把 `offset_time` 的纳秒偏移除以一百万：

@snippet{lineno} core/lio/pointcloud_preprocess.cc sensor-livox-custom-time

例如 `offset_time=25000000` 对应 25 ms。片段前还存在抽样、扫描线、tag 等筛选，后面还有高度、重复点和盲区检查；时间换算正确不意味着该点一定保留。尤其 CustomMsg 从索引 1 开始做相邻原始点比较，不能按全部输入点估计输出数。

<b>停止位置：</b>能说明“先换成秒还是直接从纳秒转毫秒、是否减消息头”即可。其他型号应继续读自己的 handler，不能照搬 Livox 判断。完整实现见 @ref lightning::PointCloudPreprocess "PointCloudPreprocess API"、@ref pointcloud_preprocess.cc "各型号预处理源码"；去畸变如何将毫秒换回秒见 @ref imu_point_reading "逐点时间与补偿"。

@anchor sensor_frame_ready_reading
# 2. 已收齐的后一帧，为什么仍然不能排出？

组帧器先用第一次有效输入的时间建立 `anchor_time_`，再按周期四舍五入到时间桶。桶号、桶时刻、完整性和过期条件集中在下面几个短函数中：

@snippet{lineno} core/lio/multi_lidar_fusion.cc sensor-frame-bin-and-expiry

`IsComplete()` 要求所有配置雷达都已出现；`IsExpired()` 比较<b>已接收输入的最大传感器时间</b>与桶时刻。它不读取系统时钟，所以输入停止后，单纯等待几秒不会自动让 `PopReady()` 排出旧桶。

落入同一个桶也不保证接收。`AddCloud()` 对已有时间戳加上新时间戳取最大、最小值，检查的是整个桶内的跨度：

@snippet{lineno} core/lio/multi_lidar_fusion.cc sensor-frame-tolerance

例如周期 0.1 s、锚点 100.000 s、容差 0.002 s，100.000 与 100.0015 属于同桶且可共存；再来另一雷达的 100.0025，跨度变为 2.5 ms，即使桶号没变也会拒绝。容差并非“每条消息离桶中心都不超过 2 ms”。重复 ID 和已经排出桶的迟到输入，会在此前由独立分支拒绝。

再看何时移动到 ready 队列：

@snippet{lineno} core/lio/multi_lidar_fusion.cc sensor-frame-promote

<b>时序例子：</b>假设配置四雷达、周期 0.1 s、乱序窗口 0.5 s，关闭自适应降载且 `min_lidars=3`。100.0 s 的旧桶只到三雷达，100.1 s 的新桶已经到齐四雷达。循环总从 `frames_.begin()` 读最早桶，旧桶尚未过期时执行 `break`，新桶也不能越过它。若缺失雷达随后补到旧桶，旧桶先排出，再排新桶；若新接受的传感器时间推进到 100.6 s，旧桶已过期且三雷达达标，就以 partial 帧排出。若只有两雷达，则丢弃旧桶再继续。

启用自适应降载时，组帧最少雷达数实际取 `tracking_min_lidars`；后续定位选择还会检查自己的主雷达和最少来源条件，不能把组帧通过当作定位输入已合格。`Flush()` 以 `force=true` 越过等待，仍执行数量门槛；它用于离线尾帧，不宜在在线每帧调用来消除等待。

<b>读到这里：</b>能解释 `begin → break → 最少数量 → erase` 即可。完整入口见 @ref lightning::MultiLidarFrameAssembler::AddCloud "AddCloud()"、@ref lightning::MultiLidarFrameAssembler::PopReady "PopReady()"、@ref lightning::MultiLidarFrameAssembler::Flush "Flush()"。各类丢弃计数含义不同；实际缺失来源应读 `present_lidar_ids/missing_lidar_ids`，不能合成一个算法丢帧率。

@anchor sensor_fusion_coordinates_reading
# 3. 副雷达点进入融合帧后，坐标和时间怎样一起变化？

`Assemble()` 取本桶各消息头的最早时间作为 `stats.begin_time`。点坐标按固定外参变到主雷达系，点时间则加上本雷达消息头相对这个共同起点的偏移：

@snippet{lineno} core/lio/multi_lidar_fusion.cc sensor-fusion-transform-time

只追 `p_primary`、`header_offset_ms` 和 `transformed.time`。前者实现
\f$p_{primary}=R_{primary,lidar}p_{lidar}+t_{primary,lidar}\f$；后两者统一相对时间基准。`PointType transformed = point` 保留其他字段，随后明确重写来源 `lidar_id`；不能通过点坐标反推它来自哪台雷达。

<b>手算例子：</b>假设桶内最早消息头为 100.000 s，副雷达头为 100.001 s，副雷达点的原始 time=20 ms。融合后 time=21 ms，绝对采样时刻仍为 100.021 s。再设旋转为单位阵、外参平移 `(1,0,0)` m，原始坐标 `(2,0,0)` m 变为 `(3,0,0)` m；时间换算不会额外改变空间位置。若这是本桶最晚的点，`stats.end_time` 就是 100.021 s。

`max(0, ...)` 会将负的合成偏移截为零，不能据此证明原输入时间合法。融合后按点时间排序；这一步只做固定外参与时间重基准，运动补偿仍在下游 `ImuProcess` 中完成。

<b>停止位置：</b>能分别算出坐标、相对毫秒和保持不变的绝对采样时刻即可。完整实现见 @ref lightning::MultiLidarFrameAssembler::Assemble "Assemble()"、@ref multi_lidar_fusion.cc "组帧源码"；继续阅读 @ref imu_point_reading "扫描末坐标系中的逐点补偿"。`DownsamplePreservingSource()` 保留真实输入点而不平均来源 ID；车体自身点过滤则使用车辆前左上轴和主雷达原点，详见 @ref lightning::FilterSelfPoints "FilterSelfPoints()"。

@anchor sensor_imu_sync_reading
# 4. 最新 IMU 已覆盖扫描末，为什么那条 IMU 仍留在队列？

`SyncPackages()` 先为队首扫描确定 `lidar_end_time_`，并用 `lidar_pushed_` 保留等待中的扫描。正常点时间用于估计扫描末；点数太少、末点时间偏短或跨度异常时另有回退，不能无条件把它写成“消息头加固定周期”。下面从扫描末已经确定的位置开始：

@snippet{lineno} core/lio/laser_mapping.cc lio-wait-for-imu

这只是最新 IMU 时间的覆盖检查。真正取出的样本由后面的循环决定：

@snippet{lineno} core/lio/laser_mapping.cc sensor-sync-consume-imu

<b>时序例子：</b>扫描末为 10.10 s，IMU 队列依次为 `[10.08, 10.10, 10.12]` s。覆盖检查通过；循环取走 10.08 和 10.10，保留 10.12。若最新消息仅为 10.09，前面的 `return false` 发生在消费之前，等待同一帧。若队首一开始就是 10.10，`imu_time < lidar_end_time_` 为 false，本次循环一条也不取；不能简单概括为“所有小于等于末时刻的 IMU 都会取走”。这是当前严格比较的边界，本页不改变实现。

走到同步末尾才移除这一帧 LiDAR 及对应元数据：

@snippet{lineno} core/lio/laser_mapping.cc sensor-sync-consume-lidar

因此，覆盖检查通过不等于 IMU 批次一定非空，更不等于初始化完成或激光更新成功。相邻帧积分还依赖 `ImuProcess` 保存的 `last_imu_`，扫描末剩余时间由 `PredictTo()` 补齐，见 @ref imu_time_reading "跨帧衔接和正时间差"。

<b>停止位置：</b>能在上述例子中指出哪些消息被 pop、哪些留给后续处理即可。扫描末估计和完整同步实现见 @ref lightning::LaserMapping::SyncPackages "SyncPackages()"；消费后是否有输出继续 @ref lio_drain_reading "RunDetailed 与排空状态"。

@anchor sensor_budget_reading
# 5. 点数降载裁剪哪份云，跟踪变差后为何又取消上限？

下面只讨论启用多雷达、自适应降载且配置了非空 `lio_point_budgets` 的预算模式。`Select()` 同时看控制器已接收的地图定位与 LIO 跟踪健康状态，分别由 `SetLocalizationGood()` 和 `Observe()` 更新：

@snippet{lineno} core/lio/multi_lidar_fusion.cc sensor-budget-health

@snippet{lineno} core/lio/multi_lidar_fusion.cc sensor-budget-cap

预算模式把采集顺序步长设为 1，健康时才根据降载级别设置空间采样上限。`AdaptiveLidarSelection::max_points` 默认是 0，表示取消这一点数上限；不是要求输出零个点。质量下降后，来源数量也恢复为全部已配置来源的目标值，但只能选当前帧实际到达的来源：

@snippet{lineno} core/lio/multi_lidar_fusion.cc sensor-budget-sources

<b>手算例子：</b>假设预算为 `[12000,8000,4000]`，关闭 predictive，当前级别 1，三台雷达均到达且满足主雷达及最少数量门槛。两种健康状态都为 true 时，`max_points=8000`；如果控制器收到任一种健康状态为 false，下一次选择使用步长 1、`max_points=0`，并尝试保留三台实际来源。它恢复的是预算前的候选点，不能恢复预处理或体素过滤已经丢掉的点。若来源门槛不满足，前面的提前返回仍可能产生空选择。

这份选择怎样作用到点云？`RunDetailed()` 从完整去畸变云构造观测用子集，随后在体素过滤（以及必要的 0.1 m 重试）之后施加点数上限：

@snippet{lineno} core/lio/laser_mapping.cc sensor-lio-observation-subset

@snippet{lineno} core/lio/laser_mapping.cc sensor-lio-point-cap

`scan_undistort_full_ → scan_undistort_ → scan_down_lidar_` 是这段的阅读顺序。`GetPublicationCloud()` 返回第一份完整云，所以不能把观测预算数当作发布云点数；这里的“完整”仍已经过预处理、去畸变和可能的车体自身点过滤。发布是否允许还要经过独立来源条件与下游门控，见 @ref lightning::AdaptiveLidarLoadController::CanPublishCloud "CanPublishCloud()" 和 @ref output_contracts "输出契约"。

空间采样按来源、8 个方位扇区、z 正负和 10 m 距离界线分箱、分配配额。最后按原输入索引复制选中的真实点：

@snippet{lineno} core/lio/multi_lidar_fusion.cc sensor-sampling-real-points

`push_back(cloud->points[i])` 保留点时间和来源，索引递增保留输入顺序；不会生成体素中心或平均雷达 ID。例如两个非空箱分别有 6、2 个有效点，均衡模式上限为 4 时各得 2 个；density 模式按密度给出 3、1 个。这里只验证配额语义，分箱均衡不构成几何可观测性证明。

<b>按需继续：</b>降载级别由 `Observe()` 的连续过载/恢复计数改变；其处理时间阈值使用 `target_latency_sec`，不是 `frame_period`。启用 predictive 时，`Select()` 还可能提前增加级别，不能只按连续帧数解释所有变化。完整入口见 @ref lightning::AdaptiveLidarLoadController "控制器 API"、@ref lightning::SampleSpatiallyBalanced "SampleSpatiallyBalanced()"、@ref lightning::LaserMapping::GetPublicationCloud "GetPublicationCloud()"；运行配置与回归方法见 @ref compute_budget "点数预算与回归"。

# 调用与生命周期边界

assembler 由 `LaserMapping` 持有，当前输入和 flush 入口在调用者的缓冲区互斥锁内使用它；组帧器自身不加锁，不能脱离该约束并发修改。对应代码在 @ref lightning::LaserMapping::ProcessPointCloud2 "ProcessPointCloud2()" 和 @ref lightning::LaserMapping::FlushMultiLidar "FlushMultiLidar()"。

后续算法阅读顺序：@ref imu_deskew "IMU 初始化与去畸变" → @ref eskf_theory "ESKF 预测与更新" → @ref lidar_residuals "LiDAR 残差与信息矩阵"。
