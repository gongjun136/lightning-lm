@page relocalization_theory 全局重定位：检索、几何验证与时间确认

# 外层问题

局部配准需要初值，全局检索提供初值。检索分数只能说明描述子相似，不能独立证明机器人就在该处。@ref lightning::loc::GlobalRelocalizer "GlobalRelocalizer" 抽象候选生成，SOLiD/BTC 提供实现，@ref lightning::loc::LidarLoc "LidarLoc" 决定是否接受和切换定位状态。

下面按四个问题核对候选生成。它们止于“给 `LidarLoc` 一个地图中的 IMU 位姿”；最终提交继续读 @ref loc_read_confirmation "候选连续确认导读"。本页数值例子均为手算，不是实测结果。

@anchor relocalization_descriptor_reading
# 1. SOLiD 相似度为 1，究竟确认了什么？

@ref lightning::loc::SolidDescriptorEngine::Compute "Compute()" 先按水平距离 r、方位角 a、仰角 e 对有效点分箱，得到距离–仰角矩阵 \f$D_{re}\f$ 与方位–仰角矩阵 \f$D_{ae}\f$。先看同一组仰角权重怎样生成两种描述子：

@snippet{lineno} core/localization/solid_descriptor.cc relocalization-solid-descriptor-projection

`colwise().sum()` 统计各仰角箱点数；计数有差异时做 min-max 归一化，全部相等时给非空箱权重 1，避开除零。最后两行分别是 \f$d_r=D_{re}w\f$、\f$d_a=D_{ae}w\f$。地点检索的打分只用前者：

@snippet{lineno} core/localization/solid_descriptor.cc relocalization-solid-range-score

<b>手算：</b>假设 \f$D_{re}=\begin{bmatrix}1&1\\0&2\end{bmatrix}\f$，两列计数是 `[1,3]`，所以 `w=[0,1]`、`range=[1,2]`。候选 `range=[2,4]` 的余弦相似度仍为 1；它既不证明点数相同，也不证明 `angle` 相同。保持范围过滤条件不变的纯绕 z 轴旋转可以保持 `range`，航向需要下一节另算。

<b>边界与停止点：</b>片段之前检查了尺寸和有限性；`Compute()` 还会拒绝没有有效点、非有限或零范数的输出。读到 `Similarity()` 中只出现 `range` 已能确认分数的含义，无需追进 Eigen 的点积。完整实现见 @ref lightning::loc::SolidDescriptorEngine::Similarity "Similarity()" 和 @ref solid_descriptor.cc "solid_descriptor.cc"；召回数量与 `min_similarity` 的筛选在 @ref lightning::loc::SolidRelocalizer::AddFrame "AddFrame()"。

@anchor relocalization_direction_reading
# 2. 航向和配准矩阵把谁变换到谁？

这里约定 \f$T_{A,B}\f$ 把 B 坐标变到 A 坐标。查询子图由 @ref lightning::loc::SolidRelocalizer::BuildQuerySubmap "BuildQuerySubmap()" 用 `T_current_lidar_odom * frame->T_odom_lidar` 合到<b>最新查询 LiDAR 坐标系</b>。下面的循环比较 `query.angle[index]` 与候选的 `index+shift`，因此返回 query → candidate 的 yaw：

@snippet{lineno} core/localization/solid_descriptor.cc relocalization-solid-yaw-direction

<b>先手算符号：</b>设 4 个角度箱，query 为 `[1,0,0,0]`、candidate 为 `[0,1,0,0]`，`shift=1` 才使距离为 0，返回 `+π/2`。query 的 x 正向旋到 candidate 的 y 正向，正好符合这个方向。`shift=2` 则折回 `-π`；严格 `<` 使并列最优保留先遍历到的 shift。对称场景仍可能有多个等价航向，此函数没有消除这种歧义。

`solid_kiss` 分支中，KISS 的 source 也是该查询子图，target 是数据库候选 LiDAR 子图。`RefineCandidateWithKissMatcher()` 只有在 `match.success` 为真后才执行下面的组合：

@snippet{lineno} core/localization/solid_relocalizer.cc relocalization-kiss-pose-composition

第一行乘法是 \f$T_{map,current\ lidar}=T_{map,candidate\ lidar}T_{candidate\ lidar,current\ lidar}\f$；末尾的 `T_imu_lidar_.inverse()` 再把 LiDAR 位姿换成 IMU 位姿。中间还调用了 `MakePlanarMapOdom()`，不能略过：

@snippet{lineno} core/localization/solid_relocalizer.cc relocalization-planar-map-odom

这里只保留原始 map ← odom 的 yaw，并重算平移，使当前 LiDAR 的地图位置保持不变。重新乘回 LIO 位姿后，roll/pitch 来自 LIO；这不是把粗配准的完整三维旋转直接交出去。

<b>再手算平移：</b>假设所有旋转为单位阵，候选 LiDAR 在地图 x=10 m，配准给出 query → candidate 的 x 平移 +2 m，`T_imu_lidar` 的 x 平移为 +0.3 m。则查询 LiDAR 在 x=12 m，返回 IMU 在 x=11.7 m。把配准矩阵求逆会得到完全不同的地点。

BTC 的 @ref lightning::loc::BtcRelocalizer::AddFrame "AddFrame()" 同样先乘 `entries_[id].T_world_lidar * T_candidate_lidar_current_lidar`，再乘外参逆；它没有这里的 `MakePlanarMapOdom()` 步骤。BTC 从平面结构、二进制特征和三角形描述生成候选，在多个查询尺度调用 `SearchLoopTopK`，再按分数及空间/航向去重；与建图 `BtcLoopDetector` 复用描述算法，但历史与接受策略不同。

<b>读到哪里停：</b>确认乘法顺序、平面化和外参方向后即可回到下一节；完整调用者见 @ref lightning::loc::SolidRelocalizer::RefineCandidateWithKissMatcher "RefineCandidateWithKissMatcher()"、@ref solid_relocalizer.cc "SOLiD 源码" 和 @ref btc_relocalizer.cc "BTC 源码"。这一步只给候选地图位姿，还没有切换跟踪状态。

@anchor relocalization_coarse_acceptance_reading
# 3. KISS 已有变换或 converged，就能接受粗配准吗？

本节适用于 `solid_kiss` 且 `refine_with_icp=true` 的粗配准分支。`MatchFeatures()` 先对 FPFH 做比例筛选和互匹配；`SolveCorrespondences()` 再按点对距离差不超过 `compatibility_tolerance`、点对距离不小于 `minimum_pair_distance` 建图，保留 maximum k-core，通过 GNC 权重估计刚体变换。到此仍可能被下面的门控拒绝。

@snippet{lineno} core/localization/kiss_matcher_registration.cc relocalization-kiss-inlier-gate

`SquaredResidual()` 是 \f$\|Tq-c\|^2\f$，与 `inlier_threshold²` 比较。注意 `inlier_ratio` 的分母是 `core_indices.size()`，不是全部 FPFH 对应数；硬门槛是内点<b>数量</b>。接着只对已选内点再拟合，统计 target 的 XY 网格覆盖，最后执行：

@snippet{lineno} core/localization/kiss_matcher_registration.cc relocalization-kiss-quality-gate

<b>逐句核对：</b>`spatial_coverage` 是“内点占据网格数 / core 占据网格数”，但拒绝条件使用的是 `spatial_cells < minimum_spatial_cells`。`score` 是质量统计，此处没有额外分数阈值。`T_target_source` 已经赋值、`converged` 已被 `iterations>0` 置真，后面仍可能因 RMSE 或网格数返回；必须读到 `success=true` 才算本函数接受。

<b>手算：</b>假设前面所有阶段都通过，core 有 20 对、内点 12 对，最终这些内点残差都为 0.3 m，则 `inlier_ratio=0.6`、RMSE=0.3 m。若 `minimum_inliers=12`、`maximum_rmse=0.5 m`、`minimum_spatial_cells=4`，但内点只落在 3 个 target XY 网格，仍返回 `insufficient_spatial_coverage`。即使 core 也只覆盖这 3 格、`spatial_coverage=1`，结果不变。

<b>适用边界：</b>最终重拟合后，这里没有重新选一遍内点，RMSE/覆盖针对先前选中的集合；不要把 `inliers` 解读成新变换下重新分类后的数量。关闭 `refine_with_icp` 会跳过这层粗配准，但后续 `LidarLoc` 的地图验证仍有自己的接受路径。读到 `success` 与调用者的失败分支即可停止，算法内部细节按需读 @ref lightning::loc::KissMatcherRegistration::SolveCorrespondences "SolveCorrespondences()"、@ref kiss_matcher_registration.cc "完整 KISS 实现"；最终状态提交转到 @ref loc_read_confirmation "连续确认"。

@anchor relocalization_database_frame_reading
# 4. 旧数据库能加载，是否就说明它和地图同坐标系？

数据库中的候选位姿参与上一节的矩阵乘法；如果地图归一化后仍配旧库，即使相对配准正确，也会带着旧世界坐标进入组合。加载检查的条件需要完整看，下面是 BTC 的实际分支：

@snippet{lineno} core/localization/btc_relocalizer.cc relocalization-btc-map-frame-gate

先看最外层 `normalize_start_ground_z`：启用时，BTC 不接受缺少坐标引用的 schema 1；schema 2/3 还要读出 `map_frame` 并通过 `SameTransform()`。在它之前已读取地图包元数据，失败同样退出。@ref lightning::loc::SolidRelocalizer::Init "SOLiD 的 Init()" 也在相同配置条件下加载地图元数据并检查数据库引用，但 SOLiD 自己的 schema 当前是 1，不能套用 BTC 的版本门槛。

<b>手算：</b>假设两个引用都标为 normalized、ID 相同、旋转相同，但地图包的 \f$T_{export,slam}\f$ z 平移是 -1.0 m、库中是 0 m。@ref lightning::map_frame::SameTransform "SameTransform()" 会得到 1 m 的平移差，超过默认 `1e-9` 容差，加载失败；只把 ID 改成相同仍无效。它也会核对 normalized 标志、ID 和旋转差。

<b>边界与停止点：</b>关闭 `normalize_start_ground_z` 会跳过这里的引用一致性检查，成功加载不能作为坐标一致的证明；`SameTransform()` 比较的是元数据，也没有逐点证明 PCD 内容正确。确认配置、引用和失败返回就足够解释这一层，地图包如何统一变换继续看 @ref map_contracts "地图坐标契约"，无需在这里追栅格导出。完整入口为 @ref lightning::loc::BtcRelocalizer::Init "BTC Init()"、@ref btc_relocalizer.cc "btc_relocalizer.cc"、@ref map_frame.cc "map_frame.cc"。

# 候选到可用定位

希望先从实现验证接受条件时，读 @ref loc_read_confirmation "候选连续确认导读"：真实片段展示 `T_map_odom` 的比较、计数复位与最终提交；本节继续解释几何含义。

`LidarLoc` 加载候选附近地图，按配置使用 NDT、点面 ICP 或可选注册后端细化，再检查重力方向、修正范围、静态地图一致性和连续帧确认。细化推导见 @ref lidar_residuals "LiDAR 残差、信息矩阵与配准参数化"。
连续确认比较候选与相对里程计之间的时空一致性；不能把重复查询同一帧算作多次独立证据。重定位通过后，下游 PGO/外推与保持状态需要开启新的时间段，不能继续融合失败前遗留的相对约束。

重力对齐可看 `T_map_odom.rotationMatrix()(2,2)`，它表示两坐标系 z 轴夹角余弦。该门控隐含地图与里程计均有可信重力方向；不是完整的三维旋转误差指标。

# 调参和故障定位

先区分“描述子无候选”“候选几何失败”“待时间确认”“确认后跟踪失败”。扩大候选数量只影响第一层，放宽匹配阈值可能绕过第二层质量控制。`lidar_loc_valid_`、融合 `valid_` 与真正 ROS 发布是不同事实，见 @ref output_contracts "定位状态、车体参考点与 ROS 输出"。

<b>源码导航：</b>@ref solid_descriptor.cc "solid_descriptor.cc"、@ref solid_relocalizer.cc "solid_relocalizer.cc"、@ref btc_relocalizer.cc "btc_relocalizer.cc"、@ref btc_loop_detector.cc "btc_loop_detector.cc"、@ref lidar_loc.cc "lidar_loc.cc"（全局重定位与 `ProcessCloud`）；算法第三方实现与许可证位于 `src/core/backend/third_party/voxel_slam_btc/`。
