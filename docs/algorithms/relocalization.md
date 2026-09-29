@page relocalization_theory 全局重定位：检索、几何验证与时间确认

# 外层问题

局部配准需要初值，全局检索提供初值。检索分数只能说明描述子相似，不能独立证明机器人就在该处。@ref lightning::loc::GlobalRelocalizer "GlobalRelocalizer" 抽象候选生成，SOLiD/BTC 提供实现，@ref lightning::loc::LidarLoc "LidarLoc" 决定是否接受和切换定位状态。

# SOLiD 描述子

`SolidDescriptorEngine::Compute` 将有效点按水平距离 r、方位角 a 和仰角 e 离散，累计距离–仰角矩阵 \f$D_{re}\f$ 与方位–仰角矩阵 \f$D_{ae}\f$。用仰角占用量的 min-max 归一化得到权重 w：
\f[
d_r=D_{re}w,\quad d_a=D_{ae}w,\quad
s(q,c)=\frac{d_r(q)^Td_r(c)}{\|d_r(q)\|\|d_r(c)\|}.
\f]
等计数退化情形对非空 bin 赋权 1，空/非有限/零范数描述子返回无结果。距离描述子用于检索，不编码绝对 yaw；航向由方位描述子的循环移位估计：
\f[
k^*=\arg\min_k\sum_i|d_a(c)_{(i+k)\bmod M}-d_a(q)_i|,
\quad\psi=2\pi k^*/M.
\f]
最终 shift 折回有符号半周区间。这个 yaw 是 `CandidateFromQuery` 方向，不要求逆两次。

# BTC 与位姿组合

BTC 从平面结构、二进制特征和三角形描述建立候选，依据描述投票和几何验证搜索数据库。`BtcRelocalizer` 在多个查询子图尺度上调用 `SearchLoopTopK`，按 score 过滤并作空间/航向去重。它与建图的 `BtcLoopDetector` 复用描述算法，但查询历史、排除近邻和最终接受策略不同。

若检索返回 \f$T_{candidate,current}\f$，数据库记录 \f$T_{world,candidate}\f$，则
\f$T_{world,current}=T_{world,candidate}T_{candidate,current}\f$。
`database.yaml` 的地图变换引用必须与地图包一致；地图归一化后只复制旧描述库会造成坐标不一致。

# 候选到可用定位

`LidarLoc` 加载候选附近地图，按配置使用 NDT、点面 ICP 或可选注册后端细化，再检查重力方向、修正范围、静态地图一致性和连续帧确认。细化推导见 @ref lidar_residuals "LiDAR 残差、信息矩阵与配准参数化"。
连续确认比较候选与相对里程计之间的时空一致性；不能把重复查询同一帧算作多次独立证据。重定位通过后，下游 PGO/外推与保持状态需要开启新的时间段，不能继续融合失败前遗留的相对约束。

重力对齐可看 `T_map_odom.rotationMatrix()(2,2)`，它表示两坐标系 z 轴夹角余弦。该门控隐含地图与里程计均有可信重力方向；不是完整的三维旋转误差指标。

# 调参和故障定位

先区分“描述子无候选”“候选几何失败”“待时间确认”“确认后跟踪失败”。扩大候选数量只影响第一层，放宽匹配阈值可能绕过第二层质量控制。`lidar_loc_valid_`、融合 `valid_` 与真正 ROS 发布是不同事实，见 @ref output_contracts "定位状态、车体参考点与 ROS 输出"。

**源码导航：** @ref solid_descriptor.cc "solid_descriptor.cc"、@ref solid_relocalizer.cc "solid_relocalizer.cc"、@ref btc_relocalizer.cc "btc_relocalizer.cc"、@ref btc_loop_detector.cc "btc_loop_detector.cc"、@ref lidar_loc.cc "lidar_loc.cc"（全局重定位与 `ProcessCloud`）；算法第三方实现与许可证位于 `src/core/backend/third_party/voxel_slam_btc/`。
