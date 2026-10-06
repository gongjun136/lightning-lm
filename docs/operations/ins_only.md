@page ins_only_operation CGI-430 纯组合导航模式

# 适用范围与入口

`run_loc_online` 启动时读取 `system.localization_mode`：缺省或 `lidar` 沿用现有激光定位；`ins_only` 创建 `InsLocSystem`。切换必须重启。纯组合导航模式直接使用 CGI-430 的导航解，不创建 LIO、地图匹配、重定位或 PGO，不做传感器融合，也不依赖地图启动。

代码：`src/app/run_loc_online.cc`、`src/core/system/ins_loc_system.cc`、`src/core/navigation/ins_navigation.cc`。坐标转换使用 GeographicLib，未引入 IC-GVINS/GLIO 的优化器或复制 GPL 实现；这些框架涉及的原始 GNSS/IMU 融合不属于本模式。

# 坐标和参考点契约

- CGI WebUI 的定位输出点必须是**后轴中心**，姿态必须是**后车体姿态**；本层不再补偿 GNSS/IMU 杆臂。
- 输入纬度/经度单位为度。固定 `origin_llh` 顺序为纬度、经度、WGS84 椭球高（米）；每次启动使用同一个现场配置，不能取本次首帧作为原点。
- 位置经过 WGS84 LLH → ECEF → 固定 ENU。`map` 的 x/y/z 是东/北/上，米。远离原点时，用 GeographicLib 的局部 ENU 到固定 ENU 旋转同步转换姿态和速度；默认拒绝离原点超过 20 km 的点。
- 高程必须核实：`ellipsoidal` 直接使用；`orthometric` 按 `h=H+N` 加入已知 `geoid_separation_m`。常数 N 只适用于经验证的现场范围，不是通用大地水准面模型。CGCS2000/ITRF 等输入不能仅改标签为 WGS84。
- CHC 航向从北顺时针，俯仰抬头为正、横滚右倾为正。后体 FLU 到测量点 ENU 的旋转为 `Rz(90°-heading) Ry(-pitch) Rx(roll)`，再左乘测量点 ENU 到固定 ENU 的旋转。业务消息仍由原有 `sany_output` 构造：`f8pry` 顺序是 **roll,pitch,yaw（度）**，yaw 从东逆时针。
- 有符号速度为导航速度在后体前向的投影；倒车为负，不取速度模长。
- 点云另需实测的 `T_rear_primary`，作用是 `p_rear=R*p_primary+t`。这与设备内部的天线杆臂是两件事。禁止把原 LIO 初始化时的重力姿态当成已知固定安装旋转。

原始手册依据：知识库 `_附件/SANY无人装载机/CGI-430/CHC®+CGI-430厘米级组合导航系统用户手册-V2.2.3-20250110.pdf`，PDF 第 18–20、53、56、63 页。状态码、卫导模式参考点变化、WebUI 安装姿态与正负号需要结合阅读。CAN 手册没有充分明确当前设备高程配置；上线前必须核实实际设置。

# 时间与严格有效性

订阅 SDK 的原始类型：`position/latitude`、`position/longitude`、`position/altitude`、`attitude`、`velocity`、三类 `sigma` 和 `ins/status`，默认前缀 `/cgi430`。不依赖简化后的 NavSatFix 状态，因为它不能区分所有设备状态和航向是否有效。

对应 SDK 缺省 CAN ID（**十进制**）为 814、813、805、810、807、806/811/808、803。设备须启用这些帧并提供与门限相容的频率；缺少任一类就不能进入 GOOD。ID 映射由 SDK 管理，本节点消费解码后的类型消息。

每类字段保存自己的 CAN 接收时间戳；所有字段必须有效、严格递增，且最新值满足配置的时间跨度和新鲜度。位置、姿态、速度五类字段各自推进后才能产生一条新的解。CAN 未提供公共测量历元，本实现是**有界接收时间配对**，不声称硬件同步。`receive_time_offset_sec` 仅用于已经测量确认的时钟偏移；不能用放宽阈值代替 CGI、主机与 LiDAR 时钟标定。诊断记录每个字段的时间戳。

默认门限如下，均为初始工程值，需要现场验证：

| 条件 | 默认值 |
|---|---|
| 工作状态 / 卫星状态 | 必须为 2 / 4，即组合导航且 RTK 固定、定向有效 |
| 接收时间新鲜度 / 字段跨度 / 未来容差 | 150 / 20 / 20 ms |
| 差分龄期 | ≤ 2 s |
| 东、北各轴位置 sigma / 高程 sigma | ≤ 0.10 / 0.20 m |
| 三轴姿态 sigma / 速度 sigma | ≤ 0.5° / 0.20 m/s |
| 恢复 | 连续 10 条新解通过全部条件 |

浮点状态 5、固定但未定向状态 8、纯惯导状态 3、卫导模式均不发布业务位姿。任何异常立即复位恢复计数并清空姿态插值历史；没有 DR 降级。业务位姿在回调检查中停止；50 ms 健康定时器检测静默并发布状态。主机单调时钟另外监视断流，因此暂停 `/clock` 也不会永久保持 GOOD。时钟回退后不会倒序发布；应恢复正确时钟或重启本模式。

`/localization/loc_status` 初始化为 INITIALIZING，通过恢复后 GOOD，首次成功后的失效为 FAIL；`fault_status` 不可用时 P0，`fault_type=2`。心跳继续，未发布新位姿时 `work_seq` 不推进。日志写失败也停止业务输出。该节点使用单线程 ROS 回调，点云处理耗时可能增加门控通知延迟；Orin 上需测量实际调度和最大延迟。

# 业务与点云输出

| Topic | 约定 |
|---|---|
| `/PosRes` | 后轴 ENU 位置、FLU 欧拉角、带符号速度 |
| `/slamPoseRaw_topic` | 同一后轴位姿 PoseStamped |
| `/localization/pose_vel` | 原 VehiclePose 类型与 comm_header |
| `/localization/loc_status`、`fault_status` | 原类型、独立于业务位姿发送 |
| `/diagnostics/heartbeat/lightning_slam` | node_id=6，独立心跳线程 |
| `/localization/ins_diagnostics` | `mode`、拒绝原因、点云丢弃数、待处理数 |
| `/localization/path` | 最多 500 个历史有效位姿，GOOD 时发布 |
| `/tf` | 可选 `map -> rear_axle`，只在有效发布时更新 |
| `/LidarDataInv` | 扫描结束时刻后轴 FLU 点云，frame=`rear_axle` |
| `/LidarDataInL` | 相同扫描在固定 ENU 中，frame=`map` |

旧 TF/位姿可能仍在接收端缓存，下游必须同时检查状态和时间戳。该模式不生成 LIO/地图匹配统计，详细诊断使用 `ins_diagnostics`。

点云复用现有 SANY Livox PointCloud2 预处理、多雷达外参和帧组装。合并到主雷达坐标后执行自车滤波，以导航位姿进行平移线性插值和四元数 slerp，逐点变换到扫描结束时刻后轴坐标。无外推；缺雷达、不完整位姿覆盖、插值间隔过大或超过等待时间就丢弃该帧。待发布队列上限 8 帧，默认等候 300 ms。当前严格策略在导航失效时**两种点云都停发**，恢复后要重新积累完整扫描区间；不以失效全球位置继续去畸变。

# 构建和启动

Ubuntu 22.04 / ROS 2 Humble，系统需安装 `libgeographic-dev`。`cgi430_interfaces` 位于 `common_msgs`，CGI 与 Livox 驱动均在同一工作区。CMake 同时兼容 GeographicLib 库名 `Geographic` 和 `GeographicLib`。

```bash
sudo apt-get install libgeographic-dev
# 从工作区根目录构建全部包（默认 2 线程）
bash src/lightning-lm/scripts/build_workspace.sh
```

可用 `LIGHTNING_LM_BUILD_WS`、`LIGHTNING_BUILD_JOBS` 覆盖工作区和并发数。旧 `build_ins_only.sh` 入口保留，但已转为统一构建。`CGI430_SDK_WS` 不再使用。安装目标包含 `ins_project_llh` 与 `georeference_map`。

源码目录迁移引起缓存冲突时，可加 `--cmake-clean-cache`；该参数重置全部包的 CMake 缓存选项，详见 @ref build_and_run "构建与缓存修复"。脚本还会去掉旧外部 `cgi430_interfaces_DIR` 缓存。

复制 `config/ins_only/sany_cgi430.yaml` 为现场配置，填入并确认原点、高程类型和主雷达到后轴的刚体外参。模板故意包含 null/confirmed:false，不能直接用于车辆。外部 `cloud.lidar_config` 相对**现场 YAML 所在目录**解析，移动配置后应改为正确绝对路径。纯位姿台架测试可设 `cloud.enabled:false`。

```bash
bash scripts/run_ins_only.sh /absolute/site.yaml
```

启动脚本只加载 ROS 与统一工作区安装环境，以仓库为 cwd；`ROS_DOMAIN_ID` 缺省 42，可显式覆盖。SDK 驱动和三雷达驱动分别启动。`LIGHTNING_LM_INSTALL_SETUP` 可覆盖 Lightning 安装环境脚本。支持 `--output_published_tum /new/file.tum`；拒绝 LIO/PGO 轨迹参数及内嵌 `--bag` 回放。

# 记录、回放与精度评估

节点每次创建 `audit_root/ins_only_<启动微秒时间>/`，保存配置快照（包括外部雷达配置）、`can_frames.csv`、`solutions.csv`、`events.csv` 和实际发布的 `published_rear_axle.tum`。raw CSV 保存九类 CAN 帧时间和字节；solutions 保存经纬高、姿态、速度、sigma、各字段时间、转换后 ENU 位姿及是否发布/拒绝原因。未通过质量条件的转换值仅供分析。TUM 拒绝覆盖已有文件。

录制/回放脚本要求调用终端已加载 ROS 和消息环境，并设置与节点相同的 domain。MCAP 存储插件应已安装：

```bash
bash scripts/record_ins_only.sh /absolute/site.yaml /new/bag_directory
# 在单独复制的回放配置中设置 ins_only.use_sim_time: true，先启动定位节点
bash scripts/run_ins_only.sh /absolute/site_replay.yaml
# 另一终端，仅回放输入 topic；不回放旧业务结果、TF、旧 /clock
bash scripts/replay_ins_only.sh /absolute/site_replay.yaml /recorded/bag --rate 1
```

记录脚本录制 CGI 输入、配置中的雷达输入、业务位姿、状态和心跳；节点 CSV 不代替原始 bag。回放暂以正常速度验收，慢速播放可能触发单调时钟断流门限；循环或回退播放应重启节点。停止/恢复原因同样保留。

实机精度评估需要独立真值，统一参考点、时钟、平面基准和高程基准，分别统计水平/高程/姿态/速度误差、P95、最大误差、可用率、断流检测与恢复时间。地图对齐使用的轨迹段不能再当作独立精度验收数据。

# 旧地图与任务迁移

`scripts/ins_georeference.py` 需要 numpy、scipy、PyYAML。fit 的两条 TUM 都必须是**后轴中心和后车体姿态**；source 属于旧地图，target 是本模式发布的固定 ENU。必须实测时间偏移，不自动以轨迹误差拟合时延。首 70% 配对用于刚体拟合，末 30% 留出验证；不估计尺度，拒绝近共线轨迹。阈值由项目验收要求显式给定，例如：

```bash
python3 scripts/ins_georeference.py fit --source old_rear.tum --target enu_rear.tum \
  --site /absolute/site.yaml --output /new/transform.yaml \
  --time-offset-sec 0 --max-rmse-m 0.10 --max-error-m 0.30 --max-angle-deg 1
python3 scripts/ins_georeference.py tum --source old_rear.tum \
  --transform /new/transform.yaml --output /new/enu_rear.tum
python3 scripts/ins_georeference.py tasks --source tasks.csv --angles degrees \
  --transform /new/transform.yaml --output /new/tasks_enu.csv
ros2 run lightning_lm georeference_map old_map.pcd /new/transform.yaml old_primary.tum /new/map_enu
```

任务 CSV 列为 `x,y,z,roll,pitch,yaw`，角度单位必须显式指定，其他列原样保留。TUM 的坐标迁移保留时间戳。变换清单包含固定地理原点、输入 SHA256、时间偏移、训练/留出误差和 accepted 标志；未通过验收的清单不能用于迁移。示例阈值不是实车精度承诺。

`georeference_map` 转换完整 PCD 并生成 `map.pcd`、`tiled/`、`georeference.yaml`；第四个输入轨迹的首位姿作为 tiled 起点，须是旧地图中**主雷达/原定位内部参考点**，不要误传后轴轨迹。起点可从原地图的建图轨迹取得。输出目录必须不存在。SOLiD/BTC 数据库须在新坐标系重新构建，不能复制旧数据库；不自动部署地图或修改任务系统。

源地图、源轨迹和源任务必须属于相同的旧地图导出坐标系；如果旧地图做过地面高度归一化，不能混入归一化前的 SLAM 轨迹。物理迁移后的激光地图路径指向新 `tiled/`，其业务输出不能再次叠加相同的 `output.fixed_map_transform`。纯组合导航分支本身直接输出固定 ENU，不应用旧地图的输出变换。新建地图也须先与该固定 ENU 原点建立明确关联，不能只复用 `frame_id: map`。

已有 PGM 在纯平移/航向变换时可用 `grid --source old/map.yaml --transform ... --output /new/grid --source-ground-z 0`，保留栅格值并转换原点/yaw。涉及倾斜的三维变换会拒绝，须利用转换后的扫描、轨迹与射线原点重新生成或重新建图。若下游只支持 yaw=0 的地图，也应重新栅格化。没有现场配对数据时，不能替现场填写变换。

# 验证入口与边界

```bash
# 在完成构建并加载工作区环境后
ctest --test-dir ../lightning_lm_ws/build/lightning_lm -R ins_navigation_test --output-on-failure
RMW_IMPLEMENTATION=rmw_cyclonedds_cpp PYTHONNOUSERSITE=1 python3 scripts/test_ins_only_ros.py \
  --binary ../lightning_lm_ws/build/lightning_lm/bin/run_loc_online --record
RMW_IMPLEMENTATION=rmw_cyclonedds_cpp PYTHONNOUSERSITE=1 python3 scripts/test_ins_only_replay.py \
  --binary ../lightning_lm_ws/build/lightning_lm/bin/run_loc_online --recorded-test runs/ins_ros_test_<生成的目录>
python3 scripts/test_ins_georeference.py --map-binary ../lightning_lm_ws/build/lightning_lm/bin/georeference_map
python3 scripts/validate_ins_m3dgr.py \
  --bag /mnt/f/datasets/M3DGR/Dataset/Standard/Outdoor01/Outdoor01.bag \
  --ground-truth /mnt/f/datasets/M3DGR/GT/Standard/Outdoor01.txt \
  --projector ../lightning_lm_ws/build/lightning_lm/bin/ins_project_llh --output /new/m3dgr_validation
```

测试脚本固定使用 domain 88 和 `ROS_LOCALHOST_ONLY=1`，只运行本机合成输入。`--record` 同时验证录包脚本；回放测试检查输出确实由新节点生成，并验证 `/clock` 停止后单调时钟仍触发 FAIL。此次 WSL 的 Fast DDS 2.6.10 连官方 talker/listener 也只能发现、不能接收业务消息；换用 Cyclone DDS 后独立对照正常，因此 ROS 验收使用 Cyclone DDS。这没有修改生产启动脚本的中间件选择。需要复现实验时安装 `ros-humble-rmw-cyclonedds-cpp`，并让测试所有进程使用相同 `RMW_IMPLEMENTATION`；Orin 上仍需验证实际生产中间件和网络。

若用户目录的 numpy 2 与 Ubuntu scipy 发生 ABI 冲突，可在坐标迁移脚本前设 `PYTHONNOUSERSITE=1` 使用配套的系统包；不要改动整个 ROS Python 环境。M3DGR 读取依赖 rosbags，可独立安装到数据分析环境。

M3DGR Outdoor01 的原始 PVT 提供经纬度和椭球高，可对照独立 WGS84 ECEF 公式验证生产 GeographicLib 转换；其真值 TUM 的重复/倒序时间及全单位四元数会单独审计，不能用来验收完整姿态。当地 GVINS school_scooter 文件只有 LiDAR，不能作为本模式完整输入。公共数据不能检验 CGI CAN 解码、WebUI 配置、真实时延、杆臂或实机绝对精度。

上线前仍需现场完成：WGS84/高程确认、固定原点、雷达到后轴外参、CAN 与 LiDAR 时间差及扫描内运动检验、静止/四向/转弯/倒车验证、RTK/定向故障注入和下游停用响应、Orin 的 CPU/延迟验收。现有工程阈值和合成测试不能替代这些测量。

参考实现：[GeographicLib LocalCartesian](https://github.com/geographiclib/geographiclib/blob/main/include/GeographicLib/LocalCartesian.hpp)、[robot_localization navsat_transform](https://github.com/cra-ros-pkg/robot_localization/blob/ros2/src/navsat_transform.cpp)、[M3DGR 官方仓库](https://github.com/sjtuyinjie/M3DGR)。
