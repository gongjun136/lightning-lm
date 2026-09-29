@page script_contracts Shell 启动与实验脚本契约

# 读取规则

表中的“参数”指脚本自身消费的参数，不包含透传给二进制/子脚本的全部 flags；“环境”只列影响控制流、路径或资源的主要变量；cwd 指脚本是否主动改变工作目录；“副作用”包括文件写入、进程/ROS 图变化和包安装。精确默认值以脚本 `usage()` 与赋值语句为准。

## 运维与薄包装入口

| 脚本 | 参数 / 主要环境 | cwd 与调用链 | 副作用 |
|---|---|---|---|
| `run.sh` | 无 flags；`LIGHTNING_LM_WS`, `ROS_DOMAIN_ID`, `SANY_*`, profile 变量 | `cd $REPO` → `run_sany_online_diagnostics.sh --config ... --map ...` | 启动长期在线定位/诊断进程；输出 runs；可选录 bag |
| `run_frontend_online.sh` | `--config` 透传；`LIGHTNING_LM_{REPO_DIR,ROS_SETUP,INSTALL_SETUP,CONFIG}` | 不改 cwd；source 两个 setup → `ros2 run ... run_frontend_online` | 取决于节点配置：发布 ROS topic/UI/log |
| `run_slam_online.sh` | `--config` + 透传；同上，另有 OUT_ROOT/RUN_NAME | 创建并 `cd` 唯一 run dir → `ros2 run ... run_slam_online` | 写 logs/run metadata；在线建图服务可能写地图 |
| `run_loc_online.sh` | `--config` + 透传；同上 | 创建并 `cd` 唯一 run dir → `ros2 run ... run_loc_online` | 写 logs；加载地图并发布定位结果 |
| `align_maps_offline.sh` | 全部参数透传 | 不改 cwd；解析安装 binary → `exec align_maps_offline` | 由二进制读写所给地图/结果路径 |
| `save_default_map.sh` | 无 | 当前 cwd → `ros2 service call /lightning/save_map` | 请求在线节点写 `new_map` |
| `install_dep.sh` | 无 | 当前 cwd → `sudo apt install ...` | **系统级**安装 OpenCV/PCL/YAML/glog/gflags/ROS PCL 包 |

**代码依据：** 上述各脚本（`source/cd/exec/ros2/sudo` 调用与变量默认值）

## 标准离线运行器

| 脚本 | 必需参数；重要可选参数 | 主要环境 | cwd / 调用链 / 副作用 |
|---|---|---|---|
| `run_frontend_offline.sh` | `--bag --config --output-dir`; sequence/repeat/playback/max frames/CPU/watchdog/UI/各输出 | ROS/INSTALL setup，repo，部分 legacy 输入/输出变量 | 不改主 shell cwd；检查 bag → inspector → `setsid taskset run_frontend_offline` + resource monitor → timing extractor；创建 results/logs、轨迹/PCD/CSV/JSON/metadata；拒绝覆盖受管产物 |
| `run_slam_offline.sh` | `--bag --config --output-dir`; 另有 map/global/TUM/frame stats/backend-only | `LIGHTNING_LM_{INPUT_BAG,CONFIG,OUT_ROOT,RUN_NAME,WAIT_UI,OUTPUT_TUM,...}` | 同上，调用 `run_slam_offline`；额外写 tiled/global map、backend diagnostics/relocalization DB，并验证 map contract |
| `run_loc_offline.sh` | `--bag --config --map --output-dir`; initial pose/reference/physical bounds 等 | 同类 setup/repo/bag/config/map/output 变量 | inspector → `run_loc_offline` + monitor → timing + localization analyzer；写轨迹、误差、结果 bag、CSV/JSON；拒绝覆盖受管产物 |
| `run_frontend_offline_batch.sh` | `--bag --config-dir|--config ... --output-root`; repeats/dry-run | repo | 枚举配置并逐个调用 `run_frontend_offline.sh`；为每组建立 run dir 和汇总，已有同名/失败计数影响退出码 |

这些脚本会用 `setsid` 建立进程组，trap/看门狗可能向整个算法进程组发送 INT、TERM、KILL；会导出 BLAS/OpenMP 线程数。输出目录不是临时目录，属于实验证据。

**代码依据：** `scripts/run_*_offline*.sh`（usage、resolve_output、cleanup_group、monitor 与最终契约判断）

## 在线诊断与回归

| 脚本 | 参数 / 环境 | cwd 与调用链 | 副作用 |
|---|---|---|---|
| `run_sany_online_diagnostics.sh` | config/map/CPU/采样；大量 `SANY_*`, `LIGHTNING_LM_*` 资源与功能开关 | source 环境；创建 run dir；启动 recorder、`run_loc_online`、资源监控、位姿静默 watchdog | 写 config 快照、topic/进程快照、logs/results/bag；管理多个进程组；可能长期运行 |
| `run_online_bag_regression.sh` | bag/config/map/mode 与 playback 时序 | `cd run_dir`；启动在线 node + bag recorder/player；可调用 save_map service 与轨迹提取器 | 写 ROS bag、地图、轨迹和日志；结束时向子进程发信号 |

**代码依据：** 对应脚本（进程启动、trap/cleanup、目录创建与 ROS 命令）

## 复现与对比实验

| 脚本 | 参数/主要环境 | 调用链 | 主要副作用 |
|---|---|---|---|
| `reproduction/.../run_sany_offline_localization.sh` | bag、mapping/localization config、map、输出、CPU、watchdog | 先 `run_slam_offline.sh` 生成/选地图，再 `run_loc_offline.sh` | 创建 mapping + localization 两阶段产物 |
| `.../run_phase_a_relocalization_matrix.sh` | bag/config/map/output matrix | 多数据集/offset 调用 `run_loc_offline.sh` | 批量 run dir、stats、失败汇总 |
| `.../create_sany_formal_fault_bags.sh` | 输入/root；scenario/mode/topic/start/duration | source ROS → Python fault-bag creator | 创建裁剪/删改 topic 的 ROS 2 故障 bag |
| `.../run_sany_formal_fault_matrix.sh` | fault bags/config/output | 遍历场景调用在线 contract runner | 批量日志/结果，聚合失败码 |
| `.../run_sany_online_contract.sh` | bag/config/output/cpuset/fault contract | source ROS+workspace；record → online node → bag play → analyzers | 修改 ROS 图、录制结果 bag、写环境/commit identity、JSON/CSV；清理进程组 |
| `formal_report/rerun_m3dgr_fastlivo2_timestamp_fix.sh` | 脚本内实验矩阵及已有 runner | 对多个序列/repeat 调用外部/仓库 runner并核验指纹/轨迹 | 大批量报告实验目录与 inventory；依赖本机数据布局 |
| `single_lidar/m3dgr/run_voxel_slam_full_backend.sh` | ROS1 数据/资源环境 | 启动 ROS1 roscore + Voxel-SLAM + rosbag + monitor + extractor | 使用独立 ROS master 端口；写轨迹/PCD/资源与日志；信号清理 |
| `run_voxel_slam_114_reference.sh` | bag/output、ROS1 workspace、topic/frame/CPU | ROS1 roscore → roslaunch → recorder/monitor → rosbag play | 建 reference 轨迹与 contract；管理 ROS1 进程组 |

**代码依据：** `scripts/reproduction/**/*.sh`、`scripts/run_voxel_slam_114_reference.sh`（循环矩阵、子运行器和输出/清理逻辑）

# 脚本修改检查清单

修改任一 `.sh` 时至少检查：

1. 参数：usage、解析 case、默认值、透传边界是否一致。
2. 环境：新变量是否有默认值、是否 export、是否进入 metadata/config snapshot。
3. cwd：相对路径在哪个 `cd` 前后求值；安装 binary 与 repo 是否可能混用。
4. 调用链：source 顺序、ROS domain/master、子进程组与信号是否正确。
5. 副作用：覆盖保护、目录所有权、临时/持久产物、ROS service/topic、系统安装。
6. 退出：trap 是否保留真实算法退出码；超时/验证失败是否区分。

文档同步判定见 @subpage documentation_rules "基于代码 diff 的文档更新规则"。
