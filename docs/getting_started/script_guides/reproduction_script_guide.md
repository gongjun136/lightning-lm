@page reproduction_script_guide 复现实验脚本导读

这组脚本固定了一部分数据布局、场景与实验条件。阅读时先确认它调度哪个 runner，再看矩阵与结果判定；通用算法入口由对应子脚本继续进入。

| 要理解的实验 | 导读 |
|---|---|
| SANY 先建图再定位 | @ref sany_offline_pipeline_guide "离线两阶段流水线" |
| 从不同传感器时间开始重定位 | @ref phase_a_matrix_guide "Phase A 起始时间矩阵" |
| 制造缺失/中断输入 | @ref fault_bags_guide "正式故障 bag 生成" |
| 遍历四个故障场景 | @ref fault_matrix_guide "正式故障矩阵" |
| 检查在线前端公开接口与恢复 | @ref sany_online_contract_guide "SANY 在线合同" |
| 重跑时间戳修复影响的 FAST-LIVO2 单元 | @ref fastlivo_rerun_guide "冻结矩阵重跑" |
| M3DGR 上的 Voxel-SLAM 完整后端对比 | @ref voxel_backend_guide "ROS 1 后端基线" |

@anchor sany_offline_pipeline_guide
# run_sany_offline_localization.sh：先建图，再用该地图定位

## 先看核心：两个顺序执行的子任务

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_sany_offline_localization.sh sany-offline-mapping

第一阶段调用 @ref slam_offline_guide "run_slam_offline.sh"，得到分块地图及 SLAM 参考轨迹。索引或轨迹缺失就结束；成功后才进入第二阶段：

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_sany_offline_localization.sh sany-offline-localization

第二阶段调用 @ref loc_offline_guide "run_loc_offline.sh"，使用刚生成的地图，并把建图轨迹交给结果分析器。两条命令都是前台执行；严格模式下第一阶段失败不会继续定位。

## 按执行顺序阅读

1. 从环境与参数选 bag、两份 YAML、输出根和 run prefix；检查文件并规范路径。
2. 生成 `<prefix>_slam_map` 与 `<prefix>_localization`，任一目录已存在就拒绝启动。
3. 构造 `common_args`，把相同的 bag、播放节奏、CPU 和停止条件交给两阶段。
4. 依次执行建图、检查产物、定位，再写两阶段路径与定位分析结果的 summary。

默认值与共享参数分别位于：

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_sany_offline_localization.sh sany-offline-settings

<details>
<summary>共享参数与结果解释</summary>

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_sany_offline_localization.sh sany-offline-shared

`--mapping-config` 与 `--localization-config` 可分别指定；兼容 `--config` 会同时设置两者。默认数据路径属于 SANY 四雷达场景，换数据时应覆盖。这个入口每次先生成地图；直接使用已有地图应调用离线定位 runner。

参考轨迹来自同一建图阶段，用于两阶段结果对照，不是独立真值。进入 C++ 的位置分别在两个子脚本的核心调用。完整源码：@ref reproduction/multi_lidar/sany_4livox/run_sany_offline_localization.sh "run_sany_offline_localization.sh"。

</details>

@anchor phase_a_matrix_guide
# run_phase_a_relocalization_matrix.sh：改变传感器起始时间

## 先看核心调用

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_phase_a_relocalization_matrix.sh phase-a-launch

每个任务调用 @ref loc_offline_guide "run_loc_offline.sh"，最多处理 90 个融合帧，并关闭 YAML 初始位姿。`--start-sensor-time` 改变数据读取起点；这里的 offset 单位是秒，不是初始位置或姿态的扰动。

## 按执行顺序阅读

| 阶段 | 目的 | 阅读入口 |
|---|---|---|
| 1. 确定矩阵 | bag root、配置、地图、输出根、数据集与 offset | @ref phase_a_inputs "矩阵输入" |
| 2. 每个数据集先跑 offset=0 | 从结果统计取得实际首个定位时间 | @ref phase_a_first_time "基准时间" |
| 3. 后续 offset 加到基准时间上 | 变成算法需要的绝对传感器时间 | @ref phase_a_start_time "时间换算" |
| 4. 顺序运行并统计失败 | 不覆盖旧目录，runner 失败后继续矩阵 | 核心调用 |

@anchor phase_a_inputs
### 1. 矩阵输入

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_phase_a_relocalization_matrix.sh phase-a-inputs

@anchor phase_a_first_time
### 2. offset=0 的结果提供基准时间

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_phase_a_relocalization_matrix.sh phase-a-first-time

首个任务结束后，从 `localization_stats.csv` 第一条数据取时间。无法读取就退出；其他 offset 依赖这个结果。

@anchor phase_a_start_time
### 3. 相对 offset 转为绝对时间

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_phase_a_relocalization_matrix.sh phase-a-start-time

尽管基准时间的读取在首个任务之后，后续循环会先完成这段时间换算，再进入核心调用。输出目录按数据集和 `offset_NNs` 分开。

<details>
<summary>失败与汇总行为</summary>

脚本使用 `set -uo pipefail`，runner 放在 `if` 条件中；单次失败增加 `runner_failures`，矩阵继续。缺 bag、已有输出或无法取得基准时间会立即退出；正常走到末尾只打印失败总数，没有据该计数返回非零，也没有自动生成跨任务分析报告。

进入 C++ 继续跟踪离线定位程序的起始时间过滤与重定位初始化。完整源码：@ref reproduction/multi_lidar/sany_4livox/run_phase_a_relocalization_matrix.sh "run_phase_a_relocalization_matrix.sh"。

</details>

@anchor fault_bags_guide
# create_sany_formal_fault_bags.sh：生成四个故障输入

## 先看核心调用

@snippet{lineno} reproduction/multi_lidar/sany_4livox/create_sany_formal_fault_bags.sh fault-bags-creator

Shell 中的 `run()` 把场景名、目标话题、模式、相对开始时间和持续时间传给 `create_sany_fault_bag.py`。定义函数不会生成 bag，后面的四次调用才执行：

@snippet{lineno} reproduction/multi_lidar/sany_4livox/create_sany_formal_fault_bags.sh fault-bags-scenarios

## 执行顺序与职责交接

加载 Humble → 消费 `INPUT OUTPUT_ROOT` → 建立输出根 → 按上面顺序生成四个场景。Python 工具负责读原始 bag、执行话题丢弃/中断并保存故障 bag 和合同；Shell 不连接现场传感器，也不启动算法。

后续阅读进入同目录的 `create_sany_fault_bag.py`；生成的数据由 @ref fault_matrix_guide "故障矩阵" 交给在线合同运行器。`drop_all` 与 `interrupt` 的实际处理规则以 Python 实现和生成合同为准，不能仅用函数参数中的 duration 推断丢弃范围。

<details>
<summary>输入、输出与退出</summary>

输入是完整 SANY 原始 bag，输出根下按四个场景名分目录。严格模式下任一次 Python 生成失败就停止后续场景。话题与时间条件是这组正式实验的固定定义。完整源码：@ref reproduction/multi_lidar/sany_4livox/create_sany_formal_fault_bags.sh "create_sany_formal_fault_bags.sh"。

</details>

@anchor fault_matrix_guide
# run_sany_formal_fault_matrix.sh：遍历四个故障场景

## 先看核心调用

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_sany_formal_fault_matrix.sh fault-matrix-run

每个场景把 bag、YAML、独立输出目录和 `fault_contract.json` 交给 @ref sany_online_contract_guide "在线合同 runner"。单次失败会计数并继续；最后 `(( failures == 0 ))` 决定矩阵整体退出状态。

## 执行顺序与矩阵来源

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_sany_formal_fault_matrix.sh fault-matrix-scenarios

消费 `FAULT_ROOT CONFIG OUTPUT_ROOT` → 建立输出根 → 遍历固定场景 → 调子脚本 → 在终端汇总成功/失败。算法进程、回放、录制与报告由子 runner 管理，矩阵脚本只负责调度。

<details>
<summary>数据与结果目录</summary>

四个场景对应 @ref fault_bags_guide "故障 bag 生成器"。每个输出目录保存在线合同 runner 的原始日志、录制 bag、恢复合同和 metadata；矩阵层的总计写到终端。完整源码：@ref reproduction/multi_lidar/sany_4livox/run_sany_formal_fault_matrix.sh "run_sany_formal_fault_matrix.sh"。

</details>

@anchor sany_online_contract_guide
# run_sany_online_contract.sh：在线前端接口与故障恢复实验

## 先看核心：录制 → 前端 → 回放

这是在线**前端**的实验入口，`template` 指向 `run_frontend_online.sh`。先启动输出录制，再启动节点：

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_sany_online_contract.sh sany-contract-record

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_sany_online_contract.sh sany-contract-node

节点和 recorder 的 ROS 图满足门槛后，才启动输入 bag：

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_sany_online_contract.sh sany-contract-playback

这三段是主线。`setsid/taskset/&` 管理独立进程组、CPU 与后台运行，三个 PID 分别用于停止及退出状态记录。前端的 C++ 入口从 @ref frontend_online_guide "在线前端子脚本" 继续进入。

## 按执行顺序阅读

| 阶段 | 目的 | 阅读入口 |
|---|---|---|
| 1. 检查输入与来源 | 确认 bag/YAML/二进制、Git 元数据和新输出目录 | 完整源码开头的检查 |
| 2. 隔离环境与冻结指纹 | 确认安装 prefix、clean worktree、独立 domain；绑定输入和工具 | @ref sany_contract_environment "环境与分析器" |
| 3. 录制并启动前端 | 先保留输出，避免漏掉启动消息 | 核心 record/node |
| 4. 等待 ROS 图再回放 | 输入订阅与输出录制均已就绪 | @ref sany_contract_readiness "回放门槛" |
| 5. 结束输入、节点、录制 | 依次排空和停止进程，保存退出码 | @ref sany_contract_stop "停止顺序" |
| 6. 分析并归档 | 检查公开 topic 或故障恢复，写产物清单 | @ref sany_contract_analysis "分析与结果" |

@anchor sany_contract_environment
### 1–2. 环境与分析器

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_sany_online_contract.sh sany-contract-selector

有故障合同就选择故障恢复分析器，否则选择公开接口分析器。安装程序比入口源码旧时拒绝运行；输出目录必须新建，且位于仓库与输入 bag 之外。

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_sany_online_contract.sh sany-contract-environment

这个正式 runner 要求工作区干净、所选 domain 没有既存节点，随后把 bag、配置、二进制和工具指纹写入 metadata。提供故障合同时还核对合同与实际 bag 的文件哈希。

@anchor sany_contract_readiness
### 4. 回放前的 ROS 图门槛

等待四路 LiDAR 和主 IMU 的订阅，以及公开输出与 recorder 的连接；30 秒内未就绪则结束。详细条件按需展开：

<details>
<summary>查看就绪条件的真实源码</summary>

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_sany_online_contract.sh sany-contract-readiness

</details>

通过后执行核心 playback，1 倍播放并发布 100 Hz 时钟，播放器有时长加余量的 timeout。`wait` 获取 `player_rc` 后额外等待 3 秒。

@anchor sany_contract_stop
### 5. 先停节点，再停录制

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_sany_online_contract.sh sany-contract-stop

录制在节点停止之后才结束，使节点最后发布的状态有机会被保存。`stop_group()` 逐步发送 INT/TERM/KILL；EXIT trap 也用它处理异常，并将未完成 metadata 标为 aborted。

@anchor sany_contract_analysis
### 6. 分析与结果

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_sany_online_contract.sh sany-contract-analysis

分析器完成后，Python 段更新 metadata、记录四类退出码，并生成带文件哈希的 artifact manifest；最后要求四类状态全部成功：

@snippet{lineno} reproduction/multi_lidar/sany_4livox/run_sany_online_contract.sh sany-contract-finish

<details>
<summary>输入、固定条件与产物</summary>

调用方式是 `BAG CONFIG OUTPUT_DIR [FAULT_CONTRACT]`。默认 CPU `0-7`、`SANY_ROS_DOMAIN_ID=217`，固定 localhost 通信。ready 条件使用此四雷达场景的固定话题，不能当作任意 YAML 的通用入口。

录制位于 topic_record/，包含五个公开前端 topic 与时钟。普通任务输出 topic contract/timing，故障任务输出 fault recovery contract；两者都保存环境身份、原始日志、运行指纹和产物清单。完整源码：@ref reproduction/multi_lidar/sany_4livox/run_sany_online_contract.sh "run_sany_online_contract.sh"。

</details>

@anchor fastlivo_rerun_guide
# rerun_m3dgr_fastlivo2_timestamp_fix.sh：重跑冻结实验中的受影响单元

## 先看核心调用

@snippet{lineno} reproduction/formal_report/rerun_m3dgr_fastlivo2_timestamp_fix.sh fastlivo-rerun-launch

这里调用的是外部 FAST-LIVO2 工作区中的 runner，`BENCH_*` 固定播放、CPU、inventory、实验指纹与 ROS 端口。本仓库脚本负责恢复实验计划，算法与单次运行管理继续在外部子脚本中阅读。

## 按执行顺序阅读

1. 读取可选输出 base，使用脚本中冻结的 12 个 sequence/repeat 单元。
2. 每个单元计算独立端口、bag 与 run dir，先查看已有 metadata。
3. 已有有效结果则跳过；已有目录但结果不满足条件则退出，不覆盖旧结果。
4. 调外部 runner，读取完成度、时间合法性和合同结果，不符合要求就停止；符合则继续下一单元。

实验计划位于：

@snippet{lineno} reproduction/formal_report/rerun_m3dgr_fastlivo2_timestamp_fix.sh fastlivo-rerun-matrix

<details>
<summary>已有结果与失败判定</summary>

@snippet{lineno} reproduction/formal_report/rerun_m3dgr_fastlivo2_timestamp_fix.sh fastlivo-rerun-existing

@snippet{lineno} reproduction/formal_report/rerun_m3dgr_fastlivo2_timestamp_fix.sh fastlivo-rerun-finish

runner 返回码会显示；是否继续的条件主要来自 metadata 合同，不直接以 `runner_rc` 判定。脚本的固定路径、指纹和单元顺序属于对应报告的恢复计划。完整源码：@ref reproduction/formal_report/rerun_m3dgr_fastlivo2_timestamp_fix.sh "rerun_m3dgr_fastlivo2_timestamp_fix.sh"。

</details>

@anchor voxel_backend_guide
# run_voxel_slam_full_backend.sh：M3DGR 完整后端基线

## 先看核心：直接启动 Voxel 二进制

@snippet{lineno} reproduction/single_lidar/m3dgr/run_voxel_slam_full_backend.sh voxel-backend-launch

启动前已经在独立 ROS master 中加载 YAML 并开启后端与地图保存。节点订阅就绪后开始资源采样，再播放固定的 Mid360 输入：

@snippet{lineno} reproduction/single_lidar/m3dgr/run_voxel_slam_full_backend.sh voxel-backend-playback

播放器结束后设置 `/finish=true`。这条入口直接运行外部 `voxelslam`，不是前面 114 参考入口的 roslaunch 调用；数据、配置与输出合同也不同。

## 按执行顺序阅读

| 阶段 | 目的 | 阅读入口 |
|---|---|---|
| 1. 解析路径与检查依赖 | bag、序列、空输出目录、配置、inventory、工具 | @ref voxel_backend_settings "输入与条件" |
| 2. 配置 ROS 1 环境、注册收尾 | 设置私有 master、仿真时间和后端参数 | @ref voxel_backend_master "master 与参数" |
| 3. 启动算法、等订阅、采样 | 保证回放输入有接收者 | 核心 launch |
| 4. 回放、请求最终优化 | 让后端处理完整数据并写结果 | 核心 playback |
| 5. 等保存、转换与统计 | 等 alidarState，提取 TUM 与回环事件 | @ref voxel_backend_results "结果处理" |
| 6. 停进程、归档并判定 | 停采样/算法/master，写 metadata | @ref voxel_backend_finish "收尾与退出" |

@anchor voxel_backend_settings
### 1. 输入与条件

@snippet{lineno} reproduction/single_lidar/m3dgr/run_voxel_slam_full_backend.sh voxel-backend-inputs

参数为 `BAG SEQUENCE OUTPUT_DIR [REPEAT]`，非空旧目录拒绝覆盖。`VOXEL_SLAM_WS` 指定外部 ROS 1 工作区，`BENCH_*` 提供实验条件。

@anchor voxel_backend_master
### 2. master 与后端参数

加载 Noetic 与 Voxel overlay、设置独立 master URI 后：

@snippet{lineno} reproduction/single_lidar/m3dgr/run_voxel_slam_full_backend.sh voxel-backend-master

这段把后端、地图保存和 save path 明确写入参数服务器；可选的 `VOXEL_LOOP_ICP_EIGVAL` 覆盖回环阈值，然后进入核心启动。

@anchor voxel_backend_results
### 5. 结果处理

请求 finish 后，脚本寻找非空 `alidarState.txt`，随后等文件大小连续稳定。等待逻辑按需展开：

<details>
<summary>查看优化结果等待过程</summary>

@snippet{lineno} reproduction/single_lidar/m3dgr/run_voxel_slam_full_backend.sh voxel-backend-optimization

</details>

取得结果后导出优化 TUM，并从算法日志提取回环候选：

@snippet{lineno} reproduction/single_lidar/m3dgr/run_voxel_slam_full_backend.sh voxel-backend-extract

随后用冻结 inventory 检查末帧和帧数，统计轨迹时间合法性、间隔和 PCD 数量。C++ 算法入口位于外部 Voxel 工作区的 `voxelslam`。

@anchor voxel_backend_finish
### 6. 收尾与退出

@snippet{lineno} reproduction/single_lidar/m3dgr/run_voxel_slam_full_backend.sh voxel-backend-stop

取消 EXIT trap 后写 metadata，再判定是否完成：

@snippet{lineno} reproduction/single_lidar/m3dgr/run_voxel_slam_full_backend.sh voxel-backend-finish

<details>
<summary>产物与固定话题</summary>

输入固定 `/livox/mid360/lidar` 与 `/livox/mid360/imu`。优化轨迹在 results/，回环事件在 data/new_map/backend_diagnostics/，Voxel 原始输出在 voxel_output/，资源记录与 metadata 在运行根目录。完整源码：@ref reproduction/single_lidar/m3dgr/run_voxel_slam_full_backend.sh "run_voxel_slam_full_backend.sh"。

</details>
