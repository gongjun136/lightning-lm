# SANY 三雷达在线定位：第一批算力裁剪

基线：`82eed51`，主 Orin 为 MAXN、12 核在线、动态调频。目标是在感知、规划、控制同时运行时持续处理 10 Hz LiDAR 观测，连续跟踪 P95 ≤80 ms、P99 <100 ms，不以固定跳帧换实时性。

## 已实现的计算路径

1. LIO 观测模型用固定 64 点块累加 Hessian/gradient，复用缓存；残差上中位数改为 `nth_element`，最大值单独统计。块顺序固定，不受 OpenMP 实际线程数影响。残差、权重和 ESKF 接受条件保持原公式，但浮点累加次序改变，需要轨迹回归。
2. LIO 在原体素滤波及低点数保护之后执行硬点数上限。按来源、8 个方位扇区、高低和远近分区分配配额，选取现有点、不额外生成质心，保留时间戳、来源及输入顺序。空间配额用于保持几何覆盖，不是可观测性证明。
3. NDT 使用独立点数上限，按原始点数比例分配空间配额，避免均分配额改变场景权重。初始定位与匹配失败后的恢复帧使用完整输入。NDT 分辨率、迭代上限和置信门槛保持原配置。
4. 新预算模式在定位与 LIO 同时健康时启用；LIO 不健康时清除降级档，地图定位不健康时恢复全来源及无上限观测。匹配前依据帧龄和近期耗时预测降档，恢复经过 75 帧健康窗口。
5. 实验用 `balanced` 档始终使用前雷达，并轮换左右雷达。三雷达正常到达且跟踪健康时，侧雷达各 5 Hz；单雷达缺失或恢复阶段不能承诺该间隔。完整发布点云保留三雷达。**此档已在回放中出现航向/RPE 超门槛，本次不放行部署。**

2026-09-08 18:06 域控包复核后，正式三雷达 SOLiD YAML 已对齐 `conservative` 候选：LIO `[2200, 1804, 1500]`、NDT 3500，保留全部三雷达且不轮换。此前提交只包含实现和候选生成器，正式 YAML 没有启用点数上限；10:46 包使用候选路径，而18:06包使用正式路径，不能混作相同配置的前后对照。四雷达和其他配置不改变。此调整仍需主 Orin 全栈实时性验收，不代表生产放行。

外层脚本现在可以继续指向 `config/reproduction/multi_lidar/sany_3livox/sany_3lidar_localization_solid.yaml`，不用改成候选路径。启动输出及 `launch_compute_config.json` 会显示所选预算和配置 SHA256；真实逐帧档位以 `LIO_BENCH_FRAME.point_budget` 为准。新增诊断说明见 [sany_resource_diagnostics.md](sany_resource_diagnostics.md)。下面的候选生成器继续保留，供显式消融或独立配置使用。

## 生成与运行

先在独立候选部署目录应用源码补丁并重新编译。补丁基于 `82eed51`，不是 ARM 二进制；不要复制本机 WSL 的可执行文件到 Orin。保留原部署目录用于回退。

```bash
git rev-parse --short HEAD
git status --short
git apply --check /path/to/sany_compute_candidate.patch
git apply /path/to/sany_compute_candidate.patch
source /opt/ros/humble/setup.bash
# 继续 source 现场原有的驱动/依赖工作空间 setup.bash。
colcon build --packages-select lightning --cmake-args -DCMAKE_BUILD_TYPE=Release
```

如 `git apply --check` 报冲突，停止应用，不要强制覆盖现场修改。

在部署仓库中生成候选，`--base` 与 `--output` 必须不同：

```bash
python3 scripts/prepare_compute_pruning_config.py \
  --base config/reproduction/multi_lidar/sany_3livox/sany_3lidar_localization_solid.yaml \
  --output config/sany_compute_candidate.yaml \
  --profile conservative --lio-points 2200 --minimum-lio-points 1500 --ndt-points 3500
```

`conservative` 保留三雷达，默认 LIO 预算为 `[2200, 1804, 1500]`。中间档为上限的 82% 且不低于指定下限；质量不健康时恢复完整点云。`balanced` 必须额外指定 `--allow-experimental-rotation` 才能生成，仅用于后续实验。

外层 `gj_ws/run.sh` 可写为：

```bash
#!/usr/bin/env bash
set -e
export LIGHTNING_LM_RUN_MODE=diagnostic
export LIGHTNING_LM_CONFIG="$PWD/lightning-lm/config/sany_compute_candidate.yaml"
export LIGHTNING_LM_LIO_THREADS=6
export LIGHTNING_LM_NDT_THREADS=4
export LIGHTNING_LM_NDT_MAX_POINTS=3500
export SANY_MAP_PATH=/home/nvidia/project/gj_ws/maps/sany_3lidar_baoma_sim_20260901
export SANY_LIDAR_LAYOUT=3
export SANY_RECORD_BAG=0
bash lightning-lm/scripts/run.sh
```

配置消融对照可把 `LIGHTNING_LM_CONFIG` 指回正式三雷达 YAML；**严格的算法基线必须使用 `82eed51` 的独立构建**，因为新二进制包含观测热循环及健康反馈修复，仅切换 YAML 不能还原旧算法。NDT 上限也可以用 `LIGHTNING_LM_NDT_MAX_POINTS` 覆盖；设为 `0` 关闭，非零合法范围为 100–1000000。线程数仍由既有 `LIGHTNING_LM_LIO_THREADS`、`LIGHTNING_LM_NDT_THREADS` 覆盖。

首轮保持 6 个 LIO 线程、4 个 NDT 线程，以便单独判断算法收益。不要同时改地图、外参、噪声、匹配接受门槛或锁定时钟。

## 回归与验收

诊断导出现在同时保留 `LIO_BENCH_FRAME` 和 `COMPUTE_BENCH_*`。LIO、组帧、NDT 和在线定位使用纳秒整数 `frame_id` 关联；逐帧 LIO 日志增加实际观测调用数、有效特征数及点数预算。

生产模式仍保留每 10 秒一条 `module=lidar_end_to_end`：范围为主雷达回调入口（预处理之前）到 NDT、PGO 处理完成，包含两级排队和 LIO。它不包含驱动/DDS 在回调之前的延迟，也不代表下游订阅者实际收到消息的时刻。诊断模式额外记录逐帧 `primary_to_pgo_ms`；值为 -1 表示未关联到入口记录，不能纳入有效延迟样本。完整指标需联合观测率、累计超时数和输入/定位队列丢弃数判断。

在线与离线均使用 NDT 的 `lidar_loc_valid_` 判定健康（`valid_` 是融合/DR 输出标志），离线回放也反馈定位健康状态，因此新的预算/轮换路径可以实际执行。`evaluate_compute_pruning.py` 使用同一 map 坐标系、1 ms 时间关联容差，统计位置差、航向差和 1/10 秒相对位姿误差，不进行拟合对齐。

```bash
python3 scripts/evaluate_compute_pruning.py \
  --baseline /path/to/baseline_run \
  --candidate /path/to/candidate_run \
  --output /path/to/comparison.json \
  --expect-point-budget 2200 --expect-lidar-count 3 --enforce-relative-gates
```

每个输入目录含 `stderr.log`、`frames.csv`、`ndt.tum`。这些结果只能证明相对基线的回归表现；主 Orin 资源、实时性及全栈相互影响需要目标机采集。

门槛：位置差 P95/P99 ≤0.15/0.30 m，航向差 P95/P99 ≤0.3/0.8°；1 秒 RPE P95 ≤0.05 m/0.2°，10 秒 RPE P95 ≤0.15 m/0.5°。要求时间关联覆盖率 ≥99.9%、初始化后不新增无效定位或匹配失败、预算和三雷达路径至少在 95% 跟踪帧实际生效。时间关联覆盖率不是点云重叠率；不设置统一的场景重叠率 0.98 门槛。

主 Orin 的基线与候选应交替运行，每个候选至少三次，每次分别检查 10 Hz、延迟长尾、跟踪丢帧、定位失效和资源占用。100 Hz 位姿发布队列的 latest-only 覆盖单独统计，不等同于 LiDAR 观测丢帧。冷启动和真正定位丢失后的全局重定位独立验收。

任何精度回归先退回上一档点数或完整三雷达。当前版本不包含减少 ESKF 迭代、对应关系跨迭代复用、PGO 调度重构或 CUDA 后端；这些变更需在第一批目标机实测后根据剩余热点决定。

## 本机回放结果（2026-09-08）

以下为 WSL/x86 离线结果，不是主 Orin 验收。基线使用修改前冻结的可执行文件及共享库，记录二进制、库、配置 SHA256。对照保持同一地图、6/4 线程、`OMP_WAIT_POLICY=PASSIVE`、诊断开关、外参及噪声配置。没有拟合对齐参考轨迹，也没有把基线当成真值。

157 秒数据：`rosbag2_2026_08_31-16_59_08`，每次完整组帧 1573 次、轨迹 1550 帧；三种候选均未新增初始化后失效或组帧丢弃，但轮换档精度不通过。

| 配置 | LIO 核心 P95/P99（ms） | NDT 匹配均值（ms） | NDT 位置差 P95（m） | NDT 航向差 P95（°） | 相对精度门槛 |
|---|---:|---:|---:|---:|---|
| 原基线 | 14.781 / 17.451 | 5.429 | — | — | 参考 |
| 三雷达，2200 点上限 | 8.718 / 10.099 | 2.822 | 0.0179 | 0.0761 | 通过 |
| 三雷达，固定 1500 点压力档 | 6.718 / 7.331 | 2.841 | 0.0199 | 0.0914 | 通过 |
| 前+侧轮换，2200 点上限 | 8.049 / 8.699 | 3.088 | 0.0759 | 0.4430 | 不通过，不部署 |

这里 LIO 核心为既有 `LIO_BENCH_FRAME.total_ms`：预处理加 IMU/滤波/匹配/地图更新分项之和，不含全部组帧、队列及日志耗时；NDT 仅匹配核。两个数字不能相加充当在线端到端时延。新增加的 `primary_to_pgo_ms` 专门用于目标机端到端验证。

三雷达 2200 点档在 99.87% 跟踪帧中实际启用点数预算，剩余为初始化过渡；全部保留三个来源。1 秒 RPE P95 为 0.0252 m / 0.1200°，10 秒为 0.0284 m / 0.1220°。LIO 核心 P95 降低约 41%，NDT 匹配均值降低约 48%。

按参考轨迹未来 1 秒位移 ≥0.1 m 或航向变化 ≥0.5° 判为运动，否则为静止，末尾无法形成完整 1 秒窗口的帧排除。运动 1313 帧、静止 227 帧；2200 点档两类 LIO P95 分别从 14.576/15.420 ms 降至 8.593/9.371 ms，位置差 P95 分别为 0.0186/0.0128 m。此分类是定位推断，不是外部静止真值。

轮换档的 1 秒旋转 RPE P95 为 0.272°，10 秒 RPE 为 0.173 m / 0.880°，均超过相应门槛，因此未通过并非只因单个峰值。早期 NDT 空间均分配额也出现偏差，已改为保留原点数比例配额后重新回放；这些失败版本没有被当作成功结果。

734 秒数据：`rosbag2_2026_09_04-13_48_09`，基线和 2200 点候选均完成 7313 次三雷达完整组帧，轨迹 7290 帧、关联覆盖率 100%，初始化后无新增定位失效或匹配失败。

| 配置 | LIO 核心 P95/P99（ms） | NDT 匹配均值（ms） | NDT 位置差 P95/P99（m） | NDT 航向差 P95/P99（°） | 相对精度门槛 |
|---|---:|---:|---:|---:|---|
| 原基线 | 16.760 / 19.563 | 5.191 | — | — | 参考 |
| 三雷达，2200 点上限 | 10.281 / 11.560 | 2.752 | 0.0140 / 0.0189 | 0.0545 / 0.0800 | 通过 |
| 三雷达，固定 1500 点压力档 | 7.924 / 9.308 | 2.765 | 0.0159 / 0.0222 | 0.0666 / 0.0929 | 通过 |

长轨迹 2200 点档的 1 秒 RPE P95 为 0.0195 m / 0.0791°，10 秒为 0.0211 m / 0.0932°；预算生效率 99.97%。运动 3309 帧、静止 3971 帧（末尾 10 帧不分类），两类位置差 P95 分别为 0.0157/0.0117 m，航向差 P95 为 0.0642/0.0447°，没有用静止段掩盖运动段误差。

固定 1500 点档在同一长轨迹也完成 7313 次完整组帧、关联 7290 帧，无新增初始化后失效或组帧丢弃。1 秒 RPE P95 为 0.0227 m / 0.0967°，10 秒为 0.0242 m / 0.1142°。正常档和最低预算均通过；在线自动降档/恢复仍需主 Orin 在真实争用下验证。

除表中 NDT 轨迹，还独立检查了 `fused.tum` 融合输出，使用完全相同的时间关联和误差门槛。2200 点档的短/长轨迹位置差 P95 分别为 0.0201/0.0155 m、航向差为 0.0820/0.0578°；1500 点档分别为 0.0224/0.0178 m、0.0987/0.0707°。四组融合输出的覆盖率、P99 和 1/10 秒 RPE 也全部通过。

构建已覆盖在线、离线入口；7 项 C++ 测试通过，涵盖 ESKF 速度门控、多雷达采样及回退、PGO、输入锁、计时、预算配置和坐标投影。Python 评估器另测了同轨迹、恒定平移、重复时间戳、失效计数及配置生成边界。

### 生产模式在线链路检查

使用同一 157 秒 bag、内嵌在线 1 倍速回放、关闭逐帧 profiling 和非必要诊断 I/O，完整运行至退出。15 个已完成的 10 秒统计窗口中，健康 LiDAR 更新率中位数 9.9977 Hz、最低 9.9869 Hz；各窗口端到端 P95 为 22.55–37.25 ms，窗口 P99 最大 53.63 ms，累计 ≥100 ms 的健康跟踪样本为 0。它们是窗口统计，不能拼成一次全程 P95/P99。

传感器队列累计丢弃为 0；定位队列在首个统计窗口已经累计丢弃 7 帧，之后所有窗口保持 7。日志同时记录了冷启动全局检索和确认，不将这次回放描述为“全程零丢帧”。全局轨迹导出 1541 帧、跨度约 154 秒，高频内部轨迹导出 26978 帧。生产模式按既有逻辑关闭发布轨迹逐帧落盘，因此 `published.tum` 无位姿行，不据此宣称下游接收链路已通过验收。

此检查验证了异步队列与端到端计时的接线，以及关闭详细诊断后的可观测性。它不包含感知/规划/控制并发，不等同于主 Orin 实车运行，也没有验证下游实际收到消息的延迟。最终是否放行仍以主 Orin 全栈 A/B、动态降档/恢复、跨场景和恢复定位结果为准。
