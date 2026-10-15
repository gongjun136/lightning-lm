@page lidar_ins_accuracy 激光建图与定位对 CGI-430 的精度评估

本页定义旁路采集与离线工具的操作契约。完整 SANY 项目方案保存在 AlgoResearchVault 的“CGI-430参考的建图与纯雷达定位精度评估方案”；现场数值须来自实际车辆标定和录制。工具不接入生产启动，不向 CAN 发送数据，不发布业务位姿或 TF。

## 评价对象与实验分工

1. A：记录建图输入和 CGI 原始参考，在离线运行 @ref offline_slam_flow "建图流程"，导出前端、前端关键帧、最终优化关键帧和辅助稠密修正轨迹。
2. B：用 A 地图独立行驶标定路线，记录纯 LiDAR 地图定位和 CGI。只用 B 拟合固定尺度 map → ENU 变换并冻结。
3. 用冻结变换评价 A；用相同变换评价另外录制的 C 纯 LiDAR 定位路线。不得分别最佳配准每路测试轨迹，不能用测试误差反调时间偏移/外参/质量门限。

CGI 不进入被评价算法的融合、初始化或接管。原有 LiDAR 模式是否包含轮速/IMU约束，应按实际配置记录，不能把现有多传感器算法描述成仅使用点云。

## 环境与路径

在 ROS 2 Humble 的源码仓库根目录执行。采集和 rosbag 导出需要 `rosbag2_py`、`rclpy`、`rosidl_runtime_py`、`cgi430_interfaces`、`geometry_msgs`；参考采集还需要 `can-utils`。MCAP 需要 `rosbag2_storage_mcap`；离线核心只需要 NumPy 和 PyYAML；绘图另需 Matplotlib。

```bash
source /opt/ros/humble/setup.bash
# 指向实际已构建工作区，勿直接沿用别的机器的集中路径配置。
source "$LIGHTNING_LM_WS/install/setup.bash"
cd "$LIGHTNING_LM_REPO_DIR"
python3 -m pip install -r scripts/accuracy_eval/requirements.txt
# 可选绘图：使用系统已有 matplotlib，或在独立环境安装。
```

新增脚本从源码 checkout 运行，不需要修改 C++ 算法或生产入口。`config/accuracy_capture.yaml`、`config/accuracy_reference.yaml`、`config/accuracy_lever_sweep.yaml` 和 `config/accuracy_normalize.yaml` 是需要核验的模板，空值不代表零。

## 采集配置与命令

先复制采集模板到实验目录，为两台机器填写相同 `source_recording_id` 与具名字符串 `epoch`。A/B/C必须不同ID。确认 `capture.confirmed/evidence`、当前实际输入Topic和从机CAN速率。原始采集不要求预先证明测量时延、ENU或参考质量；这些在正式导出时核验，避免无法录制标定数据。

```bash
# 从 Orin：CGI 接收驱动应先独立启动，接口应已按实际设备速率配置。
ip -details -statistics link show can1
bash scripts/accuracy_capture.sh record --config "$EXP/configs/capture_A.yaml" \
  --role reference --output "$EXP/A_slave" --dry-run
bash scripts/accuracy_capture.sh record --config "$EXP/configs/capture_A.yaml" \
  --role reference --output "$EXP/A_slave"

# 主 Orin：录制算法原始输入和已有业务/诊断输出，不再启动第二个定位器。
bash scripts/accuracy_capture.sh record --config "$EXP/configs/capture_A.yaml" \
  --role algorithm --output "$EXP/A_master"
```

把完整主从会话目录复制到分析机后再导出。录制器保留 typed CGI bag、所有 CAN 帧 `can_raw.log`、Topic类型/QoS、配置快照、开始/结束区间和时钟审计。`--duration-seconds N` 可设定录制时长；默认 Ctrl-C 正常结束。输出目录必须不存在。退出 0 表示正常完成；2 是参数/合同/预检问题；3 是子进程或录制失败。结束后检查 bag 的 `metadata.yaml` 和日志，不只检查目录是否存在。

`reference/algorithm` 各自本机录制为首选，`combined` 只在已验证跨机DDS实际接收和时钟后使用。采集器只检查当前CAN的UP状态/实际bitrate是否匹配，拒绝自行改变接口。不能并行启动 `ins_only` 作为参考：它与纯 LiDAR 发布相同业务Topic。CGI驱动自身可独立运行，因为只发布输入参考Topic。

## 导出原始参考与业务轨迹

正式导出前补齐 `georeference`、`height`、`reference_point`、`quality` 和 `estimate` 的实际数值/证据。固定 ENU 原点为 WGS84 纬度/经度/椭球高；设备若输出正常高，填写已确认的 `h=H+N` 转换。质量门限必须在查看算法误差以前确定；它们是参考准入条件，不是算法合格阈值。

推荐保留 CGI 的一个固定原始输出点 O，后轴转换全部离线进行。O 必须按 CGI 手册和实际输出配置确定，不能猜测它是天线或IMU中心。O 须与后轴刚性连接在同一后车体上，姿态仍统一为后车体 FLU。采集配置使用下列合同，并填写实际 `point_id` 和证据：

```yaml
reference_point:
  confirmed: true
  position: fixed_point
  point_id: FILL_VERIFIED_CGI_OUTPUT_POINT_ID
  attitude_axes: FLU
  device_lever_arm_already_applied: false
  rigid_to_rear_axle_confirmed: true
  evidence: FILL_POINT_ID_OUTPUT_CONFIGURATION_AXES_AND_RIGID_ATTACHMENT_EVIDENCE
```

这里的 `device_lever_arm_already_applied:false` 只声明设备尚未把 O 转到后轴。不要清零或改动设备内部已确认的 GNSS–IMU 杆臂和安装角，也不修改设备配置或业务发布。已有已补偿后轴的录制仍可使用旧 `position:rear_axle_center`、`device_lever_arm_already_applied:true` 合同；必须保留其真实来源，不能改标成固定原始点。

时钟和分帧历元尚未验证时只能显式接受条件性假设。`clock/assembly.confirmed:false` 必须配 `assumption_acknowledged:true`、显式整数偏移和说明。CAN Header 是内核接收时间；字段接收时间靠近不能证明它们来自同一测量历元。导出保留条件字段；核心不能用 `known offset=0` 消除这些缺口。

```bash
python3 scripts/accuracy_export_bag.py check-config \
  --config "$EXP/configs/capture_A.yaml" --purpose export
python3 scripts/accuracy_export_bag.py export --config "$EXP/configs/capture_A.yaml" \
  --reference-bag "$EXP/A_slave/bag" --reference-capture "$EXP/A_slave/capture.metadata.json" \
  --estimate-bag "$EXP/A_master/bag" --estimate-capture "$EXP/A_master/capture.metadata.json" \
  --output "$EXP/A"
```

主要产物：`reference.csv/reference.metadata.json`、`reference.raw_frames.csv`、`reference.timeline.csv`、`reference.status.jsonl`、`estimated.tum/estimated.metadata.json` 与摘要。原始 `can_raw.log` 另留在从机会话目录，typed raw CSV 不能替代未知/未成功解码CAN帧。所有拒绝参考行保留 `reference_valid=false` 和原因；标准 NavSatFix 的 GBAS 状态不能替代定向固定解、姿态、sigma和差分龄期检查。

参考列：`timestamp_ns,x,y,z,qx,qy,qz,qw,reference_valid,reason,epoch,time_basis`。固定点模式的 XYZ 为 O 的固定 ENU 位置，元数据为 `body_frame_id:cgi_output_flu`；xyzw 始终为 ENU←后车体 FLU 旋转。导出器不施加 O→后轴杆臂，这份 `reference.csv` 不能直接交给评估核心。旧后轴模式保持后轴 XYZ 和 `rear_axle_flu` 元数据。额外列保留九字段接收戳、质量和原始字节。元数据保留完整采集区间；不得用第一/最后有效位姿缩短分母。时钟回退须分成新的具名 epoch 和物理文件。

## CGI 参考点离线归一化

复制 `config/accuracy_reference.yaml` 到 `$EXP/configs/reference.yaml`，填写已测量并冻结的杆臂合同。`lever_arm_m:[x,y,z]` 是从固定原始点 O 指向目标后轴中心的完整向量，以后车体 FLU 表达：x向前、y向左、z向上。它不是 ENU 中的常向量。固定点转换为 `p_rear_ENU = p_O_ENU + R_ENU_from_rear_FLU * lever_arm_m`，姿态不变；不能只把车体向量直接加到 ENU XYZ。

| `mode` | 输入与离线位置变化 |
|---|---|
| `fixed_point` | 原始固定点 O；有效位移为 `lever_arm_m`。`source_point_id` 必须匹配原始元数据。 |
| `rear_axle_delta` | 设备已做末端后轴外推的旧录制；有效位移为 `lever_arm_m - previous_device_lever_arm_m`。两项均是 O→后轴的完整向量，并非在旧后轴位置上再加完整新杆臂。 |
| `already_rear` | 已确认的旧后轴数据恒等通过；只允许零 `lever_arm_m`，不再补偿。 |

各模式均需填写 `confirmed/evidence` 和非空 `source_point_id`（`already_rear` 可填 `rear_axle_center`），以及 `rigid_to_rear_axle_confirmed:true`。`rear_axle_delta` 额外需要已知 `previous_device_lever_arm_m`、`previous_compensation_is_rigid_translation_confirmed:true` 和 `previous_compensation_evidence`，证明设备原处理只是对同一 O 的末端刚性平移。新旧杆臂均按FLU表达，不能未经坐标转换直接照抄WebUI向量。无法证明时不能据此撤销设备内部校准。跨铰接点的动态外参不能用静态向量代替，本工具暂不支持。

```bash
python3 scripts/accuracy_reference.py normalize \
  --source "$EXP/A/reference.csv" --metadata "$EXP/A/reference.metadata.json" \
  --config "$EXP/configs/reference.yaml" --output "$EXP/A/reference_rear.csv"
```

输出 `reference_rear.csv/reference_rear.metadata.json`，归一化后 `body_frame_id:rear_axle_flu`，并记录来源、配置哈希和 `reference_normalization` 链。B、C分别从各自的原始 `reference.csv` 运行同样步骤；物理输出点和安装不变时复用同一冻结杆臂合同。已经有离线归一化链的输入会被拒绝，不能通过再跑一次补偿旧结果。

归一化只改位置，不改时间戳、姿态、参考有效性和完整录制区间；所有拒绝行和源点质量证据保留，源点准入及配对掩码不因候选杆臂改变。非零位移没有完整联合协方差可用于传播：保留的 sigma 只描述源点质量，不能称为后轴 sigma。元数据标记 `target_covariance_status:not_propagated`，报告携带 `REFERENCE_LEVER_ARM_COVARIANCE_NOT_PROPAGATED` 条件，不得把位置平移等同于已证明目标点不确定度。

## 建图轨迹与后轴转换

为新评价地图使用独立配置副本，声明 `map_export.alignment_status: not_required`，即本次原生地图无需转换到历史业务地图；这不表示它已经与ENU配准。默认 `unknown` 的地图会被定位启动检查拒绝。不要在生产配置中批量改该字段；在副本中将 `sany_vehicle_config.calibrations` 改为实际标定文件绝对路径，避免把YAML移到实验目录后丢失相对引用。记录实际VIN、雷达配置和配置哈希。

```bash
bash scripts/run_slam_offline.sh --bag "$EXP/A_master/bag" --config "$MAPPING_CONFIG" \
  --output-dir "$EXP/mapping_run_A" --output-map-dir "$EXP/map_A" \
  --output-tum "$EXP/mapping_run_A/results/trajectory_slam.tum"
```

`run_slam_offline.cc::SaveFrameTrajectoryTum/SaveKeyframeTrajectoryTum/SaveCorrectedFrameTrajectoryTum` 已导出以下同地图世界系的IMU轨迹，已应用相同地图归一化：

| 文件 | 使用方式 |
|---|---|
| `trajectory_slam.tum` | 密集建图前端 |
| `trajectory_slam_keyframes_lio.tum` | 与后端同一关键帧时刻的前端 |
| `trajectory_slam_keyframes_opt.tum` | 最终优化关键帧，后端主结论 |
| `trajectory_slam_opt.tum` | 插值关键帧最终修正重建密集轨迹；不是逐帧独立优化 |

`accuracy_normalize_tum.py` 只应用事先确认的刚体合同：IMU输入右乘 `T_IL*T_BL^-1`，主雷达输入右乘 `T_BL^-1`，已有后轴输出不再次补偿。不会用第一帧重置原点或估计安装角。需要完整三维安装外参；副雷达到主雷达的相对外参不能代替车体外参。

```bash
python3 scripts/accuracy_normalize_tum.py \
  --source "$EXP/mapping_run_A/results/trajectory_slam.tum" \
  --config "$EXP/configs/normalize_frontend.yaml" \
  --recording-metadata "$EXP/A/estimated.recording.metadata.json" --output "$EXP/A/frontend_rear.tum"
```

对前端关键帧、优化关键帧和辅助密集轨迹分别配置 `trajectory_kind` 并转换。模板的 `source_world_frame_id/world_frame_id` 要写真实坐标身份（建议 `map:<地图资产清单SHA256>`）。原始建图时间戳的 `timestamp_offset_ns` 应与采集配置的 `clock.estimate_header_offset_ns` 一致；已由导出器校正的业务TUM设为0，不能重复加偏移。录制边界用 frame-free `estimated.recording.metadata.json`，不要借用旧地图业务位姿元数据冒充新建图帧。如果业务后轴输出已经过旧地图发布变换，需要事先确认的 `world_transform` 左乘，将它转换到A导出地图系。它不能用这批测试数据拟合，且不能对已在目标地图系的轨迹重复施加。世界系变化和右侧参考点变化是两件事。

`/slamPoseRaw_topic` 或 `trajectory_published_rear_axle.tum` 是业务后轴位姿；内部 `trajectory_global.tum/trajectory_high_frequency.tum` 不可直接与CGI后轴比较。实时发布时延需另外使用接收审计，离线TUM仅能证明采样覆盖。

## 冻结配准与评估

复制 `scripts/accuracy_eval/templates/` 中B、A、C模板，它们分别读取 B/A/C 的 `reference_rear.csv/reference_rear.metadata.json`。补齐元数据路径、完整参考时间区间、具名epoch、已测量的时间偏移及最大插值间隔。声明 `t_reference=t_estimate+offset_ns`。不允许从测试轨迹寻找最佳时间偏移。B的路线应包含足够转弯、平面范围和姿态变化；近共线或静止路线不能唯一标定完整刚体变换。

```bash
python3 scripts/accuracy_eval/accuracy_eval.py validate --config "$EXP/configs/calibration_b.yaml"
python3 scripts/accuracy_eval/accuracy_eval.py fit --config "$EXP/configs/calibration_b.yaml" \
  --estimate calibration_localization --output "$EXP/frozen_map_to_enu.json"
python3 scripts/accuracy_eval/accuracy_eval.py evaluate --config "$EXP/configs/evaluation_a_mapping.yaml" \
  --output "$EXP/report_a"
python3 scripts/accuracy_eval/accuracy_eval.py evaluate --config "$EXP/configs/evaluation_c_localization.yaml" \
  --output "$EXP/report_c"
python3 scripts/accuracy_plot_report.py --report "$EXP/report_a/report.json"
python3 scripts/accuracy_plot_report.py --report "$EXP/report_c/report.json"
```

B标定残差仅是训练诊断；另留B验证路段或独立D路线评价其稳定性。A和C必须复用B同一冻结文件SHA256；更换地图/坐标/外参/授时后不能无审计复用。工具按输入内容去重、录制ID和区间拒绝校准/评价数据重叠；更改文件格式不能绕开检查。

`common_sample_groups` 给前端关键帧与最终优化关键帧取完全相同参考时刻交集，额外生成前后统计与CSV。各路总体统计仍保留，不用不同密度总体RMSE直接声称优化收益。

输出：中文 `report.md`、结构化 `report.json`、逐样本 `*_samples.csv`、完整时段 `*_coverage.csv` 和共同样本对比。图表命令生成PNG/PDF、HTML图册和输入哈希；图只画实际配对样本，不外推停发区间。

位置指标包括有符号XYZ/FLU、水平和3D范数；姿态主指标为SO(3)测地角，roll/pitch/yaw差值wrap至[-180°,180°)作为诊断。近万向节锁时保留标记。有符号误差的主分位读取 `p95_abs/max_abs`，不能把有符号负P95当精度上限。CGI是质量门控参考，差值不是雷达误差的无偏独立估计。

参考有效比例、估计采样覆盖、共同覆盖均以完整声明区间计时，停发不能变成零误差。稀疏关键帧的采样覆盖不是实时在线可用率；在线位姿的 `availability_hold_sec` 按真实发布合同确定。相对运动RPE、按实车事件自动工况分类和发布实时延迟是可扩展分析，本版没有把它们伪称为已计算指标。

`thresholds:null` 输出 `NOT_ASSESSED_NO_THRESHOLDS`，从不PASS；条件证据不能被配置清空。CLI退出0只说明处理完成。参考无效或时间/参考点/高度尚未成立时，报告条件或合同错误，不能强行给绝对精度结论。

## 杆臂敏感性扫描

扫描用于回答“保持其它条件不变时，不同候选杆臂使差值改变多少”。它不是正式精度验收，也不从测试误差中选择“最优杆臂”。正式杆臂须由实测或独立 B 标定确定并冻结；不得用 A/C 的最小误差反调杆臂。

复制 `config/accuracy_lever_sweep.yaml` 到 `$EXP/configs/lever_sweep.yaml`。`normalization_config` 指向冻结参考归一化合同，`evaluation_config` 指向 C 的评估配置，两条路径相对扫描配置所在目录解析；`candidates` 是具名 `{name, lever_arm_m}` 列表，每项仍是 O→后轴的完整 FLU 向量。候选必须有明确的工程假设，不能把未测量值标为正式标定值。

```bash
python3 scripts/accuracy_reference.py sweep \
  --source "$EXP/C/reference.csv" --metadata "$EXP/C/reference.metadata.json" \
  --config "$EXP/configs/lever_sweep.yaml" --output "$EXP/lever_sensitivity_C"
```

每项都从同一份原始参考重新归一化，复用 B 的同一冻结 map→ENU 变换、同一时间合同、参考有效性和配对样本；扫描不调用 `fit`，不针对每个候选重新配准，不排序挑选最优值，也不产生 PASS。输出是候选差值的敏感性证据，不能替代使用冻结正式杆臂生成的 A/C 精度报告。各候选子报告的实际 `report.json` 可交给 `accuracy_plot_report.py --report <子报告路径>` 绘图；协方差未传播等条件继续保留。

## 只回放输入的隔离流程

```bash
# 在另一终端先运行实验算法；同一隔离域、ROS_LOCALHOST_ONLY=1、use_sim_time=true。
bash scripts/accuracy_capture.sh playback --config "$EXP/configs/capture_C.yaml" \
  --bag "$EXP/C_master/bag" --output "$EXP/C_replay" --domain-id 221
```

回放只允许 `capture.input_topics`，不回放CGI参考、历史位姿、TF、状态或心跳。算法的use_sim_time不由采集器自动修改。重放元数据保留原bag时间区间，而不是电脑当前日期；未经完整结束的bag需审阅，不能将未完成录制标为成功。

## 验证入口

```bash
python3 -m unittest discover -s scripts/accuracy_eval/tests -v
python3 scripts/test_accuracy_normalize_tum.py
python3 scripts/test_accuracy_export_bag.py
python3 scripts/test_accuracy_reference.py -v
bash -n scripts/accuracy_capture.sh
# ROS storage smoke只写新临时目录，不连接硬件或发布生产话题。
python3 scripts/test_accuracy_export_bag.py --ros-smoke /tmp/accuracy_sqlite_unique --storage sqlite3
python3 scripts/test_accuracy_export_bag.py --ros-smoke /tmp/accuracy_mcap_unique --storage mcap
python3 scripts/accuracy_eval/tests/generate_demo.py --output "$DEMO_NEW_DIR" --run-core
```

合成演示显式标记SYNTHETIC，验证已知位置/姿态误差和停发分母；这些数字不能作为实车结果。配置未核验、未取得独立B与C和真实CGI时延证据以前，本版软件测试成功不等于车辆精度验收成功。
