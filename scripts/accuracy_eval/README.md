# 离线轨迹精度评估工具

正式数据合同、A/B/C 实验、录制与操作说明见 [激光建图与定位对 CGI-430 的精度评估](../../docs/operations/lidar_ins_accuracy.md)。此目录只包含独立离线 Python 工具、配置模板和回归测试，不接入生产启动。

| 文件 | 用途 |
|---|---|
| `accuracy_eval.py` | 验证 canonical 输入、拟合并冻结一个 SE(3)、应用冻结结果评估。 |
| `../accuracy_reference.py` | 将原始 CGI 固定输出点离线归一化到后轴；在同一冻结变换和配对合同下扫描候选杆臂。 |
| `../../config/accuracy_reference.yaml` | 参考点归一化合同；`fixed_point`、旧设备补偿差量和旧后轴恒等模式。 |
| `../../config/accuracy_lever_sweep.yaml` | 引用归一化与评估配置，列出具名候选杆臂；仅作敏感性分析。 |
| `templates/calibration_b.yaml` | B 独立地图定位标定数据，生成 map → ENU 冻结变换。 |
| `templates/evaluation_a_mapping.yaml` | A 建图前端与最终优化关键帧；重建稠密轨迹作为补充。 |
| `templates/evaluation_c_localization.yaml` | C 独立纯雷达定位；与 A 共用 B 的冻结变换。 |
| `templates/estimate.metadata.example.yaml` | 已转换至后轴 FLU 的 TUM 元数据示例。 |
| `tests/test_accuracy_eval.py` | 时间、质量、SE(3)、杆臂合同、采样缺失、去重和输出的语义测试。 |
| `tests/generate_demo.py` | 生成带明确标记的 A/B/C 数学合成数据，运行 CLI 并核对已知误差；不是实车或 CAN/ROS 验证。 |

在仓库根目录执行（模板中的空值必须先补齐）：

```bash
python3 -m pip install -r scripts/accuracy_eval/requirements.txt
# 先从每次录制的原始 CGI 参考产生后轴参考；B、C分别重复此步骤。
python3 scripts/accuracy_reference.py normalize --source EXPERIMENT/A/reference.csv --metadata EXPERIMENT/A/reference.metadata.json --config EXPERIMENT/configs/reference.yaml --output EXPERIMENT/A/reference_rear.csv
python3 scripts/accuracy_eval/accuracy_eval.py validate --config EXPERIMENT/configs/calibration_b.yaml
python3 scripts/accuracy_eval/accuracy_eval.py fit --config EXPERIMENT/configs/calibration_b.yaml --estimate calibration_localization --output EXPERIMENT/frozen_map_to_enu.json
python3 scripts/accuracy_eval/accuracy_eval.py evaluate --config EXPERIMENT/configs/evaluation_a_mapping.yaml --output EXPERIMENT/report_a
python3 scripts/accuracy_eval/accuracy_eval.py evaluate --config EXPERIMENT/configs/evaluation_c_localization.yaml --output EXPERIMENT/report_c
# 可选：对同一 C 原始参考做候选杆臂敏感性扫描。
python3 scripts/accuracy_reference.py sweep --source EXPERIMENT/C/reference.csv --metadata EXPERIMENT/C/reference.metadata.json --config EXPERIMENT/configs/lever_sweep.yaml --output EXPERIMENT/lever_sensitivity_C
python3 -m unittest discover -s scripts/accuracy_eval/tests -v
python3 scripts/test_accuracy_reference.py -v
python3 scripts/accuracy_eval/tests/generate_demo.py --output NEW_SYNTHETIC_DEMO_DIR --run-core
```

成功返回 0；合同/数据/文件错误返回 2。`fit` 以独占创建写入冻结文件，`evaluate` 要求新的输出目录；不会覆盖既有证据。返回 0 表示命令完成，验收结论必须读取 `report.json` 的 `acceptance`。没有阈值时不会给出 PASS。

输出 `report.json`、中文 `report.md`、每路逐样本 `*_samples.csv` 与按完整声明录制时长分段的 `*_coverage.csv`。参考无效样本、未配对估计和估计停发均保留证据。工具不拟合测试段时间偏移、比例尺、杆臂或地图变换。

杆臂扫描另输出 `sweep.json`、`sweep.md` 和 `summary.csv`；候选子报告位于 `<扫描目录>/<候选名>/report/report.json`。

固定原始点模式要求导出配置中的 `reference_point.position:fixed_point`、非空 `point_id`、`attitude_axes:FLU`、`device_lever_arm_already_applied:false`、`rigid_to_rear_axle_confirmed:true` 及已确认的证据。原始 `reference.csv/reference.metadata.json` 的位置保留 CGI 输出点，元数据为 `body_frame_id:cgi_output_flu`；这份固定点参考不能直接送入核心。输出点 O 必须由 CGI 手册和实际配置确认，不能假定它就是天线或IMU中心；保留设备内部已确认的 GNSS–IMU 杆臂和安装角。

归一化配置的顶层 `mode` 为 `fixed_point`、`rear_axle_delta` 或 `already_rear`，并填写 `confirmed/evidence`、非空 `source_point_id`、刚性连接证据和 `lever_arm_m`。`lever_arm_m:[x,y,z]` 表示 O→目标后轴的完整向量，以后车体 FLU 表达；固定点按当前姿态旋转该向量后加到 ENU 位置，姿态不改。旧已补偿数据仅在旧杆臂已知且原补偿确为末端刚性平移时，可用 `rear_axle_delta` 施加 `lever_arm_m - previous_device_lever_arm_m`；需额外确认 `previous_compensation_is_rigid_translation_confirmed:true` 并填写 `previous_compensation_evidence`。旧后轴恒等模式只允许零杆臂，`source_point_id` 可填 `rear_axle_center`。跨铰接点动态外参暂不支持静态向量。

归一化输出 `reference_rear.csv/reference_rear.metadata.json`，B/A/C模板均使用对应路径。时间、姿态、有效性、源点质量和完整录制分母保留，元数据记录来源哈希与 `reference_normalization` 链；已离线转换的输入不能再补偿。非零位移后的保留 sigma 只作源点质量证据，未传播为后轴协方差，标记 `target_covariance_status:not_propagated` 并携带 `REFERENCE_LEVER_ARM_COVARIANCE_NOT_PROPAGATED` 条件。

扫描配置的 `normalization_config/evaluation_config` 相对配置文件所在目录解析，`candidates` 为 `{name, lever_arm_m}` 列表。每个候选从同一原始参考开始，使用同一冻结 map→ENU、同一时间与配对样本，不调用 `fit`、不逐项重配准、不选最优值、不产生 PASS。扫描结果明确是敏感性分析；正式杆臂须实测或由独立 B 标定后冻结，不能根据 A/C 误差反调。候选子报告可用 `python3 scripts/accuracy_plot_report.py --report <候选report.json路径>` 绘图。合成测试和命令返回成功都不能代替实车验证。
