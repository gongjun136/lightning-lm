# 离线轨迹精度评估工具

正式数据合同、A/B/C 实验、录制与操作说明见 [激光建图与定位对 CGI-430 的精度评估](../../docs/operations/lidar_ins_accuracy.md)。此目录只包含独立离线 Python 工具、配置模板和回归测试，不接入生产启动。

| 文件 | 用途 |
|---|---|
| `accuracy_eval.py` | 验证 canonical 输入、拟合并冻结一个 SE(3)、应用冻结结果评估。 |
| `templates/calibration_b.yaml` | B 独立地图定位标定数据，生成 map → ENU 冻结变换。 |
| `templates/evaluation_a_mapping.yaml` | A 建图前端与最终优化关键帧；重建稠密轨迹作为补充。 |
| `templates/evaluation_c_localization.yaml` | C 独立纯雷达定位；与 A 共用 B 的冻结变换。 |
| `templates/estimate.metadata.example.yaml` | 已转换至后轴 FLU 的 TUM 元数据示例。 |
| `tests/test_accuracy_eval.py` | 时间、质量、SE(3)、杆臂合同、采样缺失、去重和输出的语义测试。 |
| `tests/generate_demo.py` | 生成带明确标记的 A/B/C 数学合成数据，运行 CLI 并核对已知误差；不是实车或 CAN/ROS 验证。 |

在仓库根目录执行（模板中的空值必须先补齐）：

```bash
python3 -m pip install -r scripts/accuracy_eval/requirements.txt
python3 scripts/accuracy_eval/accuracy_eval.py validate --config EXPERIMENT/configs/calibration_b.yaml
python3 scripts/accuracy_eval/accuracy_eval.py fit --config EXPERIMENT/configs/calibration_b.yaml --estimate calibration_localization --output EXPERIMENT/frozen_map_to_enu.json
python3 scripts/accuracy_eval/accuracy_eval.py evaluate --config EXPERIMENT/configs/evaluation_a_mapping.yaml --output EXPERIMENT/report_a
python3 scripts/accuracy_eval/accuracy_eval.py evaluate --config EXPERIMENT/configs/evaluation_c_localization.yaml --output EXPERIMENT/report_c
python3 -m unittest discover -s scripts/accuracy_eval/tests -v
python3 scripts/accuracy_eval/tests/generate_demo.py --output NEW_SYNTHETIC_DEMO_DIR --run-core
```

成功返回 0；合同/数据/文件错误返回 2。`fit` 以独占创建写入冻结文件，`evaluate` 要求新的输出目录；不会覆盖既有证据。返回 0 表示命令完成，验收结论必须读取 `report.json` 的 `acceptance`。没有阈值时不会给出 PASS。

输出 `report.json`、中文 `report.md`、每路逐样本 `*_samples.csv` 与按完整声明录制时长分段的 `*_coverage.csv`。参考无效样本、未配对估计和估计停发均保留证据。工具不拟合测试段时间偏移、比例尺、杆臂或地图变换。
