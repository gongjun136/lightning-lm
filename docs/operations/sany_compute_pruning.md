@page compute_budget 点数预算与回归

# 当前计算路径

前端固定 64 点块归约，降采样后按空间配额限制点数；定位 NDT 使用独立上限并按密度分配。预算仅在 LIO 和地图定位健康时适用，恢复阶段保留完整输入。源码和数据结构见 @ref sensor_pipeline "传感器预处理、多雷达组帧与算力预算"。

`conservative` 候选保持全部来源；来源轮换须显式启用实验选项。历史 balanced 回放出现航向/RPE 回归，不能当作已通过的部署配置。默认值以最终 YAML 和启动输出为准，不把旧批次点数当作永久默认。

# 创建候选与核验

在仓库根目录运行，base 与 output 必须不同：
```bash
python3 scripts/prepare_compute_pruning_config.py \
  --base config/reproduction/multi_lidar/sany_3livox/sany_3lidar_localization_solid.yaml \
  --output /tmp/sany_compute_candidate.yaml \
  --profile conservative --lio-points 2200 --minimum-lio-points 1500 --ndt-points 3500
```
这是一组实验参数，不是对所有地图有效的推荐阈值。脚本生成 YAML，不自动部署。

对照保持地图、外参、传感器、线程数、二进制和诊断开关一致：
```bash
python3 scripts/evaluate_compute_pruning.py \
  --baseline /path/to/baseline_run --candidate /path/to/candidate_run \
  --output /path/to/comparison.json \
  --expect-point-budget 2200 --expect-lidar-count 3 --enforce-relative-gates
```
输入包含 `stderr.log`、`frames.csv`、`ndt.tum`。该比较在同一地图系关联时间，不通过拟合对齐隐藏偏差；基线轨迹不是绝对真值。观察位置/航向、相对位姿、有效匹配、来源覆盖和实际预算生效率，再在目标机测延迟与资源。
修改配置无法还原旧二进制算法；严格基线须保存对应构建。现场采集见 @ref resource_diagnostics "资源与端到端延迟"。
