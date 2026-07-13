# 资源报告数据与图表说明

## 数据源

- 正式矩阵：`/mnt/f/datasets/SANY/mid360/data_20260701/lightning_lm_4lidar_frontend_20260710/formal_experiments_20260711_v2`
- 四雷达完整运行资源采样：`/mnt/f/datasets/SANY/mid360/data_20260701/lightning_lm_4lidar_frontend_20260710/manual_4lidar_repro_01/resource_samples.csv`
- 四雷达完整运行资源摘要：`/mnt/f/datasets/SANY/mid360/data_20260701/lightning_lm_4lidar_frontend_20260710/manual_4lidar_repro_01/resource_summary.json`

## 指标口径

- 墙钟时间：GNU time 的 `Elapsed (wall clock) time`，包含算法、ROS 2 bag 读取、结果写盘和脚本开销。
- 平均 CPU 逻辑核：GNU time 的 CPU 百分比除以 100；占 8 核配额比例再除以 8。
- 峰值 RSS：GNU time 的 `Maximum resident set size`，由 KiB 换算为 MiB。
- 时序 CPU/RSS：统一离线脚本对子进程树约每 0.2 秒采样一次。
- 阶段耗时：算法 Timer 的最近样本统计；不同阶段可能嵌套，不能相加为端到端耗时。

## 图表地图

| 报告段落 | 分析问题 | 图表 | 字段 | 支持的结论 |
|---|---|---|---|---|
| 总体资源 | 七种配置的资源代价有多大 | 横向分组比较 | 墙钟时间、峰值 RSS、平均 CPU 核 | 四雷达主要增加内存与总耗时，CPU 平均并未占满 8 核 |
| 运行时序 | C1 的资源峰值在何时出现 | 双面板时序折线 | elapsed_s、cpu_cores、rss_mb | 内存随地图增长累积，CPU 存在阶段性波动 |
| 阶段耗时 | 核心计算热点在哪里 | 分组柱状图 | average_ms、config、stage | 激光观测匹配是已埋点阶段中最重的单次调用 |

## 限制

- 正式矩阵的三次重复为确定性技术重复，不等同于三段独立数据或三台机器。
- 正式矩阵采自历史冻结安装树；当前统一分支增加了耗时导出，但未重新执行全部 21 个完整单元。
- WSL 2 + Windows 挂载盘的 I/O 与调度不代表目标车载原生 Linux 平台。
- GNU time 的 CPU 百分比是整个任务平均值，无法定位线程级并行效率。
