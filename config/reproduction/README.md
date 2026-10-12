# Lightning-LM 复现配置

本目录保存单雷达、多雷达实验参数；SANY 正式建图和定位入口统一位于 [config/sany](../sany/README.md)，只有 mapping.yaml 和 localization.yaml，两组 VIN 外参集中存放在 calibrations.yaml。

- `single_lidar/m3dgr/`：M3DGR 基准和消融配置。
- `single_lidar/sany/`：历史 114/127 单雷达配置。
- `multi_lidar/sany/sany_4lidar_frontend_final.yaml`：114/127/187/195 旧数据的历史回归。
- `multi_lidar/sany/matrix/`：C0—C6 来源缺失和噪声模型消融实验；run_frontend_offline_batch.sh 默认读取这里。

这些历史实验 YAML 不作为新老车正式部署参数。其他对照算法资产保存在各自仓库，不在本仓库复制。
