# Lightning-LM 复现配置

本目录只保存 Lightning-LM 单雷达和多雷达实验的配置。四个对照算法的复现资产分别位于其代码仓库的 `reproduction/m3dgr/`：

- `ws_fastlio/src/FAST_LIO/reproduction/m3dgr/`
- `ws_fastlivo2/src/FAST-LIVO2/reproduction/m3dgr/`
- `ws_voxel_slam/src/Voxel-SLAM/reproduction/m3dgr/`
- `ws_ground_fusion2/src/Ground-Fusion2/reproduction/m3dgr/`

Lightning-LM 不再复制 FAST-LIO、FAST-LIVO2、Voxel-SLAM 或 Ground-Fusion2 的 YAML/launch。

## 目录

```text
reproduction/
├── single_lidar/
│   ├── m3dgr/
│   │   ├── lightning_m3dgr_mid360_benchmark.yaml
│   │   ├── lightning_m3dgr_variant_acc01.yaml
│   │   ├── lightning_m3dgr_variant_classic.yaml
│   │   └── lightning_m3dgr_variant_prop_off.yaml
│   └── sany/
│       ├── sany_livox_114_20260701_velprop_no_lba.yaml
│       └── sany_livox_127.yaml
└── multi_lidar/
    └── sany_4livox/
        ├── sany_4lidar_frontend_final.yaml
        └── matrix/
            ├── C0_single114_noise.yaml
            ├── C1_four_noise.yaml
            ├── C2_four_no_noise.yaml
            ├── C3_drop114_lidar.yaml
            ├── C4_drop127.yaml
            ├── C5_drop187.yaml
            └── C6_drop195.yaml
```

## 使用约定

- M3DGR 的 `benchmark` 配置是已完成单雷达精度实验的正式配置；其余三个 `variant` 文件是 Lightning-LM 消融/调参配置。
- SANY 单雷达配置保留原文件名，便于追溯早期 114 和 127 数据处理。
- `sany_4lidar_frontend_final.yaml` 用于四雷达正常在线合同和故障注入。
- `matrix/` 中的 C0—C6 用于单 114 基线、四雷达、噪声模型消融和单路 LiDAR 缺失实验。
- `scripts/run_frontend_offline_batch.sh` 默认读取 `matrix/`，不需要再传 `--config-dir`。
