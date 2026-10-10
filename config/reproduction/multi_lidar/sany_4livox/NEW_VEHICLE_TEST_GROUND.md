# 新车测试场地配置与部署

本目录新增一对专用文件，仓库原有三/四雷达默认 YAML 保持原样：

- `sany_4lidar_new_vehicle_test_ground_mapping.yaml`：两圈建图的已验证配置。
- `sany_4lidar_new_vehicle_test_ground_localization_solid.yaml`：五包独立冷启动定位的已验证配置。

两份完整 YAML 与数据目录中的对应验证配置在解析后完全相同，只增加来源注释。四雷达矩阵为用户提供的新车标定 `user_provided_20261010_v1`，无额外地面外参修正。车辆暂记 `new_vehicle_20260922`；车辆编号、VIN、雷达序列号待确认。前/主 184、左 108、右 133、后 143 全部配置，主 IMU 为 184。

## 手动选择 YAML

现有 `scripts/run_sany_lidar_loc.sh` 使用：

```bash
config_path="${LIGHTNING_LM_CONFIG:-${default_config_path}}"
```

没有设置 `LIGHTNING_LM_CONFIG` 时，`SANY_LIDAR_LAYOUT=3/4` 仍选择原有默认文件。设置绝对路径后，该文件优先；实际传感器话题从它的 `multi_lidar.topics` 读取。布局变量不会自动修改或裁剪手动指定文件。

用户使用的《lightning-lm定位与地图数据部署》模板当前固定三雷达和旧宝马展地图；新车测试场地应把其 `run_loc.sh` 的 LiDAR 分支改为以下预设，仍调用相同仓库入口：

```bash
export LIGHTNING_LM_CONFIG="$LIGHTNING_LM_REPO_DIR/config/reproduction/multi_lidar/sany_4livox/sany_4lidar_new_vehicle_test_ground_localization_solid.yaml"
export SANY_MAP_PATH="$SANY_WS/maps/new_vehicle_4lidar_20260922_201743_calib_user_20261010_v1"
export SANY_LIDAR_LAYOUT=4
export LIGHTNING_LM_RUN_MODE=production
export SANY_RECORD_BAG=0
export SANY_ENABLE_CAN_OBSERVATION=0
exec bash "$LIGHTNING_LM_REPO_DIR/scripts/run_sany_lidar_loc.sh"
```

可直接使用仓库中的 [新车 run_loc.sh 模板](../../../../scripts/field/new_vehicle_test_ground_run_loc.sh)，替换前保留现场原入口供回退。先按同一部署流程更新/构建代码、安装依赖并传输完整新地图包，校验地图 `artifacts.sha256`。已有 SOLiD 数据库与地图成套，不能使用文档中的旧宝马展地图路径或固定坐标变换。

地图来源：`F:/datasets/SANY/sites/new_vehicle_test_ground/derived/lightning-lm/maps/new_vehicle_4lidar_20260922_201743_calib_user_20261010_v1/`。运行 YAML 的 `calibration_reference.registry_revision` 是开发端来源记录，程序不需要在域控读取这个开发端目录；域控直接读取完整 YAML。程序通过 `SANY_MAP_PATH` 选择实际地图。

## 四雷达输入

部署文档中的 `master_lidar.launch.py` 仅启动前、左、右三路，还需按新车实际接线启动后雷达 143，并把它的点云送到同一 ROS domain。运行定位前，四个话题都应有真实 publisher 和约 10 Hz 消息：

```text
/livox/lidar_192_168_3_184
/livox/lidar_192_168_1_108
/livox/lidar_192_168_2_133
/livox/lidar_192_168_4_143
```

同时确认 `/livox/imu_192_168_3_184` 有实时消息。新车网卡、主机 IP 和跨机通信应以该车配置为准，不能照搬旧车的网卡编号。启动器会按这份 YAML 预检四路点云；仅启动原三路驱动会停在输入预检。

## 已验证范围

新图完整处理 3398 个四雷达融合帧，五个独立静止录包冷启动定位 5/5 通过。仓库 YAML 的有效参数与这些实验相同：点时间纳秒缩放 `1e-9`、融合容差 `0.02 s`、新车外参、关闭旧场地固定地图变换；建图要求四雷达齐全，5 cm 栅格、Z≥0.3 m。定位阶段仍保留原有三路最低参与数及自适应负载策略，并配置全部四路输入。

这些是开发机离线功能和静止重复性验证，不能当作目标域控实时或动态精度验收。新车 LiDAR/IMU、车体几何和 CAN 系数尚未单独提供，完整 YAML 中仍保留来源参数；模板先关闭 CAN 观测。确认新车的 CAN 参数后再按现场流程启用。

启动后核对运行目录的 `config.yaml` 和 `run_metadata.txt`：配置路径应是上述专用 YAML，`raw_lidar_topics` 应包含四路，`map_path` 应是本次新图；配置中的标定版本应为 `user_provided_20261010_v1`。停止和回退沿用既有两步部署流程。
