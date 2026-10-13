# 新车测试场地

VIN 已确认是 `SANLU77E8TH000006`，LIDAR 为 `MASTER`，对应 `user_provided_20261010_v1` 标定。现在使用同目录共用 `mapping.yaml` / `localization.yaml`，详见 [统一配置说明](README.md)。无需新车专用算法 YAML。

现场统一使用 `${SANY_WS}/run_loc.sh`，由其调用 [`scripts/run_sany_lidar_loc.sh`](../../scripts/run_sany_lidar_loc.sh)。地图路径通过 `SANY_MAP_PATH` 指定。

地图版本为 `new_vehicle_4lidar_20260922_201743_calib_user_20261010_v1`，开发端位于 `F:/datasets/SANY/sites/new_vehicle_test_ground/derived/lightning-lm/maps/`。部署完整地图、分块、BTC 子图和 SOLiD 描述库，校验 `artifacts.sha256`。

历史验证基于提交 448744d 的专用配置：3398 个四雷达融合帧建图，五包静止冷启动定位 5/5 通过。历史配置使用 20 ms 帧配对和纯 solid 重定位。当前共用方案采用用户新确认的旧车 solid_kiss/plane_icp 及自适应单雷达观测策略，历史结果不能代替本次方案的实车实时与动态精度验收。

启动后核对运行目录 config.yaml 中的 VIN、LIDAR、四路话题、calibration_id；地图路径仍独立选择。新车 IMU/车体/CAN 参数尚未独立提供，现有值沿用来源配置；新车 CAN 速度换算确认前，在统一入口设置 `SANY_ENABLE_CAN_OBSERVATION=0`。
