@page md_config_sany_README SANY 共用建图与定位配置

正式算法配置只有 `mapping.yaml` 和 `localization.yaml`。`calibrations.yaml` 保存两组 VIN 的相对雷达标定和轮速换算系数，不复制共用算法参数。旧 `sany_3livox/sany_4livox` 正式配置已合并；早期 114/127/187/195 实验保存在 `config/reproduction/multi_lidar/sany/`。

## 车辆配置

车辆默认读取 `/opt/conf/vehicle_config.yaml`，离线/WSL 可通过 `SANY_VEHICLE_CONFIG` 指定另一个文件。

```yaml
VIN: "SANLV77EXSH000858"
DEVICE_CODE: "SW956E9CF7868"
LIDAR: "MASTER AND SLAVE"
```

- `SANLV77EXSH000858`：旧车标定 `legacy_vehicle`。
- `SANLU77E8TH000006`：新车标定 `user_provided_20261010_v1`。
- `MASTER AND SLAVE`：master 定位/建图使用前、左、右三台，保留旧 slave 独立驱动入口。
- `MASTER`：master 使用前、左、右、后四台；后雷达接收 IP 为 `192.168.4.100`。

VIN 与接线方式分别解析；文件缺失、未知 VIN、非法 LIDAR 或无效矩阵时拒绝启动。相对外参约定为 `T_primary_from_lidar_row_major`，平移单位 m，前雷达为 ID 0，主 IMU 同为前雷达。

2026-10-13 用户确认新车真实 VIN 为 `SANLU77E8TH000006`，此前登记的 `SW956ELCG0858` 对应同一台车，已更正 VIN 匹配值；相对雷达标定保持不变。车辆识别只读取 `VIN`，`DEVICE_CODE` 不替代 VIN，也不增加第三组标定。

`scripts/resolve_sany_config.py` 中的 `resolve()` 同时将所选 VIN 的 `wheel_speed_scale_mps_per_rpm` 写入运行快照的 `system`，与 LIDAR 接线方式无关。旧车 `SANLV77EXSH000858` 使用 `0.00120639253574024`，新车 `SANLU77E8TH000006` 使用 2026-10-18 确认的 `0.0011687`，单位均为 `(m/s)/rpm`。系数缺失、非数值、非有限或不为正时拒绝生成配置；共用模板不再保存统一系数。直接调用 C++ 前也必须先解析配置。已有运行快照不自动修改，更新后需重新启动生成快照。

## 启动

```bash
ros2 launch livox_ros_driver2_core sany_lidar.launch.py
SANY_MAP_PATH=/absolute/released/map bash scripts/run_sany_lidar_loc.sh
```

车辆布局由 LIDAR 决定，不再用 `SANY_LIDAR_LAYOUT` 选择 YAML。原 `master_lidar.launch.py`、`lidar_rear.launch.py` 保持可用，显式四台使用 `master_4lidar.launch.py`。

离线建图示例：

```bash
export SANY_VEHICLE_CONFIG=/absolute/vehicle_config.yaml
bash scripts/run_slam_offline.sh --bag /absolute/bag \
  --config "$PWD/config/sany/mapping.yaml" --output-dir /absolute/new_run
```

离线定位使用 `scripts/run_loc_offline.sh --config "$PWD/config/sany/localization.yaml"`，并提供 bag、map 和 output-dir。四个通用 online/offline Shell 入口会识别 SANY 模板并解析；其他数据集配置继续原样运行。直接调用 C++ 工具时先生成完整运行 YAML：

```bash
python3 scripts/resolve_sany_config.py \
  --base-config config/sany/localization.yaml \
  --vehicle-config /absolute/vehicle_config.yaml --output /absolute/new_run/config.yaml
```

运行快照记录 VIN、布局、标定版本以及三个源文件的 SHA-256，已有快照禁止覆盖。纯 INS 的 SANY Shell 入口也解析其 `cloud.lidar_config`，保持对共用标定的引用。

## 观测与点云发布

定位最低允许前主雷达一台提供观测。健康时使用全部可用雷达，先按 2200/1804/1500 点预算降载，再逐级减少辅助雷达，始终保留前雷达；恢复负载后逐级恢复。质量下降时恢复全部可用来源和未限额点云。观测能否被接受仍由实际匹配质量决定。

发布点云使用完整去畸变点云缓存，与 LIO 选择的来源/点数独立。只有当前配置的全部三台或四台雷达在同一融合帧中齐全，才发布拼接点云；缺任何一路时停止该帧点云发布，内部定位可继续使用满足条件的部分点云。

重定位和其他算力参数采用旧车方案：solid_kiss + plane_icp、8/64 候选/检索池、55%/75% 预检查/最终重合率，LIO 6 线程、NDT 4 线程、NDT 跟踪限额 3500 点。点时间单位按原始 bag 核实为纳秒，明确使用 `livox_point_time_scale: 1e-9`。

建图和定位共用 `match_tolerance: 0.02`（20 ms），兼容新车历史四雷达约 10 ms 的帧起点跨度；仅放宽同周期帧配对条件，不改写点的物理时间。定位重排窗口保持 0.1 s。

旧地图固定变换的数值保留，但在线启动从所选地图的 `map_manifest.yaml` 加载实际输出变换。车辆/VIN 配置不选择地图坐标。缺少清单、状态未确定、哈希或矩阵矛盾均拒绝启动；旧 `SANY_ENABLE_FIXED_MAP_TRANSFORM` 已弃用，不覆盖地图清单。历史地图补齐和新图 `map_export.alignment_status` 声明见 @ref map_contracts "地图包、坐标归一化与栅格"。
