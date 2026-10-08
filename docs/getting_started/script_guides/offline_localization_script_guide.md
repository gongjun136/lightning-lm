@page offline_localization_script_guide 离线定位脚本导读

@anchor loc_offline_guide
# run_loc_offline.sh：在已有地图上处理离线 bag

## 先看核心调用

@snippet{lineno} run_loc_offline.sh loc-offline-launch

`binary` 是本仓库安装空间中的 `run_loc_offline`。bag、YAML 与地图作为算法输入；轨迹、定位 CSV、帧统计与可选结果 bag 是算法输出。C++ 自己读取 bag，Shell 继续管理后台算法的 PID、资源采样和超时。

| 核心参数 | 作用 |
|---|---|
| `--input_bag / --config / --map_path` | 离线数据、定位 YAML、包含 `index.txt` 的分块地图 |
| `--output_tum / --output_lidar_loc_tum / --output_csv` | 全局定位轨迹、地图匹配轨迹、定位状态统计 |
| `--start_sensor_time / --max_lidar_frames` | 指定绝对传感器起点与融合帧上限 |
| `--use_config_initial_pose` | 是否使用 YAML 初始位姿 |
| `--publish_topics / --output_bag` | 可选 ROS 发布、由程序直接写入业务输出 bag |

```text
bag + YAML + 已有地图
  → C++ run_loc_offline：LIO → 地图定位 → 轨迹与统计
  → Shell：等待 → 计时提取、定位分析 → 运行结果
```

## 按执行顺序阅读

| 阶段 | 目的 | 阅读入口 |
|---|---|---|
| 1. 解析和检查输入 | 确认 bag、YAML、地图、可选参考轨迹 | @ref loc_offline_prepare "输入与运行目录" |
| 2. 准备环境和输入合同 | 核对安装树，取得时长与末帧信息 | @ref loc_offline_prepare "环境与输出" |
| 3. 注册收尾、启动算法 | 传递实际路径和停止条件 | 核心调用 |
| 4. 采样、超时和等待 | 记录资源，取回算法退出码 | @ref loc_offline_wait "等待与退出处理" |
| 5. 汇总和分析结果 | 提取计时、比较可选参考轨迹 | @ref loc_offline_analysis "定位分析" |
| 6. 决定脚本退出码 | 检查完成程度与必需产物 | @ref loc_offline_finish "结果判定" |

@anchor loc_offline_prepare
### 1–2. 输入、环境与输出

@snippet{lineno} run_loc_offline.sh loc-offline-inputs

参数检查后规范路径，拒绝覆盖目录内已有的 runner 产物，加载 ROS/overlay 并核对包 prefix。`inspect_rosbag2_sqlite.py` 生成输入合同；脚本据数据时长和播放倍率计算 watchdog。

算法产物在启动前统一命名：

@snippet{lineno} run_loc_offline.sh loc-offline-outputs

`cleanup_group()` 与 `cleanup()` 的定义出现在这附近，只在调用或 EXIT trap 触发时执行。第一次读先找到 `trap cleanup EXIT`，再顺着核心调用往下。

@anchor loc_offline_wait
### 4. 等待与退出处理

算法启动后先按 PID 启动资源采样，再循环检查进程和超时：

@snippet{lineno} run_loc_offline.sh loc-offline-watchdog

`algorithm_rc` 保留真实退出码，`watchdog_status` 说明是否超时。停止监控后再检查轨迹的时间顺序、间隔和处理尾部，并提取热点计时。

@anchor loc_offline_analysis
### 5. 定位分析与参考轨迹

@snippet{lineno} run_loc_offline.sh loc-offline-analysis

`--reference-tum` 被交给分析器，不是定位程序的输入。未提供参考时仍分析定位状态和轨迹运动；提供后增加按时间对应的误差分析。定位摘要和误差 CSV 放在 results/。

@anchor loc_offline_finish
### 6. 结果判定与退出

@snippet{lineno} run_loc_offline.sh loc-offline-finish

runner 根据算法状态、超时、完成程度、轨迹/计时和必需产物决定是否返回 4。`physical_diagnostic_pass` 会写出并显示，但这个最终条件没有把它列为硬性失败项。

进入 @ref run_loc_offline.cc "main()" 后，重点读 `LaserMapping`、`LidarLoc` 的初始化、bag 回调连接与 `rosbag.Go()`，然后跟踪最后一帧 flush 和结果保存。它有自己的离线编排；理解算法约束时可对照 @ref online_localization_flow "在线定位数据流"。

@htmlonly[block]
<details>
<summary>参数、资源采样与输出</summary>
@endhtmlonly

@snippet{lineno} run_loc_offline.sh loc-offline-usage

默认不按传感器时间节拍限速，不发布 ROS topics；程序支持按 flags 直接保存结果 bag。`--` 后的 flags 原样追加到 C++ 命令。

@snippet{lineno} run_loc_offline.sh loc-offline-monitor

算法原始日志在 logs/；轨迹、CSV、计时与分析摘要在 results/；bag 合同、watchdog、资源和 metadata 在运行根目录。完整源码：@ref run_loc_offline.sh "run_loc_offline.sh"；精确参数见 @ref script_contracts "脚本契约"。

@htmlonly[block]
</details>
@endhtmlonly
