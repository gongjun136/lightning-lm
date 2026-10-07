@page regression_script_guide 在线回归与参考实验导读

这两个入口都用已有 bag 进行实验。在线回归调用 Lightning 的 ROS 2 程序；114 参考实验调用另一个工作区的 ROS 1 Voxel-SLAM。

@anchor online_regression_guide
# run_online_bag_regression.sh：用 bag 复现在线行为

## 先看核心：默认由 C++ 内嵌播放

@snippet{lineno} run_online_bag_regression.sh online-regression-embedded

默认 `LIGHTNING_LM_EMBEDDED_PLAYBACK=1`，脚本前台调用安装空间中的在线二进制，传入 `--bag`、1 倍播放和结束等待时间。`slam` 分支保存地图，`loc` 分支保存三类轨迹。这个分支完成后直接退出，不走后面的外部录包/播放器流程。

```text
默认：Shell → run_slam_online 或 run_loc_online（内部 bag player）→ 结果
关闭内嵌播放：Shell → 在线节点 + 可选输出录制 → ros2 bag play → 收尾
```

## 按执行顺序阅读

| 阶段 | 目的 | 阅读入口 |
|---|---|---|
| 1. 解析模式和输入 | `slam/loc`、名称、YAML、bag；loc 还需要地图 | 完整源码开头的位置参数检查 |
| 2. 准备目录与环境 | 建日志目录，设置独立实验 domain，注册 cleanup | @ref online_regression_settings "运行设置" |
| 3. 选择播放方式 | 内嵌分支前台运行并结束 | 上面的核心调用 |
| 4. 外部模式启动节点和回放 | 节点先启动，loc 录制位姿，再播放输入 | @ref online_regression_external "外部模式" |
| 5. 外部模式收尾 | 等队列处理，保存地图或提取轨迹 | @ref online_regression_finish "收尾" |

@anchor online_regression_settings
### 2. 运行设置

@snippet{lineno} run_online_bag_regression.sh online-regression-settings

建立目录后加载 ROS/overlay，`ROS_DOMAIN_ID` 来自 `LIGHTNING_LM_ROS_DOMAIN_ID`，缺省 83。`stop_process()` 定义停止方式，`trap cleanup EXIT INT TERM` 注册异常退出处理；函数定义处不执行清理。

@anchor online_regression_external
### 4. 关闭内嵌播放时，先节点后输入

@snippet{lineno} run_online_bag_regression.sh online-regression-node

节点后台启动后做存活检查。loc 分支随后录制 `/slamPoseRaw_topic` 与 `/PosRes`，再调用外部播放器：

@snippet{lineno} run_online_bag_regression.sh online-regression-playback

这里的播放在前台等待，默认 discovery delay 为 10 秒。`/usr/bin/time` 分别记录默认分支的节点资源或外部模式的播放器资源，两者含义不同。

@anchor online_regression_finish
### 5. 外部模式收尾

@snippet{lineno} run_online_bag_regression.sh online-regression-finish

播放结束后默认再等 15 秒。slam 调保存服务；loc 停止录制并从结果 bag 提取 TUM，之后停止节点。这个外部流程不等同于默认内嵌分支的三轨迹直接输出。

进入 C++ 读 @ref run_slam_online.cc "在线 SLAM main()" 或 @ref run_loc_online.cc "在线定位 main()"，重点区分 ROS 事件循环与可选 bag player。

<details>
<summary>调用参数与结果目录</summary>

调用方式：`<slam|loc> <run-name> <config.yaml> <bag-dir> [map-dir]`，loc 必须提供地图。输出根缺省 `$repo_dir/runs/online_regression`，脚本在 run dir 内执行。基本目录可复用，日志/metadata 可能被截断，rosbag 自身另有输出目录约束。

`LIGHTNING_LM_EMBEDDED_PLAYBACK=0` 选择外部流程；`LIGHTNING_LM_POST_WAIT_SECONDS` 和 `LIGHTNING_LM_DISCOVERY_DELAY_SECONDS` 控制收尾/发现等待。完整源码：@ref run_online_bag_regression.sh "run_online_bag_regression.sh"；契约见 @ref script_contracts "脚本契约"。

</details>

@anchor voxel_reference_guide
# run_voxel_slam_114_reference.sh：ROS 1 参考轨迹实验

## 先看核心：launch 启动算法，bag 提供输入

@snippet{lineno} run_voxel_slam_114_reference.sh voxel-reference-launch

算法来自指定 Voxel-SLAM 工作区的 launch，保存路径和序列名由参数指定。节点就绪后启动 TF 轨迹录制和资源监控，再播放两个输入话题：

@snippet{lineno} run_voxel_slam_114_reference.sh voxel-reference-playback

`wait` 等待播放器，回放结束并让队列排空后设置 `/finish=true`，要求算法完成全局优化。接下来仍需等待最终轨迹写完。

## 按执行顺序阅读

| 阶段 | 目的 | 阅读入口 |
|---|---|---|
| 1. 检查输入与工具 | ROS 1 bag、launch/config、二进制、空输出目录 | 完整源码的参数检查 |
| 2. 检查 bag、隔离 ROS 环境 | 得到预期 LiDAR 数量与尾部时间 | @ref voxel_reference_environment "环境与 master" |
| 3. 启动算法和采集 | 等 `/voxelslam` 就绪，再开 TF 录制与采样 | 核心 launch；折叠参考中的采集 |
| 4. 播放并请求最终优化 | 只回放选定 LiDAR/IMU，使用仿真时钟 | 核心 playback |
| 5. 等保存并停止节点 | 避免在优化轨迹写入时提前停止 | @ref voxel_reference_results "优化与结果" |
| 6. 转换、汇总、退出 | 导出 TUM，检查尾部/数量/资源证据 | @ref voxel_reference_finish "结果判定" |

@anchor voxel_reference_environment
### 2. 环境与独立 master

加载 Noetic 后读取 bag 的消息数量、时间范围和 LiDAR header 时间，再加载 Voxel 工作区：

@snippet{lineno} run_voxel_slam_114_reference.sh voxel-reference-environment

@snippet{lineno} run_voxel_slam_114_reference.sh voxel-reference-master

master 使用指定端口，ROS 日志与 home 归入输出目录，`/use_sim_time=true`。`cleanup()`/trap 管理本次启动的播放、录制、launch、采样和 master 进程。

@anchor voxel_reference_results
### 5. 等优化轨迹稳定后停止节点

@snippet{lineno} run_voxel_slam_114_reference.sh voxel-reference-optimization

`alidarState.txt` 连续五次大小不变才视为稳定，超出 shutdown wait 则失败。随后停止 `/voxelslam` 和 TF 录制，再把优化状态转换为 TUM：

@snippet{lineno} run_voxel_slam_114_reference.sh voxel-reference-tum

结果中分别保留实时 TF 前端轨迹和最终优化轨迹。算法源码在外部 Voxel-SLAM 工作区，入口要从所选 launch 的节点定义继续追踪。

@anchor voxel_reference_finish
### 6. 结果判定与退出

@snippet{lineno} run_voxel_slam_114_reference.sh voxel-reference-finish

汇总检查有效轨迹、时间单调性、末帧接近程度、输出比例和资源摘要；成功后关闭 trap 并显式 cleanup。

<details>
<summary>参数、TF 录制与资源采样</summary>

@snippet{lineno} run_voxel_slam_114_reference.sh voxel-reference-usage

TF 录制的 frame 与平移修正是此参考实验的固定条件：

@snippet{lineno} run_voxel_slam_114_reference.sh voxel-reference-record

默认话题为 Livox 114 LiDAR/IMU，private master 端口 11331，CPU `0-7`，播放 1 倍。结果在 results/，算法原始文件在 data/，指纹与运行信息写 metadata。完整源码：@ref run_voxel_slam_114_reference.sh "run_voxel_slam_114_reference.sh"。

</details>
