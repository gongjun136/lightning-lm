@page mapping_script_guide 建图与地图工具导读

这组入口覆盖在线建图、离线导出，以及已有地图的保存与对齐。保存服务和对齐工具各有独立作用；建图算法从 SLAM 程序进入。

| 当前任务 | 导读 |
|---|---|
| 用实时 ROS 输入建图 | @ref slam_online_guide "run_slam_online.sh" |
| 从 bag 建图并导出 | @ref slam_offline_guide "run_slam_offline.sh" |
| 请求运行中的节点保存地图 | @ref save_map_guide "save_default_map.sh" |
| 求两份地图之间的刚体变换 | @ref align_maps_guide "align_maps_offline.sh" |

@anchor slam_online_guide
# run_slam_online.sh：运行在线 SLAM

## 先看核心调用

@snippet{lineno} run_slam_online.sh slam-online-launch

进入 run dir 后启动 `run_slam_online`，配置路径已转为绝对路径，其余 flags 由 `"$@"` 透传。调用在前台等待，日志重定向到 logs/；相对输出路径以 run dir 为基准。

## 按执行顺序阅读

| 阶段 | 目的 | 阅读入口 |
|---|---|---|
| 1. 加载环境 | 找到 ROS、消息包和 SLAM 程序 | @ref slam_online_environment "环境" |
| 2. 确定参数 | 选配置、输出根和 run name | @ref slam_online_settings "运行设置" |
| 3. 建立运行目录 | 归档日志及基本输入信息 | @ref slam_online_directory "目录" |
| 4. 进入在线 SLAM | 等待 ROS 输入，运行前端与配置选定的后端 | 核心调用 |

@anchor slam_online_environment
### 1. 加载环境

@snippet{lineno} run_slam_online.sh slam-online-environment

`-h / --help` 分支位于环境加载之后，会调用安装二进制帮助并退出。

@anchor slam_online_settings
### 2. 确定配置与名称

@snippet{lineno} run_slam_online.sh slam-online-settings

第一个非 `-` 开头参数消费为 run name，后续 flags 才传给 C++。YAML 必须存在，再转为绝对路径。

@anchor slam_online_directory
### 3. 建立运行目录

@snippet{lineno} run_slam_online.sh slam-online-directory

同名目录可以复用，基本日志与 metadata 会被截断。这里没有辅助进程管理；算法调用结束后，脚本沿用其退出状态。

## 进入 C++ 后继续读哪里

从 @ref run_slam_online.cc "main()" 沿 `SlamSystem::Init → StartSLAM → 可选内嵌 bag 播放 → Spin` 继续。前端和后端之间的交接见 @ref guide_mapping "建图与后端"。传感器或外部回放需另行启动；显式 `--bag` 的播放由 C++ 负责。

@htmlonly[block]
<details>
<summary>参数、地图保存与源码</summary>
@endhtmlonly

默认 YAML 是 `$repo_dir/config/default.yaml`，输出根是 `$repo_dir/runs`。配置、仓库、ROS/overlay、输出根和 run name 可分别通过 `LIGHTNING_LM_*` 变量覆盖。轨迹、地图等由 YAML、程序 flags 和保存服务决定；脚本自身保存 metadata 与 stdout/stderr。

服务请求见 @ref save_map_guide "保存地图脚本"。完整源码：@ref run_slam_online.sh "run_slam_online.sh"；参数边界见 @ref script_contracts "脚本契约"。

@htmlonly[block]
</details>
@endhtmlonly

@anchor slam_offline_guide
# run_slam_offline.sh：离线 SLAM 与地图导出

## 先看核心调用

@snippet{lineno} run_slam_offline.sh slam-offline-launch

`binary` 指向已安装的 `run_slam_offline`，它读取 bag，运行前端及配置选定的后端，最后导出结果。`setsid/taskset/&/algorithm_pid` 的组合让主脚本同时管理资源采样和超时。

| 核心参数 | 交给 C++ 的职责 |
|---|---|
| `--input_bag / --config` | 数据输入和算法配置 |
| `--output_tum / --output_frame_stats_csv` | SLAM 轨迹、融合帧统计 |
| `--output_map_dir / --output_global_map` | 分块地图目录、全局 PCD |
| `--backend_evaluation_only` | 仅评估后端时跳过地图/重定位数据导出 |

## 按执行顺序阅读

| 阶段 | 目的 | 阅读入口 |
|---|---|---|
| 1. 解析输入、建立目录 | 得到 bag/YAML/输出，拒绝覆盖旧产物 | @ref slam_offline_prepare "运行准备" |
| 2. 加载环境和输入合同 | 确认本仓库二进制，取得时长与末帧 | @ref slam_offline_prepare "环境" |
| 3. 注册收尾、启动算法 | 允许 EXIT 时清理进程组 | 核心调用 |
| 4. 采样并等待 | 保存资源记录，处理超时 | @ref slam_offline_wait "等待与超时" |
| 5. 统计并判断结果 | 检查处理完成程度和地图产物 | @ref slam_offline_finish "结果与退出" |

@anchor slam_offline_prepare
### 1–2. 准备输入、输出与安装环境

@snippet{lineno} run_slam_offline.sh slam-offline-inputs

输出目录中若已有 results/logs/data 等脚本产物就拒绝覆盖。未显式给 `--output-dir` 时，兼容环境变量中的 output root/run name。

@snippet{lineno} run_slam_offline.sh slam-offline-outputs

`resolve_output()` 把相对输出参数放到 run dir 下；默认地图目录为 `data/new_map`。

@snippet{lineno} run_slam_offline.sh slam-offline-environment

随后 `inspect_rosbag2_sqlite.py` 写入输入合同，得到预期结束时间与时长，再计算 watchdog。`cleanup()` 的定义负责进程组清理，`trap cleanup EXIT` 注册它；主线在后面的核心调用启动算法。

@anchor slam_offline_wait
### 4. 等待与超时

@snippet{lineno} run_slam_offline.sh slam-offline-watchdog

算法启动后先启动资源采样，再进入这个循环。结束或超时后，用 `wait` 得到 `algorithm_rc`，停止监控并取消 trap，然后进行统计。bag 播放结束后的 flush、后端最终优化与导出仍是 C++ 的运行阶段。

@anchor slam_offline_finish
### 5. 结果与退出

@snippet{lineno} run_slam_offline.sh slam-offline-finish

普通建图还要求分块索引、全局地图和一定数量的地图点；`backend_evaluation_only=true` 时不要求这部分地图产物。所有任务都要求轨迹、帧统计、计时和资源记录满足 runner 条件。

从 @ref run_slam_offline.cc "main()" 继续：初始化 LIO 与后端 → `rosbag.Go()` → 最后一帧 flush → 等待后端 → 地图导出。配置可选择新 `BackendPipeline`、legacy 或关闭后端，不能仅凭脚本名判断用了哪一种。

@htmlonly[block]
<details>
<summary>参数与资源采样</summary>
@endhtmlonly

@snippet{lineno} run_slam_offline.sh slam-offline-usage

资源采样依赖算法 PID，发生在核心调用之后：

@snippet{lineno} run_slam_offline.sh slam-offline-monitor

轨迹和计时在 results/，地图在 data/new_map/，原始日志在 logs/，资源与 metadata 在运行根目录。完整源码：@ref run_slam_offline.sh "run_slam_offline.sh"。

@htmlonly[block]
</details>
@endhtmlonly

@anchor save_map_guide
# save_default_map.sh：请求保存当前地图

## 核心调用与执行顺序

@snippet{lineno} save_default_map.sh save-map-request

整个脚本只有这个服务请求：当前终端的 ROS 环境 → 找到 `/lightning/save_map` → 发送 `map_id: new_map` → 等待服务响应。实际地图写入由已运行的节点处理，脚本不启动 SLAM，也不加载环境。

请求中的 `new_map` 是传给服务的地图标识；实际落盘位置取决于服务实现及节点运行目录。需要改变地图标识时直接阅读并调整服务调用参数。后续看 @ref run_slam_online.cc "在线 SLAM 入口" 与地图保存实现，完整源码见 @ref save_default_map.sh "save_default_map.sh"。

@anchor align_maps_guide
# align_maps_offline.sh：调用已有地图对齐程序

## 核心调用与执行顺序

@snippet{lineno} align_maps_offline.sh align-maps-launch

执行顺序是：确定仓库 → 优先找仓库内 install 中的程序 → 找不到时尝试 bin/ → 两处都不可执行则退出 2 → `exec` 启动程序。此脚本没有 `source`、`cd`、录包或在线节点管理。

`"$@"` 把全部 flags 原样交给对齐程序，Shell 不解析地图或输出参数；相对路径沿用调用终端的 cwd。进入 @ref align_maps_offline.cc "align_maps_offline.cc 的 main()"，再读候选搜索、匹配和刚体变换输出。

@htmlonly[block]
<details>
<summary>路径与参数边界</summary>
@endhtmlonly

这个入口按仓库本地 install/bin 布局找二进制，没有通过 `ros2 pkg prefix` 解析当前 overlay。具体地图、搜索配置和结果路径由 C++ flags 指定。完整源码：@ref align_maps_offline.sh "align_maps_offline.sh"；参数与副作用见 @ref script_contracts "脚本契约"。

@htmlonly[block]
</details>
@endhtmlonly
