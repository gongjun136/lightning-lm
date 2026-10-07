@page build_navigation_script_guide 构建与导航数据工具导读

构建工具产出工作区安装环境；导航数据工具负责录制或回放，定位节点由 @ref sany_ins_launcher_guide "SANY INS 入口" 单独启动。

| 任务 | 导读 |
|---|---|
| 统一构建 CGI、雷达与 Lightning | @ref build_workspace_guide "build_workspace.sh" |
| 理解旧 INS 构建入口 | @ref build_ins_guide "build_ins_only.sh" |
| 安装开发依赖 | @ref install_dep_guide "install_dep.sh" |
| 录制 INS 输入与业务输出 | @ref ins_record_guide "record_ins_only.sh" |
| 用录制输入重跑 INS | @ref ins_replay_guide "replay_ins_only.sh" |

@anchor build_workspace_guide
# build_workspace.sh：统一构建工作区

## 先看核心调用

@snippet{lineno} build_workspace.sh build-workspace-launch

`cd "$ws"` 后执行 `colcon build`。`--executor sequential` 控制包之间顺序构建，`MAKEFLAGS` 与 `CMAKE_BUILD_PARALLEL_LEVEL` 控制单个构建的并行度。Release 是构建类型，`BUILD_TESTING=ON` 启用测试目标的构建配置；此命令没有执行测试。

## 按执行顺序阅读

| 阶段 | 目的 | 阅读入口 |
|---|---|---|
| 1. 确定工作区、jobs 与 cache 开关 | 找到正确的 src/ 与构建并发 | @ref build_workspace_settings "工作区与参数" |
| 2. 检查包布局 | 缺包时在构建前结束 | @ref build_workspace_packages "所需包" |
| 3. 加载 ROS、进入工作区 | 让 colcon 解析包及依赖 | 核心调用 |
| 4. 构建与退出 | 将构建结果交回调用者 | 核心调用中的 colcon；严格模式保留失败 |

@anchor build_workspace_settings
### 1. 工作区与参数

@snippet{lineno} build_workspace.sh build-workspace-settings

仓库在工作区 src/ 内时向上推导；兼容布局使用相邻 `lightning_lm_ws`。可用 `LIGHTNING_LM_BUILD_WS` 覆盖。jobs 默认 2，只接受正整数；脚本自身只消费可选的 `--cmake-clean-cache`。

@anchor build_workspace_packages
### 2. 所需包

@snippet{lineno} build_workspace.sh build-workspace-packages

这份清单要求 CGI 消息、CGI 驱动、Livox 驱动/SDK 与 Lightning 共用工作区。`-Ucgi430_interfaces_DIR` 清除旧外部 CGI 消息包的 CMake 缓存指向，让 colcon 使用本工作区依赖。

<details>
<summary>产物与副作用</summary>

构建在所选工作区写 build/install/log，`--cmake-clean-cache` 要求清理 CMake 缓存。安装环境之后由启动脚本加载；构建工具本身不启动导航或定位节点。布局见 @ref workspace_layout "工作区说明"，完整源码见 @ref build_workspace.sh "build_workspace.sh"。

</details>

@anchor build_ins_guide
# build_ins_only.sh：兼容入口转交统一构建

## 核心调用与执行顺序

@snippet{lineno} build_ins_only.sh build-ins-delegate

确定自身目录后，用 `exec bash` 进入同目录的 `build_workspace.sh`，全部参数由 `"$@"` 转交。当前 Shell 被替换，因此没有第二套 INS 构建流程，也没有独立收尾阶段。接着读 @ref build_workspace_guide "统一构建导读" 即可。

工作区和并发仍使用子脚本的环境变量；允许的参数仍为可选 `--cmake-clean-cache`。完整源码：@ref build_ins_only.sh "build_ins_only.sh"。

@anchor install_dep_guide
# install_dep.sh：安装基础开发软件包

## 核心调用与执行顺序

@snippet{lineno} install_dep.sh install-dependencies

整个脚本执行一次 `sudo apt install`：由 apt 解析软件包、按需请求权限与安装确认，再安装 OpenCV、PCL、YAML、glog、gflags 和 ROS Humble PCL 转换包。它会修改系统软件包，不创建运行目录，不构建或启动算法。

这里没有额外的参数解析、环境加载和后台进程。后续构建从 @ref build_workspace_guide "build_workspace.sh" 进入。完整源码：@ref install_dep.sh "install_dep.sh"。

@anchor ins_record_guide
# record_ins_only.sh：录制配置输入与导航输出

## 先看核心调用

@snippet{lineno} record_ins_only.sh ins-record-launch

`topics` 是从现场 YAML 得到的输入清单，后面追加业务位姿、定位/故障状态、INS 诊断与心跳。`-s mcap` 指定存储格式，`-o "$2"` 指定输出 bag；`exec` 替换 Shell，后续由 rosbag 持续录制，停止时由它完成落盘。

## 按执行顺序阅读

@snippet{lineno} record_ins_only.sh ins-record-inputs

1. 要求两个参数：`site.yaml output_bag`。
2. `ins_topics.py` 读取 YAML，返回实际配置需要的输入话题；`mapfile` 按行转成 Bash 数组。
3. 进入上面的录制命令，等待输入与输出消息，直到录制进程停止。

输入列表使用数组展开，每个话题保持独立参数。算法入口仍由 @ref sany_ins_launcher_guide "run_sany_ins_only.sh" 启动；此处继续深入应读 `ins_topics.py` 与 rosbag 命令，而不是 C++ 入口。

<details>
<summary>环境、cwd 与结果</summary>

调用终端先加载 ROS/消息工作区环境，并设置与节点一致的 domain；本脚本没有 `source`。它不改变 cwd，相对 YAML 与 output_bag 路径相对调用终端。结果是包含配置输入及业务输出的 MCAP bag。

完整源码：@ref record_ins_only.sh "record_ins_only.sh"；导航运行边界见 @ref ins_only_operation "纯组合导航操作契约"。

</details>

@anchor ins_replay_guide
# replay_ins_only.sh：只回放输入，重新生成结果

## 先看核心调用

@snippet{lineno} replay_ins_only.sh ins-replay-launch

`--topics` 使用配置的输入白名单，`--clock` 生成本次回放时钟；不回放旧业务结果、TF 或旧 `/clock`。定位节点在另一终端先启动，重新计算本次输出。

## 按执行顺序阅读

| 阶段 | 目的 | 阅读入口 |
|---|---|---|
| 1. 消费 YAML 与 bag 参数 | 剩余参数留给 rosbag play | @ref ins_replay_settings "回放配置" |
| 2. 核对仿真时间 | 保证节点使用本次回放时钟 | 同一段的 Python 检查 |
| 3. 得到输入白名单 | 只注入导航所需的原始输入 | @ref ins_replay_topics "输入清单" |
| 4. 开始回放 | 等待播放结束，以命令状态退出 | 核心调用 |

@anchor ins_replay_settings
### 1–2. 回放配置

@snippet{lineno} replay_ins_only.sh ins-replay-config

YAML 必须开启 `ins_only.use_sim_time: true`。`shift 2` 后，其余 flags 由 `"$@"` 交给 rosbag。

@anchor ins_replay_topics
### 3. 输入清单

@snippet{lineno} replay_ins_only.sh ins-replay-inputs

得到数组后进入核心调用。对照录制工具可以看到：录制包含输入和输出，回放只选择输入；这是重新运行算法的关键边界。

<details>
<summary>环境、节点与参数</summary>

调用方式为 `site_replay.yaml bag [ros2 bag play options]`。先加载 ROS/消息工作区环境，保证 domain 一致，再单独启动 @ref sany_ins_launcher_guide "INS 定位节点"。脚本不改变 cwd、也不自动建立结果目录；导航输出按节点 YAML/flags 保存。

完整源码：@ref replay_ins_only.sh "replay_ins_only.sh"。

</details>
