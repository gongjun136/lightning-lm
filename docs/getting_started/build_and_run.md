@page build_and_run 构建、目标与启动

# 构建边界

CGI-430 模式新增 `cgi430_interfaces` 和 GeographicLib 依赖，以及 `ins_project_llh`、`georeference_map` 安装目标。消息工作区加载、构建脚本和配置见 @ref ins_only_operation "纯组合导航构建与启动"。

工程是 ROS 2 `ament_cmake` 包，包名 `lightning_lm`。顶层 CMake 载入 `cmake/packages.cmake` 的三方依赖，进入 `src/` 生成库和程序，最后调用 `ament_package()`。

**代码依据：** `CMakeLists.txt`、`package.xml`（`project(lightning_lm)`、`add_subdirectory(src)` 与 `build_type=ament_cmake`）

## 生产构建目标

| CMake 目标 | 类型 | 来源/职责 | 是否安装 |
|---|---|---|---|
| `lightning_lm.libs` | shared library | common、I/O、LIO、backend、localization、maps、system、UI、wrapper | `lib/` |
| `miao.core` | shared library | 图优化核心、求解器、鲁棒核 | `lib/` |
| `miao.utils` | shared library | `miao` 采样辅助 | `lib/` |
| `run_frontend_offline`, `run_frontend_online` | executable | 前端离线/在线入口 | `lib/lightning_lm/` |
| `run_slam_offline`, `run_slam_online` | executable | 完整建图入口 | `lib/lightning_lm/` |
| `run_loc_offline`, `run_loc_online` | executable | 定位入口 | `lib/lightning_lm/` |
| `convert_voxel_slam_map`, `build_solid_database`, `align_maps_offline` | executable | 地图/重定位数据准备 | `lib/lightning_lm/` |
| `pos_res_recorder` | executable | ROS 位姿结果记录 | `lib/lightning_lm/` |
| `test_ui`, `test_ui_init_quit` | executable | UI 手工测试 | **未安装** |
| `docs` | custom target | Doxygen + Graphviz 阅读站点 | 仅 `BUILD_DOCS=ON` 时存在 |

测试目标由 `BUILD_TESTING` 控制；`src/test/*.cc` 里多数目标同时注册为 CTest。`miao/examples` 在当前构建中被注释，不属于默认目标。

**代码依据：** `src/CMakeLists.txt`、`src/app/CMakeLists.txt`、`src/core/miao/CMakeLists.txt`（库、程序、测试和安装声明）

# 标准构建

定位、CGI 驱动、Livox 驱动和消息包现已放在同一个工作区。源码布局及导入清单见 @ref workspace_layout "统一工作区布局"。无需加载外部 CGI SDK 安装空间。

在工作区根目录执行：

```bash
bash src/lightning-lm/scripts/build_workspace.sh
source /opt/ros/humble/setup.bash
source install/local_setup.bash
```

脚本自动加载 ROS，按依赖顺序编译全部 15 个包，缺省最多 2 个编译任务。`LIGHTNING_BUILD_JOBS=4` 可用于内存足够的 WSL；服务器可设为 16。`LIGHTNING_LM_BUILD_WS` 可指定工作区绝对路径。旧入口 `build_ins_only.sh` 转调同一脚本。

手动构建也只需要基础 ROS 环境：

```bash
source /opt/ros/humble/setup.bash
MAKEFLAGS="-j2 -l2" CMAKE_BUILD_PARALLEL_LEVEL=2 colcon build --executor sequential
source install/local_setup.bash
```

迁移前缓存若显式设置过外部 `cgi430_interfaces_DIR`，第一次请使用统一脚本；它会取消这项旧缓存，使 CMake 通过 colcon 找到本工作区接口。无需 source `cgi430_build_check_20261006` 或旧 `sdk/cgi430_sdk/install`。

`--packages-up-to lightning_lm` 只构建定位及其接口依赖，不包含独立的传感器驱动。完整交付应使用全工作区构建。

## 新终端与服务器构建环境

`source` 只改变当前 shell 及其子进程的环境。在另一条 SSH 会话中编译成功，不代表新登录终端已加载 ROS。若接口包提示找不到 `ament_cmakeConfig.cmake`，先加载 `/opt/ros/humble/setup.bash`；无需重新安装已有的 `ament_cmake`。

`cloud-gongjun` 的统一构建入口：

```bash
cd /data1/gongjun/code/lightning_lm_ws
LIGHTNING_BUILD_JOBS=16 bash src/lightning-lm/scripts/build_workspace.sh
```

不要在同一工作区同时运行两次构建。服务器的二进制和 build/install 缓存不能直接复制到 WSL 复用。

## 源码目录改名后的缓存修复

若 `src/lightning_lm_msgs` 改为 `src/common_msgs`，或驱动源码从旧 SDK 工作区迁入当前工作区，旧 CMake cache 可能保留原绝对路径。出现 `The source ... does not match the source ... used to generate cache` 时，在完整统一工作区执行一次：

```bash
bash src/lightning-lm/scripts/build_workspace.sh --cmake-clean-cache
```

该参数重建全部包的 CMake cache，保留源码及 build/install 目录；自定义缓存选项需重新指定。完成后普通增量构建不再需要该参数。

## 仅构建文档

Doxygen 与 Graphviz 是可选开发依赖。只更新文档时，在仓库根目录独立构建，无需 ROS、C++ 编译器或已有业务构建：

```bash
cmake -S docs -B build_docs
cmake --build build_docs --target docs
python3 docs/maintenance/check_links.py --html build_docs/docs/html
```

入口为 `build_docs/docs/html/index.html`。公式使用 MathJax，首次显示需要浏览器能够加载配置的 MathJax 资源；API 与源码导航来自本地 HTML。

若需与业务程序一同构建，在工作区根目录执行：

```bash
sudo apt-get install doxygen graphviz
source /opt/ros/humble/setup.bash
colcon build --packages-up-to lightning_lm \
  --cmake-args -DBUILD_DOCS=ON -DBUILD_TESTING=OFF
```

`BUILD_DOCS=ON` 会把 `docs` 加入默认 `ALL` 目标，因此上述 `colcon build` 完成时文档已经生成，入口页为 `build/lightning_lm/docs/html/index.html`。如只需重新生成文档，可执行：

`BUILD_TESTING=OFF` 只关闭 CTest/测试目标，用来缩短这次文档验证；它不控制 Doxygen。需要同时编译测试时可以省略该参数，`BUILD_DOCS` 的行为不变。

```bash
cmake --build build/lightning_lm --target docs
```

`--packages-up-to` 会先准备工作区内的消息/接口依赖；若改用裸 `cmake` 或 `--packages-select`，需要先 `source install/setup.bash`。普通 `colcon build` 的初始默认值是 `BUILD_DOCS=OFF`；在同一 build 目录设置过 ON 后，CMake cache 会保留它，需要显式传 `-DBUILD_DOCS=OFF` 才能关闭。

# 启动约定

优先从 `scripts/` 的运行器启动，而不是直接使用 `build/.../bin`：标准脚本会 source ROS 与安装空间、验证 `ros2 pkg prefix lightning_lm` 指向当前仓库、规范输出目录，并执行输入/输出契约检查。

离线 SLAM 示例：

```bash
scripts/run_slam_offline.sh \
  --bag /absolute/path/to/bag \
  --config "$PWD/config/default.yaml" \
  --output-dir "$PWD/runs/example" \
  --wait-ui false
```

在线入口示例：

```bash
LIGHTNING_LM_CONFIG="$PWD/config/default.yaml" scripts/run_slam_online.sh
```

完整脚本契约见 @ref script_contracts "Shell 脚本契约"；在线定位的逐步调用路径见 @ref online_localization_flow "在线定位端到端流程"；离线建图路径见 @ref offline_slam_flow "离线 SLAM 端到端流程"。

# 常见失败边界

- `ros2 pkg prefix lightning_lm` 指到别的 overlay：标准离线脚本会拒绝运行。
- bag 缺少 `metadata.yaml`、配置不存在或目标输出已存在：脚本返回 2，防止覆盖实验产物。
- 无关键帧、地图/轨迹/诊断导出失败：`run_slam_offline` 返回 4。
- 文档配置阶段提示找不到 Doxygen：只对显式打开 `BUILD_DOCS` 的构建生效；安装可选依赖后重新 configure。
- 裸 CMake 提示找不到 `diagnostic_monitor_interfaces` 等工作区包：先 source 基础 ROS 与当前工作区 `install/setup.bash`，或改用 `colcon build --packages-up-to lightning_lm`。

**代码依据：** `scripts/run_slam_offline.sh`、`src/app/run_slam_offline.cc`（前缀校验、覆盖保护以及返回码分支）
