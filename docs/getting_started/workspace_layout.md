@page workspace_layout 工作区与源码布局

# 统一定位与驱动工作区

Livox SN 使用各雷达节点的 `device_info.*` ROS 2 参数按需查询，不再发布专用身份话题。
旧 `LidarDeviceInfo` 消息已删除；更新后使用干净的 build/install 目录，避免残留旧生成文件。
接口检查器会拒绝旧消息与旧发布器组合，依赖版本以 `workspace.repos` 为准。
可运行 `ros2 param get /livox_front_publisher device_info.serial_number`，或用
`ros2 run livox_ros_driver2_core collect_lidar_identity.py --output /tmp/lidar_identity.json` 导出清单。
主从 DDS 不互通时分别指定 `--roles front left right` 和 `--roles rear`。
参数为最近确认身份，不表示设备在线。现场升级仍需驱动重启窗口与实际收数核验。

主、从 Orin 使用同一份已安装 `run_lidar.sh`，调用 `run_lidar_bundle.py --role auto`，
按 hostname 区分 `ubuntuM`（主）和 `ubuntuS`（从），未知 hostname 拒绝启动。
两端须分别更新和构建；从机仅需 `colcon build --packages-up-to livox_ros_driver2_core`。
每次启动保存配置、完整日志、短时数据/SN/HMS 检查和节点退出码，
通过 `runs/latest_lidar_master` 或 `latest_lidar_slave` 查找最近运行。需按驱动 `docs/operations/run_bundle.md`
更新根脚本；仅更新依赖并不会替换库外现场入口。主从需分别升级，从机旧版不支持 SN 参数时检查会明确失败。

2026-10-27 主从 Orin 已部署统一入口，分别原生构建并通过 8 项功能测试，四台雷达数据、SN 和设备状态复查通过。实测二进制构建版本为 `9a33d6b`，清单后续 `1718003` 仅补充部署文档。Ctrl-C 仍复现既有退出码 -11，运行包已正确捕获；一次 HMS 超时和后续独立复查均保留。定位进程恢复后仍因全局重定位匹配失败无有效位姿，不能宣称整车健康通过。完整证据与回滚备份位置见驱动手册 `docs/operations/run_bundle.md`。

同一工作区包含 16 个 ROS/CMake 包：定位、CGI-430 驱动、Livox SDK2、Livox 驱动，以及 `common_msgs` 内的 12 个接口包。`cgi430_interfaces` 的消息契约保持不变，已从独立 CGI SDK 仓库迁入 `common_msgs`。

2026-10-06 同步公司 `main` 后新增 `diagnostic_runtime_interfaces` 与 `driving_evaluation_interfaces`；`CurrentGear` 迁入 `geosun_msgs`，旧 `domain_vcu_can_bridge` 接口包退出源码。部分现有 CAN 消息新增 `comm_header`，部署时需同步重编发布端和订阅端。旧录包兼容性不由本次编译验证保证。

```text
lightning_lm_ws/src/
├── lightning-lm/                  # 包 lightning_lm
├── common_msgs/
│   ├── cgi430_interfaces/
│   ├── livox_ros_driver/          # 消息包 livox_ros_driver2
│   └── ...                       # 其他共享接口
├── cgi430_driver/                 # 仓库 gongjun136/cgi430_driver
├── livox_ros_driver2/             # 实际驱动包 livox_ros_driver2_core
└── Livox-SDK2/                    # 包 livox_sdk2
```

先检出 Lightning 仓库，再从工作区根目录导入配套仓库。清单使用 GitHub SSH 远端，开发机需已配置对应仓库的 SSH 访问：

```bash
vcs import src < src/lightning-lm/workspace.repos
bash src/lightning-lm/scripts/build_workspace.sh
```

`workspace.repos` 是唯一依赖清单，固定消息、驱动和 SDK 的配套源码版本。不要把旧 `livox_sdk/src/common_msgs` 再复制进来，每个包名在工作区中只能出现一次。

更新定位分支后，依赖应重新按该分支的清单检出。`git fetch` 会获取远端提交，
清单中的固定 SHA 决定实际使用版本；重复导入同一份清单不会自动改用依赖的 `main`。
2026-10-13 的 DTC 配套版本为 `common_msgs@3f918b3` 和 `livox_ros_driver2@713e92e`，
分别提供 DTC 消息与四雷达 DTC 聚合节点；`/safety/throstle` 复用已有的
`geosun_msgs/msg/Safety` 的 `runstate` 字段。
发布新的跨仓库接口时，须同时更新配套清单并核对导入后的实际 HEAD。
源码检出完成后，还须重新构建并加载当前工作区的安装环境。

2026-10-20 修正：当前消息依赖固定为 `common_msgs@b7c005f`，包含心跳代码
使用的 `NodeHeartbeat.active_fault`。旧版本 `3f918b3` 不包含该字段。
依赖变更必须同步更新 `workspace.repos`；仅更新定位仓库或复用本机旧安装目录不构成验证。
构建脚本会先检查依赖 HEAD、未提交修改及其与 `workspace.repos` 的一致性。手工构建前也应执行：

```bash
python3 src/lightning-lm/scripts/check_workspace_repos.py --workspace "$PWD"
```

雷达依赖以 workspace.repos 的精确哈希为准。四台驱动各自在自身心跳 active_fault
上报 DTC，定位直接消费四路 diagnostic_state；独立聚合节点及 ID 13 已移除。
检查器同时验证心跳字段、同步证据结构/版本、四路心跳 ID 和上报路径。
用 `python3 src/lightning-lm/scripts/test_workspace_interfaces.py` 运行接口回归测试。
部署时停止旧聚合进程，使用新的安装目录，实测四路诊断、五路心跳、Safety 与位姿；
静态接口检查不能替代现场消息和故障恢复验收。

四雷达部署还必须同步从 Orin 的后雷达。节点启动稳定后，在与生产节点相同的
ROS 域和 FastDDS 配置环境中执行：

```bash
python3 src/lightning-lm/scripts/check_runtime_interfaces.py --workspace "$PWD" --duration 25 --require-healthy
```

该只读探针不订阅点云；要求四路实际诊断的角色、源 ID 和 `config_version` 一致，
各设备 active_fault 仅含自身当前码，Safety、五路心跳和位姿持续收到。健康验收要求采样开始 5 秒后 Safety 为 0、
四路可用性为 1。失败返回非零，输出 JSON 原始快照和错误。当前用于四雷达拓扑。
图发现不等于消息可接收；尤其要与生产节点使用相同的 DDS 网卡白名单。

检查仅验证源码版本，不验证安装空间。消息更新后用 `colcon build --packages-up-to
lightning_lm` 重编消息与消费者，并重新加载安装环境；不要仅构建消费者。

`cgi430_driver` 仓库保留原 Git 历史，根目录直接是 `cgi430_driver` 包，协议库与运行脚本随驱动维护。`common_msgs` 只存接口；定位与驱动分别声明对接口的依赖，由 colcon 确定构建顺序。不再加载外部 CGI 安装空间。

## Local compatibility mapping

On the maintained Windows workspace, the physical historical-data checkout
remains at:

```text
F:\SLAM_AI_KnowledgeBase\code\WSL_Ubuntu_22.04\lightning-lm
```

The workspace entry below is a directory junction to that existing checkout:

```text
F:\SLAM_AI_KnowledgeBase\code\WSL_Ubuntu_22.04\lightning_lm_ws\src\lightning-lm
```

Consequently, historical `data/`, `runs/`, `.tmp_sany_analysis/`, build,
install, and log paths remain available through the old path and are also
visible through the workspace path. Do not delete the old path while the
junction is in use.

Livox SDK、Livox 驱动和 CGI 驱动在本地也通过目录联接纳入 `src`，复用原 Git 检出；服务器采用上述普通目录布局。目录联接目标仍是源码唯一副本，不能删除。
