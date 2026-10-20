@page workspace_layout 工作区与源码布局

# 统一定位与驱动工作区

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
