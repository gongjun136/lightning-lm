@page workspace_layout 工作区与源码布局

# 统一定位与驱动工作区

同一工作区包含 15 个 ROS/CMake 包：定位、CGI-430 驱动、Livox SDK2、Livox 驱动，以及 `common_msgs` 内的 11 个接口包。`cgi430_interfaces` 的消息契约保持不变，已从独立 CGI SDK 仓库迁入 `common_msgs`。

```text
lightning_lm_ws/src/
├── lightning-lm/                  # 包 lightning_lm
├── common_msgs/
│   ├── cgi430_interfaces/
│   ├── livox_ros_driver/          # 消息包 livox_ros_driver2
│   └── ...                       # 其他共享接口
├── cgi430_driver/                 # 仓库 URL 仍为 gongjun136/cgi430_sdk
├── livox_ros_driver2/             # 实际驱动包 livox_ros_driver2_core
└── Livox-SDK2/                    # 包 livox_sdk2
```

先检出 Lightning 仓库，再从工作区根目录导入配套仓库。清单使用 GitHub SSH 远端，开发机需已配置对应仓库的 SSH 访问：

```bash
vcs import src < src/lightning-lm/workspace.repos
bash src/lightning-lm/scripts/build_workspace.sh
```

`workspace.repos` 固定配套源码版本；`common_msgs.repos` 仍保留仅导入消息仓库的用途。不要把旧 `livox_sdk/src/common_msgs` 再复制进来，每个包名在工作区中只能出现一次。

`cgi430_sdk` 仓库保留原 Git 历史，根目录现在直接是 `cgi430_driver` 包，协议库与运行脚本随驱动维护。`common_msgs` 只存接口；定位与驱动分别声明对接口的依赖，由 colcon 确定构建顺序。不再加载外部 CGI 安装空间。

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
