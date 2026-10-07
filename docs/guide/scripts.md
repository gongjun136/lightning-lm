@page guide_scripts 脚本启动主线

先回答“这个脚本最终启动谁”，再沿执行阶段理解它如何准备输入、管理进程和保存结果。已熟悉 SLAM 时，可从在线定位入口直接进入 C++；其他模式按任务查阅。

| 当前任务 | 阅读入口 |
|---|---|
| 在线地图定位、CGI 纯组合导航 | @subpage shell_script_guide "脚本目录与在线定位入口"：SANY LiDAR、SANY INS、通用在线定位；也提供全部脚本的导读与源码目录 |
| 在线前端、单次与批量离线前端 | @subpage frontend_script_guide "前端脚本导读" |
| 在线/离线建图、保存地图、地图对齐 | @subpage mapping_script_guide "建图与地图工具导读" |
| 使用已有地图处理离线 bag | @subpage offline_localization_script_guide "离线定位脚本导读" |
| 构建、依赖、INS 录制与回放 | @subpage build_navigation_script_guide "构建与导航数据工具导读" |
| bag 复现在线行为、Voxel-SLAM 参考轨迹 | @subpage regression_script_guide "在线回归与参考实验导读" |
| SANY 故障/重定位矩阵、M3DGR 对比实验 | @subpage reproduction_script_guide "复现实验脚本导读" |

## 怎样对照源码阅读

每份导读按“核心调用 → 执行阶段 → 按需参考”组织。核心调用中的 `--config`、bag、地图和输出参数，是 Shell 与算法程序之间的交接点；参数解析、统计和资源记录围绕这个交接点展开。

源码片段由 `@snippet` 在构建时从真实脚本读取，用成对的语义注释标识区段。完整源码页自动显示行号，正文链接指向稳定的阅读区段。代码位置移动时，无需维护正文中的行号；某个区段的职责改变时，应同步调整对应说明。

首次阅读先沿调用点走。`usage()`、`cleanup()` 等函数的定义不会立即执行，`trap` 是注册退出处理，`&` 才表示后台启动，`wait` 负责取回对应进程的结果。较长的 helper 和参数清单放在折叠参考里。

需要精确调整 flags、环境、cwd 或输出副作用时，查 @ref script_contracts "脚本契约"。所有 Shell 完整源码仍可从顶部“脚本源码”进入。
