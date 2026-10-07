@page guide_overview 1. 系统全貌与上手

先建立“输入、算法模块、运行入口、结果”的对应关系。读完这一章，应能区分 LIO 里程计、建图、地图定位与 CGI 纯组合导航，并找到当前工作区真正使用的源码和配置。

**阅读顺序**

1. @subpage architecture "系统架构与模块地图"：先看模块输入输出和入口；对象管理细节可在修改相关代码时回看。
2. @subpage workspace_layout "工作区与源码布局"：区分项目根、colcon 工作区、源码目录与 ROS 包名。
3. @subpage build_and_run "构建与运行入口"：确认依赖、构建目标、环境加载和各模式的启动方式。
4. @subpage configuration "配置入口与运行证据"：从实际 YAML 读取位置查参数，明确单位和运行产物。
5. @subpage script_contracts "启动与实验脚本"：需要运行或修改脚本时，查参数、环境变量、工作目录及副作用。

**抓住三个问题**

- 在线定位的 `LaserMapping → LidarLoc → PGO` 各提供什么，地图绝对约束在哪一步进入？
- 建图的 `BackendPipeline` 与定位的 PGO 分别优化什么，在哪些运行模式中使用？
- `run_loc_online` 的 `lidar` 和 `ins_only` 分支分别由谁编排？

答案从 @ref architecture "架构页" 与 @ref index "首页运行模式表" 对照。已具备构建环境时，可直接进入下一章，再按需回查工作区与脚本。

下一章：@ref guide_online "在线定位主线"。
