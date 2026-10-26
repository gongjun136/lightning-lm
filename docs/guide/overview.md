@page guide_overview 1. 系统全貌与上手

先建立“输入、算法模块、运行入口、结果”的对应关系。读完这一章，应能区分 LIO 里程计、建图、地图定位与 CGI 纯组合导航，并找到当前工作区真正使用的源码和配置。

**阅读顺序**

1. @subpage architecture "系统架构与模块地图"：先看模块输入输出和入口；对象管理细节可在修改相关代码时回看。
2. @subpage workspace_layout "工作区与源码布局"：区分项目根、colcon 工作区、源码目录与 ROS 包名。
3. @subpage build_and_run "构建与运行入口"：确认依赖、构建目标、环境加载和各模式的启动方式。
4. @subpage configuration "配置入口与运行证据"：从实际 YAML 读取位置查参数，明确单位和运行产物。
5. @subpage md_config_sany_README "SANY 共用建图与定位配置"：查正式模板、VIN/接线方式、配置解析与运行快照。
6. @subpage guide_scripts "脚本启动主线"：先按任务选择入口，沿真实启动命令和阶段导读进入 C++。
7. @subpage script_contracts "启动与实验脚本契约"：需要调整运行设置时，查参数、环境变量、工作目录及副作用。

**抓住三个问题**

- 在线定位的 `LaserMapping → LidarLoc → PGO` 各提供什么，地图绝对约束在哪一步进入？
- 建图的 `BackendPipeline` 与定位的 PGO 分别优化什么，在哪些运行模式中使用？
- `run_loc_online` 的 `lidar` 和 `ins_only` 分支分别由谁编排？

答案从 @ref architecture "架构页" 与 @ref index "首页运行模式表" 对照。已具备构建环境时，可直接进入下一章，再按需回查工作区与脚本。

改配置之前，先用 @ref config_multilidar_enable_reading "显式关闭多雷达时的解析分支" 和 @ref config_extrinsic_matrix_reading "外参矩阵读取的实际校验范围" 核对“写进 YAML”与“通过验证并生效”的区别。

下一章：@ref guide_online "在线定位主线"。
