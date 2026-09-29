@page configuration 配置入口、单位与运行边界

# 从实际加载点查参数

配置采用 YAML + 程序 flags + 脚本环境变量。字段在不同入口可能有不同默认值和覆盖关系；不要维护一份脱离源码的全部默认值副本。先保留运行时配置与哈希，再按下表找到解释层。

| 参数域 | 解释与加载代码 | 修改时必须核对 |
|---|---|---|
| 传感器、时间、外参、多雷达 | `LaserMapping::Init`、`LoadMultiLidarConfig` | topic 类型、秒/毫秒、主雷达、外参方向 |
| IMU 初始化与滤波 | `ImuProcess`、`ESKF::Options` | 静止条件、噪声尺度、状态传播与门控 |
| LIO 体素、点数、关键帧 | `LaserMapping`、自适应控制器 | 点数保护、健康反馈、几何覆盖 |
| `backend.mode/local_ba/btc/hba` | `ReadBackendMode`、`BackendPipeline::Init` | system 总开关、提交条件与线程 |
| `relocalization` | SOLiD/BTC、`LidarLoc::Init` | 数据库与地图一致、细化后端、确认阈值 |
| 定位队列/输出 | `Localization::Init`、`LocSystem::Init` | 队列容量、最大帧龄、失效停止发布 |
| `map_export` | `map_frame::ReadExportOptions` | 地面归一化、PGM 高度带、覆盖行为 |
| `output.fixed_map_transform` | `LoadFixedMapTransform` | 仅输出边界左乘，不修改内部地图 |

算法参数深入见 @ref eskf_theory "18 维迭代 ESKF：从预测到观测注入"、@ref backend_optimization "体素 BA、回环位姿图与分层优化"、@ref relocalization_theory "全局重定位：检索、几何验证与时间确认"。部署参数见 @ref script_contracts "Shell 启动与实验脚本契约"、@ref resource_diagnostics "资源与端到端延迟"。

# 重现一轮运行的最小信息

记录源码 commit/未提交修改、二进制与共享库、最终 YAML 和覆盖环境、传感器 topic 与时间范围、地图 metadata/数据库、线程数、运行模式和原始日志。固定这些条件后才讨论算法回归；“相同 YAML”不能证明二进制相同。
可用 `scripts/summarize_run_config.py`、`summarize_rosbag2_topics.py`、标准运行脚本生成清单；参数解析以各脚本 `--help` 为准。

# 输出目录

运行数据保存在 `runs/` 或外部数据目录，生成的网站放在 build 树。`docs/assets/` 只放仍被说明引用的资源；不要把完整 bag、一次实验的大量中间图或构建产物放进文档目录。
