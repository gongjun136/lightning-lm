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
| `map_export` | `map_frame::ReadExportOptions`、`map_manifest::ReadExportOptions` | 地面归一化、PGM 高度带、`site_id/map_id/alignment_status` |
| 地图包 `map_manifest.yaml` 的 `output.fixed_map_transform` | `map_manifest::LoadOutputConfig` → `LoadFixedMapTransform` | 所选地图的清单优先；启动核对几何哈希；仅输出边界左乘 |

算法参数深入见 @ref eskf_theory "18 维迭代 ESKF：从预测到观测注入"、@ref backend_optimization "体素 BA、回环位姿图与分层优化"、@ref relocalization_theory "全局重定位：检索、几何验证与时间确认"。部署参数见 @ref script_contracts "Shell 启动与实验脚本契约"、@ref resource_diagnostics "资源与端到端延迟"。

SANY 正式模板、车辆选择和配置解析见 @ref md_config_sany_README "SANY 共用建图与定位配置"；标定来源、版本冻结与验收见 @ref vehicle_calibration_management "车辆标定档案与验收"。

@anchor config_multilidar_enable_reading
# 没写 enabled，就一定关闭多雷达吗？

先区分结构体默认值与 YAML 加载规则。`MultiLidarConfig` 的默认开关是：

@snippet{lineno} core/lio/multi_lidar_fusion.h config-multi-default-disabled

而 @ref lightning::LoadMultiLidarConfig "LoadMultiLidarConfig()" 每次先重置配置，再按整个节点是否存在决定分支：

@snippet{lineno} core/lio/multi_lidar_fusion.cc config-multi-enable

<b>只追三个返回位置：</b>完全没有 `multi_lidar` 块时保持默认关闭并返回 true；块存在但漏写 `enabled` 时，三目运算取 true；显式 `enabled:false` 则立即返回，后面的 topics、外参和时间参数都不读取。因此“读取成功”与“功能开启”是两个条件，关闭分支的成功也不证明其中暂存的参数有效。

<b>手工推演：</b>旧对象即使先前启用过，重新加载不含该块的 YAML 后也会被第一行重置为关闭。`multi_lidar: {}` 则先被视为开启，随后因缺少 topics/extrinsics 失败。若显式关闭，块内一个错误类型的 `frame_period` 不会到达转换语句；以后再开启仍需重新检查它。

开启后还要通过实际的时间关系：

@snippet{lineno} core/lio/multi_lidar_fusion.cc config-multi-timing

在其他配置均有效、数值有限的假设下，周期 0.1 s、容差 0.02 s、重排窗口 0.1 s 能通过这段比较；容差改成 0.05 s 会因 `>= 0.5*frame_period` 失败，重排窗口改成 0.099 s 也失败。这里展示的是比较条件，不是对所有时间字段的完整有限性校验。

<b>停止位置：</b>确认默认构造、节点存在性与提前返回即可。这个 true 只属于多雷达配置读取器，不代表整个系统已初始化或收到传感器。组帧运行时如何使用这些时间继续 @ref sensor_frame_ready_reading "时间桶与排出条件"；完整加载与校验见 @ref multi_lidar_fusion.cc "multi_lidar_fusion.cc"，默认字段见 @ref lightning::MultiLidarConfig "MultiLidarConfig"。

@anchor config_extrinsic_matrix_reading
# 外参能够加载，是否就说明方向和车辆标定正确？

这里验证 C++ 读取器如何解释 `multi_lidar.extrinsics.lidarN.T`。`ReadTransform()` 接受行优先的 3×4 或 4×4 数列；先检查长度、有限性和可选齐次末行：

@snippet{lineno} core/lio/multi_lidar_fusion.cc config-extrinsic-shape

再读矩阵元素的索引和旋转检查，而不是只看报错中的 SO(3)：

@snippet{lineno} core/lio/multi_lidar_fusion.cc config-extrinsic-rotation

平移取第 4、8、12 个元素，旋转取每行前三个。代码比较 `R.transpose()*R-I` 的最大元素误差和行列式误差，门限均为 `1e-3`；它不会在这里把矩阵自动正交化。4×4 末行另用 `1e-9` 门限。

<b>手算：</b>数组 `[1,0,0,1, 0,1,0,2, 0,0,1,3]` 得到 R=I、t=(1,2,3) m。按当前“雷达到主雷达”的方向，点 `(4,0,0)` 会成为 `(5,2,3)`。若误填逆变换，t=(-1,-2,-3) 仍能通过这些数值检查，却得到 `(3,-2,-3)`；同样，误把毫米数值当米也不一定违反矩阵检查。实际点变换证据见 @ref sensor_fusion_coordinates_reading "融合坐标与点时间"。

<b>边界与停止点：</b>“数值合法”不证明变换方向、平移单位、传感器身份或标定精度正确。档案生成脚本 `scripts/resolve_vehicle_config.py` 的 `resolve()` 还要求 16 个元素、主雷达单位变换，并以 `1e-6` 检查旋转；不能把直接调用 C++ 的接受条件等同于档案验收。读到 R/t 索引和拒绝条件即可，完整读取器见 @ref multi_lidar_fusion.cc "ReadTransform / LoadMultiLidarConfig 源码"，物理来源和验收见 @ref vehicle_calibration_management "车辆标定档案"。

# 重现一轮运行的最小信息

记录源码 commit/未提交修改、二进制与共享库、最终 YAML 和覆盖环境、传感器 topic 与时间范围、地图 metadata/数据库、线程数、运行模式和原始日志。固定这些条件后才讨论算法回归；“相同 YAML”不能证明二进制相同。
可用 `scripts/summarize_run_config.py`、`summarize_rosbag2_topics.py`、标准运行脚本生成清单；参数解析以各脚本 `--help` 为准。

# 输出目录

运行数据保存在 `runs/` 或外部数据目录，生成的网站放在 build 树。`docs/assets/` 只放仍被说明引用的资源；不要把完整 bag、一次实验的大量中间图或构建产物放进文档目录。
