@page guide_mapping 5. 后端优化与建图

这一章完整覆盖离线和在线建图。先理解关键帧怎样进入后端、优化位姿怎样写回，再对照 BA、回环位姿图与 HBA 的目标和源码，最后沿运行流程看到地图导出。

**阅读顺序**

1. @subpage backend_module "后端模块：关键帧、回环与优化"：区分局部 BA、BTC 检测、PGO 和 HBA 的职责与提交条件。
2. @subpage backend_optimization "后端公式：体素 BA 与分层优化"：对照最小特征值目标、充分统计、位姿增量和 LM 接受条件。
3. @subpage offline_slam_flow "离线建图：bag 到地图包"：完整理解同步回放、前后端交接、尾帧处理与导出。
4. @subpage online_slam_flow "在线建图：ROS 到地图服务"：对比 ROS 回调、后端 worker、保存服务与退出行为。
5. @subpage backend_configuration "后端配置与评测方法"：需要运行对照实验时查参数、评测口径与 BTC 来源。

**关键公式怎样落地**

| 问题 | 代码入口 | 阅读重点 |
|---|---|---|
| 为什么平面体素的目标是最小特征值？ | `voxel_bundle_adjustment.cc` 的 `EvaluateFactors` | 法向被消去，当前目标按有效体素累加 |
| 怎样避免每轮重新变换所有点？ | `PointCluster::Transformed/Covariance` | 点数、一阶和二阶统计在刚体变换下的组合 |
| 优化成功为何还可能不提交？ | `BackendPipeline`、`HbaLoop` | 几何求解与系统接受/提交门控的区别 |

图优化残差与信息权重在 @ref pose_graph_theory "位姿图专题" 中统一解释。建图 PGO 调整历史关键帧，定位 PGO 融合地图观测与相对运动；阅读同一求解器时仍需确认调用场景。

地图导出后的格式与坐标：@ref guide_map_output "地图、坐标与输出"。前端基础：@ref guide_lio "传感器与 LIO"。
