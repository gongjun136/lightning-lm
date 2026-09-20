# M3DGR 后端评测结果

评测日期：2026-07-15。所有 Lightning-LM 运行使用同一 Livox Mid360 + 内置 IMU
前端输入；RTK 只在运行结束后参与评测，不进入任何优化器。轨迹采用固定尺度 SE(3)
对齐，禁止 Sim(3)。真值回环定义为时间间隔至少 30 s 且 RTK 距离不超过 5 m。

## 主结果

| 序列 | 方法 | ATE RMSE (m) | RPE@10m RMSE (m) | Loop P/R | Wall (s) | Peak RSS (MB) |
|---|---|---:|---:|---:|---:|---:|
| Grass02 | Lightning 旧后端 | 0.4364 | 0.7382 | N/A | 276.36 | 1346.9 |
| Grass02 | BA + BTC + HBA | **0.4297** | 0.7363 | **1.000 / 1.000** | **109.06** | 1834.9 |
| Grass02 | 完整 Voxel-SLAM | 0.4313 | **0.7305** | N/A / 0.000 | 192.94 | 1746.7 |
| Outdoor04 | Lightning 旧后端 | 5.4570 | 1.7755 | N/A | 929.24 | 4693.3 |
| Outdoor04 | BA + BTC + HBA | 0.4808 | 0.1510 | **1.000 / 0.200** | **537.38** | 6490.7 |
| Outdoor04 | 完整 Voxel-SLAM | **0.2806** | **0.1087** | N/A / 0.000 | 823.54 | **3680.3** |

Grass02 上新后端 ATE 优于旧后端和完整 Voxel-SLAM。Outdoor04 上新后端消除了旧后端
错误距离回环造成的约 10.9 m 终点校正，ATE 相对旧后端下降 91.2%，但仍未达到完整
Voxel-SLAM 的前后端精度。Voxel-SLAM 在两条序列都检索到了 BTC 候选，但默认阈值和
ICP 没有接受任何 `addedge`；其精度优势主要来自自身前端和局部/全局 BA，而非闭环边。

耗时不是完全同构的实时性对比：Lightning 离线读取不做传感器时间节拍，Voxel-SLAM
通过 ROS1 rosbag 以 1× 播放。表中 Wall 是运行元数据的完整墙钟时间，Peak RSS 是同一
8 核分配下的进程树采样。Outdoor04 当前 HBA 与大点云同时驻留，内存仍是后续优化点。

## 回环核验

- Grass02：接受描述子 152→2，时间间隔 150 s，RTK 距离 3.157 m；该回环残差仅
  0.097 m / 0.458°，因此未强制改图。
- Outdoor04：接受描述子 753→0 和 769→2，RTK 距离分别为 0.929 m 和 0.935 m；
  两者残差都低于图优化门限，未施加不必要校正。
- Outdoor04 的严格 RTK 协议包含 5 个连续回环事件，当前覆盖终点事件，因此 episode
  Recall 为 0.2；候选级统计和全部拒绝原因保存在评测 JSON。
- 完整 Voxel-SLAM 的 Grass02/Outdoor04 接受边均为 0，所以 Precision 无定义、Recall
  为 0，而不是将 Precision 人为记为 0。

## HBA 提交消融

`require_applied_loop_for_commit=false` 会允许没有图修正时的 HBA 结果直接提交：

| 序列 | 事务回滚 ATE (m) | 强制提交 ATE (m) | 结论 |
|---|---:|---:|---|
| Grass02 | 0.4297 | **0.4277** | 小幅改善 |
| Outdoor04 | **0.4808** | 0.4843 | ATE/RPE 均退化 |

因此默认保持事务策略：有效 BTC 回环仍立即触发后台 HBA，数据结束/ROS 服务仍强制
执行；只有闭环约束确实改变过全局图时才提交点云层优化，否则完成计算后回滚。这一
默认值优先保证两条验证序列都不退化。

完整机器可读报告位于工作区外部实验目录：

- `_m3dgr_work/bench_v2/results/Grass02_all_backends.json`
- `_m3dgr_work/bench_v2/results/Outdoor04_all_backends.json`
