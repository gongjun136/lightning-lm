@page resource_diagnostics 资源与端到端延迟

# 主 Orin 资源诊断与配置交付

本页是 `process_resources.jsonl` 与计算、排队耗时的口径参考，适合核对字段、单位和采样限制。需要把资源窗口与业务事件、CAN融合或静止保持联合分析时，从 @ref causal_diagnostics "因果链与时间关联" 开始；启动、采集和结果路径见 @ref shell_script_guide "脚本使用指南"。分析已有run时应保留代码版本、配置与启动快照；本机采集测试不能替代目标域控的实时性验收。

## 默认行为与开销控制

`LIGHTNING_LM_RUN_MODE=diagnostic` 默认开启现有逐帧/分阶段计时，并默认每1秒进行一次独立Python资源采样。只使用Linux `/proc` 和标准库，不需sudo、psutil、pidstat或新增ROS消息。采样器根据本次启动的进程组和完整可执行文件名寻找真实 `run_loc_online` PID，不将 `ros2 run` Python父进程或另一套定位混入。进程退出后不跟踪重启实例。

`production` 默认关闭此资源采样和现有高频诊断，保留已有10秒端到端汇总。环境显式覆盖仍优先：

```bash
export LIGHTNING_LM_RUN_MODE=diagnostic
export LIGHTNING_LM_RESOURCE_PROFILE=1       # 0可单独关闭资源采样
export LIGHTNING_LM_RESOURCE_INTERVAL_SEC=1  # 0.5–60秒；建议首轮1秒
```

RSS等轻量信息每次读取，PSS/smaps_rollup默认30秒一次。记录采样器读取CPU/墙钟耗时及累计自身CPU时间；这不是零开销工具，1秒采样不能证明毫秒级调度事件的完整因果链。共享库哈希、全机所有进程命令行和完整环境变量不收集；只记录目标可执行文件SHA256与OpenMP/线程预算相关白名单环境变量。

## 数据包与资源字段

| 文件/字段 | 意义与单位 |
|---|---|
| `launch_compute_config.json` | YAML摘要、SHA256、线程/NDT环境覆盖、声明的LIO预算；不是运行占用 |
| `logs/launch_compute_config.log` | 启动时直接打印预算，未启用会警告 |
| `results/process_resources.jsonl` | metadata、逐秒sample、end/error；保留墙钟纳秒和单调时钟 |
| `cpu_core_equivalents` | 进程全部线程CPU时间增量 / 单调墙钟间隔；3.2表示该窗口平均消耗约3.2个逻辑核的计算时间 |
| `cpu_percent_one_core` | 上述值×100；320%不代表越界 |
| `cpu_percent_allowed_capacity` | 上述核当量 / 主线程允许CPU数量×100；线程独立绑核另见各线程字段 |
| `allowed_cpus/allowed_cpu_count` | 允许运行的逻辑核集合/数量；不能当成实际占满的核数 |
| `threads[].cpu_core_equivalents` | 各TID窗口CPU消耗；首次或新线程为null；短寿命线程可能无法匹配但其消耗仍计入进程CPU |
| `active_thread_count` | 窗口中有CPU tick增量的已匹配线程数量，不等于占用核数，也不是精确并行度 |
| `threads[].last_cpu/state` | 读取时最后CPU编号/状态快照，不是窗口内运行过的完整CPU集合 |
| `threads[].sched_wait_ns_delta` | 内核runqueue累计等待增量；不是互斥锁等待、sleep或业务消息排队。kernel_sched_schedstats=0时零值不能解读为没有等待 |
| `threads[].voluntary/involuntary_switches_delta` | 自愿/非自愿上下文切换计数增量 |
| `memory` | VmRSS、VmHWM、VmSize、VmSwap、RssAnon/RssFile及低频PSS，单位kB；缺权限/不可得为null或缺字段 |
| `io_delta` | 进程读写字节/调用次数增量；与系统缓存/设备I/O口径不同 |
| `threads[].wchan` | 读取瞬间的内核等待通道，不是窗口内完整等待栈 |
| `system_process_cpu.top` | 区间CPU最多的20个已配对进程；结合PID和start_ticks识别，短命进程可能遗漏 |
| `host_cpus` | 整机每逻辑核busy和iowait百分比，与定位进程占用分开；GPU/温度仍由tegrastats记录 |
| `run_metadata.txt: resource_monitor_exit_code` | 采样器退出状态；非零要查看stderr，不能把缺失文件当成零占用 |

每个sample记录`elapsed_s`，CPU分母使用实际间隔而不是假定1秒。同一进程和run中线程名和TID可与glog线程ID关联；结合 `start_ticks` 区分TID复用。OpenMP线程通常共用名字，不能仅凭名称精确归属LIO/NDT。统计采用内核tick分辨率，跨线程文件读取不是原子全进程快照。低频PSS读取失败记录错误，不阻止定位运行。

## 处理时长与排队时长

已有逐帧/窗口计时继续记录预处理、组帧、LIO、NDT、PGO、发布回调的wall/thread CPU，以及同帧`primary_to_pgo_ms`。进程CPU计时字段包含并发其他线程，不归为某个函数独占CPU；不能相加各模块P95/P99。

队列诊断字段：

- `COMPUTE_BENCH_SUMMARY module=sensor_queue source=imu|lidar`：分别按消息类型汇总进入传感器消息队列前到消费回调入口的墙钟时间，含入队锁开销，含三雷达消息但不是融合帧数；每秒输出P95/P99/max。
- `COMPUTE_BENCH_FRAME module=lidar_loc_pipeline queue_wait_ms`：对应整帧进入地图定位队列前到消费回调入口的单调墙钟时间；离线直调或无入队记录时为−1。后续处理锁等待不含在此字段内。

这些字段使用steady_clock，不与传感器Header/ROS时钟相减。不开启compute profile时不启用这些队列窗口记录；低频端到端汇总仍可能写入日志，`compute_profile.log` 的存在不能证明该开关已开启。队列丢弃数、陈旧观测丢弃数、高频latest-only队列覆盖数和相同历元候选拒绝数必须分别报告，不能混加成LiDAR丢帧。

## 域控运行

源码更新后，在原ROS及驱动依赖环境重新编译lightning；代码含C++队列计时，不只更新脚本。现场通过项目根入口 `bash "${SANY_WS}/run.sh"` 启动，LiDAR 分支进入 @ref run_sany_lidar_loc.sh "scripts/run_sany_lidar_loc.sh"；外层脚本选择正式三雷达YAML。例如一组候选配置的启动打印如下；实际参数以所选配置为准：

```text
LIO point budgets: [2200, 1804, 1500]
```

将完整run目录导出，包括`launch_compute_config.json`、`process_resources.jsonl`、stderr和profile。下次记录感知/规划/控制启停时刻，分别保留运动、停车、起步段；先不改变线程数、点数和门槛。无原始bag/独立真值时，仍不能做严格同输入精度验收或证明所有速度台阶的真实原因。

测试入口：`python3 scripts/test_resource_diagnostics.py -v`。包括CPU核当量单位、TID复用、缺失计数、正式/候选配置一致性及Linux真实进程组采样；Windows会明确跳过Linux集成测试。
