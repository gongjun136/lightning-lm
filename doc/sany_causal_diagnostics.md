# 主 Orin 长尾与 CAN 速度因果诊断

## 部署范围

本轮只新增诊断，不调整CAN融合权重、速度步长门槛、静止检测条件、线程数、300/500 ms保护或雷达点数。正式三雷达配置的 `output.fixed_map_transform.enabled=false` 延续用户确认的新场地坐标约定；其他外参仍保留。

需重新编译C++并部署脚本。外层脚本继续使用 `LIGHTNING_LM_RUN_MODE=diagnostic` 即自动开启。新增可选覆盖：

```bash
export LIGHTNING_LM_CAUSAL_TRACE=1   # diagnostic默认1，production默认0
export LIGHTNING_LM_RESOURCE_PROFILE=1
export LIGHTNING_LM_RESOURCE_INTERVAL_SEC=1
```

运行脚本自动将异步追踪限制在本次run的 `results/causal_trace.csv`，清除继承的旧路径；不复用旧文件。直接调用二进制的离线回放可显式设置 `LIGHTNING_LM_CAUSAL_TRACE_PATH`，必须是不存在的新文件。没有设置路径时不创建后台线程或文件。

## 新增证据

| 记录 | 可回答的问题 |
|---|---|
| `process_resources.jsonl` schema_version=2 / `system_process_cpu.top` | 每个采样区间CPU最多的20个进程：PID、start_ticks、进程名、进程组、CPU核等效量。定位是否在其他进程繁忙时失去CPU时间？ |
| `threads[].wchan` | 读取瞬间线程位于哪个内核等待通道；与已有线程CPU、切换次数、schedstat联合检查。不是区间内完整等待栈。 |
| `can_input` | 实际收到的CAN换算速度、扭矩、Header历元及本机接收追踪时刻；即使关闭速度融合仍记录。 |
| `can_update` | 实际送入ESKF的CAN样本时刻、滤波器历元、更新前后纵向速度、创新、NIS、观测标准差、接受/拒绝原因。source区分lio_predict、imu_predict、imu_replay。 |
| `lidar_update` | 同一LIO观测更新前后纵向速度、位置变化范数、yaw变化（rad）、是否接受及frame_id。 |
| `hf_rebuild` | 从LIO状态重建高频状态的前后速度/历元、位置与yaw变化；包含内部预测和CAN重放的总变化，不可与内部事件重复相加。 |
| `static_decision` | 静止进入/退出所选的CAN样本时刻和速度；退出reason是最后一个防抖样本触发的CAN/LIO/IMU条件组合，不代表三个样本都由同一条件触发。 |
| `static_transition` | 高速滤波器实际保持切换前后速度、状态时刻。进入时归零，退出不直接恢复历史速度。 |
| `sensor_process` | 每个LiDAR消息及慢IMU消息的Header、frame_id、入队等待、消费回调墙钟/调用线程CPU。IMU详细记录条件为排队≥5ms或处理≥20ms，不能用这些详细记录计算全量IMU分位数；原全量窗口摘要仍保留。 |
| `lock_wait` | IMU/LiDAR入口的生命周期锁与处理锁等待≥0.2ms的记录；未涵盖所有第三方库或OpenMP内部锁。 |

追踪每行包含wall_ns、steady_ns、TID、唯一sequence。不同事件用stamp/input_stamp/prior_stamp区分输出状态、所选输入和更新前状态历元。未测字段为NaN、未适用布尔为−1，不是零。两个高频重建状态可能处于不同历元，先比较历元才能解释速度差。重放滤波器分支中的CAN接受次数不等于独立CAN消息数。

CAN速度更新只改变速度状态；reason区分stale、invalid_measurement、nonfinite_prediction、innovation_gate、invalid_variance、nis_gate、velocity_step_gate、accepted。当前0.35m/s步长门槛用于拒绝异常更新，不是平滑器，也不是下游发布相邻速度步长的绝对上限。

## 开销与完整性

- C++追踪用8192条固定容量的缓冲、try-lock入队，独立线程约100ms一批写出并flush；估计线程不等待文件I/O或缓冲空间。满缓冲/争用会丢诊断记录，不丢传感器消息。字段数值格式化在写线程执行。
- 正常结束末尾输出 `# received=... written=... dropped=... io_errors=...`，运行元数据会摘录它。分析必须核对footer、行数、丢弃和I/O错误。SIGKILL、断电或磁盘错误可能缺末尾，不能将缺失事件当成未发生。
- 进程采样只读 `/proc/*/stat`、定位线程wchan及已有资源字段；不读取全机命令行、其他进程环境、内核栈。新/退出进程没有配对CPU差分会记录未配对/消失数；短命进程可能遗漏。没有权限的记录不等于零CPU。
- 采样器自身耗时仍保留。Top20是截断榜单，不代表所有进程CPU之和。进程名可能被截断或重复，结合PID及start_ticks判断。
- 不自动改调度优先级、CPU亲和性或内核参数。`kernel_sched_schedstats=0`时sched等待零值无解释力；wchan也仅是瞬时线索。若仍不能确定原因，再在域控获准进行短时调度跟踪，不能承诺仅靠应用日志区分所有抢占、锁和OpenMP等待。

## CAN A/B 回放约定

同一输入/地图/参数/二进制，顺序回放两组，仅切换 `system.enable_wheel_speed_dr_observation`。保持 `enable_wheel_speed_observation=true`，让两组静止检测仍用同样的CAN输入。不能用 `SANY_ENABLE_CAN_OBSERVATION=0` 替代此开关，否则同时改变静止检测。

对照中的CAN换算速度来自电机转速，不是独立车速真值。时间分析先画未对齐原始曲线，再计算动态段偏移；正偏移定义为estimate(t+lag)最接近CAN(t)，即估计曲线相对CAN滞后。不能把相关性偏移直接当成CAN网络传输时延，也不能由此断言应修改消息时间戳。

本机回放使用localhost隔离ROS域，不能将回放定位输出发送到车辆控制域；本机处理时间不作为Orin实时性验收。

## 测试

`python3 scripts/test_causal_trace.py -v`：Linux独立编译并发写入测试，核对footer/丢弃计数、不覆盖已有文件以及关闭路径。

`python3 scripts/test_resource_diagnostics.py -v`：CPU单位、PID/TID复用、系统进程排序、正式配置及真实/proc采样。

## 2026-09-09 本机回放验证

使用 `rosbag2_2026_09_04-13_48_09`，两组完整1倍速、同地图/二进制，只切换CAN速度融合；CAN静止检测保留，LIO预算固定1500。两组36593条原始CAN输入；实际诊断写入95708/152498条、丢诊断0/1条、I/O错误0，传感器队列丢弃均0。定位队列累计丢弃7/5条，健康跟踪帧7282/7284，不能当位级确定性实验。

共同20 ms网格33399个样本，关闭/开启融合：对CAN速度RMSE 0.09072/0.08868 m/s，20 ms绝对速度变化P99 0.07388/0.07157 m/s；动态段P99则为0.09857/0.09906 m/s，不能说所有工况均改善。CAN是换算电机速度，非独立精度真值。

开启融合组最大原始发布单步0.267091 m/s，事件记录指向LIO校正后高频重建0.226732 m/s，再叠加下一IMU历元0.040359 m/s变化；对应高频历元之间无直接CAN更新。已记录单次CAN更新最大0.026966 m/s。关闭融合仍有0.384971 m/s台阶，因此不能把速度不平滑全部归于直接CAN融合。上述只归因本次回放事件，不直接归因9月8日域控另一包的0.272 m/s事件。

动态段最佳相对偏移两组均为−60 ms（估计略领先CAN）；关闭融合的21个20秒动态窗口范围−360至+100 ms，且零偏移落在整体RMSE最小值5%以内。没有证据支持统一修改一个固定时间偏置。

本机已编译在线程序，通过7项相关C++测试及上述Linux资源/追踪测试；CAN接受、创新/NIS/步长拒绝原因额外有断言。Orin全栈并发实时性和新增诊断开销仍需目标机复测。此次没有修改融合权重、静止门槛或输出平滑策略。
