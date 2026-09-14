# 速度平滑影子诊断

在外层运行脚本中设置（编译本次代码后生效）：

```bash
export LIGHTNING_LM_RUN_MODE=diagnostic
export LIGHTNING_LM_SPEED_SMOOTHING_SHADOW=1
bash lightning-lm/scripts/run.sh
```

默认关闭，只有显式diagnostic模式、开关为1且非精简诊断时启用；production始终禁用。取消开关或设为0即可关闭。运行清单记录请求值，启动日志的`SPEED_SMOOTHING_SHADOW enabled=`记录实际生效值。

固定30ms一阶因果平滑只在在线程序的诊断路径运行。`/PosRes`及派生位姿消息仍使用原速度，**不新增话题、不改变估计器、不改变控制输入**。每次PosRes发布后记录`SPEED_SMOOTHING_SHADOW stamp_sec=... stamp_nanosec=... raw_mps=... filtered_mps=... tau_ms=30 parking=... reset=... output=raw`。

已判定停车时影子值立即为0；恢复运动后继续因果平滑。首次、超过100ms间隔或发布门控关闭后重置；无效输入不掩盖并清空历史。一般加减速仍有滤波滞后，不能把停车优先等同于制动全过程没有延迟。

下一包需包含起步、正常减速、停车和倒车。按相同时间戳比较raw/filtered的台阶、延迟与停车切换，同时检查新增日志对处理时间长尾的影响。诊断逐帧日志有I/O成本，production无此影子计算或逐帧日志。尚未授权将filtered替换控制输出；需要影子验证和控制闭环验收。
