# 2026-10-08 SANY CAN 消息版本冲突与启动静默退出

## 现场证据

域控 `sany-03-master` 部署提交 `5597e25`，配置 SHA256 为
`acfe837565193c978809f94190f30cc3e3936bb647807ff65777175556321e4b`。
运行目录 `sany_3lidar_20261008_111251` 的日志显示：

- 11:13:05，CAN Header 比 IMU 落后约 190,095,821 秒，轮速融合因无新鲜样本退化。
- CAN 桥的 `geosun_msgs/SpeThrCAN4` 定义是 `Header,x,y`；新定位工作空间的同名定义是 `TopicCommHeader,Header,x,y`。后者在开头增加了 `uint16,uint32,uint64` 三个字段。
- 用 CAN 桥自身的接口读取，Header 时间正常，frame_id 为 `spe_thr_can4`；用新接口读取旧数据，Header 和 comm_header 都发生错位，部分 frame_id 字节被解释成时间戳。
- 11:24:19 起出现地图匹配失败和停止发布，此前已有持续的雷达速度更新降级，内部速度估计达到约 9 m/s。接口冲突确定存在；缺少本次原始传感器 bag，不能仅凭日志断言它解释了全部运动工况下的失效。
- 控制工作空间与定位工作空间的 `lightning/VehiclePose` 字段布局一致，只有 speed 注释不同。

## 修复

定位使用序列化订阅，在 CAN 入口显式识别上述两种 CDR1 布局。必须完整消费消息，检查时间字段、字符串长度和终止符、有限数值及消息定义的转速/转矩范围；无法唯一识别的消息直接拒绝。保留原始 Header 时间和既有轮速缩放，不用接收时刻替代传感器时刻。启动日志打印实际识别的 CAN 布局。

启动脚本使用 `check_localization_inputs` 接收真实输入，要求每路至少三条时间递增且最近一秒收到的消息，并检查 CAN 与主 IMU 时间差小于 0.25 秒。只存在于 ROS 图中的话题不再被视为就绪；超时会逐项显示接收计数、接收年龄和时间差，保存在 `logs/input_check.log`。

旧脚本的 `actual_type="$(ros2 topic type ...)"` 在子命令返回非零时会被 `set -e` 终止，而子命令的 stdout 被收进变量，没有打印。用“topic list 成功、topic type 输出 Unknown topic 后返回 1”的故障注入，能够复现只打印 Waiting 行就返回终端的现象。现场失败实例没有保存退出码和这一步的输出，无法追溯当次 CLI 失败的底层原因。新流程移除这组 CLI 查询，并为其他意外错误打印命令、行号、退出码和运行目录。

## 部署与验证

需一起更新算法库、`run_loc_online`、`check_localization_inputs` 和启动脚本，不能只复制脚本。CAN 桥、控制程序和它们的消息包无需改变。

现场传感器使用 `$HOME/Documents/ros2/fastdds.xml` 的 UDP 配置。非交互 SSH 不会自动继承交互 shell 的该变量；启动时必须加载相同的 `FASTRTPS_DEFAULT_PROFILES_FILE`。缺少配置时现场可发现话题却收不到消息，新检查会明确失败。

验证命令：

```bash
ctest --test-dir "$LIGHTNING_LM_WS/build/lightning_lm" \
  -R 'sany_wheel_speed_wire_test|sany_localization_output_test|localization_input_locking_test|localization_pgo_test' \
  --output-on-failure
python3 "$LIGHTNING_LM_REPO_DIR/scripts/test_sany_lidar_launcher.py"
```

解码测试包含旧布局固定字节样本、当前 ROS 生成器的序列化数据、不同字符串对齐、大小端、截断、额外尾部、非法时间和 NaN。启动回归测试覆盖输入检查失败时不启动算法，以及算法错误退出时保留退出码与日志。

### 本次域控验证结果

4 项 C++ 回归与 2 项脚本回归全部通过。修复版本通过正式 `~/Sany/run_loc.sh` 启动，运行目录为 `sany_3lidar_20261008_123449`；约 190 秒后由诊断用 timeout 发送 SIGINT 正常停止，不是定位故障。

| 指标 | 修复前短时对照 | 修复后 |
| --- | ---: | ---: |
| CAN 输入累计数 | 3,116 | 9,075 |
| CAN 时间戳拒绝累计数 | 3,115 | 0 |
| 雷达滤波器接受轮速累计数 | 0 | 8,962 |
| IMU 滤波器接受轮速累计数 | 0 | 9,189 |
| 独立订阅 pose_vel 条数 | 9,264 | 21,658 |
| 独立订阅连续输出跨度 | 57.24 秒 | 139.71 秒 |
| 最大接收间隔 | 47.8 毫秒 | 50.8 毫秒 |
| 大于 0.5 秒的接收间隔 | 0 | 0 |

累计计数的统计时长不同，不能用于比较性能。两次都是当前现场短时观测，修复后观测期间 CAN 速度为零、输出静止，没有复现此前的运动工况。本次没有修改算法门限。

后续现场反馈：用户已使用修复版本完成行驶测试，反馈没有问题，并确认提交 GitHub。这是用户提供的实车验证结果；本记录未采集该次行驶测试的原始 bag，也不将上述静止观测指标视为行驶指标。

旧程序、旧启动脚本和原始环境配置备份于域控
`~/Sany/lightning_lm_ws/deployment_backups/20261008_can_wire`。
现场 `~/.config/sany/workspace.bash` 已补充 DDS 配置变量，使 SSH 与交互终端使用相同传输配置。
