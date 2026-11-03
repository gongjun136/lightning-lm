@page record_localization_sensors 主 Orin 多传感器录制

脚本在主 Orin 录制定位输入、CGI 解码消息和可见的业务/诊断输出，保留设备原输出点，供后续离线建图与精度评价。脚本不启动或停止定位及驱动，不配置 CAN，不在线补杆臂。

## 1. 从 Orin 准备 CGI 驱动

当前确认 CGI 使用 `can1`、500000 bit/s。若接口为 DOWN，且现场已确认接线和速率，由从机用户在交互终端执行；`sudo` 需要该用户输入本机密码。接口已正常 UP 时先检查现有状态与速率，不重复修改：

```bash
sudo ip link set can1 type can bitrate 500000
sudo ip link set can1 up
```

在从机新终端加载 ROS 与统一工作区，确认驱动包可发现，再启动 CGI 驱动并保持运行。当前现场 Domain 为 42；若实际配置改变，主从使用经确认的一致设置。

```bash
source "$HOME/.config/sany/workspace.bash"
source /opt/ros/humble/setup.bash
source "$LIGHTNING_LM_WS/install/setup.bash"
export ROS_DOMAIN_ID=42
ros2 pkg prefix cgi430_driver
ros2 pkg prefix cgi430_interfaces
ros2 launch cgi430_driver cgi430.launch.py config:="$SANY_WS/config/cgi430.yaml" namespace:=cgi430
```

如果 `ros2 pkg prefix` 找不到包，先确认使用的是从机运行用户并重新加载上述环境；需要的话在统一工作区构建对应驱动和消息接口。Topic 可见不等于主机实际收到消息。

## 2. 主 Orin 录制

```bash
source "$HOME/.config/sany/workspace.bash"
bash "$LIGHTNING_LM_REPO_DIR/scripts/record_localization_sensors.sh"
```

脚本优先读取正在运行定位进程的有效配置与通信环境；没有运行配置时按车辆配置解析 SANY 模板。必需输入实际收到后才启动 MCAP 录制；结束时按一次 `Ctrl-C` 并等待写盘，屏幕会显示输出目录和退出码。

默认必须实际收到配置所需雷达和 IMU、启用的轮速、CGI 14 类解码字段与 `driver_status`；缺少必需项返回 2，不开始录包。默认不要求固定解，录制中参考失效的区间也保留。需要起录前严格固定解门禁时：

```bash
bash "$LIGHTNING_LM_REPO_DIR/scripts/record_localization_sensors.sh" --require-fixed
```

该选项要求九类质量时序、2/4状态、差分龄期及 sigma 门限，并有连续10条五个主字段全部推进的新解；不等于实车精度合格。

## 3. 检查与常用参数

```bash
bash "$LIGHTNING_LM_REPO_DIR/scripts/record_localization_sensors.sh" --check-only
```

`--check-only` 实际订阅并报告输入和参考质量，不录 bag；同时明确提供 `--output` 时保存审计报告。其他参数：

| 参数 | 用途 |
| --- | --- |
| `--output DIR` | 指定尚不存在的输出目录。 |
| `--config FILE` | 指定已审核定位配置，默认使用运行中配置。 |
| `--duration-seconds N` | 从录制器就绪后计时；默认0，按 `Ctrl-C` 结束。 |
| `--wait-seconds N` | 等待实际输入，默认20秒。 |
| `--cgi-prefix /NAME` | 修改默认 `/cgi430` 命名空间。 |
| `--extra-topic /TOPIC` | 附加必需话题，可重复指定。 |
| `--require-fixed` | 启用上述固定解与质量门禁。 |

脚本动态匹配 `/cgi430/`、`/livox/`、`/sensors/livox/`、`/localization/`、`/diagnostics/`，并录制轮速、转角、车体惯性、Safety 状态、TF、clock、rosout 和 parameter_events 等清单。实际类型以现场为准；图中发现但无消息或缺少消息类型依赖的话题会在报告中列出，不等于已录。录包默认4 GiB分包、256 MiB缓存，剩余空间低于2 GiB时停止并保留已有数据。

输出包括 `bag/`（MCAP）、`capture.metadata.json`、`quality_and_inputs.jsonl`、`recorder.log` 和 `snapshots/`；达到周期采样时间时另有 `clock.jsonl`。用屏幕显示的目录检查 `bag/metadata.yaml`、`ros2 bag info`、实际消息计数和退出状态；缺失、零消息及进程异常需如实保留。脚本自动记录主机开始、结束和周期时钟信息；CGI所在从机需另行保存 chrony 状态和同步配置证据。NTP跟踪误差不能证明跨机或雷达/CGI来自同一硬件测量历元。

退出码 `0` 表示按合同结束，`2` 表示参数、配置或录前必需输入/固定解门禁未通过，`3` 表示录制器异常、强制结束、空间不足或录包中缺少必需消息。退出码为0不代表参考或定位精度达标。

CGI解码消息保留对应CAN ID、载荷、Header和有效位，但并不覆盖未知CAN帧和错误帧；需要完整总线审计时，在CGI所在从机另行被动捕获并保存原始CAN。脚本不修改WebUI，不施加设备点到后轴的杆臂。离线评估前须确认真实输出点、后体姿态、固定ENU、高程、时间与杆臂；参见 @ref lidar_ins_accuracy "激光与CGI精度评估"。

目标环境为 Orin 的 ROS 2 Humble/Fast DDS；在 WSL 或其他中间件上使用时，须重新验证发现、录制与退出行为。
