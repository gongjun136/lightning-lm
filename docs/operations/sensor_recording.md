# 域控传感器录制

在主域控启动实际启用的雷达驱动；需要 CGI 参考时，还需启动 CGI 驱动并确保 ROS Domain 一致。录制不要求在线定位进程运行。先加载该机工作区配置：

```bash
source "$HOME/.config/sany/workspace.bash"
```

基础采集默认要求所有配置的雷达和主 IMU 实收，其他 IMU、CGI 和轮速可选：

```bash
bash "$LIGHTNING_LM_REPO_DIR/scripts/record_localization_sensors.sh" --mode basic --check-only
```

预检通过后录制：

```bash
bash "$LIGHTNING_LM_REPO_DIR/scripts/record_localization_sensors.sh" --mode basic
```

需要 CGI 精度参考时使用 `reference`，它额外要求纬度、经度、高程、姿态、速度、位置／姿态／速度 sigma、导航状态、时间和驱动状态话题实收：

```bash
bash "$LIGHTNING_LM_REPO_DIR/scripts/record_localization_sensors.sh" --mode reference --check-only
```

预检通过后录制：

```bash
bash "$LIGHTNING_LM_REPO_DIR/scripts/record_localization_sensors.sh" --mode reference
```

预检不创建会话目录，输出的 `missing_required` 应为空，订阅不得报错。两种模式都尝试录制配置的轮速话题，包括录制过程中才出现的话题；配置未指定时使用 `/SpeThrCAN4_topic`。加 `--require-wheel` 将该话题设为必需。仅 `reference` 可加 `--require-fixed`，要求开始前连续 10 条参考解通过工程质量门限；这不保证后续全程固定。预检和正式录制使用相同的门槛参数。

录制器会加载 ROS、工作区及 DDS 环境，自动创建 `sensors_时间戳` 会话，输出根目录优先使用 `LIGHTNING_LM_OUT_ROOT`，以终端打印的路径为准。可选 `--output <新目录>` 指定位置；它不能与 `--check-only` 同用。开始预检或录制前应确认配置对应现场车辆和实际传感器；默认优先读取运行中定位器配置，没有定位器时解析仓库模板。如需指定配置，在预检和录制均使用 `--config <文件>`。

看到“录制器已就绪”后采集；结束时按一次 Ctrl-C，等待写盘。也可加 `--duration-seconds <秒数>`，从就绪时计时并自动停止。录制失败可能留下未完成会话，检查原因后使用新目录重试。

保留整个会话：`bag/` 包含 MCAP 和 `metadata.yaml`，`config.yaml` 是解析后的采集配置，`snapshots/` 保存可用输入配置，`capture.metadata.json`、`quality_and_inputs.jsonl` 和 `recorder.log` 记录采集证据。检查退出码、`missing_required_in_bag`、各话题 bag 消息数和 `wheel_summary`；轮速监测的到达间隔与 bag 消息数应结合判断，不能将发现话题或监测断流直接当作源端丢包率。需要轮速换算时核对会话中的 Topic 和比例。快照缺失或其来源、哈希与现场配置不一致时，先找回可信配置；复制数据时应连同整个会话一起校验。

bag 是 ROS 消息记录，不是未知 CAN ID 或错误帧的完整原始总线。设备状态、时间同步、配置正确性以及参考质量需另行核验。
