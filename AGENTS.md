# Dependency and deployment changes

- Before changing interfaces or consuming changes from another repository, inspect
  `workspace.repos` and the actual dependency Git HEADs. A successful build against
  an existing install overlay does not prove the pinned dependencies are compatible.
- When a dependency changes, update its exact commit in `workspace.repos` in the
  same delivery. Use `workspace.repos` as the only dependency manifest. Ensure the
  dependency commit is available on the declared remote before pushing the consumer.
- Before build validation and before committing/pushing, run
  `python3 scripts/check_workspace_repos.py --workspace <workspace-root>`.
  Resolve mismatched or dirty dependency checkouts; never reset user changes.
- Rebuild changed message packages and their consumers. Report the actual tested
  dependency commits and distinguish source checks, compilation, tests, and deployment.
- Matching Git hashes alone does not establish compatibility. Run the interface
  checks in `check_workspace_repos.py` and `test_workspace_interfaces.py`; verify
  live producer/consumer topics and received data after deployment. Never treat
  graph discovery alone as proof that messages are received or Safety is healthy.
- For dependency fixes, validate from the pinned source in a fresh build/install
  directory so cached generated headers cannot conceal missing message fields.

The 2026-10-20 Orin failure was caused by `common_msgs` remaining pinned to
`3f918b3` although the consumer needed `NodeHeartbeat.active_fault` from `b7c005f`.
The local newer dependency hid this mistake. Keep the checks above in future work.

# Diagnostic contract maintenance

- Changes to Safety, lidar diagnostics, incident capture or relevant configuration must update the diagnostic parser and regression cases, then review the contract before running `python3 scripts/check_diagnostic_contract.py --update`. Do not merely refresh hashes. CMake configuration and CTest check the reviewed contract.
- Keep the incident guide and diagnostic output aligned with actual evidence, build markers and coverage limits; distinguish source checks from validation of the running binary.

# 代码修改前同步

- 每次开始新的代码修改任务前，确认仓库、分支、目标上游及工作区和暂存区状态，先 fetch 目标远端；干净且可快进时使用 `git pull --ff-only`，否则通过安全整合或基于最新目标上游的隔离工作树完成同步后再修改。同步失败时不得在未同步的旧代码上继续新代码编辑。
- 保留已有本地和暂存改动；不得擅自丢弃、覆盖、stash、强推或切换无关分支。目标上游不明确时先确定正确目标，不猜测分支。
- 修改固定依赖前读取工作区根目录的 `../../AGENTS.md` 和 `workspace.repos`，并继续遵守本文件的依赖提交、接口验证与交付约束。

# 定位日志与现场轨迹约定

- 定位日志和 `/localization/debug_message` 的人类可读正文已中文化；后续定位程序、诊断工具和相关运行入口新增或修改面向人的日志、状态、警告和错误提示时使用中文。稳定的 `event=` 事件码、机器字段名、单位、枚举、ROS 接口以及第三方原始错误详情保持原有契约；诊断解析和历史日志兼容性必须同步维护。
- SANY LiDAR 现场默认只保存一份全程 TUM 轨迹 `results/trajectory_published_rear_axle.tum`：按每次 `/localization/pose_vel` 发布调用完成记录对应消息的后轴中心业务地图位姿，时间列由 ROS 消息头的整数秒和纳秒写成九位小数，避免纪元时间转 double 的低位损失。标准 TUM 八列不含速度，文件不能证明下游订阅者已收到消息；使用缓冲写入与周期、正常退出时 flush，强杀可能丢失尚未 flush 的尾段。故障采集另有分析用途，应继续保留。
- PGO 全局校正和定位器内部高频轨迹仅供显式调试时导出，现场默认关闭相应记录和无界内存累积。修改轨迹行为时同步检查使用说明，不要把默认现场轨迹改成内部位姿或降频轨迹。
