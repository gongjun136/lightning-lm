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
