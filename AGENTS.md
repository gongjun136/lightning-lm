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
- For dependency fixes, validate from the pinned source in a fresh build/install
  directory so cached generated headers cannot conceal missing message fields.

The 2026-10-20 Orin failure was caused by `common_msgs` remaining pinned to
`3f918b3` although the consumer needed `NodeHeartbeat.active_fault` from `b7c005f`.
The local newer dependency hid this mistake. Keep the checks above in future work.
