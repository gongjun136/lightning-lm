# Lightning-LM workspace layout

`lightning_lm` depends on two ROS 2 interface packages maintained in the
separate `common_msgs` repository:

- `diagnostic_monitor_interfaces`
- `lightning` (source directory: `lightning_interfaces`)

The supported source layout is:

```text
lightning_lm_ws/
└── src/
    ├── lightning-lm/
    └── common_msgs/
```

Import the public message repository from the workspace root with:

```bash
vcs import src < src/lightning-lm/common_msgs.repos
```

The company checkout uses these remotes for `common_msgs`:

```text
origin   git@github.com:gongjun136/common_msgs.git
company  git@10.233.88.16:common/message.git
```

## Local compatibility mapping

On the maintained Windows workspace, the physical historical-data checkout
remains at:

```text
F:\SLAM_AI_KnowledgeBase\code\WSL_Ubuntu_22.04\lightning-lm
```

The workspace entry below is a directory junction to that existing checkout:

```text
F:\SLAM_AI_KnowledgeBase\code\WSL_Ubuntu_22.04\lightning_lm_ws\src\lightning-lm
```

Consequently, historical `data/`, `runs/`, `.tmp_sany_analysis/`, build,
install, and log paths remain available through the old path and are also
visible through the workspace path. Do not delete the old path while the
junction is in use.
