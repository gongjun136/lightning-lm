# BA + BTC + HBA backend

This backend keeps Lightning-LM's Livox/IMU LIO frontend and replaces the
odometry-distance loop trigger with a descriptor-first mapping pipeline:

1. plane-voxel local bundle adjustment runs on every sliding keyframe window;
2. BTC builds a descriptor from a non-overlapping local submap and searches the
   complete descriptor database without an odometry-distance retrieval gate;
3. a matched loop is checked by plane registration, observability, drift,
   temporal confirmation, and a post-registration safety gate;
4. accepted constraints enter a robust pose graph;
5. HBA runs in the background after a graph-changing loop, on manual request,
   and at data end. Its bottom-up levels form submaps and its global top layer
   jointly optimizes spatially overlapping coarse submaps before corrections
   are propagated to keyframes.

RTK is never read by the SLAM executable. It is used only by the evaluation
script after a run has finished.

## Configuration

The M3DGR reference configuration is
`config/reproduction/single_lidar/m3dgr/lightning_m3dgr_mid360_benchmark.yaml`.
The main switch is:

```yaml
backend:
  mode: ba_btc_hba  # ba_btc_hba | legacy | disabled
```

`backend.local_ba.enabled`, `backend.btc.enabled`, and `backend.hba.enabled`
provide module-level ablations. `backend.hba.require_applied_loop_for_commit`
defaults to `true`: end/manual requests still execute HBA, but their point-level
updates are rolled back unless a validated BTC constraint changed the pose
graph. This prevents a forced end-of-data pass from degrading an already
consistent LIO trajectory. Set it to `false` for an HBA-only experiment.

Generate self-contained ablation files with:

```bash
python3 scripts/reproduction/single_lidar/m3dgr/generate_backend_ablation_configs.py \
  --base config/reproduction/single_lidar/m3dgr/lightning_m3dgr_mid360_benchmark.yaml \
  --output-dir /tmp/lightning_m3dgr_ablations
```

## Offline and online operation

The normal offline runner accepts any generated YAML:

```bash
bash scripts/run_slam_offline.sh \
  --bag /path/to/ros2_bag \
  --config /path/to/lightning_m3dgr_ba_btc_hba.yaml \
  --output-dir /path/to/run
```

It writes raw and optimized frame trajectories, keyframe trajectories, the
optimized map, `backend_summary.yaml`, and `btc_loop_candidates.csv`. Online
SLAM publishes the optimized `map -> odom` transform without rewriting the
high-rate odometry stream. A global pass can be requested with:

```bash
ros2 service call /lightning/optimize_backend std_srvs/srv/Trigger '{}'
```

Map save and shutdown force a final HBA request and wait for the background
optimizer to finish.

## M3DGR evaluation protocol

Only `/livox/mid360/lidar` and `/livox/mid360/imu` are consumed. Outdoor RTK is
used as ground truth with fixed-scale SE(3) alignment; Sim(3) alignment is not
allowed. A truth loop has at least 30 s temporal separation and at most 5 m RTK
distance. Both values are configurable evaluation arguments.

```bash
python3 scripts/reproduction/single_lidar/m3dgr/evaluate_backend_runs.py \
  --ground-truth /path/to/rtk.txt \
  --run legacy=/path/to/legacy/results/trajectory_slam_opt.tum \
  --run ba_btc_hba=/path/to/new/results/trajectory_slam_opt.tum \
  --output /path/to/report.json
```

The report contains translation ATE, distance-based RPE, loop precision/recall,
CPU time, and peak RSS. M3DGR outdoor quaternion fields are identity
placeholders, so the script explicitly marks orientation error unavailable.

## BTC provenance and license

The descriptor implementation under
`src/core/backend/third_party/voxel_slam_btc` is adapted from the BTC module in
the referenced Voxel-SLAM repository. It remains licensed under GPL-2.0; the
corresponding license text is included in that directory. Distributing a binary
that contains this source must comply with GPL-2.0 obligations.
