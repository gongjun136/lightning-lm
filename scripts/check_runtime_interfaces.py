#!/usr/bin/env python3
"""Read-only deployed interface/health check; use the production ROS/DDS environment."""
import argparse
import json
from pathlib import Path
import time
from collections import Counter
import yaml
import rclpy
from rclpy.qos import QoSProfile, ReliabilityPolicy
from rosidl_runtime_py.utilities import get_message
from rosidl_runtime_py.convert import message_to_ordereddict

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--workspace', type=Path, required=True)
parser.add_argument('--duration', type=float, default=25)
parser.add_argument('--require-healthy', action='store_true')
args = parser.parse_args()
if args.duration < 10:
    parser.error('--duration must be at least 10 seconds')
params = yaml.safe_load((args.workspace / 'install/livox_ros_driver2_core/share/livox_ros_driver2_core/config/dtc_aggregator.yaml').read_text())['lidar_dtc_aggregator']['ros__parameters']
radar_topic = params['internal_output_topic']
rclpy.init()
node = rclpy.create_node('safety_repair_probe')
topics = {
 '/safety/throstle': 'geosun_msgs/msg/Safety',
 radar_topic: 'diagnostic_monitor_interfaces/msg/DtcObservationArray',
 '/localization/pose_vel': 'lightning/msg/VehiclePose',
}
for role in ['front', 'rear', 'left', 'right']:
    topics[f'/sensors/livox/{role}/diagnostic_state'] = 'livox_ros_driver2/msg/LidarDiagnosticState'
counts, last, states, gaps, previous = Counter(), {}, Counter(), {}, {}
bounded_counts, sync_fault_counts, offset_ranges = Counter(), Counter(), {}
errors = set()
roles = ['front', 'rear', 'left', 'right']
def receive(topic, msg):
    now = time.monotonic()
    if topic in previous:
        gaps[topic] = max(gaps.get(topic, 0), now - previous[topic])
    previous[topic] = now
    counts[topic] += 1
    last[topic] = message_to_ordereddict(msg)
    if topic == '/safety/throstle':
        states[msg.runstate] += 1
        if args.require_healthy and now - start >= 5 and msg.runstate != 0:
            errors.add(f'Safety.runstate={msg.runstate} after probe warmup')
    if topic.endswith('/diagnostic_state'):
        role = roles.index(topic.split('/')[3])
        if msg.config_version != params['config_version']:
            errors.add(f'{topic}: config_version={msg.config_version}, expected={params["config_version"]}')
        if msg.role != role or msg.comm_header.source_id != params['expected_source_ids'][role]:
            errors.add(f'{topic}: role/source_id mismatch')
        if msg.bounded_time_sync_warning and (
                abs(msg.time_sync_offset_ns) > 1000000 or msg.time_sync_evidence_age_us > 2500000):
            errors.add(f'{topic}: invalid bounded time-sync evidence')
        bounded_counts[roles[role]] += int(msg.bounded_time_sync_warning)
        if len(msg.observations) != 11:
            errors.add(f'{topic}: expected 11 source observations')
        else:
            sync_fault_counts[roles[role]] += int(msg.observations[1].physical_condition == 2)
        bounds = offset_ranges.setdefault(roles[role], [msg.time_sync_offset_ns, msg.time_sync_offset_ns])
        bounds[0] = min(bounds[0], msg.time_sync_offset_ns)
        bounds[1] = max(bounds[1], msg.time_sync_offset_ns)
        if args.require_healthy and now - start >= 5 and msg.usability != 1:
            errors.add(f'{topic}: usability={msg.usability}')
    if topic == radar_topic:
        if msg.comm_header.source_id != params['functional_safety.node_id'] or msg.interface_version != 1 or len(msg.observations) != 46:
            errors.add('aggregator source/version/observation-count mismatch')
        if msg.config_version != params['config_version']:
            errors.add('aggregator config_version mismatch')
subs = [node.create_subscription(get_message(kind), topic,
        lambda msg, topic=topic: receive(topic, msg),
        QoSProfile(depth=20, reliability=ReliabilityPolicy.BEST_EFFORT))
        for topic, kind in topics.items()]
start = time.monotonic()
while time.monotonic() - start < args.duration:
    rclpy.spin_once(node, timeout_sec=0.1)
for topic in topics:
    if counts[topic] < 2:
        errors.add(f'{topic}: fewer than two messages received; check publisher, QoS and ROS/DDS environment')
    elif time.monotonic() - previous[topic] > 1:
        errors.add(f'{topic}: no recent message')
result = dict(duration=time.monotonic()-start, counts=dict(counts),
              bounded_warning_counts=dict(bounded_counts), sync_fault_counts=dict(sync_fault_counts),
              reported_offset_range_ns=offset_ranges,
              safety_states=dict(states), max_gap_sec=gaps, errors=sorted(errors), last=last)
print(json.dumps(result, indent=2), flush=True)
node.destroy_node()
rclpy.shutdown()
raise SystemExit(1 if errors else 0)
