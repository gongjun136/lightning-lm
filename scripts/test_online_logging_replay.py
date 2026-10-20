#!/usr/bin/env python3
"""Run from sourced ROS Humble environment; never connect to vehicle DDS."""
import argparse
import json
import os
import subprocess
import threading
import time
import re
from pathlib import Path
os.environ['ROS_DOMAIN_ID'] = '173'
os.environ['RMW_IMPLEMENTATION'] = 'rmw_cyclonedds_cpp'
os.environ['ROS_LOCALHOST_ONLY'] = '1'
os.environ['LIGHTNING_LM_COMPUTE_PROFILE'] = '1'
os.environ['LIGHTNING_LM_REDUCE_NONESSENTIAL_OVERHEAD'] = '0'
os.environ['SANY_ENABLE_CAN_OBSERVATION'] = '0'
import rclpy
from rclpy.qos import qos_profile_sensor_data
from lightning.msg import VehiclePose
parser = argparse.ArgumentParser(description=__doc__)
for key in ('binary', 'config', 'map', 'bag', 'output'):
    parser.add_argument('--' + key, required=True)
parser.add_argument('--slow', action='store_true')
options = parser.parse_args()
out = Path(options.output).resolve()
out.mkdir(parents=True, exist_ok=True)
slow = options.slow
name = 'slow' if slow else 'normal'
args = [str(Path(options.binary).resolve()), '--config=' + str(Path(options.config).resolve()), '--map=' + str(Path(options.map).resolve()), '--bag=' + str(Path(options.bag).resolve()), '--post_wait_seconds=2']
rclpy.init()
node = rclpy.create_node('async_logging_validation')
poses = []

def receive(msg):
    poses.append((time.monotonic(), msg.header.stamp.sec + msg.header.stamp.nanosec * 1e-09))
sub = node.create_subscription(VehiclePose, '/localization/pose_vel', receive, qos_profile_sensor_data)
proc = subprocess.Popen(args, cwd=out, stdout=open(out / (name + '.stdout'), 'wb'), stderr=subprocess.PIPE)
pause = threading.Event()
blocks = []
blocked_writers = []
deprecated_publishers = 0
graph = []

def drain():
    with (out / (name + '.stderr')).open('wb') as f:
        while True:
            if pause.is_set():
                blocks.append(('start', time.monotonic()))
                time.sleep(5)
                blocks.append(('end', time.monotonic()))
                pause.clear()
            buf = os.read(proc.stderr.fileno(), 4096)
            if not buf:
                break
            f.write(buf)
thread = threading.Thread(target=drain)
thread.start()
start = time.monotonic()
injected = False
while proc.poll() is None and time.monotonic() - start < 180:
    rclpy.spin_once(node, timeout_sec=0.02)
    if not graph and poses:
        graph = node.get_topic_names_and_types()
        deprecated_publishers = node.count_publishers('/PosRes')
    if blocks and len(blocks) == 1 and (not blocked_writers) and (time.monotonic() - blocks[0][1] > 1):
        for task in Path('/proc', str(proc.pid), 'task').glob('*'):
            try:
                wchan = (task / 'wchan').read_text().strip()
                if 'pipe_write' in wchan:
                    blocked_writers.append({'tid': int(task.name), 'wchan': wchan, 'observed_at': time.monotonic()})
            except OSError:
                pass
    if slow and poses and (not injected) and (time.monotonic() - poses[0][0] > 8):
        pause.set()
        injected = True
if proc.poll() is None:
    proc.terminate()
try:
    proc.wait(timeout=20)
except subprocess.TimeoutExpired:
    proc.kill()
    proc.wait()
thread.join(timeout=10)
node.destroy_node()
rclpy.shutdown()
gaps = [poses[i][0] - poses[i - 1][0] for i in range(1, len(poses))]
result = {'command': args, 'deprecated_publishers': deprecated_publishers, 'graph': graph, 'returncode': proc.returncode, 'poses': len(poses), 'injected': injected, 'blocked_writers': blocked_writers, 'blocks': blocks, 'max_wall_gap': max(gaps, default=None), 'sensor_time_monotonic': all((poses[i][1] > poses[i - 1][1] for i in range(1, len(poses))))}
if blocks:
    a, b = (blocks[0][1], blocks[-1][1])
    inside = [p for p in poses if a <= p[0] <= b]
    result.update(poses_during_block=len(inside), max_gap_during_block=max([poses[i][0] - poses[i - 1][0] for i in range(1, len(poses)) if poses[i][0] >= a and poses[i - 1][0] <= b], default=None))
text = (out / (name + '.stderr')).read_text(errors='replace')
seq = [int(x) for x in re.findall('async_seq=(\\d+)', text)]
result['log_sequence_ordered'] = seq == list(range(len(seq)))
result['summary'] = re.findall('ASYNC_LOG_SUMMARY.*', text)
(out / (name + '.json')).write_text(json.dumps(result, indent=2))
(out / (name + '_poses.json')).write_text(json.dumps(poses))
print(json.dumps(result, indent=2))
assert result['returncode'] == 0 and len(poses) > 100, 'online localization did not produce enough poses'
assert result['log_sequence_ordered'] and result['sensor_time_monotonic'] and (deprecated_publishers == 0)
if slow:
    assert blocked_writers and result.get('poses_during_block', 0) > 100 and (result['max_gap_during_block'] < 0.5)
