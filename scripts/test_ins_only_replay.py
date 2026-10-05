#!/usr/bin/env python3
"""Replay a bag created by test_ins_only_ros.py --record; check paused-clock fail closure."""
import argparse
from pathlib import Path
import os
import signal
import subprocess
import tempfile
import time
import yaml
import rclpy
from rclpy.qos import qos_profile_sensor_data
from geometry_msgs.msg import PoseStamped
from lightning.msg import LocalizationStatus
from diagnostic_monitor_interfaces.msg import NodeHeartbeat


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', required=True)
    p.add_argument('--recorded-test', required=True)
    a = p.parse_args()
    recorded = Path(a.recorded_test).resolve()
    out = Path(tempfile.mkdtemp(prefix='replay_', dir=recorded))
    config = yaml.safe_load((recorded/'test.yaml').read_text())
    config['ins_only'].update(use_sim_time=True, audit_root=str(out))
    cfg = out/'replay.yaml'
    cfg.write_text(yaml.safe_dump(config))
    os.environ['ROS_DOMAIN_ID'], os.environ['ROS_LOCALHOST_ONLY'] = '88', '1'
    rclpy.init()
    n = rclpy.create_node('ins_replay_observer')
    poses, states, heartbeats = [], [], []
    subs = [n.create_subscription(PoseStamped,'/slamPoseRaw_topic',poses.append,qos_profile_sensor_data),
            n.create_subscription(LocalizationStatus,'/localization/loc_status',states.append,qos_profile_sensor_data),
            n.create_subscription(NodeHeartbeat,'/diagnostics/heartbeat/lightning_slam',heartbeats.append,qos_profile_sensor_data)]
    node_log, bag_log = (out/'node.log').open('w'), (out/'player.log').open('w')
    node = subprocess.Popen([a.binary,'--config',str(cfg)],stdout=node_log,stderr=subprocess.STDOUT)
    player = None
    def spin(duration):
        end = time.monotonic()+duration
        while time.monotonic()<end:
            rclpy.spin_once(n,timeout_sec=.01)
            if node.poll() is not None:
                raise AssertionError('node exited; inspect '+str(out/'node.log'))
    try:
        spin(1)
        player = subprocess.Popen(['bash',str(Path(__file__).with_name('replay_ins_only.sh')),str(cfg),str(recorded/'input_bag'),
                                   '--delay','2','--rate','1'],stdout=bag_log,stderr=subprocess.STDOUT)
        deadline = time.monotonic()+45
        while player.poll() is None:
            if time.monotonic()>deadline:
                raise AssertionError('replay did not finish')
            spin(.05)
        assert player.returncode == 0, 'replay_ins_only.sh failed'
        spin(.7)
        assert len(poses)>10, 'no real recomputed replay poses'
        stamps = [m.header.stamp.sec+m.header.stamp.nanosec*1e-9 for m in poses]
        audit = list(out.glob('ins_only_*/published_rear_axle.tum'))
        assert len(audit)==1
        published = [float(line.split()[0]) for line in audit[0].read_text().splitlines()]
        assert len(published)>10 and all(any(abs(t-p)<1e-6 for p in published) for t in stamps), 'observed poses were not produced by the replay node'
        assert all(y>x for x,y in zip(stamps,stamps[1:])), 'replayed old business messages or nonmonotonic output'
        assert all(abs(m.pose.position.x)+abs(m.pose.position.y)+abs(m.pose.position.z)<1e-6 for m in poses)
        assert any(s.status==LocalizationStatus.STATUS_GOOD for s in states)
        assert states[-1].status==LocalizationStatus.STATUS_FAIL, 'frozen /clock held GOOD'
        count, beats, work = len(poses), len(heartbeats), heartbeats[-1].work_seq
        spin(.4)
        assert len(poses)==count and len(heartbeats)>beats and heartbeats[-1].work_seq==work
        report=dict(passed=True, poses=len(poses), check='input-only clocked replay, monotonic recomputed poses, wall-time timeout after /clock freezes')
        (out/'result.yaml').write_text(yaml.safe_dump(report))
        print('PASS',out,report,sep='\n')
    finally:
        for process in [player,node]:
            if process and process.poll() is None:
                process.send_signal(signal.SIGINT)
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.terminate(); process.wait(timeout=5)
        node_log.close();bag_log.close()
        n.destroy_node();rclpy.shutdown()


if __name__=='__main__':
    main()
