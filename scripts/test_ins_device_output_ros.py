#!/usr/bin/env python3
"""Real-node regression: unknown device output must never become rear-axle business data."""
import argparse
from collections import defaultdict
import math
import os
from pathlib import Path
import signal
import subprocess
import tempfile
import time

import rclpy
from rclpy.qos import qos_profile_sensor_data
from geometry_msgs.msg import PoseStamped
from lightning.msg import LocalizationStatus
from std_msgs.msg import String
from tf2_msgs.msg import TFMessage
from cgi430_interfaces import msg as ci
import yaml


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', required=True)
    parser.add_argument('--output-root', required=True)
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    root = Path(args.output_root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    out = Path(tempfile.mkdtemp(prefix='ins_device_output_', dir=root))
    repo = Path(__file__).resolve().parents[1]
    config = yaml.safe_load((repo/'config/ins_only/sany_cgi430.yaml').read_text())
    ins = config['ins_only']
    ins['georeference'].update(confirmed=True, origin_llh=[0, 0, 0], input_height='ellipsoidal')
    ins['quality'].update(recovery_samples=3, max_age_sec=.3, max_skew_sec=.08)
    ins['audit_root'] = str(out)
    config['system']['pub_tf'] = True
    cfg = out/'site.yaml'
    cfg.write_text(yaml.safe_dump(config))
    os.environ.update(ROS_DOMAIN_ID='89', ROS_LOCALHOST_ONLY='1')
    rclpy.init()
    node = rclpy.create_node('ins_device_output_regression')
    data = defaultdict(list)
    subscriptions = [node.create_subscription(typ, topic, lambda m, t=topic: data[t].append(m),
                                              qos_profile_sensor_data)
                     for typ, topic in [(PoseStamped, '/localization/ins_reference/pose'),
                                        (LocalizationStatus, '/localization/ins_reference/status'),
                                        (String, '/localization/ins_reference/diagnostics'), (TFMessage, '/tf')]]
    definitions = [
        (ci.Latitude, 'position/latitude', dict(latitude_deg=0.)),
        (ci.Longitude, 'position/longitude', dict(longitude_deg=0.)),
        (ci.Altitude, 'position/altitude', dict(altitude_m=3.25)),
        (ci.Attitude, 'attitude', dict(heading_deg=0., pitch_deg=0., roll_deg=0.)),
        (ci.EarthVelocity, 'velocity', dict(east_mps=1., north_mps=0., up_mps=0.)),
        (ci.PositionSigma, 'position/sigma', dict(east_m=.01, north_m=.01, up_m=.02)),
        (ci.AttitudeSigma, 'attitude/sigma', dict(heading_deg=.1, pitch_deg=.1, roll_deg=.1)),
        (ci.VelocitySigma, 'velocity/sigma', dict(east_mps=.01, north_mps=.01, up_mps=.01)),
        (ci.InsStatus, 'ins/status', dict(system_state=2, satellite_status=4, differential_age_s=.1)),
    ]
    pubs = [(node.create_publisher(typ, '/cgi430/'+name, qos_profile_sensor_data), typ, values)
            for typ, name, values in definitions]
    log = (out/'node.log').open('w')
    process = subprocess.Popen([args.binary, '--config', str(cfg)], stdout=log, stderr=subprocess.STDOUT)
    recorder = record_log = None

    def spin(seconds):
        end = time.monotonic()+seconds
        while time.monotonic() < end:
            rclpy.spin_once(node, timeout_sec=.002)
            assert process.poll() is None, 'node exited; inspect '+str(out/'node.log')

    def cycle(status=4):
        stamp = node.get_clock().now().to_msg()
        for pub, typ, values in pubs[-1:]+pubs[:-1]:
            message = typ(**values)
            message.frame.header.stamp = stamp
            message.frame.valid, message.frame.dlc = True, 8
            if isinstance(message, ci.InsStatus):
                message.satellite_status = status
            pub.publish(message)
        spin(.012)

    forbidden = ['/slamPoseRaw_topic', '/localization/pose_vel', '/localization/loc_status',
                 '/localization/fault_status', '/localization/path', '/localization/ins_diagnostics',
                 '/diagnostics/heartbeat/lightning_slam', '/LidarDataInv', '/LidarDataInL']
    try:
        deadline = time.monotonic()+30
        while not all(pub.get_subscription_count() for pub, *_ in pubs):
            assert time.monotonic() < deadline, 'CGI discovery timeout'
            spin(.05)
        spin(.5)
        if args.record:
            record_log = (out/'recorder.log').open('w')
            recorder = subprocess.Popen(['bash', str(repo/'scripts/record_ins_only.sh'), str(cfg), str(out/'bag')],
                                        stdout=record_log, stderr=subprocess.STDOUT)
            deadline = time.monotonic()+30
            while node.count_subscribers('/localization/ins_reference/pose') < 2:
                assert recorder.poll() is None and time.monotonic() < deadline, 'recorder discovery timeout'
                cycle()
            # Count discovery is only readiness; the finalized bag must still
            # contain received poses below. Keep valid navigation flowing.
            for _ in range(150):
                cycle()
        for _ in range(100):
            cycle()
        poses = data['/localization/ins_reference/pose']
        assert len(poses) >= 10, 'reference pose not received'
        pose = poses[-1]
        assert pose.header.frame_id == 'ins_enu'
        assert abs(pose.pose.position.x) < 1e-7 and abs(pose.pose.position.y) < 1e-7
        assert abs(pose.pose.position.z-3.25) < 1e-7, 'device output point was shifted'
        q = pose.pose.orientation
        assert abs(q.x)+abs(q.y) < 1e-7 and abs(q.z-math.sqrt(.5)) < 1e-7 and abs(q.w-math.sqrt(.5)) < 1e-7
        assert all(node.count_publishers(t) == 0 for t in forbidden), 'reference mode created a business publisher'
        transforms = [t for message in data['/tf'] for t in message.transforms]
        assert transforms and all(t.child_frame_id == 'cgi_output_flu' for t in transforms)
        assert any('attitude_reference=device_output' in d.data and 'output_reference_confirmed=0' in d.data
                   for d in data['/localization/ins_reference/diagnostics'])
        for _ in range(10):
            cycle(status=5)
        count = len(poses)
        for _ in range(10):
            cycle(status=5)
        assert len(poses) == count, 'invalid RTK reference published'
        for _ in range(30):
            cycle()
        assert len(poses) > count, 'reference failed to recover'
        spin(.7)
        assert data['/localization/ins_reference/status'][-1].status == LocalizationStatus.STATUS_FAIL
    finally:
        process.send_signal(signal.SIGINT) if process.poll() is None else None
        process.wait(timeout=15)
        log.close()
        if recorder is not None:
            recorder.send_signal(signal.SIGINT) if recorder.poll() is None else None
            recorder.wait(timeout=20)
            record_log.close()
            assert recorder.returncode == 0, 'reference recording failed'
        node.destroy_node()
        rclpy.shutdown()
    audits = list(out.glob('ins_only_*'))
    assert len(audits) == 1
    metadata = yaml.safe_load((audits[0]/'reference_metadata.yaml').read_text())
    assert metadata['output_reference'] == 'device_output' and metadata['output_reference_confirmed'] is False
    assert metadata['output_point_id'] == 'unverified_cgi_output' and metadata['lever_arm_applied'] is False
    assert metadata['attitude_reference'] == 'device_output' and metadata['output_frame'] == 'cgi_output_flu'
    assert (audits[0]/'published_device_output.tum').stat().st_size > 0
    assert not (audits[0]/'published_rear_axle.tum').exists()
    if args.record:
        bag = yaml.safe_load((out/'bag/metadata.yaml').read_text())['rosbag2_bagfile_information']
        counts = {t['topic_metadata']['name']: t['message_count'] for t in bag['topics_with_message_count']}
        assert counts.get('/localization/ins_reference/pose', 0) > 0
        assert not any(t in counts for t in forbidden), 'recorder selected business outputs'
    for name in ['unconfirmed_rear', 'unknown_reference', 'raw_cloud', 'raw_rear_frame', 'unknown_attitude', 'confirmed_unknown_point']:
        bad = yaml.safe_load(cfg.read_text())
        c = bad['ins_only']
        if name == 'unconfirmed_rear':
            c.update(output_reference='rear_axle', attitude_reference='rear_body',
                     output_point_id='rear_axle_center', output_frame='rear_axle')
        elif name == 'unknown_reference':
            c['output_reference'] = 'unknown'
        elif name == 'raw_cloud':
            bad['cloud']['enabled'] = True
        elif name == 'raw_rear_frame':
            c['output_frame'] = 'rear_axle'
        elif name == 'unknown_attitude':
            c['attitude_reference'] = 'unknown'
        else:
            c['output_reference_confirmed'] = True
        path = out/(name+'.yaml')
        path.write_text(yaml.safe_dump(bad))
        rejected = subprocess.run([args.binary, '--config', str(path)], capture_output=True, text=True, timeout=20)
        (out/(name+'.log')).write_text(rejected.stdout+rejected.stderr)
        assert rejected.returncode in [1, 2], name+' did not refuse startup'
    (out/'result.yaml').write_text(yaml.safe_dump(dict(passed=True, recorded=args.record,
        checks=['unknown physical point retained', 'device attitude preserved', 'truthful TF and metadata',
                'no business publishers or Safety heartbeat', 'quality loss and recovery', 'invalid contracts rejected'])))
    print('PASS', out)


if __name__ == '__main__':
    main()
