#!/usr/bin/env python3
"""Synthetic end-to-end ROS test; run in an isolated ROS_DOMAIN_ID.

Starts the real run_loc_online binary, no map/IMU/LIO or CGI hardware required.
Artifacts are retained in a new test directory for diagnosis.
"""
import argparse
from collections import defaultdict
from pathlib import Path
import os
import signal
import struct
import subprocess
import tempfile
import time

import rclpy
from rclpy.qos import qos_profile_sensor_data
from rclpy.utilities import get_rmw_implementation_identifier
from sensor_msgs.msg import PointCloud2, PointField
from geometry_msgs.msg import PoseStamped
from geosun_msgs.msg import PosRes
from lightning.msg import VehiclePose, LocalizationStatus, FaultStatus
from diagnostic_monitor_interfaces.msg import NodeHeartbeat
from cgi430_interfaces import msg as ci
import yaml


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', required=True)
    p.add_argument('--output-root', default='runs')
    p.add_argument('--record', action='store_true', help='also exercise record_ins_only.sh and retain an input bag')
    args = p.parse_args()
    root = Path(args.output_root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    out = Path(tempfile.mkdtemp(prefix='ins_ros_test_', dir=root))
    config = yaml.safe_load((Path(__file__).resolve().parents[1] / 'config/ins_only/sany_cgi430.yaml').read_text())
    ins = config['ins_only']
    ins['georeference'].update(confirmed=True, origin_llh=[0, 0, 0], input_height='ellipsoidal')
    ins['audit_root'] = str(out)
    ins['quality'].update(recovery_samples=3, max_age_sec=.3, max_skew_sec=.08)
    config['cloud'].pop('lidar_config')
    config['cloud'].update(voxel_leaf_m=0, rear_from_primary=dict(confirmed=True, translation_m=[1, 0, 0], quaternion_xyzw=[0, 0, 0, 1]))
    config['multi_lidar'] = dict(enabled=True, primary_lidar_id=0, frame_period=.1,
        match_tolerance=.002, reorder_window=.1, min_lidars=2,
        topics=dict(lidar_0='/test/lidar0', lidar_1='/test/lidar1'),
        extrinsics=dict(lidar0=dict(T=[1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]),
                        lidar1=dict(T=[1,0,0,0, 0,1,0,1, 0,0,1,0, 0,0,0,1])))
    config['fasterlio'] = dict(lidar_type=1, blind=.1, point_filter_num=1, scan_line=4, time_scale=.001, livox_point_time_scale=1)
    config['self_point_filter'] = dict(enabled=False)
    config['roi'] = dict(height_min=-10, height_max=10)
    cfg = out / 'test.yaml'
    cfg.write_text(yaml.safe_dump(config))
    os.environ['ROS_DOMAIN_ID'] = '88'
    os.environ['ROS_LOCALHOST_ONLY'] = '1'
    rclpy.init()
    middleware = get_rmw_implementation_identifier()
    node = rclpy.create_node('ins_only_test_driver')
    data = defaultdict(list)
    subscriptions = []
    for typ, topic in [(PosRes, '/PosRes'), (VehiclePose, '/localization/pose_vel'),
                       (PoseStamped, '/slamPoseRaw_topic'), (LocalizationStatus, '/localization/loc_status'),
                       (FaultStatus, '/localization/fault_status'), (NodeHeartbeat, '/diagnostics/heartbeat/lightning_slam'),
                       (PointCloud2, '/LidarDataInv'), (PointCloud2, '/LidarDataInL')]:
        subscriptions.append(node.create_subscription(typ, topic, lambda m, t=topic: data[t].append(m), qos_profile_sensor_data))
    definitions = [
        (ci.Latitude, 'position/latitude', dict(latitude_deg=0.)),
        (ci.Longitude, 'position/longitude', dict(longitude_deg=0.)),
        (ci.Altitude, 'position/altitude', dict(altitude_m=0.)),
        (ci.Attitude, 'attitude', dict(heading_deg=90., pitch_deg=0., roll_deg=0.)),
        (ci.EarthVelocity, 'velocity', dict(east_mps=-2., north_mps=0., up_mps=0.)),
        (ci.PositionSigma, 'position/sigma', dict(east_m=.01, north_m=.01, up_m=.02)),
        (ci.AttitudeSigma, 'attitude/sigma', dict(heading_deg=.1, pitch_deg=.1, roll_deg=.1)),
        (ci.VelocitySigma, 'velocity/sigma', dict(east_mps=.01, north_mps=.01, up_mps=.01)),
        (ci.InsStatus, 'ins/status', dict(system_state=2, satellite_status=4, differential_age_s=.1)),
    ]
    pubs = [(node.create_publisher(typ, '/cgi430/'+name, qos_profile_sensor_data), typ, name, values)
            for typ, name, values in definitions]
    lidar_pubs = [node.create_publisher(PointCloud2, '/test/lidar'+str(i), qos_profile_sensor_data) for i in range(2)]

    def spin(duration):
        end = time.monotonic()+duration
        while time.monotonic() < end:
            rclpy.spin_once(node, timeout_sec=min(.002, max(0, end-time.monotonic())))
            if process.poll() is not None:
                raise AssertionError('localization exited; inspect '+str(out/'node.log'))

    def stamp(sec):
        from builtin_interfaces.msg import Time
        ns = round(sec*1e9)
        return Time(sec=ns//10**9, nanosec=ns%10**9)

    def cycle(status=4, skip=None, cloud=False):
        now = node.get_clock().now().nanoseconds*1e-9-.01
        # Put status first so the invalid state is observed before new motion fields.
        for pub, typ, name, values in pubs[-1:]+pubs[:-1]:
            if name == skip:
                continue
            m = typ(**values)
            m.frame.header.stamp = stamp(now)
            m.frame.valid, m.frame.dlc = True, 8
            if name == 'ins/status':
                m.satellite_status = status
            pub.publish(m)
        if cloud:
            for i, pub in enumerate(lidar_pubs):
                m = PointCloud2()
                m.header.stamp = stamp(now-.04)
                m.header.frame_id = 'test_lidar'+str(i)
                m.height, m.width, m.point_step, m.row_step = 1, 3, 26, 78
                m.fields = [PointField(name=n, offset=o, datatype=d, count=1) for n,o,d in
                            [('x',0,7),('y',4,7),('z',8,7),('intensity',12,7),('tag',16,2),('line',17,2),('timestamp',18,8)]]
                m.data = b''.join(struct.pack('<ffffBBd', 10., 2.-i, 0., 30., 0, 0, dt) for dt in [0.,.01,.02])
                pub.publish(m)
        spin(.01)

    report = []
    recorder, record_log = None, None
    if args.record:
        record_log = (out/'recorder.log').open('w')
        recorder = subprocess.Popen(['bash',str(Path(__file__).with_name('record_ins_only.sh')),str(cfg),str(out/'input_bag')],
                                    stdout=record_log,stderr=subprocess.STDOUT)
    log = (out/'node.log').open('w')
    process = subprocess.Popen([args.binary, '--config', str(cfg)], stdout=log, stderr=subprocess.STDOUT)
    try:
        deadline = time.monotonic()+30
        while not all(pub.get_subscription_count() for pub, *_ in pubs):
            if time.monotonic()>deadline:
                raise AssertionError('CGI subscription discovery timeout')
            spin(.05)
        spin(1)
        for i in range(300):
            cycle(cloud=i > 10 and i % 10 == 0)
            if i>=50 and len(data['/PosRes'])>=30 and data['/LidarDataInv'] and data['/LidarDataInL']:
                break
        (out/'observed.yaml').write_text(yaml.safe_dump({k:len(v) for k,v in data.items()}))
        assert len(data['/PosRes']) >= 10, 'no accepted navigation'
        pose = data['/localization/pose_vel'][-1]
        assert max(abs(pose.x), abs(pose.y), abs(pose.z), abs(pose.yaw)) < 1e-6
        assert abs(pose.speed+2) < 1e-6, 'reverse speed sign'
        assert data['/localization/loc_status'][-1].status == LocalizationStatus.STATUS_GOOD
        for topic, frame in [('/LidarDataInv','rear_axle'),('/LidarDataInL','map')]:
            assert data[topic], 'cloud output missing: '+topic
            cloud = data[topic][-1]
            assert cloud.header.frame_id == frame and cloud.width == 6
            for j in range(cloud.width):
                x,y,z = struct.unpack_from('<fff', cloud.data, j*cloud.point_step)
                assert abs(x-11)<1e-5 and abs(y-2)<1e-5 and abs(z)<1e-5, 'cloud/extrinsic mismatch'
        report.append('valid rear-axle pose, reverse speed, two-lidar body/map cloud outputs')
        for bad in [5, 8]:
            for _ in range(12):
                cycle(status=bad)
            n = len(data['/PosRes'])
            clouds = (len(data['/LidarDataInv']), len(data['/LidarDataInL']))
            for _ in range(12):
                cycle(status=bad, cloud=True)
            assert len(data['/PosRes']) == n, f'published in invalid status {bad}'
            assert clouds == (len(data['/LidarDataInv']), len(data['/LidarDataInL'])), 'published cloud across quality loss'
            assert data['/localization/loc_status'][-1].status == LocalizationStatus.STATUS_FAIL
            assert data['/localization/fault_status'][-1].level == FaultStatus.LEVEL_P0
            for _ in range(2):
                cycle()
            assert len(data['/PosRes']) == n, 'recovered before three new valid solutions'
            for _ in range(23):
                cycle()
            assert len(data['/PosRes']) > n, 'failed to recover'
        report.append('RTK float and fixed-without-heading stop publication; recovery works')
        for _ in range(45):
            cycle(skip='attitude')
        n = len(data['/PosRes'])
        for _ in range(10):
            cycle(skip='attitude')
        assert len(data['/PosRes']) == n, 'stale attitude published'
        for _ in range(30):
            cycle()
        assert len(data['/PosRes']) > n
        spin(.45)
        n, heartbeats = len(data['/PosRes']), len(data['/diagnostics/heartbeat/lightning_slam'])
        work = data['/diagnostics/heartbeat/lightning_slam'][-1].work_seq
        spin(.3)
        assert len(data['/PosRes']) == n
        assert len(data['/diagnostics/heartbeat/lightning_slam']) > heartbeats
        assert data['/diagnostics/heartbeat/lightning_slam'][-1].work_seq == work
        assert data['/localization/loc_status'][-1].status == LocalizationStatus.STATUS_FAIL
        report.append('stale attitude and complete dropout stop output; heartbeat stays alive without advancing work')
    finally:
        if process.poll() is None:
            process.send_signal(signal.SIGINT)
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                process.terminate()
                process.wait(timeout=5)
        log.close()
        if recorder is not None:
            if recorder.poll() is None:
                recorder.send_signal(signal.SIGINT)
                recorder.wait(timeout=15)
            record_log.close()
            assert recorder.returncode == 0, 'record_ins_only.sh failed; inspect recorder.log'
        node.destroy_node()
        rclpy.shutdown()
    for name in ['unconfirmed_origin', 'unconfirmed_extrinsic', 'legacy_transform', 'unknown_mode']:
        bad = yaml.safe_load(cfg.read_text())
        if name == 'unconfirmed_origin':
            bad['ins_only']['georeference']['confirmed'] = False
        elif name == 'unconfirmed_extrinsic':
            bad['cloud']['rear_from_primary']['confirmed'] = False
        elif name == 'legacy_transform':
            bad['output']['fixed_map_transform'] = dict(enabled=True)
        else:
            bad['system']['localization_mode'] = 'unknown'
        path = out/(name+'.yaml')
        path.write_text(yaml.safe_dump(bad))
        result = subprocess.run([args.binary,'--config',str(path)],capture_output=True,text=True,timeout=15)
        (out/(name+'.log')).write_text(result.stdout+result.stderr)
        assert result.returncode in [1,2], name+' should refuse startup'
    result = subprocess.run([args.binary,'--config',str(cfg),'--map','unused'],capture_output=True,text=True,timeout=15)
    assert result.returncode == 2, 'must reject misleading map flag'
    report.append('unconfirmed calibration, legacy transform, unknown mode and map flag refuse startup')
    (out/'result.yaml').write_text(yaml.safe_dump(dict(passed=True, middleware=middleware, ros_domain_id=88,
                                                    input_bag=str(out/'input_bag') if args.record else None, checks=report)))
    print('PASS', out, *report, sep='\n')


if __name__ == '__main__':
    main()
