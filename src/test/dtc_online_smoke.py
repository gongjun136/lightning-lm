#!/usr/bin/env python3
"""Validate the real run_loc_online publishers with a temporary, empty map.

source install/setup.bash
RMW_IMPLEMENTATION=rmw_cyclonedds_cpp ROS_DOMAIN_ID=172 python3 \
  src/lightning-lm/src/test/dtc_online_smoke.py build/lightning_lm/bin/run_loc_online
"""
import argparse
import csv
import json
import os
from pathlib import Path
import statistics
import struct
import re
import shutil
import subprocess
import tempfile

import rclpy
import yaml
from rclpy.clock import Clock, ClockType
from rclpy.qos import QoSProfile, ReliabilityPolicy, DurabilityPolicy
from diagnostic_monitor_interfaces.msg import DtcObservationArray, NodeHeartbeat
from geosun_msgs.msg import Safety
from sensor_msgs.msg import Imu, PointCloud2, PointField
from rclpy.serialization import deserialize_message


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('--safety-fault-reporting', choices=('default', 'true', 'false'), default='default')
    parser.add_argument('--radar-time-sync-shadow', choices=('true', 'false'), default='true')
    parser.add_argument('--log-output', type=Path, help='Preserve process log after the test')
    parser.add_argument('--incident-output', type=Path, help='Enable and validate in-process MCAP capture in this empty directory')
    parser.add_argument('--startup-timeout-sec', type=float, default=15.0,
                        help='Process loading/discovery allowance; does not change DTC startup grace')
    parser.add_argument('--sany-lidars', type=int, choices=(3, 4), help='Test a resolved SANY deployment topology')
    args = parser.parse_args()
    if args.startup_timeout_sec <= 0:
        parser.error('--startup-timeout-sec must be positive')
    os.environ.setdefault('ROS_DOMAIN_ID', '172')
    os.environ['ROS_LOCALHOST_ONLY'] = '1'
    executable = args.executable.resolve()
    source_root = Path(__file__).resolve().parents[2]
    config = yaml.safe_load((source_root / 'config/default.yaml').read_text())
    config['system']['with_ui'] = False
    config['relocalization'] = {'enable_map_consistency': False, 'backend': 'btc'}
    config['map_export'] = {'normalize_start_ground_z': False}
    config['maps']['save_dyn_when_quit'] = False
    config['diagnostics'] = {'startup_grace_sec': 2.0, 'config_version': 'smoke-test'}
    config['diagnostics']['radar_time_sync_shadow_enabled'] = args.radar_time_sync_shadow == 'true'
    if args.sany_lidars:
        sany = yaml.safe_load((source_root / 'config/sany/localization.yaml').read_text())
        for key in ('config_version', 'consistency_unknown_grace_sec', 'major_sigma_limit_m',
                    'lateral_sigma_limit_m', 'yaw_sigma_limit_deg', 'covariance_variance_scale',
                    'position_sigma_floor_m', 'angle_sigma_floor_deg'):
            config['diagnostics'][key] = sany['diagnostics'][key]
        config['multi_lidar'] = sany['multi_lidar']
        config['multi_lidar']['topics'] = {
            key: config['multi_lidar']['topics'][key] for i in range(args.sany_lidars)
            for key in (f'lidar_{i}', f'imu_{i}')}
        config['multi_lidar']['extrinsics'] = {
            key: value for key, value in yaml.safe_load(
                (source_root / 'config/sany/calibrations.yaml').read_text())['vehicles']['SANLV77EXSH000858']['extrinsics'].items()
            if key in {f'lidar{i}' for i in range(args.sany_lidars)}}
        config['multi_lidar']['adaptive_load']['cloud_publish_min_lidars'] = args.sany_lidars
        config['vehicle_configuration'] = {
            'LIDAR': 'MASTER' if args.sany_lidars == 4 else 'MASTER AND SLAVE',
            'lidar_count': args.sany_lidars}
    if args.safety_fault_reporting != 'default':
        config['diagnostics']['safety_fault_reporting_enabled'] = args.safety_fault_reporting == 'true'
    if args.incident_output:
        args.incident_output = args.incident_output.resolve()
        if args.incident_output.exists() and any(args.incident_output.iterdir()):
            parser.error('--incident-output must be empty to avoid reusing old evidence')
        config['diagnostics']['incident_capture'] = {
            'enabled': True, 'pre_seconds': 1, 'post_seconds': 3, 'max_event_seconds': 5,
            'startup_timeout_seconds': 2, 'output_timeout_seconds': .3, 'cooldown_seconds': 0,
            'memory_mb': 16, 'max_message_mb': 4, 'event_mb': 4, 'disk_mb': 16,
            'min_free_mb': 0, 'segment_mb': 1, 'write_mb_per_second': 20,
        }
    with tempfile.TemporaryDirectory(prefix='dtc_online_', dir=executable.parent.parent) as directory:
        fixture = Path(directory)
        (fixture / 'index.txt').write_text('0 0 0\n')
        # Localization validates map identity before starting its publishers.
        # Build the temporary fixture through the same manifest exporter used
        # for real maps, including the geometry checksum and native frame.
        subprocess.run([str(executable.with_name('export_map_manifest')),
                        '--map=' + str(fixture), '--status=not_required',
                        '--site_id=dtc_smoke_fixture'], check=True, capture_output=True, text=True)
        path = fixture / 'config.yaml'
        path.write_text(yaml.safe_dump(config))
        rclpy.init()
        node = rclpy.create_node('localization_dtc_smoke_probe')
        steady_clock = Clock(clock_type=ClockType.STEADY_TIME)

        def steady_time():
            # Match the C++ wall timer's ROS steady clock; system clock slewing
            # can change Python time.monotonic() relative to this clock on WSL.
            return steady_clock.now().nanoseconds * 1e-9

        qos = QoSProfile(depth=100, reliability=ReliabilityPolicy.RELIABLE,
                         durability=DurabilityPolicy.VOLATILE)
        values, safety_messages, snapshots, heartbeats = [], [], [], []

        def observe_safety(msg):
            safety_messages.append(msg)
            values.append((steady_time(), msg.runstate))

        node.create_subscription(Safety, '/safety/throstle', observe_safety, qos)
        node.create_subscription(DtcObservationArray, '/diagnostics/source_dtc/localization', snapshots.append, qos)
        node.create_subscription(NodeHeartbeat, '/diagnostics/heartbeat/lightning_slam', heartbeats.append,
                                 QoSProfile(depth=100, reliability=ReliabilityPolicy.BEST_EFFORT))
        input_publishers = []
        inputs_sent = False
        if args.incident_output and args.sany_lidars:
            topics = config['multi_lidar']['topics']
            primary = config['multi_lidar']['primary_lidar_id']
            imu_publisher = node.create_publisher(Imu, topics[f'imu_{primary}'], 10)
            cloud_publishers = [node.create_publisher(PointCloud2, topics[f'lidar_{i}'], 10)
                                for i in range(args.sany_lidars)]
            input_publishers = [imu_publisher, *cloud_publishers]
        log_path = fixture / 'process.log'
        with log_path.open('w') as log:
            command = [str(executable), '--config=' + str(path), '--map=' + str(fixture)]
            if args.incident_output:
                command.append('--incident_dir=' + str(args.incident_output))
            process = subprocess.Popen(command,
                                       cwd=fixture, stdout=log, stderr=log)
            try:
                startup_deadline = steady_time() + args.startup_timeout_sec
                end = None
                while end is None or steady_time() < end:
                    rclpy.spin_once(node, timeout_sec=0.01)
                    if values and input_publishers and not inputs_sent and all(p.get_subscription_count() for p in input_publishers):
                        # One valid sample checks the real callback-to-MCAP path,
                        # without supplying enough IMU history to initialize LIO.
                        stamp = node.get_clock().now().to_msg()
                        imu = Imu()
                        imu.header.stamp = stamp
                        imu.linear_acceleration.z = 9.81
                        imu_publisher.publish(imu)
                        cloud = PointCloud2()
                        cloud.header.stamp = stamp
                        cloud.header.frame_id = 'incident_smoke'
                        cloud.width = cloud.height = 1
                        cloud.point_step = cloud.row_step = 32
                        cloud.fields = [PointField(name=name, offset=offset, datatype=kind, count=1)
                                        for name, offset, kind in (
                                            ('x', 0, PointField.FLOAT32), ('y', 4, PointField.FLOAT32),
                                            ('z', 8, PointField.FLOAT32), ('intensity', 12, PointField.FLOAT32),
                                            ('tag', 16, PointField.UINT8), ('line', 17, PointField.UINT8),
                                            ('timestamp', 24, PointField.FLOAT64))]
                        cloud.data = struct.pack('<ffffBB6xd', 5., 0., 0., 42., 0, 0, 0.)
                        for publisher in cloud_publishers:
                            publisher.publish(cloud)
                        inputs_sent = True
                    if process.poll() is not None:
                        raise RuntimeError(log_path.read_text()[-5000:])
                    if end is None:
                        if values:
                            end = steady_time() + 6.0
                        elif steady_time() >= startup_deadline:
                            raise RuntimeError('no Safety output within startup timeout:\n' + log_path.read_text()[-5000:])
                assert len(values) >= 150, f'too few safety samples: {len(values)}'
                if args.safety_fault_reporting == 'false':
                    assert all(value == 0 for _, value in values), 'disabled fault reporting must always give zero'
                else:
                    assert any(value == 4 for _, value in values), 'startup UNKNOWN must be 4'
                    assert values[-1][1] == -9, 'no valid output after startup grace must be -9'
                    assert all(value in (4, -9) for _, value in values), 'uninitialized localization must never give zero'
                publishers = node.get_publishers_info_by_topic('/safety/throstle')
                assert len(publishers) == 1, f'expected one safety publisher, found {len(publishers)}'
                assert publishers[0].topic_type == 'geosun_msgs/msg/Safety'
                assert all(msg.header.stamp.sec > 0 for msg in safety_messages), 'Safety must carry a publish timestamp'
                assert all((msg.timesync, msg.errorcode, msg.strategies, msg.reserved) == (0, 0, 0, 0)
                           for msg in safety_messages), 'unassigned Safety fields must retain their default values'
                intervals = [values[i][0] - values[i - 1][0] for i in range(1, len(values))]
                median_hz = 1.0 / statistics.median(intervals)
                average_hz = (len(values) - 1) / (values[-1][0] - values[0][0])
                stamps = [msg.header.stamp.sec + msg.header.stamp.nanosec * 1e-9 for msg in safety_messages]
                header_intervals = [stamps[i] - stamps[i - 1] for i in range(1, len(stamps))]
                assert all(interval > 0 for interval in header_intervals), 'Safety publish timestamps must advance'
                header_hz = 1.0 / statistics.median(header_intervals)
                print(f'Safety.runstate samples={len(values)}, ROS steady average={average_hz:.1f} Hz, '
                      f'ROS steady median={median_hz:.1f} Hz, header median={header_hz:.1f} Hz', flush=True)
                assert 45 <= median_hz <= 55 and 43 <= average_hz <= 57, (median_hz, average_hz)
                assert not snapshots, 'legacy localization DTC must not be published'
                assert node.count_publishers('/diagnostics/source_dtc/localization') == 0
                assert node.count_subscribers('/internal/diagnostics/livox_dtc') == 0, 'Safety must not depend on the retired aggregator'
                assert node.count_publishers('/diagnostics/heartbeat/livox_dtc') == 0
                for role in ('front', 'rear', 'left', 'right'):
                    expected = 0 if args.sany_lidars == 3 and role == 'rear' else 1
                    assert node.count_subscribers(f'/sensors/livox/{role}/diagnostic_state') == expected, f'incorrect required radar: {role}'
                assert node.count_subscribers('/diagnostics/source_dtc/livox') == 0, 'Safety must not use the retired external topic'
                if args.incident_output:
                    assert not input_publishers or inputs_sent, 'raw input fixture never matched localization'
                    assert not any('rosbag' in name for name in node.get_node_names()), 'collector must not launch a recorder node'
                    if args.sany_lidars:
                        topics = config['multi_lidar']['topics']
                        for i in range(args.sany_lidars):
                            assert node.count_subscribers(topics[f'lidar_{i}']) == 1, 'collector added a point-cloud subscription'
                        primary = int(config['multi_lidar'].get('primary_lidar_id', 0))
                        for i in range(args.sany_lidars):
                            assert node.count_subscribers(topics[f'imu_{i}']) == (1 if i == primary else 0), 'collector added an unused IMU subscription'
                assert len(heartbeats) >= 20, 'heartbeat must continue during UNKNOWN'
                assert all(msg.node_id == 6 and msg.boot_id and not msg.active_fault for msg in heartbeats)
                assert all(msg.work_seq == 0 for msg in heartbeats), 'heartbeat is not business work'
                assert len({msg.boot_id for msg in heartbeats}) == 1
                assert all(heartbeats[i].heartbeat_seq > heartbeats[i-1].heartbeat_seq
                           for i in range(1, len(heartbeats)))
                log_text = log_path.read_text()
                if args.sany_lidars:
                    loaded = [dict(re.findall(r'(\w+)=([^\s]+)', line))
                              for line in log_text.splitlines() if 'LOCALIZATION_DTC_CONFIG ' in line]
                    assert len(loaded) == 1, 'resolved SANY DTC config must be logged once'
                    assert loaded[0]['version'] == sany['diagnostics']['config_version']
                    for key in ('major_sigma_limit_m', 'lateral_sigma_limit_m', 'yaw_sigma_limit_deg'):
                        assert float(loaded[0][key]) == sany['diagnostics'][key], (key, loaded[0])
                audit = [dict(re.findall(r'(\w+)=([^\s]+)', line))
                         for line in log_text.splitlines() if 'SAFETY_DIAGNOSTIC ' in line]
                assert audit, 'fault audit must be present independently of reporting and shadow switches'
                assert {int(row['raw_runstate']) for row in audit} == {4, -9}, 'retain startup and timeout diagnoses'
                assert len({row['sequence'] for row in audit}) == len(audit), 'audit sequence must advance'
                assert all('steady_sec' in row and 'covariance_valid' in row and 'radar_all_failed' in row
                           for row in audit), 'audit must include timing and causal evidence'
                if args.safety_fault_reporting == 'false':
                    assert all(row['published_runstate'] == '0' and row['suppressed'] == '1'
                               and row['reporting_enabled'] == '0' for row in audit)
                else:
                    assert all(row['published_runstate'] == row['raw_runstate']
                               and row['suppressed'] == '0' for row in audit)
                if args.radar_time_sync_shadow == 'false':
                    assert 'SAFETY_TIME_SYNC_SHADOW ' not in log_text, 'optional shadow logging is off'
                print(f'dtc_online_smoke passed: Safety {average_hz:.1f} Hz; '
                      f'safety_fault_reporting={args.safety_fault_reporting}; '
                      f'shadow={args.radar_time_sync_shadow}; audit_rows={len(audit)}; heartbeat DTC only, startup UNKNOWN')
            finally:
                process.terminate()
                try:
                    process.wait(timeout=30 if args.incident_output else 5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
                if args.log_output:
                    args.log_output.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(log_path, args.log_output)
                node.destroy_node()
                rclpy.shutdown()
        if args.incident_output:
            import rosbag2_py
            summaries = list(args.incident_output.glob('event_*/summary.json'))
            assert summaries, 'no incident exported by real localization process'
            reasons, recorded_topics, raw_values = set(), set(), set()
            messages = 0
            for path in summaries:
                summary = json.loads(path.read_text())
                assert summary['status'] in ('complete', 'partial'), summary
                assert not summary['error'], summary
                reasons.update(summary['reasons'])
                assert (path.parent / 'summary.html').is_file()
                assert (path.parent / 'config.yaml').is_file()
                with (path.parent / 'timeline.csv').open() as timeline:
                    for row in csv.DictReader(timeline):
                        raw_values.add(int(row['raw']))
                        if args.safety_fault_reporting == 'false':
                            assert int(row['published']) == 0, row
                for segment in summary['closed_segments']:
                    reader = rosbag2_py.SequentialReader()
                    reader.open(rosbag2_py.StorageOptions(uri=str(path.parent / segment), storage_id='mcap'),
                                rosbag2_py.ConverterOptions('', ''))
                    while reader.has_next():
                        topic, data, stamp = reader.read_next()
                        if args.sany_lidars and topic in config['multi_lidar']['topics'].values():
                            if topic == config['multi_lidar']['topics'][f'imu_{primary}']:
                                assert deserialize_message(data, Imu).linear_acceleration.z == 9.81
                            else:
                                assert bytes(deserialize_message(data, PointCloud2).data) == struct.pack('<ffffBB6xd', 5., 0., 0., 42., 0, 0, 0.)
                        recorded_topics.add(topic)
                        messages += 1
            assert {'raw_4', 'raw_-9', 'startup_without_output'} <= reasons, reasons
            assert {4, -9} <= raw_values, raw_values
            assert {'/safety/throstle', '/diagnostics/heartbeat/lightning_slam',
                    '/localization/pipeline_diagnostics'} <= recorded_topics, recorded_topics
            if args.sany_lidars:
                topics = config['multi_lidar']['topics']
                assert {topics[f'lidar_{i}'] for i in range(args.sany_lidars)} <= recorded_topics, recorded_topics
                assert topics[f'imu_{primary}'] in recorded_topics
                assert all(topics[f'imu_{i}'] not in recorded_topics for i in range(args.sany_lidars) if i != primary)
            print(f'in-process capture passed: {len(summaries)} events, {messages} MCAP messages; '
                  f'masked raw faults retained; no additional input DDS subscriptions')


if __name__ == '__main__':
    main()
