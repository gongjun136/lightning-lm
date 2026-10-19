#!/usr/bin/env python3
"""Validate the real run_loc_online publishers with a temporary, empty map.

source install/setup.bash
RMW_IMPLEMENTATION=rmw_cyclonedds_cpp ROS_DOMAIN_ID=172 python3 \
  src/lightning-lm/src/test/dtc_online_smoke.py build/lightning_lm/bin/run_loc_online
"""
import argparse
import os
from pathlib import Path
import statistics
import subprocess
import tempfile

import rclpy
import yaml
from rclpy.clock import Clock, ClockType
from rclpy.qos import QoSProfile, ReliabilityPolicy, DurabilityPolicy
from diagnostic_monitor_interfaces.msg import DtcObservationArray, NodeHeartbeat
from geosun_msgs.msg import Safety


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('--safety-fault-reporting', choices=('default', 'true', 'false'), default='default')
    args = parser.parse_args()
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
    if args.safety_fault_reporting != 'default':
        config['diagnostics']['safety_fault_reporting_enabled'] = args.safety_fault_reporting == 'true'
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
        log_path = fixture / 'process.log'
        with log_path.open('w') as log:
            process = subprocess.Popen([str(executable), '--config=' + str(path), '--map=' + str(fixture)],
                                       cwd=fixture, stdout=log, stderr=log)
            try:
                startup_deadline = steady_time() + 15.0
                end = None
                while end is None or steady_time() < end:
                    rclpy.spin_once(node, timeout_sec=0.01)
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
                assert node.count_subscribers('/internal/diagnostics/livox_dtc') == 1, 'Safety must retain its internal radar subscription'
                assert node.count_subscribers('/diagnostics/source_dtc/livox') == 0, 'Safety must not use the retired external topic'
                assert len(heartbeats) >= 20, 'heartbeat must continue during UNKNOWN'
                assert all(msg.node_id == 6 and msg.boot_id and not msg.active_fault for msg in heartbeats)
                assert all(msg.work_seq == 0 for msg in heartbeats), 'heartbeat is not business work'
                assert len({msg.boot_id for msg in heartbeats}) == 1
                assert all(heartbeats[i].heartbeat_seq > heartbeats[i-1].heartbeat_seq
                           for i in range(1, len(heartbeats)))
                print(f'dtc_online_smoke passed: Safety {average_hz:.1f} Hz; '
                      f'safety_fault_reporting={args.safety_fault_reporting}; heartbeat DTC only, startup UNKNOWN')
            finally:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
                node.destroy_node()
                rclpy.shutdown()


if __name__ == '__main__':
    main()
