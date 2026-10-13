#!/usr/bin/env python3
"""Validate the real run_loc_online publishers with a temporary, empty map.

source install/setup.bash
RMW_IMPLEMENTATION=rmw_cyclonedds_cpp ROS_DOMAIN_ID=172 python3 \
  src/lightning-lm/src/test/dtc_online_smoke.py build/lightning_lm/bin/run_loc_online
"""
import os
from pathlib import Path
import statistics
import subprocess
import sys
import tempfile
import time

import rclpy
import yaml
from rclpy.qos import QoSProfile, ReliabilityPolicy, DurabilityPolicy
from diagnostic_monitor_interfaces.msg import DtcObservationArray
from geosun_msgs.msg import Safty


def main():
    os.environ.setdefault('ROS_DOMAIN_ID', '172')
    os.environ.setdefault('ROS_LOCALHOST_ONLY', '1')
    executable = Path(sys.argv[1]).resolve()
    source_root = Path(__file__).resolve().parents[2]
    config = yaml.safe_load((source_root / 'config/default.yaml').read_text())
    config['system']['with_ui'] = False
    config['relocalization'] = {'enable_map_consistency': False, 'backend': 'btc'}
    config['map_export'] = {'normalize_start_ground_z': False}
    config['maps']['save_dyn_when_quit'] = False
    config['diagnostics'] = {'startup_grace_sec': 2.0, 'config_version': 'smoke-test'}
    with tempfile.TemporaryDirectory(prefix='dtc_online_', dir=executable.parent.parent) as directory:
        fixture = Path(directory)
        (fixture / 'index.txt').write_text('0 0 0\n')
        path = fixture / 'config.yaml'
        path.write_text(yaml.safe_dump(config))
        rclpy.init()
        node = rclpy.create_node('localization_dtc_smoke_probe')
        qos = QoSProfile(depth=100, reliability=ReliabilityPolicy.RELIABLE,
                         durability=DurabilityPolicy.VOLATILE)
        values, snapshots = [], []
        node.create_subscription(Safty, '/safety/throstle', lambda msg: values.append((time.monotonic(), msg.data)), qos)
        node.create_subscription(DtcObservationArray, '/diagnostics/source_dtc/localization', snapshots.append, qos)
        log_path = fixture / 'process.log'
        with log_path.open('w') as log:
            process = subprocess.Popen([str(executable), '--config=' + str(path), '--map=' + str(fixture)],
                                       cwd=fixture, stdout=log, stderr=log)
            try:
                end = time.monotonic() + 6.0
                while time.monotonic() < end:
                    rclpy.spin_once(node, timeout_sec=0.01)
                    if process.poll() is not None:
                        raise RuntimeError(log_path.read_text()[-5000:])
                assert len(values) >= 150, f'too few safety samples: {len(values)}'
                assert any(value == 4 for _, value in values), 'startup UNKNOWN must be 4'
                assert values[-1][1] == -9, 'no valid output after startup grace must be -9'
                assert all(value in (4, -9) for _, value in values), 'uninitialized localization must never give zero'
                intervals = [values[i][0] - values[i - 1][0] for i in range(1, len(values))]
                median_hz = 1.0 / statistics.median(intervals)
                average_hz = (len(values) - 1) / (values[-1][0] - values[0][0])
                assert 45 <= median_hz <= 55 and 43 <= average_hz <= 57, (median_hz, average_hz)
                assert len(snapshots) >= 20
                for msg in snapshots:
                    assert msg.comm_header.source_id == 6 and msg.boot_id and msg.interface_version == 1
                    assert [obs.dtc for obs in msg.observations] == [0x30811813, 0x40450198]
                    assert all(obs.physical_condition == 0 for obs in msg.observations)
                    assert all(obs.evaluation_seq == 0 for obs in msg.observations)
                print(f'dtc_online_smoke passed: real Safty publisher {average_hz:.1f} Hz, median {median_hz:.1f} Hz; startup 4 -> timeout -9; complete localization DTC snapshots')
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
