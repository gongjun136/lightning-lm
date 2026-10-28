#!/usr/bin/env python3
"""Kill only the spawned recorder fixture; read back already closed MCAP segments."""
import argparse
import json
from pathlib import Path
import signal
import subprocess
import time

import rosbag2_py
from rclpy.serialization import deserialize_message
from sensor_msgs.msg import PointCloud2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path, help='incident_recorder_test binary')
    parser.add_argument('output', type=Path, help='Empty directory; retained for inspection')
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error('output must be empty')
    output.mkdir(parents=True, exist_ok=True)
    archive = output / 'incidents'
    with (output / 'process.log').open('w') as log:
        process = subprocess.Popen([str(args.executable.resolve()), '--crash-fixture', str(archive)],
                                   stdout=log, stderr=log)
        try:
            deadline = time.monotonic() + 20
            while time.monotonic() < deadline:
                if process.poll() is not None:
                    raise AssertionError((output / 'process.log').read_text())
                reports = list(archive.glob('event_*/summary.json'))
                if reports:
                    report = json.loads(reports[0].read_text())
                    if report['closed_segments']:
                        break
                time.sleep(.02)
            else:
                raise AssertionError('no closed segment before timeout')
            process.kill()
            assert process.wait(timeout=5) == -signal.SIGKILL
        finally:
            if process.poll() is None:
                process.kill()
                process.wait(timeout=5)
    report = json.loads(reports[0].read_text())
    assert report['status'] == 'capturing', 'SIGKILL must not be reported as clean completion'
    assert report['closed_segments']
    messages = 0
    for segment in report['closed_segments']:
        reader = rosbag2_py.SequentialReader()
        reader.open(rosbag2_py.StorageOptions(uri=str(reports[0].parent / segment), storage_id='mcap'),
                    rosbag2_py.ConverterOptions('', ''))
        while reader.has_next():
            topic, data, stamp = reader.read_next()
            cloud = deserialize_message(data, PointCloud2)
            assert topic == '/test/cloud' and len(cloud.data) == 100000
            assert all(value == 42 for value in cloud.data)
            messages += 1
    assert messages > 0
    print(f'SIGKILL recovery passed: {len(report["closed_segments"])} closed MCAP segments, '
          f'{messages} intact point clouds; report correctly remains capturing')


if __name__ == '__main__':
    main()
