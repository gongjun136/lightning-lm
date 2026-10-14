#!/usr/bin/env python3
"""Verify nontrivial sensor/body composition and preservation of recording coverage."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import yaml

import accuracy_normalize_tum as normalizer


class NormalizeTests(unittest.TestCase):
    def test_rotated_lever_arm_and_imu_lidar_rotation(self):
        s = np.sqrt(.5)
        config = dict(schema_version=1, source_frame='imu',
                      T_imu_lidar=dict(confirmed=True, evidence='bench',
                          translation_m=[1, 0, 0], quaternion_xyzw=[0, 0, s, s]),
                      T_rear_axle_lidar=dict(confirmed=True, evidence='survey',
                          translation_m=[0, 2, 0], quaternion_xyzw=[s, 0, 0, s]))
        q, t = normalizer.source_from_body(config)
        # T_IL Rz90 with t=(1,0,0); inverse of T_BL Rx90,t=(0,2,0).
        np.testing.assert_allclose(t, [1, 0, 2], atol=1e-12)
        # A vector on body Y maps into IMU -Z, proving noncommuting rotations.
        np.testing.assert_allclose(normalizer.rotate(q, np.array([0, 1, 0.])), [0, 0, -1], atol=1e-12)
        # Applying world Rz90 rotates the lever arm; must not subtract a world constant.
        np.testing.assert_allclose(normalizer.rotate(np.array([0, 0, s, s]), t), [0, 1, 2], atol=1e-12)

    def test_no_double_compensation_and_no_duration_shrink(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            source = d / 'input.tum'
            source.write_text('2 3 4 5 0 0 0 1\n3 4 5 6 0 0 0 1\n', encoding='utf-8')
            cfg = dict(schema_version=1, source_frame='rear_axle_flu', timestamp_offset_ns=0,
                       published_reference_point_confirmed=True,
                       published_reference_point_evidence='publication contract',
                       world_frame_id='map:abc', trajectory_kind='localization')
            (d / 'config.yaml').write_text(yaml.safe_dump(cfg), encoding='utf-8')
            record = dict(source_recording_id='C', time_basis='utc_sensor_header',
                          recorded_intervals=[dict(start_ns=0, end_ns=10_000_000_000, epoch='C_0')])
            (d / 'record.json').write_text(json.dumps(record), encoding='utf-8')
            normalizer.normalize(source, d / 'config.yaml', d / 'record.json', d / 'out.tum')
            np.testing.assert_allclose(np.loadtxt(d / 'out.tum'), np.loadtxt(source))
            meta = json.loads((d / 'out.metadata.json').read_text(encoding='utf-8'))
            self.assertEqual(meta['recorded_intervals'], record['recorded_intervals'])
            self.assertEqual(meta['body_frame_id'], 'rear_axle_flu')

    def test_known_output_world_change_is_left_multiplied(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            (d / 'input.tum').write_text('2.000000123 2 0 0 0 0 0 1\n', encoding='utf-8')
            s = float(np.sqrt(.5))
            cfg = dict(schema_version=1, source_frame='rear_axle_flu', timestamp_offset_ns=0,
                       published_reference_point_confirmed=True, published_reference_point_evidence='contract',
                       source_world_frame_id='legacy:1', world_frame_id='map:1', trajectory_kind='localization',
                       world_transform=dict(enabled=True, confirmed=True, evidence='inverse manifest output transform',
                                            translation_m=[10, 0, 0], quaternion_xyzw=[0, 0, s, s]))
            (d / 'config.yaml').write_text(yaml.safe_dump(cfg), encoding='utf-8')
            record = dict(source_recording_id='C', time_basis='utc_sensor_header',
                          recorded_intervals=[dict(start_ns=0, end_ns=10_000_000_000, epoch='C_0')])
            (d / 'record.json').write_text(json.dumps(record), encoding='utf-8')
            normalizer.normalize(d / 'input.tum', d / 'config.yaml', d / 'record.json', d / 'out.tum')
            row = np.loadtxt(d / 'out.tum')
            np.testing.assert_allclose(row[1:4], [10, 2, 0], atol=1e-12)
            self.assertTrue((d / 'out.tum').read_text().startswith('2.000000123 '))

    def test_unconfirmed_mount_rejected(self):
        with self.assertRaisesRegex(ValueError, 'confirmed'):
            normalizer.source_from_body(dict(schema_version=1, source_frame='imu'))

    def test_pose_frame_relabel_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            (d / 'input.tum').write_text('2 0 0 0 0 0 0 1\n', encoding='utf-8')
            cfg = dict(schema_version=1, source_frame='rear_axle_flu', timestamp_offset_ns=0,
                       published_reference_point_confirmed=True, published_reference_point_evidence='contract',
                       world_frame_id='map:1', source_world_frame_id='map:1', trajectory_kind='localization')
            (d / 'config.yaml').write_text(yaml.safe_dump(cfg), encoding='utf-8')
            record = dict(source_recording_id='C', time_basis='utc_sensor_header',
                          world_frame_id='legacy:1', body_frame_id='rear_axle_flu',
                          recorded_intervals=[dict(start_ns=0, end_ns=10_000_000_000, epoch='C_0')])
            (d / 'record.json').write_text(json.dumps(record), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'cannot be relabeled'):
                normalizer.normalize(d / 'input.tum', d / 'config.yaml', d / 'record.json', d / 'out.tum')

    def test_complete_output_outage_still_exports_recording_denominator(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            (d / 'input.tum').write_text('', encoding='utf-8')
            cfg = dict(schema_version=1, source_frame='rear_axle_flu', timestamp_offset_ns=0,
                       published_reference_point_confirmed=True, published_reference_point_evidence='contract',
                       world_frame_id='map:1', trajectory_kind='localization')
            (d / 'config.yaml').write_text(yaml.safe_dump(cfg), encoding='utf-8')
            record = dict(source_recording_id='C', time_basis='utc_sensor_header',
                          conditions=['SYNTHETIC'],
                          recorded_intervals=[dict(start_ns=0, end_ns=10_000_000_000, epoch='C_0')])
            (d / 'record.json').write_text(json.dumps(record), encoding='utf-8')
            result = normalizer.normalize(d / 'input.tum', d / 'config.yaml', d / 'record.json', d / 'out.tum')
            self.assertEqual(result['rows'], 0)
            self.assertEqual(json.loads((d / 'out.metadata.json').read_text())['recorded_intervals'], record['recorded_intervals'])


if __name__ == '__main__':
    unittest.main()
