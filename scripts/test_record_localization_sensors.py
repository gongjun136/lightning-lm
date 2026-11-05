#!/usr/bin/env python3
"""Acquisition contracts: raw source selection and conservative quality reports."""
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace as NS
import unittest

import record_localization_sensors as capture


class CaptureTest(unittest.TestCase):
    def test_sensor_only_session_uses_existing_site_dds_profile(self):
        with tempfile.TemporaryDirectory() as directory:
            profile = Path(directory) / 'fastdds.xml'
            environment = {}
            capture.configure_transport_environment(environment, profile)
            self.assertEqual(environment, {})
            profile.write_text('<profiles/>')
            capture.configure_transport_environment(environment, profile)
            self.assertEqual(environment['FASTRTPS_DEFAULT_PROFILES_FILE'], str(profile))
            self.assertEqual(environment['FASTRMW_DEFAULT_PROFILES_FILE'], str(profile))

    def test_explicit_transport_and_other_rmw_are_preserved(self):
        environment = {'FASTRMW_DEFAULT_PROFILES_FILE': '/operator/dds.xml'}
        capture.configure_transport_environment(environment, Path('/unused'))
        self.assertEqual(environment['FASTRTPS_DEFAULT_PROFILES_FILE'], '/operator/dds.xml')
        environment = {'RMW_IMPLEMENTATION': 'rmw_cyclonedds_cpp'}
        capture.configure_transport_environment(environment, Path(__file__))
        self.assertEqual(environment, {'RMW_IMPLEMENTATION': 'rmw_cyclonedds_cpp'})

    def samples(self):
        values = {'latitude': [30.], 'longitude': [110.], 'altitude': [100.],
                  'attitude': [20., 0., 0.], 'velocity': [0., 0., 0.],
                  'position_sigma': [.01, .01, .02], 'attitude_sigma': [.1, .1, .1],
                  'status': [2, 4, .5], 'velocity_sigma': [.01, .01, .01]}
        return {k: (NS(frame=NS(valid=True), **dict(zip(attrs, values[k]))), 10., 1_000_000_000)
                for k, (_, _, attrs) in capture.FIELDS.items()}

    def test_plan_requires_navigation_reference_and_keeps_auxiliary_inputs(self):
        cfg = {'multi_lidar': {'enabled': True, 'topics': {'lidar_0': '/cloud0', 'imu_0': '/imu0',
                 'lidar_1': '/cloud1', 'imu_1': '/imu1'}},
               'system': {'enable_wheel_speed_observation': True, 'wheel_speed_topic': '/wheel'}}
        p = capture.make_plan(cfg, '/reference', ['/steering'], mode='reference')
        for topic in ['/cloud0', '/imu0', '/steering', '/reference/attitude/sigma',
                      '/reference/time', '/reference/driver_status']:
            self.assertIn(topic, p['required_topics'])
        for topic in ['/imu1', '/wheel', '/reference/imu/acceleration_raw', '/reference/imu/angular_rate_raw',
                      '/reference/vehicle/acceleration', '/reference/vehicle/angular_rate']:
            self.assertIn(topic, p['optional_input_topics'])
            self.assertNotIn(topic, p['required_topics'])
        matcher = capture.re.compile(p['topic_regex'])
        self.assertTrue(matcher.fullmatch('/reference/new_output'))
        self.assertTrue(matcher.fullmatch('/livox/imu_new_sensor'))
        self.assertTrue(matcher.fullmatch('/wheel'))
        self.assertTrue(matcher.fullmatch('/reference/imu/angular_rate_raw'))
        self.assertFalse(matcher.fullmatch('/camera/front/image'))
        strict = capture.make_plan(cfg, '/reference', ['/wheel', '/reference/imu/acceleration_raw'])
        self.assertIn('/wheel', strict['required_topics'])
        self.assertIn('/reference/imu/acceleration_raw', strict['required_topics'])
        self.assertNotIn('/wheel', strict['optional_input_topics'])

    def test_basic_needs_primary_imu_but_not_reference_or_auxiliary_imu(self):
        cfg = {'multi_lidar': {'enabled': True, 'primary_lidar_id': 1,
               'topics': {'lidar_0': '/c0', 'lidar_1': '/c1', 'imu_0': '/i0', 'imu_1': '/i1'}}}
        plan = capture.make_plan(cfg, '/cgi430', [])
        self.assertEqual(plan['required_topics'], ['/c0', '/c1', '/i1'])
        self.assertIn('/i0', plan['optional_input_topics'])
        self.assertIn('/cgi430/time', plan['optional_input_topics'])
        self.assertIn('/SpeThrCAN4_topic', plan['optional_input_topics'])
        self.assertIsNotNone(capture.re.fullmatch(plan['topic_regex'], '/SpeThrCAN4_topic'))
        for mode in ('basic', 'reference'):
            strict = capture.make_plan(cfg, '/cgi430', [], mode, True)
            self.assertIn('/SpeThrCAN4_topic', strict['required_topics'])
            self.assertNotIn('/SpeThrCAN4_topic', strict['optional_input_topics'])

    def test_missing_primary_imu_rejected(self):
        with self.assertRaises(ValueError):
            capture.make_plan({'common': {'lidar_topic': '/cloud'}}, '/cgi430', [])

    def test_wheel_late_arrival_gap_recovery_and_tail(self):
        wheel = capture.WheelContinuity(0)
        self.assertEqual(wheel.report(5)['status'], 'not_received')
        wheel.receive(5)
        wheel.receive(5.1)
        self.assertEqual(wheel.report(5.2)['status'], 'received')
        self.assertEqual(wheel.report(5.2)['first_arrival_delay_sec'], 5)
        wheel.receive(8)
        self.assertEqual(wheel.report(8.1)['gap_count'], 1)
        self.assertEqual(wheel.report(11)['gap_count'], 2)
        wheel.receive(11)
        self.assertEqual(wheel.report(11)['gap_count'], 2)
        self.assertEqual(wheel.report(11)['status'], 'interrupted')

    def test_invalid_mode_flags_fail_before_ros_or_session_creation(self):
        for args in (['--mode', 'basic', '--require-fixed'],
                     ['--check-only', '--output', '/must-not-create'],
                     ['--duration-seconds', 'nan']):
            with self.assertRaises(SystemExit) as error:
                capture.main(args)
            self.assertEqual(error.exception.code, 2)

    def test_fixed_flag_alone_does_not_establish_quality(self):
        samples = self.samples()
        good = capture.quality_report(samples, 1_010_000_000, 10.01)
        self.assertTrue(good['device_fixed_flag'])
        self.assertTrue(good['instant_quality_ok'])
        samples['position_sigma'][0].east_m = .5
        bad = capture.quality_report(samples, 1_010_000_000, 10.01)
        self.assertTrue(bad['device_fixed_flag'])
        self.assertFalse(bad['instant_quality_ok'])
        self.assertIn('position_sigma:limit', bad['issues'])

    def test_missing_stale_future_and_skew_fail(self):
        for change, issue in [(lambda s: s.pop('latitude'), 'latitude:missing'),
                              (lambda s: s.update(latitude=(s['latitude'][0], 9., 1_000_000_000)), 'latitude:stale'),
                              (lambda s: s.update(latitude=(s['latitude'][0], 10., 1_050_000_000)), 'latitude:future'),
                              (lambda s: s.update(latitude=(s['latitude'][0], 10., 970_000_000)), 'field_skew')]:
            samples = self.samples()
            change(samples)
            result = capture.quality_report(samples, 1_010_000_000, 10.01)
            self.assertFalse(result['instant_quality_ok'])
            self.assertIn(issue, result['issues'])

    def test_nonfinite_values_remain_json_auditable(self):
        samples = self.samples()
        samples['status'][0].differential_age_s = float('nan')
        result = capture.quality_report(samples, 1_010_000_000, 10.01)
        self.assertFalse(result['instant_quality_ok'])
        json.dumps(result, allow_nan=False)


if __name__ == '__main__':
    unittest.main()
