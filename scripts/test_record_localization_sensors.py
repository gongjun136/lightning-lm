#!/usr/bin/env python3
"""Acquisition contracts: raw source selection and conservative quality reports."""
import json
from types import SimpleNamespace as NS
import unittest

import record_localization_sensors as capture


class CaptureTest(unittest.TestCase):
    def samples(self):
        values = {'latitude': [30.], 'longitude': [110.], 'altitude': [100.],
                  'attitude': [20., 0., 0.], 'velocity': [0., 0., 0.],
                  'position_sigma': [.01, .01, .02], 'attitude_sigma': [.1, .1, .1],
                  'status': [2, 4, .5], 'velocity_sigma': [.01, .01, .01]}
        return {k: (NS(frame=NS(valid=True), **dict(zip(attrs, values[k]))), 10., 1_000_000_000)
                for k, (_, _, attrs) in capture.FIELDS.items()}

    def test_plan_requires_every_configured_imu_and_raw_cgi(self):
        cfg = {'multi_lidar': {'enabled': True, 'topics': {'lidar_0': '/cloud0', 'imu_0': '/imu0',
                 'lidar_1': '/cloud1', 'imu_1': '/imu1'}},
               'system': {'enable_wheel_speed_observation': True, 'wheel_speed_topic': '/wheel'}}
        p = capture.make_plan(cfg, '/reference', ['/steering'])
        for topic in ['/cloud0', '/imu1', '/wheel', '/steering', '/reference/imu/acceleration_raw',
                      '/reference/time', '/reference/driver_status']:
            self.assertIn(topic, p['required_topics'])
        matcher = capture.re.compile(p['topic_regex'])
        self.assertTrue(matcher.fullmatch('/reference/new_output'))
        self.assertTrue(matcher.fullmatch('/livox/imu_new_sensor'))
        self.assertFalse(matcher.fullmatch('/camera/front/image'))

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
