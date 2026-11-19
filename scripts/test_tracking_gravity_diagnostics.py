"""跟踪独立重力门限的实际日志与旧证据兼容性。"""
import math
from pathlib import Path
import tempfile
import unittest
import diagnose_localization as diagnose
import diagnostic_summary as summary


class GravityDiagnosticTest(unittest.TestCase):
    def test_startup_log_preserves_independent_thresholds(self):
        with tempfile.TemporaryDirectory() as folder:
            data = diagnose.overlap_observation_report({
                'tracking_map_overlap_check_enabled': '0', 'map_consistency_enabled': '1',
                'tracking_map_max_gravity_angle_deg': '20', 'tracking_map_min_gravity_cos': str(math.cos(math.radians(20))),
                'initialization_map_min_gravity_cos': '.95', 'relocalization_map_min_gravity_cos': '.95'}, [], Path(folder), None)
        text = summary.render_summary({'overlap_observations': data})
        self.assertIn('跟踪夹角上限 20°', text)
        self.assertIn('初始化余弦下限 0.95', text)
        self.assertIn('重定位余弦下限 0.95', text)

    def test_legacy_missing_threshold_never_infers_twenty_degrees(self):
        with tempfile.TemporaryDirectory() as folder:
            data = diagnose.overlap_observation_report({}, [], Path(folder), None)
        self.assertIsNone(data['policy']['tracking_map_max_gravity_angle_deg'])
        text = summary.render_summary({'overlap_observations': data,
            'matching_thresholds': {'tracking_map_max_gravity_angle_deg': 20}})
        self.assertNotIn('跟踪夹角上限 20', text)
        view = summary.rejection_view({'reason_code': 'map_gravity', 'evidence': {
            'map_checked': 1, 'map_gravity_cos': .9498}})['text']
        self.assertIn('实际下限=未记录', view)
        self.assertNotIn('20', view)

    def test_historical_comparison_uses_its_own_actual_threshold(self):
        for threshold in (.95, math.cos(math.radians(20))):
            raw = f'rejection_reason=map_gravity map_checked=1 map_points=100 map_gravity_cos=0.9498 map_min_gravity_cos={threshold:.17g}'
            text, _ = diagnose.matching_summary({'time': 'T+1', 'text': raw}, {})
            self.assertEqual('已证实未满足项: map_gravity_cos' in text, .9498 < threshold)
            view = summary.rejection_view({'reason_code': 'map_gravity', 'evidence': {
                'map_checked': 1, 'map_gravity_cos': .9498, 'map_min_gravity_cos': threshold}})['text']
            self.assertIn('实际下限=', view)
            self.assertIn('实际上限=', view)

    def test_unchecked_map_does_not_claim_angle_measurement(self):
        text=summary.rejection_view({'evidence': {'map_checked': 0, 'map_gravity_cos': 0, 'map_min_gravity_cos': None}})['text']
        self.assertNotIn('重力余弦=', text)


if __name__ == '__main__': unittest.main()
