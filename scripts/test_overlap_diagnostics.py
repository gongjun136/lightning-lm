"""overlap_diagnostics 的合成回归测试（unittest，仅标准库，不依赖项目文件）。"""

import copy
import json
import unittest

import overlap_diagnostics as od


def make_row(**overrides):
    row = {
        'map_overlap_check_enabled': 0,
        'map_checked': 1,
        'map_overlap': 0.59,
        'map_min_overlap': 0.60,
        'accepted': 1,
    }
    row.update(overrides)
    return row


def drop(row, *keys):
    return {key: value for key, value in row.items() if key not in keys}


def joined(info):
    return '|'.join(info['lines'])


class GateAndResultTest(unittest.TestCase):
    def test_observation_event_is_accepted_low(self):
        info = od.describe_overlap(make_row(event=od.LOW_OVERLAP_EVENT))
        self.assertEqual(info['gate'], 'disabled')
        self.assertEqual(info['result'], 'bypassed_low')
        self.assertIs(info['accepted'], True)
        self.assertIs(info['map_checked'], True)
        self.assertTrue(info['is_observation'])
        self.assertIn('不是故障', joined(info))
        self.assertEqual(info['warnings'], [])
        self.assertAlmostEqual(info['ratio'], 0.59)
        self.assertAlmostEqual(info['threshold'], 0.60)

    def test_threshold_equality_passes_for_both_gates(self):
        for gate in (0, 1, 'true', 'false'):
            for ratio in (0.60, '0.60', '0.6'):
                with self.subTest(gate=gate, ratio=ratio):
                    info = od.describe_overlap(make_row(map_overlap_check_enabled=gate, map_overlap=ratio))
                    self.assertEqual(info['result'], 'pass')

    def test_point_59_bypassed_or_rejected(self):
        self.assertEqual(od.describe_overlap(make_row())['result'], 'bypassed_low')
        info = od.describe_overlap(make_row(map_overlap_check_enabled=1, accepted=0, rejection_reason='map_overlap'))
        self.assertEqual(info['gate'], 'enabled')
        self.assertEqual(info['result'], 'fail')
        self.assertIs(info['accepted'], False)
        self.assertFalse(info['is_observation'])

    def test_enforced_low_but_accepted_is_inconsistent(self):
        info = od.describe_overlap(make_row(map_overlap_check_enabled=True, accepted=True))
        self.assertEqual(info['result'], 'inconsistent')
        self.assertTrue(info['warnings'])

    def test_map_overlap_reason_with_gate_disabled_is_inconsistent(self):
        info = od.describe_overlap(make_row(accepted=0, rejection_reason='map_overlap'))
        self.assertEqual(info['result'], 'inconsistent')

    def test_map_overlap_reason_with_passing_ratio_is_inconsistent(self):
        info = od.describe_overlap(make_row(map_overlap_check_enabled=1, map_overlap=0.60, accepted=0,
                                            rejection_reason='map_overlap'))
        self.assertEqual(info['result'], 'inconsistent')

    def test_rejected_frame_bypassed_low_not_rendered_as_accepted(self):
        info = od.describe_overlap(make_row(accepted=0, rejection_reason='other_reason'))
        self.assertEqual(info['result'], 'bypassed_low')
        self.assertIs(info['accepted'], False)
        self.assertFalse(info['is_observation'])
        self.assertNotIn('本帧已被接受', joined(info))
        self.assertNotIn('不是故障', joined(info))

    def test_absent_accepted_does_not_claim_acceptance(self):
        info = od.describe_overlap(drop(make_row(), 'accepted'))
        self.assertIsNone(info['accepted'])
        self.assertEqual(info['result'], 'bypassed_low')
        self.assertFalse(info['is_observation'])
        self.assertNotIn('本帧已被接受', joined(info))

    def test_map_checked_false_is_unavailable_despite_stale_metric(self):
        for gate in (0, 1):
            for stale in (0.1, float('nan'), 0.9):
                with self.subTest(gate=gate, stale=stale):
                    info = od.describe_overlap(make_row(map_checked=0, map_overlap=stale,
                                                        map_overlap_check_enabled=gate))
                    self.assertEqual(info['result'], 'unavailable')
                    self.assertIsNone(info['ratio'])
                    self.assertIsNone(info['threshold'])
                    self.assertIs(info['map_checked'], False)

    def test_missing_map_checked_is_unknown(self):
        info = od.describe_overlap(drop(make_row(), 'map_checked'))
        self.assertEqual(info['result'], 'unknown')
        self.assertIsNone(info['map_checked'])

    def test_missing_map_checked_with_explicit_contradiction(self):
        row = drop(make_row(accepted=1, rejection_reason='map_overlap'), 'map_checked')
        self.assertEqual(od.describe_overlap(row)['result'], 'inconsistent')

    def test_invalid_metrics_are_not_healthy(self):
        bad_values = (float('nan'), float('inf'), float('-inf'), -0.01, 1.01, 'nan', 'inf', 'abc', True, [], None)
        for field in ('map_overlap', 'map_min_overlap'):
            for value in bad_values:
                with self.subTest(field=field, value=value):
                    info = od.describe_overlap(make_row(**{field: value}))
                    self.assertEqual(info['result'], 'invalid')
                    self.assertTrue(info['warnings'])
                    self.assertFalse(info['is_observation'])
        info = od.describe_overlap(drop(make_row(), 'map_overlap'))
        self.assertEqual(info['result'], 'invalid')
        nan_info = od.describe_overlap(make_row(map_overlap=float('nan')))
        self.assertIsNone(nan_info['ratio'])
        self.assertEqual(nan_info['result'], 'invalid')

    def test_old_log_without_gate_is_unknown(self):
        old = drop(make_row(), 'map_overlap_check_enabled')
        info = od.describe_overlap(old)
        self.assertEqual(info['gate'], 'unknown')
        self.assertEqual(info['result'], 'unknown')
        passing = od.describe_overlap(drop(make_row(map_overlap=0.8), 'map_overlap_check_enabled'))
        self.assertEqual(passing['gate'], 'unknown')
        self.assertEqual(passing['result'], 'pass')
        rejected = od.describe_overlap(drop(make_row(accepted=0, rejection_reason='map_overlap'),
                                            'map_overlap_check_enabled'))
        self.assertEqual(rejected['gate'], 'unknown')
        self.assertEqual(rejected['result'], 'unknown')

    def test_event_contradictions(self):
        event = od.LOW_OVERLAP_EVENT
        self.assertEqual(od.describe_overlap(make_row(event=event, accepted=0))['result'], 'inconsistent')
        self.assertEqual(od.describe_overlap(make_row(event=event, map_overlap_check_enabled=1))['result'],
                         'inconsistent')
        self.assertEqual(od.describe_overlap(make_row(event=event, map_overlap=0.7))['result'], 'inconsistent')
        self.assertEqual(od.describe_overlap(make_row(event=event, map_checked=0))['result'], 'inconsistent')
        historical = od.describe_overlap(drop(make_row(event=event), 'map_overlap_check_enabled'))
        self.assertEqual(historical['gate'], 'unknown')
        self.assertEqual(historical['result'], 'unknown')


class StrictJsonTest(unittest.TestCase):
    NONFINITE = (float('nan'), float('inf'), float('-inf'), 'nan', 'NaN', 'inf', '-Infinity', '1e400')
    FIELDS = (('map_overlap', 'ratio'), ('map_min_overlap', 'threshold'))

    def assert_strict_json(self, info):
        json.dumps(info, allow_nan=False, ensure_ascii=False)

    def test_nonfinite_metrics_serialize_strictly(self):
        for field, key in self.FIELDS:
            for value in self.NONFINITE:
                for checked_absent in (False, True):
                    with self.subTest(field=field, value=value, checked_absent=checked_absent):
                        row = make_row(**{field: value})
                        if checked_absent:
                            row = drop(row, 'map_checked')
                        info = od.describe_overlap(row)
                        self.assertIsNone(info[key])
                        self.assertEqual(info['result'], 'unknown' if checked_absent else 'invalid')
                        self.assertTrue(info['warnings'])
                        self.assert_strict_json(info)

    def test_both_metrics_nonfinite_serialize_strictly(self):
        for checked_absent in (False, True):
            with self.subTest(checked_absent=checked_absent):
                row = make_row(map_overlap=float('nan'), map_min_overlap=float('inf'))
                if checked_absent:
                    row = drop(row, 'map_checked')
                info = od.describe_overlap(row)
                self.assertIsNone(info['ratio'])
                self.assertIsNone(info['threshold'])
                self.assert_strict_json(info)

    def test_contradictory_and_missing_evidence_serialize_strictly(self):
        rows = (
            make_row(event=od.LOW_OVERLAP_EVENT, map_overlap=float('nan'), observation_streak_sec=float('nan')),
            make_row(map_overlap_check_enabled=1, accepted=1, map_min_overlap=float('-inf')),
            drop(make_row(rejection_reason='map_overlap', map_overlap=float('inf')), 'map_checked'),
            drop(make_row(), 'map_overlap', 'map_min_overlap', 'map_checked'),
        )
        for row in rows:
            with self.subTest(row=repr(row)):
                info = od.describe_overlap(row)
                self.assertIn(info['result'], ('invalid', 'inconsistent', 'unknown'))
                self.assert_strict_json(info)


class StrictBoolTest(unittest.TestCase):
    def test_bool_strings_do_not_truthify_arbitrary_words(self):
        for value in ('yes', 'no', 'on', 'off', 'enabled', 'tru', ' true', 'true ', '', '2', '01', 'tRuE',
                      2, -1, 1.0, 0.0, [], {}):
            with self.subTest(value=value):
                info = od.describe_overlap(make_row(map_overlap_check_enabled=value, accepted=value))
                self.assertEqual(info['gate'], 'unknown')
                self.assertIsNone(info['accepted'])
                self.assertTrue(info['warnings'])
                self.assertNotIn('本帧已被接受', joined(info))

    def test_exact_bool_spellings(self):
        cases = ((True, 'enabled'), (False, 'disabled'), (1, 'enabled'), (0, 'disabled'),
                 ('1', 'enabled'), ('0', 'disabled'), ('true', 'enabled'), ('True', 'enabled'),
                 ('TRUE', 'enabled'), ('false', 'disabled'), ('False', 'disabled'), ('FALSE', 'disabled'))
        for value, expected in cases:
            with self.subTest(value=value):
                self.assertEqual(od.describe_overlap(make_row(map_overlap_check_enabled=value))['gate'], expected)


class ObservationCounterTest(unittest.TestCase):
    def test_counters_are_optional(self):
        info = od.describe_overlap(make_row())
        for name in ('streak_frames', 'streak_seconds', 'total_low_frames', 'capture_suppressed'):
            self.assertIsNone(info[name])
        self.assertEqual(info['warnings'], [])

    def test_valid_counters(self):
        info = od.describe_overlap(make_row(observation_streak_frames=3, observation_streak_sec='0.2',
                                            observation_total_low_frames='12',
                                            observation_capture_suppressed=11))
        self.assertEqual(info['streak_frames'], 3)
        self.assertAlmostEqual(info['streak_seconds'], 0.2)
        self.assertEqual(info['total_low_frames'], 12)
        self.assertEqual(info['capture_suppressed'], 11)
        self.assertEqual(info['warnings'], [])

    def test_large_integer_counters_preserved_exactly(self):
        big = 2 ** 53 + 1
        info = od.describe_overlap(make_row(observation_total_low_frames=big,
                                            observation_capture_suppressed=str(big),
                                            observation_streak_frames=2 ** 64 + 7))
        self.assertEqual(info['total_low_frames'], big)
        self.assertIsInstance(info['total_low_frames'], int)
        self.assertEqual(info['capture_suppressed'], big)
        self.assertEqual(info['streak_frames'], 2 ** 64 + 7)
        self.assertIn(str(big), joined(info))

    def test_large_float_counter_rejected_not_rounded(self):
        info = od.describe_overlap(make_row(observation_total_low_frames=float(2 ** 53 + 1)))
        self.assertIsNone(info['total_low_frames'])
        self.assertTrue(any('observation_total_low_frames' in w for w in info['warnings']))
        exact = od.describe_overlap(make_row(observation_total_low_frames=5.0))
        self.assertEqual(exact['total_low_frames'], 5)

    def test_malformed_counters_warn_instead_of_zero(self):
        for value in (-1, '-1', '1.5', 1.5, 'abc', '', True, float('nan'), [], '1e3'):
            with self.subTest(value=value):
                info = od.describe_overlap(make_row(observation_streak_frames=value))
                self.assertIsNone(info['streak_frames'])
                self.assertTrue(any('observation_streak_frames' in w for w in info['warnings']))

    def test_malformed_duration_warns(self):
        for value in ('abc', -0.5, float('nan'), float('inf'), True, []):
            with self.subTest(value=value):
                info = od.describe_overlap(make_row(observation_streak_sec=value))
                self.assertIsNone(info['streak_seconds'])
                self.assertTrue(any('observation_streak_sec' in w for w in info['warnings']))


class PurityTest(unittest.TestCase):
    def test_input_not_mutated(self):
        row = make_row(event=od.LOW_OVERLAP_EVENT, observation_streak_frames='5', extra=[1, 2])
        snapshot = copy.deepcopy(row)
        od.describe_overlap(row)
        self.assertEqual(row, snapshot)

    def test_repeat_calls_are_independent(self):
        first = od.describe_overlap(make_row())
        first['lines'].append('外部修改')
        second = od.describe_overlap(make_row())
        self.assertNotIn('外部修改', second['lines'])
        self.assertEqual(od.describe_overlap(make_row()), second)

    def test_non_mapping_input(self):
        for value in (None, [], 'row', 3):
            with self.subTest(value=value):
                info = od.describe_overlap(value)
                self.assertEqual(info['gate'], 'unknown')
                self.assertEqual(info['result'], 'unknown')
                self.assertTrue(info['warnings'])

    def test_help_summary_labels_observation(self):
        text = od.help_summary()
        self.assertIn(od.LOW_OVERLAP_EVENT, text)
        self.assertIn('观测', text)
        self.assertIn('不是故障', text)
        for field in ('accepted=1', 'map_checked=1', 'map_overlap_check_enabled=0'):
            self.assertIn(field, text)


if __name__ == '__main__':
    unittest.main()
