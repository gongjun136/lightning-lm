'''translation_diagnostics 的单元测试(仅合成数据)。

这些测试只验证诊断证据的解析与说明契约,不能替代整车回放验证。
'''

import copy
import math
import unittest

import translation_diagnostics as td

T0 = 1767225600.0


def abstract_limits(distance):
    # 仅用于生成合成的“调用方”记录;被测模块本身从不计算阈值。
    return (0.5 + 0.1 * distance, 0.3 + 0.1 * distance)


def comparison(slot, frm, to, distance, predicted, measured, limits=None, raw_distinct=None,
               fail=None):
    # 模拟调用方:误差为 |实测 - 预测|,判定为严格大于阈值(相等通过)。
    if limits is None:
        limits = abstract_limits(distance)
    if raw_distinct is None:
        raw_distinct = slot.endswith('_raw')
    rx = measured[0] - predicted[0]
    ry = measured[1] - predicted[1]
    ex = abs(rx)
    ey = abs(ry)
    fx = ex > limits[0] if fail is None else fail[0]
    fy = ey > limits[1] if fail is None else fail[1]
    values = {
        'available': '1',
        'raw_distinct': '1' if raw_distinct else '0',
        'from_stamp': repr(frm),
        'to_stamp': repr(to),
        'dt_sec': repr(to - frm),
        'distance_m': repr(distance),
        'predicted_x_m': repr(predicted[0]),
        'predicted_y_m': repr(predicted[1]),
        'measured_x_m': repr(measured[0]),
        'measured_y_m': repr(measured[1]),
        'residual_x_m': repr(rx),
        'residual_y_m': repr(ry),
        'error_x_m': repr(ex),
        'error_y_m': repr(ey),
        'limit_x_m': repr(limits[0]),
        'limit_y_m': repr(limits[1]),
        'fail_x': '1' if fx else '0',
        'fail_y': '1' if fy else '0',
    }
    return slot, values


def adj(slot, predicted, measured, **kwargs):
    return comparison(slot, T0 + 0.4, T0 + 0.5, 0.1, predicted, measured, **kwargs)


def win(slot, predicted, measured, **kwargs):
    return comparison(slot, T0, T0 + 0.5, 0.5, predicted, measured, **kwargs)


def passing_set():
    return [
        adj('adjacent_final', (0.1, 0.0), (0.11, 0.0)),
        adj('adjacent_raw', (0.1, 0.0), (0.12, 0.0)),
        win('window_final', (0.5, 0.0), (0.52, 0.01)),
        win('window_raw', (0.5, 0.0), (0.53, 0.01)),
    ]


def with_slots(*overrides):
    slots = dict(passing_set())
    for slot, values in overrides:
        slots[slot] = values
    return list(slots.items())


def make_row(comparisons=(), selected='adjacent', origin='current', sequence=7,
             stamp=T0 + 0.5, condition=None, latched=None, available=None, extra=None):
    row = {
        td.VERSION_KEY: '1',
        td.SEQUENCE_KEY: str(sequence),
        td.STAMP_KEY: repr(stamp),
        td.SELECTED_KEY: selected,
        td.ORIGIN_KEY: origin,
    }
    filled = dict(comparisons)
    for slot, _interval, _branch in td.SLOTS:
        prefix = td.slot_prefix(slot)
        values = filled.get(slot)
        for suffix in td.SLOT_SUFFIXES:
            if values is not None:
                row[prefix + suffix] = values[suffix]
            elif suffix in td.FLAG_SUFFIXES:
                row[prefix + suffix] = '0'
            else:
                row[prefix + suffix] = 'nan'
    if available is None:
        available = selected in ('adjacent', 'window') and any(
            slot.startswith(selected + '_') for slot in filled)
    row[td.AVAILABLE_KEY] = '1' if available else '0'
    if condition is not None:
        row[td.CONDITION_KEY] = str(condition)
    if latched is not None:
        row[td.LATCH_KEY] = '1' if latched else '0'
    if extra:
        row.update(extra)
    return row


def slot_names(items):
    return [item['slot'] for item in items]


def find_line(lines, prefix):
    matches = [line for line in lines if line.startswith(prefix)]
    return matches[0] if matches else None


ADJ_FINAL_FAIL = adj('adjacent_final', (0.1, 0.0), (0.75, 0.0))
ADJ_RAW_FAIL = adj('adjacent_raw', (0.1, 0.0), (0.75, 0.0))
WIN_FINAL_FAIL_Y = win('window_final', (0.5, 0.0), (0.5, 0.4))


class SelectionAndTriggerTest(unittest.TestCase):

    def test_final_only_failure(self):
        row = make_row(with_slots(ADJ_FINAL_FAIL), selected='adjacent', condition=2)
        result = td.describe_translation(row)
        self.assertEqual(result['status'], 'complete')
        self.assertEqual(slot_names(result['triggers']), ['adjacent_final'])
        line = find_line(result['lines'], '[当前比较] 相邻区间/发布分支')
        self.assertIsNotNone(line)
        for text in ('实测', '阈值', '超限', '实际时长', 'D=0.1000 m', repr(T0 + 0.4),
                     repr(T0 + 0.5), 'X 实测误差 0.6500 m / 阈值 0.5100 m:超限'):
            self.assertIn(text, line)

    def test_raw_only_failure(self):
        row = make_row(with_slots(ADJ_RAW_FAIL), condition=2)
        result = td.describe_translation(row)
        self.assertEqual(result['status'], 'complete')
        self.assertEqual(slot_names(result['triggers']), ['adjacent_raw'])
        line = find_line(result['lines'], '[当前比较] 相邻区间/未平滑分支')
        self.assertIsNotNone(line)
        self.assertNotIn('回退', line)

    def test_adjacent_failure_wins_over_passing_window(self):
        row = make_row(with_slots(ADJ_FINAL_FAIL), selected='adjacent', condition=2)
        result = td.describe_translation(row)
        self.assertEqual(result['selected_window'], 'adjacent')
        roles = {c['slot']: c['role'] for c in result['comparisons']}
        self.assertEqual(roles['window_final'], 'other_interval_pass')
        self.assertEqual(roles['window_raw'], 'other_interval_pass')
        self.assertEqual(slot_names(result['triggers']), ['adjacent_final'])

    def test_window_failure_while_adjacent_passes(self):
        row = make_row(with_slots(WIN_FINAL_FAIL_Y), selected='window', condition=2)
        result = td.describe_translation(row)
        self.assertEqual(result['status'], 'complete')
        self.assertEqual(slot_names(result['triggers']), ['window_final'])
        line = find_line(result['lines'], '[当前比较] 窗口区间/发布分支')
        self.assertIsNotNone(line)
        self.assertIn('Y 实测误差 0.4000 m / 阈值 0.3500 m:超限', line)
        self.assertNotIn('0.3100', line)

    def test_differing_distance_and_thresholds_across_intervals(self):
        row = make_row(with_slots(WIN_FINAL_FAIL_Y), selected='window', condition=2)
        result = td.describe_translation(row, include_passing=True)
        comps = {c['slot']: c for c in result['comparisons']}
        self.assertEqual(comps['adjacent_final']['distance_m'], 0.1)
        self.assertEqual(comps['window_final']['distance_m'], 0.5)
        self.assertEqual(comps['adjacent_final']['limit_y_m'], 0.3 + 0.1 * 0.1)
        self.assertEqual(comps['window_final']['limit_y_m'], 0.3 + 0.1 * 0.5)
        passing = find_line(result['lines'], '[当前比较][通过,非故障] 相邻区间/发布分支')
        self.assertIsNotNone(passing)
        self.assertIn('阈值 0.3100 m', passing)

    def test_both_branches_fail_are_shown_but_correlated(self):
        row = make_row(with_slots(ADJ_FINAL_FAIL, ADJ_RAW_FAIL), condition=2)
        result = td.describe_translation(row)
        self.assertEqual(slot_names(result['triggers']), ['adjacent_final', 'adjacent_raw'])
        self.assertIsNotNone(find_line(result['lines'], '[当前比较] 相邻区间/发布分支'))
        self.assertIsNotNone(find_line(result['lines'], '[当前比较] 相邻区间/未平滑分支'))
        self.assertTrue(any('相关事件' in w for w in result['warnings']))

    def test_other_interval_failure_is_preserved_but_not_trigger(self):
        window_raw_fail = win('window_raw', (0.5, 0.0), (1.2, 0.0))
        row = make_row(with_slots(ADJ_FINAL_FAIL, window_raw_fail), selected='adjacent',
                       condition=2)
        result = td.describe_translation(row)
        self.assertEqual(slot_names(result['triggers']), ['adjacent_final'])
        roles = {c['slot']: c['role'] for c in result['comparisons']}
        self.assertEqual(roles['window_raw'], 'other_interval_failure')
        self.assertTrue(any('非判定区间' in line and '窗口区间/未平滑分支' in line
                            for line in result['lines']))


class EdgeCaseTest(unittest.TestCase):

    def test_exact_equality_passes(self):
        equal = win('window_final', (0.0, 0.0), (0.5, 0.3), limits=(0.5, 0.3))
        row = make_row(with_slots(equal), selected='window')
        result = td.describe_translation(row, include_passing=True)
        self.assertEqual(result['status'], 'complete')
        self.assertEqual(result['triggers'], [])
        comps = {c['slot']: c for c in result['comparisons']}
        self.assertIs(comps['window_final']['fail_x'], False)
        self.assertIs(comps['window_final']['fail_y'], False)
        line = find_line(result['lines'], '[当前比较][通过,非故障] 窗口区间/发布分支')
        self.assertIsNotNone(line)
        self.assertIn('X 实测误差 0.5000 m / 阈值 0.5000 m:未超限', line)
        self.assertIn('Y 实测误差 0.3000 m / 阈值 0.3000 m:未超限', line)

    def test_equality_recorded_as_failure_is_inconsistent(self):
        equal_fail = win('window_final', (0.0, 0.0), (0.5, 0.3), limits=(0.5, 0.3),
                         fail=(True, False))
        row = make_row(with_slots(equal_fail), selected='window', condition=2)
        result = td.describe_translation(row)
        self.assertEqual(result['status'], 'inconsistent')
        self.assertEqual(result['triggers'], [])

    def test_reverse_motion_signed_increments(self):
        reverse = adj('adjacent_final', (-0.1, 0.0), (-0.8, 0.02))
        row = make_row(with_slots(reverse), condition=2)
        result = td.describe_translation(row)
        self.assertEqual(slot_names(result['triggers']), ['adjacent_final'])
        trigger = result['triggers'][0]
        self.assertEqual(trigger['predicted_x_m'], -0.1)
        self.assertEqual(trigger['measured_x_m'], -0.8)
        self.assertEqual(trigger['residual_x_m'], -0.8 - -0.1)
        self.assertLess(trigger['residual_x_m'], 0.0)
        self.assertEqual(trigger['error_x_m'], abs(-0.8 - -0.1))

    def test_raw_fallback_is_labelled(self):
        fallback = adj('adjacent_raw', (0.1, 0.0), (0.75, 0.0), raw_distinct=False)
        row = make_row(with_slots(fallback), condition=2)
        result = td.describe_translation(row)
        self.assertEqual(slot_names(result['triggers']), ['adjacent_raw'])
        self.assertIs(result['triggers'][0]['raw_fallback'], True)
        line = find_line(result['lines'], '[当前比较] 相邻区间/未平滑分支')
        self.assertIn('回退', line)
        self.assertTrue(any('回退' in w for w in result['warnings']))

    def test_latched_failure_during_current_unknown(self):
        row = make_row(with_slots(ADJ_FINAL_FAIL), origin='latched_failure', sequence=41,
                       stamp=T0 + 0.5, condition=0, latched=True)
        result = td.describe_translation(row)
        self.assertEqual(result['status'], 'complete')
        self.assertEqual(result['origin'], 'latched_failure')
        self.assertEqual(result['sequence'], 41)
        self.assertEqual(result['evaluation_stamp'], T0 + 0.5)
        self.assertIn('UNKNOWN', result['lines'][0])
        self.assertIn('历史', result['lines'][0])
        self.assertIsNotNone(find_line(result['lines'], '[锁存的历史失败证据] 相邻区间/发布分支'))
        self.assertFalse(any('[当前比较]' in line for line in result['lines']))

    def test_latched_failure_contradictions(self):
        cases = {
            'current_fail': dict(origin='latched_failure', condition=2, latched=True),
            'current_pass': dict(origin='latched_failure', condition=1, latched=True),
            'invalid_condition': dict(origin='latched_failure', condition=9, latched=True),
            'latch_cleared': dict(origin='latched_failure', condition=0, latched=False),
            'latch_not_shown': dict(origin='current', condition=0, latched=True),
        }
        for name, kwargs in cases.items():
            with self.subTest(name):
                row = make_row(with_slots(ADJ_FINAL_FAIL), **kwargs)
                result = td.describe_translation(row)
                self.assertEqual(result['status'], 'inconsistent')
                self.assertEqual(result['triggers'], [])

    def test_missing_current_condition_is_not_inferred(self):
        row = make_row(with_slots(ADJ_FINAL_FAIL), origin='latched_failure', latched=True)
        result = td.describe_translation(row)
        self.assertIn('未记录', result['lines'][0])
        self.assertNotIn('UNKNOWN', result['lines'][0])

    def test_timestamp_order_cannot_be_hidden_by_epoch_roundoff(self):
        for end in (T0, math.nextafter(T0, -math.inf)):
            row = make_row(with_slots(ADJ_FINAL_FAIL), condition=2)
            prefix = td.slot_prefix('adjacent_final')
            row.update({prefix + 'from_stamp': repr(T0), prefix + 'to_stamp': repr(end),
                        prefix + 'dt_sec': '1e-8'})
            result = td.describe_translation(row)
            self.assertEqual(result['status'], 'inconsistent')
            self.assertEqual(result['triggers'], [])

    def test_inconsistent_flags_do_not_override_recorded_pass(self):
        row = make_row(with_slots(ADJ_FINAL_FAIL), condition=1)
        result = td.describe_translation(row)
        self.assertEqual(result['status'], 'inconsistent')
        self.assertNotIn('当前平移判定:失败', '\n'.join(result['lines']))

    def test_selected_none_is_unknown_not_pass(self):
        row = make_row([adj('adjacent_final', (0.1, 0.0), (0.11, 0.0)),
                        adj('adjacent_raw', (0.1, 0.0), (0.12, 0.0))], selected='none')
        result = td.describe_translation(row)
        self.assertEqual(result['status'], 'partial')
        self.assertEqual(result['triggers'], [])
        text = ''.join(result['lines'])
        self.assertIn('UNKNOWN', text)
        self.assertIn('不视为通过', text)

    def test_passing_details_only_on_request(self):
        row = make_row(passing_set(), selected='window')
        default = td.describe_translation(row)
        self.assertFalse(any('非故障' in line for line in default['lines']))
        detailed = td.describe_translation(row, include_passing=True)
        self.assertTrue(any('[通过,非故障]' in line for line in detailed['lines']))
        self.assertEqual(detailed['triggers'], [])

    def test_timestamp_round_trip_epoch_like(self):
        frm = 1767225600.1234567
        to = math.nextafter(frm + 0.1, math.inf)
        precise = comparison('adjacent_final', frm, to, 0.1, (0.1, 0.0), (0.75, 0.0))
        row = make_row(with_slots(precise), condition=2)
        result = td.describe_translation(row)
        trigger = result['triggers'][0]
        self.assertEqual(trigger['from_stamp'], frm)
        self.assertEqual(trigger['to_stamp'], to)
        self.assertEqual(trigger['dt_sec'], to - frm)
        line = find_line(result['lines'], '[当前比较] 相邻区间/发布分支')
        self.assertIn(repr(frm), line)
        self.assertIn(repr(to), line)


class EvidenceQualityTest(unittest.TestCase):

    def test_legacy_aggregate_only_is_missing(self):
        row = {'error_x': '0.62', 'error_y': '0.33', 'translation_condition': '2'}
        result = td.describe_translation(row)
        self.assertEqual(result['status'], 'missing')
        self.assertEqual(result['comparisons'], [])
        self.assertEqual(result['triggers'], [])
        text = ''.join(result['lines'])
        self.assertIn('阈值', text)
        self.assertIn('无法', text)
        self.assertTrue(any('error_x' in w for w in result['warnings']))

    def test_independent_maxima_cannot_identify_trigger(self):
        # 两组不同的精确证据产生相同的逐轴最大值,但触发项不同。
        aggregate = {'error_x': '0.52', 'error_y': '0.32'}
        row_a = make_row(with_slots(adj('adjacent_final', (0.1, 0.0), (0.15, 0.32)),
                                    win('window_raw', (0.5, 0.0), (1.02, 0.0))),
                         condition=2, extra=aggregate)
        row_b = make_row(with_slots(adj('adjacent_raw', (0.1, 0.0), (0.62, 0.01)),
                                    win('window_final', (0.5, 0.0), (0.52, 0.32))),
                         condition=2, extra=aggregate)
        result_a = td.describe_translation(row_a)
        result_b = td.describe_translation(row_b)
        self.assertEqual(slot_names(result_a['triggers']), ['adjacent_final'])
        self.assertEqual(slot_names(result_b['triggers']), ['adjacent_raw'])
        trigger_a = result_a['triggers'][0]
        self.assertIs(trigger_a['fail_y'], True)
        self.assertIs(trigger_a['fail_x'], False)
        self.assertLess(trigger_a['error_x_m'], 0.1)  # 最大 X 来自另一比较
        legacy_a = td.describe_translation(dict(aggregate, translation_condition='2'))
        legacy_b = td.describe_translation(dict(aggregate, translation_condition='2'))
        self.assertEqual(legacy_a, legacy_b)
        self.assertEqual(legacy_a['status'], 'missing')
        self.assertEqual(legacy_a['triggers'], [])

    def test_missing_field_is_partial_without_trigger(self):
        row = make_row(with_slots(ADJ_FINAL_FAIL), condition=2)
        del row['translation_adjacent_final_limit_x_m']
        result = td.describe_translation(row)
        self.assertEqual(result['status'], 'partial')
        self.assertEqual(result['triggers'], [])
        roles = {c['slot']: c['role'] for c in result['comparisons']}
        self.assertEqual(roles['adjacent_final'], 'unconfirmed')
        self.assertTrue(any('未确认' in line for line in result['lines']))

    def test_nonfinite_required_fields_are_inconsistent(self):
        cases = (
            ('translation_adjacent_final_error_x_m', 'nan'),
            ('translation_adjacent_final_limit_y_m', 'inf'),
            ('translation_adjacent_final_from_stamp', '-inf'),
        )
        for key, value in cases:
            with self.subTest(key):
                row = make_row(with_slots(ADJ_FINAL_FAIL), condition=2)
                row[key] = value
                result = td.describe_translation(row)
                self.assertEqual(result['status'], 'inconsistent')
                self.assertEqual(result['triggers'], [])

    def test_contradictory_slot_fields_are_inconsistent(self):
        p = 'translation_adjacent_final_'
        cases = {
            'fail_flag_without_exceedance': {p + 'fail_x': '1'},
            'residual_mismatch': {p + 'residual_x_m': '0.9'},
            'error_mismatch': {p + 'error_x_m': '0.3'},
            'dt_mismatch': {p + 'dt_sec': '0.2'},
            'nonpositive_limit': {p + 'limit_x_m': '0'},
            'negative_dt': {p + 'from_stamp': repr(T0 + 0.5), p + 'to_stamp': repr(T0 + 0.4),
                            p + 'dt_sec': repr((T0 + 0.4) - (T0 + 0.5))},
            'bad_flag': {p + 'available': '2'},
        }
        for name, changes in cases.items():
            with self.subTest(name):
                row = make_row(passing_set(), selected='window')
                row.update(changes)
                result = td.describe_translation(row)
                self.assertEqual(result['status'], 'inconsistent')
                self.assertEqual(result['triggers'], [])

    def test_contradictory_masks_and_selection_are_inconsistent(self):
        three = passing_set()[:3]
        unavailable_flag = make_row(three, selected='window')
        unavailable_flag['translation_window_raw_fail_x'] = '1'
        unavailable_number = make_row(three, selected='window')
        unavailable_number['translation_window_raw_limit_x_m'] = '0.55'
        cases = {
            'unavailable_with_flag': unavailable_flag,
            'unavailable_with_number': unavailable_number,
            'availability_mismatch': make_row(passing_set(), selected='window', available=False),
            'adjacent_fail_not_selected': make_row(with_slots(ADJ_FINAL_FAIL), selected='window',
                                                   condition=2),
            'fail_condition_without_record': make_row(passing_set(), selected='window',
                                                      condition=2),
        }
        for name, row in cases.items():
            with self.subTest(name):
                result = td.describe_translation(row)
                self.assertEqual(result['status'], 'inconsistent')
                self.assertEqual(result['triggers'], [])

    def test_unknown_or_malformed_version(self):
        for version in ('2', 'abc', '', '1.5', None):
            with self.subTest(version=version):
                row = make_row(passing_set(), selected='window')
                row[td.VERSION_KEY] = version
                result = td.describe_translation(row)
                self.assertEqual(result['status'], 'inconsistent')
                self.assertEqual(result['comparisons'], [])
                self.assertEqual(result['triggers'], [])

    def test_malformed_enums(self):
        for key, value in ((td.SELECTED_KEY, 'adjacnt'), (td.ORIGIN_KEY, 'bogus')):
            with self.subTest(key=key):
                row = make_row(passing_set(), selected='window')
                row[key] = value
                result = td.describe_translation(row)
                self.assertEqual(result['status'], 'inconsistent')
                self.assertEqual(result['triggers'], [])


class ParsingContractTest(unittest.TestCase):

    def test_native_types_match_string_rows(self):
        row = make_row(with_slots(ADJ_FINAL_FAIL), condition=2, latched=True)
        native = {}
        for key, value in row.items():
            if key in (td.SELECTED_KEY, td.ORIGIN_KEY):
                native[key] = value
            elif key in (td.VERSION_KEY, td.SEQUENCE_KEY, td.CONDITION_KEY):
                native[key] = int(value)
            elif key == td.LATCH_KEY or any(key.endswith('_' + s) for s in td.FLAG_SUFFIXES):
                native[key] = value == '1'
            else:
                native[key] = float(value)
        from_strings = td.describe_translation(row)
        from_native = td.describe_translation(native)
        self.assertEqual(from_native['status'], from_strings['status'])
        self.assertEqual(slot_names(from_native['triggers']), slot_names(from_strings['triggers']))
        self.assertEqual(from_native['lines'], from_strings['lines'])

    def test_does_not_mutate_input(self):
        row = make_row(with_slots(ADJ_FINAL_FAIL), condition=2, extra={'error_x': '0.65'})
        snapshot = copy.deepcopy(row)
        td.describe_translation(row, include_passing=True)
        td.extract_translation_evidence(row)
        self.assertEqual(row, snapshot)

    def test_extract_only_schema_fields(self):
        row = make_row(passing_set(), selected='window', condition=1, latched=False,
                       extra={'error_x': '0.1', 'foo': 'bar'})
        extracted = td.extract_translation_evidence(row)
        self.assertIsNot(extracted, row)
        self.assertNotIn('error_x', extracted)
        self.assertNotIn('foo', extracted)
        self.assertIn(td.CONDITION_KEY, extracted)
        self.assertIn(td.LATCH_KEY, extracted)
        self.assertEqual(len(extracted), len(td.SCHEMA_KEYS) + 2)

    def test_schema_column_count_matches_cpp_contract(self):
        # C++ 契约:6 个元数据字段 + 4 个槽位 × 18 个字段,顺序一致。
        self.assertEqual(len(td.SCHEMA_KEYS), 6 + 4 * 18)
        self.assertEqual(td.SCHEMA_KEYS[0], 'translation_evidence_version')
        self.assertEqual(td.SCHEMA_KEYS[6], 'translation_adjacent_final_available')
        self.assertEqual(td.SCHEMA_KEYS[-1], 'translation_window_raw_fail_y')

    def test_classifier_attribution(self):
        self.assertTrue(td.needs_translation_finding({'translation_condition': '2'}))
        self.assertTrue(td.needs_translation_finding({'translation_failure_latched': '1'}))
        yaw_only = {'failure_latched': '1', 'translation_failure_latched': '0',
                    'translation_condition': '0'}
        self.assertFalse(td.needs_translation_finding(yaw_only))
        self.assertIs(td.translation_latch_attribution(yaw_only), False)
        old = {'failure_latched': '1'}
        self.assertFalse(td.needs_translation_finding(old))
        self.assertIsNone(td.translation_latch_attribution(old))


class DiagnosticIntegrationTest(unittest.TestCase):
    def test_unavailable_comparisons_remain_valid_strict_json(self):
        import json
        findings = self.classify(make_row([ADJ_RAW_FAIL], condition=2, latched=True))
        encoded = json.dumps(findings, allow_nan=False)
        self.assertIn('translation_condition', encoded)
        self.assertEqual(findings[0]['translation_detail']['comparisons'][2]['error_x_m'], '非有限值(nan)')

    def classify(self, evidence):
        import diagnose_localization as diagnosis
        row = dict(raw_runstate=4, published_runstate=4, covariance_valid=1,
                   localization_uncertainty=1, localization_consistency=2,
                   yaw_condition=1, confirmed_failure=1)
        row.update(diagnosis.fields(' '.join(f'{k}={v}' for k, v in evidence.items())))
        row['localization_consistency'] = row['translation_condition']
        return diagnosis.classify(row, {})

    def test_default_summary_retains_trigger_and_actual_limits(self):
        import diagnostic_summary as summary
        evidence = make_row(with_slots(ADJ_RAW_FAIL), condition=2, latched=True)
        findings = self.classify(evidence)
        self.assertEqual([f['key'] for f in findings], ['translation_condition'])
        self.assertEqual(findings[0]['translation_detail']['status'], 'complete')
        view = summary.finding_view(findings[0], {}, 'C1')
        lines = []
        summary._finding_compact(view, lines)
        text = '\n'.join(lines)
        for required in ('相邻区间/未平滑分支', '实际时长', '0.6500 m / 阈值 0.5100 m',
                         repr(T0 + 0.4), repr(T0 + 0.5)):
            self.assertIn(required, text)

    def test_translation_latch_has_specific_cause_without_duplicate_generic_fault(self):
        findings = self.classify(make_row(with_slots(ADJ_RAW_FAIL), origin='latched_failure',
                                          condition=0, latched=True))
        self.assertEqual([f['key'] for f in findings], ['translation_condition'])
        self.assertIn('历史', findings[0]['title'])
        self.assertEqual(findings[0]['translation_detail']['origin'], 'latched_failure')

    def test_old_log_does_not_invent_dynamic_bounds(self):
        import diagnostic_summary as summary
        findings = self.classify({'translation_condition': 2, 'error_x': .713239})
        self.assertEqual(findings[0]['translation_detail']['status'], 'missing')
        text = '\n'.join(summary.measures('translation_condition', findings[0]['evidence'],
                                         {'longitudinal_tolerance_m': .5}))
        self.assertIn('均未记录', text)
        self.assertNotIn('阈值 0.5', text)


if __name__ == '__main__':
    unittest.main()
