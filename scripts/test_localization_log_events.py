"""Chinese/legacy diagnosis parity and stable event-code boundaries."""
import copy
import datetime as dt
from pathlib import Path
import tempfile
import unittest

import analyze_localization_health_log as health
import diagnose_localization as d
from localization_log_events import any_event, has_event, LEGACY_EVENTS


class LocalizedLogTests(unittest.TestCase):
    def test_every_legacy_event_and_arbitrary_chinese_prose(self):
        for code, markers in LEGACY_EVENTS.items():
            with self.subTest(code=code):
                self.assertTrue(has_event(' '.join(markers), code))
                self.assertTrue(has_event(f'event={code} 中文说明可以调整 field=1', code))
                self.assertTrue(any_event(f'event={code} 不依赖正文措辞', (code,)))

    def test_explicit_codes_take_precedence_and_do_not_match_prefixes(self):
        self.assertFalse(has_event('event=LOC_MATCH_FAILED_EXTRA localization failed!', 'LOC_MATCH_FAILED'))
        self.assertFalse(has_event('event=UNKNOWN_EVENT localization failed!', 'LOC_MATCH_FAILED'))
        self.assertFalse(has_event('not_event=LOC_MATCH_FAILED 定位失败', 'LOC_MATCH_FAILED'))
        self.assertTrue(has_event('event=periodic localization failed!', 'LOC_MATCH_FAILED'))
        self.assertFalse(has_event('GLOBAL_RELOCALIZATION[solid] rejected: reason=none', 'LOC_RELOCALIZATION_ACCEPTED'))

    def test_diagnosis_preserves_faults_boundaries_protections_and_priority(self):
        bodies = [
            'SAFETY_DIAGNOSTIC sequence=1 raw_runstate=0 published_runstate=0',
            'localization failed! score: 2.28 map_checked=1 map_overlap=0.72 map_min_overlap=0.75',
            'sensor queue lag reached 0.7 sec; keeping fresh lidar admissible',
            'drop stale lidar correction before matching: age=0.68 sec, hard_deadline=0.3',
            'adaptive lidar load step 0 -> 2, processing_ms=45 latency_ms=680 tracking_healthy=0',
            'SAFETY_DIAGNOSTIC sequence=2 raw_runstate=-9 published_runstate=-9 output_stale=1 localization_lost=1 output_age_sec=0.68',
            'localization pose and map outputs disabled after rejected localization: consecutive_lost_frames=4',
            'GLOBAL_RELOCALIZATION[solid_kiss] accepted: candidate=7 overlap=0.98',
            'localization pose and map outputs recovered after a fresh valid lidar match',
            'sensor queue lag recovered to 0 sec; resuming lidar input',
            'SAFETY_DIAGNOSTIC sequence=3 raw_runstate=0 published_runstate=0 localization_lost=0',
        ]
        translations = {
            1: 'event=LOC_MATCH_FAILED 匹配未通过 score: 2.28 map_checked=1 map_overlap=0.72 map_min_overlap=0.75',
            2: 'event=LOC_SENSOR_QUEUE_LAG 输入处理落后 lag_sec=0.7 sec',
            3: 'event=LIO_STALE_CORRECTION_DROPPED 丢弃过期校正 age=0.68 sec, hard_deadline=0.3',
            4: 'event=LIO_ADAPTIVE_LOAD_CHANGED 调整雷达负载 from_step=0 to_step=2 processing_ms=45 latency_ms=680 tracking_healthy=0',
            6: 'event=LOC_OUTPUT_DISABLED 关闭输出 consecutive_lost_frames=4',
            7: 'GLOBAL_RELOCALIZATION[solid_kiss] event=LOC_RELOCALIZATION_ACCEPTED 接受重定位 candidate=7 overlap=0.98',
            8: 'event=LOC_OUTPUT_RECOVERED 恢复输出',
            9: 'event=LOC_SENSOR_QUEUE_RECOVERED 处理进度恢复 lag_sec=0 sec',
        }
        stamps = (0, 0.98, 0.99, 0.991, 0.992, 1, 1.01, 1.4, 1.41, 1.42, 1.5)
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            (run / 'logs').mkdir()
            (run / 'config.yaml').write_text('diagnostics: {}\n', encoding='utf-8')
            (run / 'run_metadata.txt').write_text(f'run_dir={run}\nstarted_at=2026-10-31T12:00:00+00:00\n')
            log = run / 'logs/run_loc_online.stderr.log'
            reports = []
            for language in ('legacy', 'chinese', 'mixed'):
                lines = []
                for index, original in enumerate(bodies):
                    body = translations.get(index, original) if language == 'chinese' or language == 'mixed' and index % 2 else original
                    lines.append(f'I1031 12:00:{stamps[index]:09.6f} 42 fixture.cc:1] {body}')
                log.write_text('\n'.join(lines) + '\n', encoding='utf-8')
                # Keep evidence byte offsets anchored to the actual UTF-8 file.
                clues, warnings = [], []
                d.log_rows(log, d.MIB, warnings, reference=dt.datetime(2026, 10, 31, tzinfo=dt.timezone.utc), context=clues)
                raw = log.read_bytes()
                for clue in clues:
                    self.assertTrue(raw[clue['byte_offset']:].decode('utf-8').startswith(clue['text']))
                reports.append(d.analyze(run, None, '指定路径', []))

        def comparable(report):
            def clean(value):
                if isinstance(value, dict):
                    return {k: clean(v) for k, v in value.items()
                            if k not in ('text', 'source', 'onset_source', 'first_source', 'last_source',
                                         'byte_offset', '_source', '_scope')}
                if isinstance(value, list):
                    return [clean(v) for v in value]
                return value
            return clean(copy.deepcopy(report['findings'])), clean(copy.deepcopy(report['runtime_protections']))

        expected = comparable(reports[0])
        self.assertTrue(reports[0]['findings'])
        self.assertTrue(any('处理进度落后 0.7s' in s['text'] for s in reports[0]['priority_summary']))
        for report in reports[1:]:
            self.assertEqual(comparable(report), expected)
            self.assertEqual([s['text'] for s in report['priority_summary']],
                             [s['text'] for s in reports[0]['priority_summary']])
        steps = next(p for p in reports[1]['runtime_protections'] if p['key'] == 'adaptive_lidar_load_step')
        self.assertEqual((steps['evidence']['from_step'], steps['evidence']['to_step']), (0, 2))
        self.assertIsInstance(steps['evidence']['from_step'], int)
        self.assertIsInstance(steps['evidence']['to_step'], int)

    def test_invalid_or_partial_load_steps_do_not_become_numeric_evidence(self):
        for invalid in ('bad', 'nan', '-1', '1.5'):
            with self.subTest(invalid=invalid):
                protections = {}
                d.record_runtime_protection(protections, 'adaptive_lidar_load_step',
                    'event=LIO_ADAPTIVE_LOAD_CHANGED 负载变化',
                    dict(_time='t', _source='fixture', from_step=invalid, to_step=2, processing_ms=45))
                evidence = protections['adaptive_lidar_load_step']['evidence']
                self.assertNotIn('from_step', evidence)
                self.assertNotIn('to_step', evidence)
                self.assertEqual(evidence['processing_ms'], 45)

    def test_health_analyzer_accepts_chinese_events_and_preserves_audit_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory) / 'sany_20261031_fixture'
            (run / 'logs').mkdir(parents=True)
            bodies = [
                'event=LOC_MATCH_FAILED 激光匹配失败 map_overlap=0.72',
                'GLOBAL_RELOCALIZATION[solid] event=LOC_RELOCALIZATION_ACCEPTED 重定位接受 candidate=7',
                'event=LOC_TRACKING_LOST 跟踪丢失',
                'event=LOC_OUTPUT_DISABLED 停止输出',
                'event=SMOOTHER_MOTION_EXCESSIVE 平滑增量过大 distance=2',
                'LIDAR_LOC_AUDIT confidence=3 published_valid=1 guess_xyz=[0 0 0] ndt_xyz=[0.1 0 0] balanced_xyz=[0.1 0 0]',
                'PGO_OUTPUT_AUDIT output_lidar_delta_xy=0.02',
            ]
            (run / 'logs/run_loc_online.stderr.log').write_text(''.join(
                f'I1031 12:00:{i:02}.000000 42 fixture.cc:1] {body}\n' for i, body in enumerate(bodies)), encoding='utf-8')
            result = health.analyze(run)
            self.assertEqual(result['sampled_matches'], 1)
            self.assertAlmostEqual(result['max_raw_correction']['raw_correction_m'], 0.1)
            self.assertEqual(len(result['events']), 4)
            self.assertEqual(result['smoother_motion_rejections'], 1)


if __name__ == '__main__':
    unittest.main()
