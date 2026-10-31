#!/usr/bin/env python3
"""Evidence integrity and diagnosis regressions; no ROS or vehicle access."""
import json
from pathlib import Path
import tempfile
import unittest

import compute_profile_report as profile
import diagnose_localization as diagnosis


class ComputeProfileReportTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.run = Path(self.tmp.name)
        (self.run / 'results').mkdir()
        (self.run / 'logs').mkdir()
        (self.run / 'run_metadata.txt').write_text('compute_profile=1\ncompute_profile_storage=dedicated-v1\n')
        self.path = self.run / profile.PROFILE

    def tearDown(self):
        self.tmp.cleanup()

    def evidence(self, sequence=0, **counters):
        footer = dict(received=1, written=1, dropped=0, critical_dropped=0, io_errors=0, truncated=0)
        footer.update(counters)
        text = profile.HEADER.decode()
        text += f'I1031 10:00:00.000000 1 driver.cc:1] COMPUTE_BENCH_FRAME module=lidar_loc frame_id=42 profile_seq={sequence}\n'
        text += 'COMPUTE_PROFILE_SUMMARY ' + ' '.join(f'{k}={v}' for k, v in footer.items()) + '\n'
        self.path.write_bytes(text.encode())

    def test_complete_file_records_module_counts_and_sequential_integrity(self):
        self.evidence()
        result = profile.validate(self.run, True)
        self.assertEqual(result['status'], 'complete')
        self.assertEqual(result['modules'], {'lidar_loc': 1})
        self.assertEqual(result['records'], 1)

    def test_loss_errors_truncation_and_counter_mismatch_cannot_claim_completeness(self):
        for change, expected in [({'received': 2, 'dropped': 1}, 'dropped'),
                                 ({'io_errors': 1}, 'io_errors'), ({'truncated': 1}, 'truncated'),
                                 ({'written': 2}, 'counter_mismatch')]:
            with self.subTest(change=change):
                self.evidence(**change)
                result = profile.validate(self.run, True)
                self.assertEqual(result['status'], 'incomplete')
                self.assertIn(expected, result['problems'])

    def test_missing_footer_after_crash_and_sequence_gaps_are_visible(self):
        self.evidence(sequence=2)
        self.assertIn('sequence_gaps_or_partial_records', profile.validate(self.run, True)['problems'])
        self.path.write_text(self.path.read_text().split('COMPUTE_PROFILE_SUMMARY')[0])
        self.assertIn('missing_or_invalid_footer', profile.validate(self.run, True)['problems'])

    def test_post_close_write_error_is_not_hidden_by_a_complete_file_footer(self):
        self.evidence()
        (self.run / 'logs/run_loc_online.stderr.log').write_text(
            'I1031 10:00:00.000000 1 driver.cc:1] COMPUTE_PROFILE_STATUS path=test COMPUTE_PROFILE_SUMMARY '
            'received=1 written=1 dropped=0 critical_dropped=0 io_errors=1 truncated=0\n')
        report = profile.validate(self.run, True)
        self.assertEqual(report['status'], 'incomplete')
        self.assertIn('post_close_counter_mismatch', report['problems'])

    def test_missing_file_and_disabled_sampling_are_distinct(self):
        self.assertEqual(profile.validate(self.run, True)['status'], 'missing')
        self.assertEqual(profile.validate(self.run, False)['status'], 'disabled')
        self.assertEqual(profile.inspect(self.run)['status'], 'missing')

    def test_recorded_sampling_override_is_respected_in_either_mode(self):
        for mode in ('diagnostic', 'production'):
            for enabled in (False, True):
                with self.subTest(mode=mode, enabled=enabled):
                    (self.run / 'run_metadata.txt').write_text(
                        f'run_mode={mode}\ncompute_profile={int(enabled)}\ncompute_profile_storage=dedicated-v1\n')
                    result = profile.inspect(self.run)
                    self.assertEqual(result['enabled'], enabled)
                    self.assertEqual(result['status'], 'missing' if enabled else 'disabled')

    def test_cached_validation_is_invalidated_if_evidence_changes(self):
        self.evidence()
        (self.run / profile.STATUS).write_text(json.dumps(profile.validate(self.run, True)))
        self.assertEqual(profile.inspect(self.run)['status'], 'complete')
        with self.path.open('a') as stream:
            stream.write('partial')
        self.assertEqual(profile.inspect(self.run)['status'], 'unverified')
        self.assertEqual(profile.inspect(self.run, live=True)['status'], 'recording')
        for invalid in ('[]', '{"schema": 1, "enabled": true}'):
            (self.run / profile.STATUS).write_text(invalid)
            self.assertEqual(profile.inspect(self.run)['status'], 'unverified')

    def test_legacy_package_is_not_mistaken_for_new_dedicated_storage(self):
        (self.run / 'run_metadata.txt').write_text('compute_profile=1\n')
        self.assertEqual(profile.inspect(self.run)['status'], 'legacy_unverified')

    def test_diagnosis_reports_missing_performance_evidence_without_inventing_algorithm_fault(self):
        report = diagnosis.analyze(self.run, None, 'explicit', [])
        self.assertEqual(report['compute_profile']['status'], 'missing')
        self.assertEqual(report['findings'], [])
        self.assertTrue(any('逐帧性能证据' in item for item in report['warnings']))
        self.assertIn('逐帧性能采样: missing', diagnosis.render(report))


if __name__ == '__main__':
    unittest.main()
