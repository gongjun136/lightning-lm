#!/usr/bin/env python3
"""Selection and causal diagnosis regressions; never start/stop localization."""
import argparse
import csv
import datetime as dt
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import diagnose_localization as d


class DiagnosisTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name).resolve()
        self.run = self.root / 'run_new'
        (self.run / 'logs').mkdir(parents=True)
        (self.run / 'config.yaml').write_text('diagnostics:\n  major_sigma_limit_m: 0.5\n  lateral_sigma_limit_m: 0.5\n  yaw_sigma_limit_deg: 5\n  consistency_unknown_grace_sec: 0.15\n')
        (self.run / 'run_metadata.txt').write_text(f'run_dir={self.run}\nstarted_at=2026-10-28T12:00:00+00:00\n')
        self.log = self.run / 'logs/run_loc_online.stderr.log'
        self.log.write_text('')
        self.archive = self.root / 'incidents'
        self.archive.mkdir()
        self.p = dict(pid=42, comm='loc_executor', exe='/install/run_loc_online',
                      run=str(self.run), config=str(self.run / 'config.yaml'), log=str(self.log),
                      archive=str(self.archive), start_ticks=100)

    def tearDown(self):
        self.tmp.cleanup()

    def event(self, ns, pid=42, run=None):
        event = self.archive / f'event_{ns}_{pid}'
        event.mkdir()
        (event / 'summary.json').write_text(json.dumps(dict(trigger_wall_ns=ns, status='complete', reasons=['raw_4'])))
        run = self.run if run is None else run
        (event / 'run_metadata.txt').write_text(f'run_dir={run}\nstarted_at=2026-10-28T12:00:00+00:00\n')
        (event / 'config.yaml').write_text((self.run / 'config.yaml').read_text())
        return event

    def args(self, **kw):
        return argparse.Namespace(**(dict(path=None, pid=None, root=None) | kw))

    @unittest.skipUnless(os.name == 'posix', 'Linux /proc symlink fixture')
    def test_renamed_comm_and_argv_still_identified_by_exe(self):
        proc = self.root / 'proc'
        for pid, exe in [(42, '/install/run_loc_online (deleted)'), (43, '/usr/bin/python3')]:
            p = proc / str(pid)
            (p / 'fd').mkdir(parents=True)
            (p / 'exe').symlink_to(exe)
            (p / 'cwd').symlink_to(self.root)
            (p / 'fd/2').symlink_to(self.log)
            (p / 'comm').write_text('loc_executor\n')
            (p / 'cmdline').write_bytes(b'custom_name\0--config\0run_new/config.yaml\0--output_tum=run_new/results/a.tum\0')
            (p / 'stat').write_text(f'{pid} (loc_executor) ' + ' '.join(['S'] + ['0'] * 18 + ['123']))
        found = d.find_processes(proc)
        self.assertEqual([p['pid'] for p in found], [42])
        self.assertEqual(found[0]['config'], str(self.run / 'config.yaml'))
        self.assertEqual(found[0]['archive'], str(self.run / 'results/incidents'))

    def test_live_selects_own_event_not_newer_other_run(self):
        own = self.event(100)
        self.event(200, 43, self.root / 'other_run')
        got, p, _, _ = d.select(self.args(), [self.p])
        self.assertEqual(got, own)
        self.assertEqual(p['comm'], 'loc_executor')

    def test_live_without_event_uses_current_log(self):
        self.event(300, 43, self.root / 'old')
        got, _, _, warnings = d.select(self.args(), [self.p])
        self.assertEqual(got, self.run)
        self.assertIn('尚无', warnings[0])

    def test_reused_pid_does_not_select_previous_lifetime(self):
        self.event(100)
        self.p['start_wall_ns'] = 200
        self.assertEqual(d.select(self.args(), [self.p])[0], self.run)

    def test_no_process_defaults_to_newest_event(self):
        # Keep metadata run timestamp below events so the most recent event wins.
        older = self.event(1793190000000000000)
        newer = self.event(1793191000000000000)
        os.utime(older, (2000000000, 2000000000))
        got, p, _, _ = d.select(self.args(), [], [self.root])
        self.assertEqual(got, newer)  # trigger time, not copying/directory mtime
        self.assertIsNone(p)

    def test_explicit_event_is_not_replaced_by_live(self):
        event = self.event(100)
        self.event(200)
        self.assertEqual(d.select(self.args(path=str(event / 'summary.json')), [self.p])[0], event)

    def test_explicit_run_only_its_own_events(self):
        own = self.event(100)
        self.event(200, 99, self.root / 'other')
        self.assertEqual(d.choose_path(self.run), own)

    def test_copied_run_prefers_its_local_log_over_original_absolute_path(self):
        (self.run / 'run_metadata.txt').write_text('run_dir=/unavailable/original/run\n')
        self.log.write_text('I1028 21:10:58.502960 42] SAFETY_DIAGNOSTIC raw_runstate=4 published_runstate=4 localization_consistency=2 yaw_condition=2 error_yaw_rad=.04\n')
        report = d.analyze(self.run, None, '指定路径', [])
        self.assertEqual(report['findings'][0]['key'], 'yaw_condition')
        self.assertIn(str(self.log), report['findings'][0]['source'])

    def test_pid_selector_rejects_missing_instance(self):
        with self.assertRaisesRegex(ValueError, 'PID'):
            d.select(self.args(pid=99), [self.p])

    def test_default_root_parses_workspace_without_execution(self):
        config = self.root / '.config/sany'
        config.mkdir(parents=True)
        (config / 'workspace.bash').write_text('export SANY_WS="$HOME/vehicle"\nexport SANY_WS="${SANY_WS%/}"\nexport LIGHTNING_LM_OUT_ROOT="$SANY_WS/runs"\ntouch /tmp/should_not_execute\n')
        self.assertEqual(d.default_roots({}, self.root)[0], self.root / 'vehicle/runs')

    def causes(self, **kwargs):
        row = dict(raw_runstate=4, published_runstate=4, covariance_valid=1,
                   pose_unobservable=0, localization_uncertainty=1,
                   localization_consistency=1, yaw_condition=1, translation_condition=1,
                   major_sigma=.24, lateral_sigma=.24, yaw_sigma_rad=.024)
        row.update(kwargs)
        return d.classify(row, d.thresholds(self.run / 'config.yaml'))

    def test_yaw_motion_failure_not_uncertainty(self):
        out = self.causes(localization_consistency=2, yaw_condition=2, error_yaw_rad=.0428282, confirmed_failure=1)
        self.assertEqual([f['key'] for f in out], ['yaw_condition'])
        self.assertAlmostEqual(out[0]['evidence']['error_yaw_deg'], 2.45387595, places=4)

    def test_unknown_grace_expired_without_input_silence(self):
        out = self.causes(localization_consistency=0, yaw_condition=0, translation_condition=0,
                          unknown_tolerated=0, optional_can_unavailable=0, yaw_pass_age_sec=.156812,
                          imu_arrival_age_sec=.001, imu_stamp_lag_sec=-.045)
        self.assertEqual(out[0]['key'], 'consistency_unknown')
        self.assertIn('0.156812s > 0.15s', out[0]['title'])

    def test_sigma_converts_degrees(self):
        out = self.causes(localization_uncertainty=2, yaw_sigma_rad=.10)
        self.assertEqual(out[0]['key'], 'yaw_sigma_rad')
        self.assertIn('5.72958 > 5', out[0]['title'])

    def test_unobservable_not_guessed_sigma(self):
        out = self.causes(localization_uncertainty=2, pose_unobservable=1, covariance_valid=0)
        self.assertEqual([f['key'] for f in out], ['unobservable'])

    def test_missing_threshold_does_not_invent_exceedance(self):
        out = d.classify(dict(raw_runstate=4, covariance_valid=1, major_sigma=10), {})
        self.assertEqual(out[0]['confidence'], '证据不足')

    def test_masked_raw_fault_is_still_analyzed(self):
        out = self.causes(localization_uncertainty=2, major_sigma=.6, published_runstate=0)
        self.assertEqual(out[0]['published'], 0)
        self.assertEqual(out[0]['key'], 'major_sigma')

    def test_severe_and_radar_mapping(self):
        severe = d.classify(dict(raw_runstate=-9, output_stale=1, output_age_sec=.4), {})
        self.assertEqual(severe[0]['key'], 'output_stale')
        radar = d.classify(dict(raw_runstate=2, time_sync_fault_mask=12, non_sync_fault=0, non_sync_all_pass=1), {})
        self.assertIn('left, right', radar[0]['title'])

    def test_bounded_log_offsets_and_nan_safe(self):
        self.log.write_bytes(('x' * 100 + '\nI1028 21:10:58.502960 42 loc.cc:1] SAFETY_DIAGNOSTIC raw_runstate=4 published_runstate=4 major_sigma=nan\n').encode())
        warnings = []
        rows, _ = d.log_rows(self.log, 150, warnings)
        self.assertEqual(len(rows), 1)
        self.assertIn('byte=101', rows[0]['_source'])
        json.dumps(rows, allow_nan=False)

    def test_prefers_modern_audit_over_duplicate_shadow(self):
        self.log.write_text('I1028 21:10:58.502960 42] SAFETY_TIME_SYNC_SHADOW published_runstate=4\nI1028 21:10:58.502960 42] SAFETY_DIAGNOSTIC raw_runstate=4 published_runstate=0\n')
        rows, _ = d.log_rows(self.log, d.MIB, [])
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['published_runstate'], 0)

    def test_legacy_shadow_keeps_unknown_raw_warning(self):
        self.log.write_text('I1028 21:10:58.502960 42] SAFETY_TIME_SYNC_SHADOW published_runstate=4 localization_uncertainty=2\n')
        rows, _ = d.log_rows(self.log, d.MIB, [])
        self.assertIn('原始码未知', rows[0]['_scope'])

    def test_partial_csv_line_and_standalone_detailed_evidence(self):
        event = self.event(100)
        (event / 'timeline.csv').write_text('wall_ns,raw,published,consistency,uncertainty,yaw_condition,error_yaw_rad\n1793193058502960000,4,0,2,1,2,0.0428282\n1793193058,4')
        report = d.analyze(event, None, '指定路径', [])
        self.assertEqual(report['findings'][0]['key'], 'yaw_condition')
        self.assertEqual(report['findings'][0]['published'], 0)
        self.assertTrue(any('尚未写完' in w for w in report['warnings']))

    def test_explicit_event_excludes_unrelated_log_time(self):
        ns = int(dt.datetime(2026, 10, 28, 21, 10, 58, tzinfo=dt.timezone.utc).timestamp() * 1e9)
        event = self.event(ns)
        (event / 'timeline.csv').write_text(f'wall_ns,raw,published,consistency,uncertainty,yaw_condition\n{ns},0,0,1,1,1\n')
        self.log.write_text('I1028 19:44:32.114810 42] SAFETY_DIAGNOSTIC raw_runstate=4 published_runstate=4 localization_uncertainty=2\n')
        report = d.analyze(event, None, '指定路径', [])
        self.assertEqual(report['findings'], [])
        self.assertEqual(report['log_samples'], 0)

    def test_collector_startup_label_conflict_is_not_algorithm_fault(self):
        event = self.event(100)
        (event / 'summary.json').write_text(json.dumps(dict(status='partial', reasons=['startup_without_output'])))
        (event / 'timeline.csv').write_text('wall_ns,raw,published,publish_end_age\n1793193058502960000,0,0,0.02\n')
        report = d.analyze(event, None, '指定路径', [])
        self.assertEqual(report['findings'], [])
        self.assertTrue(any('并存' in w for w in report['warnings']))

    def test_no_evidence_is_not_reported_as_healthy(self):
        event = self.event(100)
        with patch.object(d, 'bag_rows', return_value=[]):
            report = d.analyze(event, None, '指定路径', [])
        self.assertFalse(report['has_safety_evidence'])
        self.assertTrue(any('不能给出故障原因' in w for w in report['warnings']))


if __name__ == '__main__':
    unittest.main()
