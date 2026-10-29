#!/usr/bin/env python3
"""Selection and causal diagnosis regressions; never start/stop localization."""
import argparse
import csv
import datetime as dt
import json
import io
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import diagnose_localization as d
import check_diagnostic_contract as contract


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
        self.assertEqual([f['device'] for f in radar], ['left', 'right'])

    def test_advisory_sync_warning_is_not_blamed_for_other_radar_fault(self):
        found = d.classify(dict(raw_runstate=2, time_sync_fault_mask=4,
                               bounded_policy_enabled=1, bounded_sync_warning_mask=4,
                               non_sync_fault=1, non_sync_all_pass=0), {})
        self.assertEqual([f['key'] for f in found], ['radar_other'])

    def test_confirmed_failure_latch_precedes_unknown_explanation(self):
        found = self.causes(localization_consistency=0, yaw_condition=0, translation_condition=0,
                            confirmed_failure=1, unknown_tolerated=0, optional_can_unavailable=0)
        self.assertEqual([f['key'] for f in found], ['latched_failure'])

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

    def log_row(self, second, code, **fields):
        return dict(raw_runstate=code, published_runstate=code, _time=f't={second}',
                    _wall_ns=second * 10**9, _source=f'log:{second}',
                    _scope='运行日志（有界尾部）', **fields)

    def test_each_radar_and_recurrence_survive_aggregation(self):
        rows = [self.log_row(1, 2, time_sync_fault_mask=4), self.log_row(2, 0),
                self.log_row(3, 2, time_sync_fault_mask=1), self.log_row(4, 0),
                self.log_row(5, 2, time_sync_fault_mask=8), self.log_row(6, 0),
                self.log_row(7, 2, time_sync_fault_mask=4)]
        found = d.fault_episodes(rows, {})
        self.assertEqual([f['device'] for f in found], ['left', 'front', 'right', 'left'])
        self.assertEqual(found[0]['state'], 'recovered')
        self.assertEqual(found[-1]['state'], 'open')
        self.assertEqual(found[0]['observed_duration_sec'], 1)

    def test_startup_and_runtime_output_loss_remain_distinct(self):
        rows = [self.log_row(1, -9, output_stale=1, output_age_sec='inf'),
                self.log_row(2, 0), self.log_row(3, -9, output_stale=1, output_age_sec=.77),
                self.log_row(4, 4, localization_uncertainty=0)]
        found = d.fault_episodes(rows, {})
        self.assertEqual([f.get('phase') for f in found[:2]], ['startup', 'runtime'])
        self.assertIn('启动', found[0]['title'])
        self.assertEqual(found[1]['state'], 'predicate_changed')
        self.assertEqual(found[1]['next_raw'], 4)

    def test_peak_sample_keeps_first_joint_input_stall(self):
        rows = [self.log_row(1, -9, output_stale=1, output_age_sec=.524,
                             imu_arrival_age_sec=.44, can_arrival_age_sec=.443),
                self.log_row(2, -9, output_stale=1, output_age_sec=.77,
                             imu_arrival_age_sec=.001, can_arrival_age_sec=.01),
                self.log_row(3, 0)]
        found = d.fault_episodes(rows, {'max_input_gap_sec': .1})[0]
        self.assertEqual(found['onset_evidence']['imu_arrival_age_sec'], .44)
        self.assertEqual(found['evidence']['output_age_sec'], .77)
        self.assertEqual(found['samples'], 2)

    def test_newer_timeline_wins_over_older_log(self):
        event = self.event(1793193058502960000)
        (event / 'timeline.csv').write_text('wall_ns,raw,published\n1793193058502960000,0,0\n')
        with patch.object(d, 'log_rows', return_value=([self.log_row(1, 4)], {})):
            report = d.analyze(event, self.p, '运行中的定位实例', [])
        self.assertEqual(report['latest']['raw_runstate'], 0)
        self.assertEqual(report['current_findings'], [])

    def test_unknown_timestamps_use_last_record_not_first(self):
        rows = [dict(raw_runstate=4, published_runstate=4, _time='unknown', _scope='log'),
                dict(raw_runstate=0, published_runstate=0, _time='unknown', _scope='log')]
        with patch.object(d, 'log_rows', return_value=(rows, {})):
            report = d.analyze(self.run, None, '指定路径', [])
        self.assertEqual(report['latest']['raw_runstate'], 0)

    def test_missing_capture_is_reported_without_claiming_never_recorded(self):
        event = self.event(100 * 10**9)
        s = d.summary_at(event)
        s['streams'] = [dict(type='geosun_msgs/msg/Safety', written=20,
                             first_relative_sec=-30, last_relative_sec=10)]
        (event / 'summary.json').write_text(json.dumps(s))
        with patch.object(d, 'log_rows', return_value=([
                self.log_row(90, -9, output_stale=1, output_age_sec=.4), self.log_row(91, 0),
                self.log_row(120, -9, output_stale=1, output_age_sec=.77)], {})):
            report = d.analyze(event, self.p, '运行中的定位实例', [])
        self.assertEqual(report['findings'][0]['capture_coverage'], 'range_overlap')
        self.assertEqual(report['findings'][1]['capture_coverage'], 'no_retained_overlap')
        self.assertTrue(any('不代表从未采集' in w for w in report['warnings']))

    def test_scope_intervals_do_not_double_count_samples(self):
        a = self.log_row(1, 4, yaw_condition=2, error_yaw_rad=.04)
        b = dict(a, _scope='事件时间线', _source='csv:2')
        found = d.fault_episodes([b, a], {})
        self.assertEqual(len(found), 2)
        self.assertEqual([f['samples'] for f in found], [1, 1])

    def test_compiled_contract_is_not_inferred_from_source_metadata(self):
        warnings = []
        unknown = d.runtime_contract({}, warnings)
        self.assertEqual(unknown['status'], 'unverified')
        self.assertTrue(any('只代表启动时源码' in w for w in warnings))
        expected = json.loads(d.CONTRACT_PATH.read_text())['contract_id']
        self.assertEqual(d.runtime_contract({'_build': {'contract_id': expected, 'schema': 2}}, [])['status'], 'matched')
        self.assertEqual(d.runtime_contract({'_build': {'contract_id': expected, 'schema': 99}}, [])['status'], 'mismatch')
        self.assertEqual(d.runtime_contract({'_build': {'contract_id': 'old'}}, [])['status'], 'mismatch')

    def test_log_reads_compiled_version_and_relocalization_clue(self):
        self.log.write_text('I1028 21:10:58.002000 42] LOCALIZATION_DIAGNOSTIC_BUILD schema=2 contract_id=abc source_revision=123\n'
                            'I1028 21:10:58.102000 42] GLOBAL_RELOCALIZATION[solid_kiss] accepted: candidate=221\n'
                            'I1028 21:10:58.502960 42] SAFETY_DIAGNOSTIC raw_runstate=0 published_runstate=0\n')
        context = []
        rows, limits = d.log_rows(self.log, d.MIB, [], dt.datetime(2026, 10, 28, tzinfo=dt.timezone.utc), context)
        self.assertEqual(limits['_build']['contract_id'], 'abc')
        self.assertEqual(len(rows), 1)
        self.assertEqual(len(context), 1)
        self.assertIn('accepted', context[0]['text'])

    def test_color_defaults_respect_terminal_and_environment(self):
        terminal = io.StringIO()
        terminal.isatty = lambda: True
        self.assertTrue(d.color_enabled('auto', terminal, {}))
        self.assertFalse(d.color_enabled('auto', terminal, {'NO_COLOR': ''}))
        self.assertFalse(d.color_enabled('auto', terminal, {'TERM': 'dumb'}))
        self.assertFalse(d.color_enabled('auto', io.StringIO(), {}))
        self.assertTrue(d.color_enabled('always', io.StringIO(), {'NO_COLOR': '1'}))
        self.assertFalse(d.color_enabled('never', terminal, {}))

    def test_json_remains_plain_even_with_forced_color(self):
        self.log.write_text('I1028 21:10:58.502960 42] SAFETY_DIAGNOSTIC raw_runstate=0 published_runstate=0\n')
        with patch.object(d, 'low_impact', return_value={}), patch('sys.stdout', new_callable=io.StringIO) as out:
            status = d.main([str(self.run), '--json', '--color', 'always'])
            raw = out.getvalue()
        self.assertEqual(status, 0)
        self.assertNotIn('\x1b', raw)
        self.assertEqual(json.loads(raw)['schema'], 2)

    def test_text_color_and_truncation_are_explicit(self):
        rows = []
        for i in range(15):
            rows += [self.log_row(2*i+1, 2, time_sync_fault_mask=1), self.log_row(2*i+2, 0)]
        with patch.object(d, 'log_rows', return_value=(rows, {})):
            report = d.analyze(self.run, None, '指定路径', [])
        plain = d.render(report)
        colored = d.render(report, color=True)
        full = d.render(report, all_events=True)
        self.assertNotIn('\x1b', plain)
        self.assertIn('\x1b[1;32m', colored)
        self.assertIn('--all', plain)
        self.assertEqual(plain.count('历史记录 [2]'), 12)
        self.assertEqual(full.count('历史记录 [2]'), 15)

    def test_reviewed_source_contract_is_current(self):
        self.assertEqual(contract.check(), json.loads(d.CONTRACT_PATH.read_text())['contract_id'])

    def test_unknown_future_code_is_not_silently_healthy(self):
        finding = d.classify(dict(raw_runstate=7, published_runstate=7), {})[0]
        self.assertEqual(finding['key'], 'unknown_code')
        self.assertEqual(finding['confidence'], '规则缺失')

    def test_log_gap_does_not_invent_continuous_fault_or_recovery(self):
        a = self.log_row(1, 2, time_sync_fault_mask=1, sequence=1)
        b = self.log_row(3, 0, sequence=5)
        found = d.fault_episodes([a, b], {})
        self.assertEqual(found[0]['state'], 'observation_gap')
        self.assertIsNone(found[0]['end_time'])
        self.log.write_text('I1028 21:10:58.002000 42] SAFETY_DIAGNOSTIC sequence=1 raw_runstate=2\n'
                            'I1028 21:10:59.002000 42] SAFETY_DIAGNOSTIC sequence=5 raw_runstate=0\n')
        warnings = []
        d.log_rows(self.log, d.MIB, warnings)
        self.assertTrue(any('不连续' in w for w in warnings))

    def test_malformed_code_is_not_usable_safety_evidence(self):
        self.log.write_text('I1028 21:10:58.002000 42] SAFETY_DIAGNOSTIC sequence=1 raw_runstate=nan\n')
        with patch.object(d, 'bag_rows', return_value=[]):
            report = d.analyze(self.run, None, '指定路径', [])
        self.assertFalse(report['has_safety_evidence'])
        self.assertTrue(any('格式变化' in w for w in report['warnings']))

    def test_contract_detects_unreviewed_safety_change(self):
        manifest = self.root / 'contract.json'
        original = contract.snapshot()
        manifest.write_text(json.dumps(original))
        changed = json.loads(json.dumps(original))
        changed['reviewed_sources'][contract.SOURCES[0]] = 'changed'
        with patch.object(contract, 'snapshot', return_value=changed):
            with self.assertRaisesRegex(ValueError, '未同步'):
                contract.check(manifest=manifest)

    @unittest.skipUnless(os.name == 'posix', 'Bash dispatch fixture')
    def test_launcher_diagnosis_dispatch_precedes_all_launch_side_effects(self):
        folder = self.root / 'scripts with spaces'
        folder.mkdir()
        launcher = folder / 'run_sany_lidar_loc.sh'
        shutil.copyfile(Path(d.__file__).with_name(launcher.name), launcher)
        helper = folder / 'diagnose_localization.sh'
        helper.write_text('#!/bin/bash\nprintf "argument=<%s>\\n" "$@"\nexit 0\n')
        result = subprocess.run(['bash', str(launcher), '--diagnose', '/path with spaces/event', '--color', 'never'],
                                text=True, capture_output=True, timeout=5,
                                env=dict(os.environ, LIGHTNING_LM_RUN_MODE='invalid',
                                         LIGHTNING_LM_OUT_ROOT=str(self.root/'must-not-exist')))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('argument=</path with spaces/event>', result.stdout)
        self.assertFalse((self.root/'must-not-exist').exists())


if __name__ == '__main__':
    unittest.main()
