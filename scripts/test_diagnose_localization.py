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

    def test_unknown_without_any_output_is_explained_as_startup(self):
        fields = dict(covariance_valid=0, localization_uncertainty=0,
                      localization_consistency=0, yaw_condition=0, translation_condition=0,
                      unknown_tolerated=0, optional_can_unavailable=0)
        startup = self.causes(**fields, output_age_sec='inf')
        self.assertEqual([f['key'] for f in startup], ['uncertainty_unknown', 'consistency_unknown'])
        for finding in startup:
            self.assertEqual(finding['phase'], 'startup')
            self.assertIn('启动后尚无有效定位输出', finding['title'])
            self.assertEqual(finding['evidence']['output_age_sec'], 'inf')
            self.assertIn('初始位姿、全局重定位候选拒绝原因及接受时间', finding['advice'])
            self.assertNotIn('IMU 时间', finding['advice'])
        runtime = self.causes(**fields, output_age_sec=.1)
        self.assertTrue(all('phase' not in f for f in runtime))

    def test_unknown_without_yaw_pass_explains_window_without_assigning_phase(self):
        row = dict(raw_runstate=4, covariance_valid=1, localization_uncertainty=1,
                   localization_consistency=0, yaw_condition=0, translation_condition=0,
                   unknown_tolerated=0, optional_can_unavailable=0, yaw_pass_age_sec='inf',
                   imu_arrival_age_sec=.003847, imu_stamp_lag_sec=-.069712)
        for output_age in (.002572, .25, 2.0):
            with self.subTest(output_age=output_age):
                out = d.classify(dict(row, output_age_sec=output_age), {'motion_window_sec': .5})
                self.assertEqual([f['key'] for f in out], ['consistency_unknown'])
                self.assertNotIn('phase', out[0])
                self.assertEqual(out[0]['title'], '运动一致性证据不足，未满足 UNKNOWN 宽容条件')
                self.assertIn('尚无可用航向 PASS', out[0]['advice'])
                self.assertIn('核查是否处于运动评估窗口建立/重建阶段', out[0]['advice'])
                self.assertIn('也可能是 IMU 历史区间覆盖不足', out[0]['advice'])
                self.assertIn('不能仅凭当前新鲜判定根因', out[0]['advice'])
                self.assertNotIn('正常预热', out[0]['advice'])
                self.assertEqual(out[0]['evidence']['motion_window_sec'], .5)

    def test_unknown_with_previous_yaw_pass_keeps_input_advice(self):
        out = self.causes(localization_consistency=0, yaw_condition=0, translation_condition=0,
                          unknown_tolerated=0, optional_can_unavailable=0, yaw_pass_age_sec=.156812,
                          output_age_sec=.002572, imu_arrival_age_sec=.003847)
        self.assertEqual([f['key'] for f in out], ['consistency_unknown'])
        self.assertNotIn('phase', out[0])
        self.assertIn('0.156812s > 0.15s', out[0]['title'])
        self.assertEqual(out[0]['advice'], '检查 IMU 时间区间覆盖、时间戳跳变与输出间隔；当前 IMU 到达新鲜也不证明历史积分区间完整。')

    def test_startup_without_yaw_pass_keeps_initialization_advice(self):
        out = self.causes(covariance_valid=0, localization_uncertainty=0,
                          localization_consistency=0, yaw_condition=0, translation_condition=0,
                          unknown_tolerated=0, optional_can_unavailable=0,
                          yaw_pass_age_sec='inf', output_age_sec='inf')
        self.assertEqual([f['key'] for f in out], ['uncertainty_unknown', 'consistency_unknown'])
        for finding in out:
            self.assertEqual(finding['phase'], 'startup')
            self.assertIn('启动后尚无有效定位输出', finding['title'])
            self.assertEqual(finding['advice'], '先核对初始位姿、全局重定位候选拒绝原因及接受时间，再核对首个有效定位输出。')

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

    def test_severe_radar_conditions_are_side_evidence_with_latest_source(self):
        event = self.event(3 * 10**9)
        (event / 'timeline.csv').write_text(
            'wall_ns,raw,published,output_stale,output_age_sec,time_sync_fault_mask,bounded_sync_warning_mask\n'
            '3000000000,-9,-9,1,inf,4,0\n')
        rows = [self.log_row(1, -9, output_stale=1, output_age_sec='inf', time_sync_fault_mask=4),
                self.log_row(2, -9, output_stale=1, output_age_sec='inf',
                             time_sync_fault_mask=4, non_sync_fault=1)]
        with patch.object(d, 'log_rows', return_value=(rows, {})):
            report = d.analyze(event, self.p, '运行中的定位实例', [])
        self.assertEqual(len(report['findings']), 2)  # One output-stale interval per source.
        self.assertTrue(all(f['code'] == -9 for f in report['findings'] + report['current_findings']))
        clues = {c['key']: c for c in report['coexisting_radar_conditions']}
        self.assertEqual(clues['radar_time_sync']['device'], 'left')
        self.assertEqual(clues['radar_time_sync']['raw_runstate'], -9)
        self.assertEqual(clues['radar_time_sync']['last_ns'], 3 * 10**9)
        self.assertEqual(clues['radar_time_sync']['source'], f'{event / "timeline.csv"}:2')
        self.assertEqual(clues['radar_other']['source'], 'log:2')
        rendered = d.render(report)
        self.assertIn('并存的较低优先级雷达诊断条件', rendered)
        self.assertIn('对应记录原始码为 -9', rendered)
        self.assertIn('不代表当前仍异常', rendered)
        self.assertIn('不能据此认定启动失败原因', rendered)
        self.assertIn('时间同步诊断未通过: left', rendered)
        self.assertNotIn('历史记录 [2]', rendered)

    def test_severe_radar_side_evidence_excludes_enabled_bounded_warnings(self):
        for policy, mask, expected in [(1, 4, []), (1, 12, ['right']), (0, 4, ['left'])]:
            with self.subTest(policy=policy, mask=mask):
                row = self.log_row(1, -9, output_stale=1, output_age_sec='inf',
                                   time_sync_fault_mask=mask, bounded_policy_enabled=policy,
                                   bounded_sync_warning_mask=4, non_sync_fault=0)
                with patch.object(d, 'log_rows', return_value=([row], {})):
                    report = d.analyze(self.run, None, '指定路径', [])
                self.assertEqual([c['device'] for c in report['coexisting_radar_conditions']], expected)
                self.assertEqual([f['code'] for f in report['findings']], [-9])
                if not expected:
                    self.assertNotIn('并存的较低优先级雷达诊断条件', d.render(report))

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

    def test_render_distinguishes_predicate_change_from_recovery(self):
        rows = [self.log_row(1, 4, localization_uncertainty=0),
                self.log_row(2, -9, output_stale=1, output_age_sec='inf'),
                self.log_row(3, -9, output_stale=1, output_age_sec=.8),
                self.log_row(4, 0)]
        with patch.object(d, 'log_rows', return_value=(rows, {})):
            report = d.analyze(self.run, None, '指定路径', [])
        rendered = d.render(report)
        self.assertIn('时段状态: 原始码切换为 -9（未恢复为 0）', rendered)
        self.assertIn('时段状态: 判据已变化，原始码仍为 -9（未恢复为 0）', rendered)
        self.assertIn('时段状态: 已观测到原始码恢复为 0', rendered)

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
        self.assertTrue(all(f['capture_coverage_basis'] == 'episode_start' for f in report['findings']))
        self.assertTrue(any('不代表从未采集' in w for w in report['warnings']))
        self.assertIn('不保证整段覆盖或中间完整', d.render(report))

    def test_capture_coverage_text_only_describes_episode_start(self):
        event = self.event(100 * 10**9)
        s = d.summary_at(event)
        s['streams'] = [dict(type='geosun_msgs/msg/Safety', written=20,
                             first_relative_sec=0, last_relative_sec=10)]
        (event / 'summary.json').write_text(json.dumps(s))
        with patch.object(d, 'log_rows', return_value=([
                self.log_row(90, -9, output_stale=1, output_age_sec='inf'),
                self.log_row(120, -9, output_stale=1, output_age_sec='inf')], {})):
            report = d.analyze(event, self.p, '运行中的定位实例', [])
        finding = report['findings'][0]
        # The event overlaps the middle of this episode, but not its onset.
        self.assertEqual(finding['capture_coverage'], 'no_retained_overlap')
        self.assertEqual(finding['capture_coverage_basis'], 'episode_start')
        self.assertIn('故障起点数据覆盖: 起始时刻未找到保留事件范围覆盖（不代表整段未采集）',
                      d.render(report))

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

    def test_runtime_protections_are_visible_with_healthy_safety(self):
        self.log.write_text(
            'I1028 12:00:00.000000 42] SAFETY_DIAGNOSTIC sequence=1 raw_runstate=0 published_runstate=0\n'
            'W1028 12:00:01.000000 42] drop stale lidar correction before matching: age=0.309895 sec, hard_deadline=0.3\n'
            'W1028 12:00:02.000000 42] drop stale lidar correction before matching: age=0.310037 sec, hard_deadline=0.3\n'
            'W1028 12:00:03.000000 42] adaptive lidar load step 0 -> 1, processing_ms=70, latency_ms=20, tracking_healthy=1\n'
            'W1028 12:00:04.000000 42] BLOCK_SOLVER_INCREMENTAL_REBUILD solver_pose_count=5 retained_pose_count=0 indexed_pose_count=2\n'
            'W1028 12:00:05.000000 42] No point, skip this scan!\n'
            'I1028 12:00:06.000000 42] SAFETY_DIAGNOSTIC sequence=2 raw_runstate=0 published_runstate=0\n')
        report = d.analyze(self.run, None, '指定路径', [])
        self.assertEqual(report['findings'], [])
        self.assertEqual(report['current_findings'], [])
        self.assertEqual(report['log_samples'], 2)
        clues = {c['key']: c for c in report['runtime_protections']}
        self.assertEqual(set(clues), {'stale_lidar_correction', 'adaptive_lidar_load_step'})
        stale = clues['stale_lidar_correction']
        self.assertEqual(stale['count'], 2)
        self.assertEqual(stale['first_time'], '10-28 12:00:01.000000')
        self.assertEqual(stale['last_time'], '10-28 12:00:02.000000')
        self.assertEqual(stale['evidence'], {'age': .310037, 'hard_deadline': .3})
        self.assertIn('age=0.310037', stale['text'])
        self.assertIn(str(self.log) + ':byte=', stale['source'])
        self.assertNotEqual(stale['first_source'], stale['source'])
        self.assertEqual(clues['adaptive_lidar_load_step']['evidence'],
                         dict(from_step=0, to_step=1, processing_ms=70, latency_ms=20, tracking_healthy=1))
        self.assertTrue(all('code' not in clue for clue in clues.values()))
        rendered = d.render(report)
        self.assertIn('独立于 Safety 故障判据，不计入故障段，不等于定位失败', rendered)
        self.assertIn('故障判据段: 共 0 段', rendered)
        self.assertIn('age=0.310037', rendered)
        self.assertNotIn('BLOCK_SOLVER_INCREMENTAL_REBUILD', rendered)
        self.assertNotIn('No point', rendered)
        json.dumps(report, allow_nan=False)

    def test_protections_do_not_split_or_add_safety_fault_episodes(self):
        first = 'I1028 12:00:00.000000 42] SAFETY_DIAGNOSTIC sequence=1 raw_runstate=4 published_runstate=4 yaw_condition=2\n'
        last = 'I1028 12:00:02.000000 42] SAFETY_DIAGNOSTIC sequence=2 raw_runstate=0 published_runstate=0\n'
        self.log.write_text(first + last)
        baseline = d.analyze(self.run, None, '指定路径', [])['findings']
        self.log.write_text(first +
            'W1028 12:00:01.000000 42] adaptive lidar load step 1 -> 2, processing_ms=80, latency_ms=20, tracking_healthy=1\n' + last)
        report = d.analyze(self.run, None, '指定路径', [])
        self.assertEqual(report['findings'], baseline)
        self.assertEqual(len(report['runtime_protections']), 1)

    def test_explicit_event_protections_use_exact_timeline_window(self):
        ns = int(dt.datetime(2026, 10, 28, 12, 0, 1, tzinfo=dt.timezone.utc).timestamp() * 1e9)
        event = self.event(ns)
        (event / 'timeline.csv').write_text(
            f'wall_ns,raw,published\n{ns},0,0\n{ns + 2_000_000_000},0,0\n')
        self.log.write_text(
            'W1028 12:00:00.999999 42] drop stale lidar correction before matching: age=0.9 sec, hard_deadline=0.3\n'
            'W1028 12:00:00.999999 42] adaptive lidar load step 8 -> 9, processing_ms=900\n'
            'W1028 12:00:01.000000 42] drop stale lidar correction before matching: age=0.31 sec, hard_deadline=0.3\n'
            'W1028 12:00:02.000000 42] adaptive lidar load step 0 -> 1, processing_ms=70, latency_ms=20, tracking_healthy=1\n'
            'W1028 12:00:03.000000 42] drop stale lidar correction before matching: age=0.32 sec, hard_deadline=0.3\n'
            'W1028 12:00:03.000001 42] drop stale lidar correction before matching: age=0.8 sec, hard_deadline=0.3\n'
            'W1028 12:00:03.000001 42] adaptive lidar load step 9 -> 10, processing_ms=800\n'
            'drop stale lidar correction before matching: age=0.7 sec, hard_deadline=0.3\n')
        for with_summary in (True, False):
            with self.subTest(with_summary=with_summary):
                if not with_summary:
                    (event / 'summary.json').unlink()
                report = d.analyze(event, None, '指定路径', [])
                clues = {c['key']: c for c in report['runtime_protections']}
                self.assertEqual(clues['stale_lidar_correction']['count'], 2)
                self.assertEqual(clues['stale_lidar_correction']['first_ns'], ns)
                self.assertEqual(clues['stale_lidar_correction']['last_ns'], ns + 2_000_000_000)
                self.assertEqual(clues['stale_lidar_correction']['evidence']['age'], .32)
                self.assertEqual(clues['adaptive_lidar_load_step']['count'], 1)
                self.assertEqual(clues['adaptive_lidar_load_step']['evidence']['to_step'], 1)
                self.assertEqual(report['findings'], [])
                self.assertIn('指定事件时间线', report['runtime_protection_scope'])
                self.assertTrue(any('缺少可解析时间戳' in warning for warning in report['warnings']))

    def test_event_without_timeline_does_not_count_run_protections(self):
        event = self.event(100)
        self.log.write_text('W1028 12:00:01.000000 42] drop stale lidar correction before matching: age=0.31 sec, hard_deadline=0.3\n')
        with patch.object(d, 'bag_rows', return_value=[]):
            report = d.analyze(event, None, '指定路径', [])
        self.assertEqual(report['runtime_protections'], [])
        self.assertTrue(any('缺少指定事件时间线范围' in warning for warning in report['warnings']))

    def test_runtime_protection_counts_respect_scan_budget(self):
        clue = 'W1028 12:00:01.000000 42] drop stale lidar correction before matching: age=0.31 sec, hard_deadline=0.3\n'
        tail = clue * 2
        self.log.write_bytes((clue + 'x' * 200 + '\n' + tail).encode())
        protections, warnings = {}, []
        rows, _ = d.log_rows(self.log, len(tail.encode()) + 1, warnings,
                            protections=protections)
        self.assertEqual(rows, [])
        self.assertEqual(protections['stale_lidar_correction']['count'], 2)
        self.assertIn(f'byte={len(clue.encode()) + 201}', protections['stale_lidar_correction']['first_source'])
        self.assertTrue(any('日志仅检查末尾' in warning for warning in warnings))
        # Streaming summaries keep the full scanned count even beyond 1000 occurrences.
        self.log.write_text(clue * 1002)
        protections = {}
        d.log_rows(self.log, d.MIB, [], protections=protections)
        self.assertEqual(protections['stale_lidar_correction']['count'], 1002)

    def test_log_context_limit_preserves_latest_relocalization_evidence(self):
        start = dt.datetime(2026, 10, 28, 12, tzinfo=dt.timezone.utc)
        lines = ['I1028 12:00:00.000000 42] SAFETY_DIAGNOSTIC raw_runstate=-9 output_stale=1 output_age_sec=inf']
        for i in range(1001):
            stamp = (start + dt.timedelta(milliseconds=100 * (i + 1))).strftime('%m%d %H:%M:%S.%f')
            lines.append(f'I{stamp} 42] GLOBAL_RELOCALIZATION[solid_kiss] rejected: candidate={i}')
        lines += ['I1028 12:01:41.000000 42] GLOBAL_RELOCALIZATION[solid_kiss] accepted: candidate=221',
                  'I1028 12:01:42.000000 42] SAFETY_DIAGNOSTIC raw_runstate=-9 output_stale=1 output_age_sec=inf']
        self.log.write_text('\n'.join(lines) + '\n')
        context, warnings = [], []
        rows, _ = d.log_rows(self.log, d.MIB, warnings, start, context)
        self.assertEqual(len(rows), 2)
        self.assertEqual(len(context), 1000)
        self.assertTrue(context[0]['text'].endswith('candidate=2'))
        self.assertIn('accepted: candidate=221', context[-1]['text'])
        self.assertEqual(sum('最新 1000 条' in w for w in warnings), 1)
        report = d.analyze(self.run, None, '指定路径', [])
        self.assertIn('accepted: candidate=221', report['findings'][0]['context'][-1]['text'])
        self.assertIn('accepted: candidate=221', d.render(report))

    def test_tracking_rejection_onset_survives_relocalization_flood_and_recovery(self):
        start = dt.datetime(2026, 10, 29, 17, 29, 21, tzinfo=dt.timezone.utc)
        def line(seconds, text):
            stamp = (start + dt.timedelta(seconds=seconds)).strftime('%m%d %H:%M:%S.%f')
            return f'I{stamp} 580916 lidar_loc.cc:1820] {text}'
        lines = [
            line(.9, 'SAFETY_DIAGNOSTIC sequence=1 raw_runstate=0 published_runstate=0'),
            line(.993317, 'localization failed! score: 2.28838, iterations=1, raw_odom_valid=1, raw_odom_delta=0.00848839, map_checked=1, map_overlap=0.728835'),
            line(.998689, 'SAFETY_DIAGNOSTIC sequence=2 raw_runstate=-9 published_runstate=-9 localization_lost=1 output_stale=0 output_age_sec=0.00526'),
            line(1, 'localization pose and map outputs disabled after rejected localization: consecutive_lost_frames=1'),
        ]
        lines += [line(2 + i * .5, f'GLOBAL_RELOCALIZATION[solid_kiss] rejected: attempt={i}')
                  for i in range(1105)]
        lines += [
            line(880, 'SOLID_KISS_REJECT candidate=235, inliers=5, reason=too_few_robust_inliers'),
            line(882, 'GLOBAL_RELOCALIZATION[solid_kiss] accepted: candidate=236'),
            line(882.1, 'localization pose and map outputs recovered after a fresh valid lidar match'),
            line(882.2, 'SAFETY_DIAGNOSTIC sequence=3 raw_runstate=0 published_runstate=0'),
        ]
        self.log.write_text('\n'.join(lines) + '\n')
        report = d.analyze(self.run, None, '指定路径', [])
        self.assertEqual(report['latest']['raw_runstate'], 0)
        self.assertEqual([f['key'] for f in report['findings']], ['localization_lost'])
        finding = report['findings'][0]
        self.assertEqual(finding['state'], 'recovered')
        self.assertEqual(len(finding['onset_context']), 2)
        self.assertIn('map_overlap=0.728835', finding['onset_context'][0]['text'])
        self.assertIn('outputs disabled', finding['onset_context'][1]['text'])
        self.assertLessEqual(len(finding['context']), 6)
        rendered = d.render(report)
        for evidence in ('map_overlap=0.728835', 'outputs disabled',
                         'candidate=235, inliers=5', 'accepted: candidate=236'):
            self.assertIn(evidence, rendered)
        self.assertIn('起点附近日志（相关线索）', rendered)
        self.assertIn('不是已证明的物理根因', rendered)
        for clue in finding['onset_context']:
            self.assertTrue(self.log.read_bytes()[clue['byte_offset']:].startswith(clue['text'].encode()))

    def test_boundary_context_cache_is_bounded_and_warns_when_exhausted(self):
        start = dt.datetime(2026, 10, 28, 12, tzinfo=dt.timezone.utc)
        lines = []
        for i in range(1200):
            stamp = (start + dt.timedelta(milliseconds=100 * i)).strftime('%m%d %H:%M:%S.%f')
            lines.append(f'I{stamp} 42] localization failed! map_overlap=0.7 frame={i}')
        self.log.write_text('\n'.join(lines) + '\n')
        context, warnings = [], []
        d.log_rows(self.log, d.MIB, warnings, start, context)
        self.assertEqual(len(context), 1128)
        self.assertEqual(len({c['source'] for c in context}), len(context))
        self.assertTrue(context[0]['text'].endswith('frame=0'))
        self.assertTrue(context[-1]['text'].endswith('frame=1199'))
        self.assertEqual(sum('关联日志超过 1000 条' in w for w in warnings), 1)
        self.assertEqual(sum('关键边界线索超过 128 条' in w for w in warnings), 1)

    def test_boundary_context_does_not_expand_log_scan_budget(self):
        tail = ('I1028 12:01:00.000000 42] GLOBAL_RELOCALIZATION[solid_kiss] accepted: candidate=236\n'
                'I1028 12:01:00.100000 42] SAFETY_DIAGNOSTIC raw_runstate=0\n')
        self.log.write_bytes(('I1028 12:00:00.000000 42] localization failed! map_overlap=0.728835\n' +
                              'x' * 200 + '\n' + tail).encode())
        context, warnings = [], []
        rows, _ = d.log_rows(self.log, len(tail.encode()) + 1, warnings,
                             dt.datetime(2026, 10, 28, tzinfo=dt.timezone.utc), context)
        self.assertEqual(len(rows), 1)
        self.assertEqual(len(context), 1)
        self.assertIn('accepted: candidate=236', context[0]['text'])
        self.assertTrue(any('日志仅检查末尾' in w for w in warnings))

    def test_merged_context_sorts_by_timestamp_then_numeric_offset(self):
        lines = [
            'I1028 12:00:02.000000 42] localization pose and map outputs disabled',
            'I1028 12:00:01.000000 42] localization failed! map_overlap=0.728835',
            'I1028 12:00:02.000000 42] GLOBAL_RELOCALIZATION[solid_kiss] accepted: candidate=236',
        ]
        self.log.write_text('\n'.join(lines) + '\n')
        context = []
        d.log_rows(self.log, d.MIB, [], dt.datetime(2026, 10, 28, tzinfo=dt.timezone.utc), context)
        self.assertEqual([c['text'] for c in context], [lines[1], lines[0], lines[2]])
        self.assertEqual(len({c['source'] for c in context}), 3)

    def test_solid_kiss_failure_clues_keep_details_and_evidence_offsets(self):
        clues = [
            'failed to load SOLiD ICP cloud: /map/candidate_221.pcd',
            'failed to build SOLiD+KISS target: candidate=221, error=empty target',
            'failed to prepare SOLiD+KISS query: query_frames=3, error=too few points',
            'SOLID_KISS_REJECT candidate=221, query_frames=3, rough_matches=44, core_matches=3, inliers=0, reason=insufficient inliers',
            'GLOBAL_RELOCALIZATION[solid_kiss] rejected: candidate=221',
        ]
        lines = ['I1028 12:00:00.000000 42] SAFETY_DIAGNOSTIC raw_runstate=-9 output_stale=1 output_age_sec=inf']
        lines += [f'I1028 12:00:00.{i + 1}00000 42] {clue}' for i, clue in enumerate(clues)]
        lines.append('I1028 12:00:01.000000 42] SAFETY_DIAGNOSTIC raw_runstate=-9 output_stale=1 output_age_sec=inf')
        self.log.write_text('\n'.join(lines) + '\n')
        report = d.analyze(self.run, None, '指定路径', [])
        self.assertEqual(report['log_samples'], 2)
        self.assertEqual([f['key'] for f in report['findings']], ['output_stale'])
        context = report['findings'][0]['context']
        self.assertEqual(len(context), len(clues))
        log_bytes = self.log.read_bytes()
        for clue, entry in zip(clues, context):
            self.assertIn(clue, entry['text'])
            self.assertIsNotNone(entry['wall_ns'])
            offset = int(entry['source'].rsplit('=', 1)[1])
            self.assertTrue(log_bytes[offset:].startswith(entry['text'].encode()))
        rendered = d.render(report)
        self.assertIn('rough_matches=44, core_matches=3, inliers=0', rendered)
        self.assertIn('不是已证明的物理根因', rendered)

    def test_latest_solid_rejection_survives_later_global_summaries(self):
        lines = [
            'I1028 12:00:00.000000 42] SAFETY_DIAGNOSTIC raw_runstate=-9 output_stale=1 output_age_sec=inf',
            'I1028 12:00:00.100000 42] SOLID_KISS_REJECT candidate=221, inliers=3, reason=too_few_robust_inliers',
            'I1028 12:00:00.200000 42] SOLID_KISS_REJECT candidate=222, inliers=5, reason=too_few_robust_inliers',
        ]
        lines += [f'I1028 12:00:0{i}.000000 42] GLOBAL_RELOCALIZATION[solid_kiss] rejected: attempt={i}'
                  for i in range(1, 9)]
        lines.append('I1028 12:00:09.000000 42] SAFETY_DIAGNOSTIC raw_runstate=-9 output_stale=1 output_age_sec=inf')
        self.log.write_text('\n'.join(lines) + '\n')
        report = d.analyze(self.run, None, '指定路径', [])
        context = report['findings'][0]['context']
        self.assertEqual(len(context), 6)
        self.assertIn('SOLID_KISS_REJECT candidate=222', context[0]['text'])
        self.assertIn('attempt=8', context[-1]['text'])
        rendered = d.render(report)
        self.assertIn('candidate=222, inliers=5, reason=too_few_robust_inliers', rendered)
        self.assertIn('attempt=8', rendered)
        self.assertNotIn('candidate=221', rendered)

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
        with patch.object(d, 'log_rows', return_value=([a, b], {})):
            report = d.analyze(self.run, None, '指定路径', [])
        rendered = d.render(report)
        self.assertIn('时段状态: 后续存在观测缺口，结束或恢复时刻未知', rendered)
        self.assertNotIn('时段状态: 已观测到原始码恢复为 0', rendered)
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
