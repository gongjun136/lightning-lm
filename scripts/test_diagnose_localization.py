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

    def test_latest_partial_launch_failure_beats_old_event_and_names_loader_layer(self):
        old_ns = int(dt.datetime(2026, 10, 29, tzinfo=dt.timezone.utc).timestamp() * 1e9)
        self.event(old_ns)
        failed = self.root / 'sany_lidar_20261030_102537'
        (failed / 'logs').mkdir(parents=True)
        (failed / 'config.yaml').write_text('diagnostics: {}\n')
        message = '/install/export_map_manifest: error while loading shared libraries: liblightning_lm.libs.so: cannot open shared object file: No such file or directory'
        (failed / 'logs/map_manifest.log').write_text(message + '\n[ros2run]: Process exited with failure 127\n')
        os.utime(failed, (1000000000, 1000000000))  # copied/changed mtime must not replace launch time
        selected, process, selection, warnings = d.select(self.args(), [], [self.root])
        self.assertEqual(selected, failed)
        self.assertIsNone(process)
        with patch.object(d, 'bag_rows', side_effect=AssertionError('preflight failure is not a rosbag')):
            report = d.analyze(selected, process, selection, warnings)
        self.assertTrue(report['has_launch_evidence'])
        self.assertFalse(report['has_safety_evidence'])
        self.assertEqual(len(report['launch_findings']), 1)
        finding = report['launch_findings'][0]
        self.assertEqual(finding['key'], 'shared_library_missing')
        self.assertEqual(finding['stage'], 'runtime_loader')
        self.assertIn('liblightning_lm.libs.so', finding['title'])
        self.assertEqual(finding['evidence'], message)
        self.assertIn('map_manifest.log:byte=0', finding['source'])
        self.assertNotEqual(finding['key'], 'map_preflight_rejected')
        self.assertIn('工具执行前', finding['advice'])
        self.assertIn('入口阻断证据', d.render(report))
        with patch.object(d, 'low_impact', return_value={}), patch.object(d, 'driver_diagnosis') as driver, \
                patch('sys.stdout', new_callable=io.StringIO) as output:
            status = d.main([str(failed), '--json'])
        self.assertEqual(status, 0)
        driver.assert_not_called()
        self.assertEqual(json.loads(output.getvalue())['launch_findings'][0]['key'], 'shared_library_missing')

    def test_old_partial_run_does_not_displace_newer_healthy_history(self):
        old = self.root / 'sany_lidar_20261027_102537'
        (old / 'logs').mkdir(parents=True)
        (old / 'logs/map_manifest.log').write_text('[ros2run]: Process exited with failure 127\n')
        expected = self.event(int(dt.datetime(2026, 10, 29, tzinfo=dt.timezone.utc).timestamp() * 1e9))
        self.assertEqual(d.select(self.args(), [], [self.root])[0], expected)

    def test_launch_preflight_reasons_are_kept_distinct(self):
        cases = [('Map manifest rejected: required map_manifest.yaml is missing: /maps/a/map_manifest.yaml\n',
                  'map_manifest.log', 'map_preflight_rejected'),
                 ('[ros2run]: Process exited with failure 127\n', 'map_manifest.log', 'preflight_tool_failed'),
                 ('ERROR: sensor delivery check timed out; check publishers, ROS_DOMAIN_ID\n',
                  'input_check.log', 'input_check_timeout'),
                 ('SANY vehicle configuration rejected: unknown VIN\n', 'config_resolution.log', 'configuration_rejected')]
        for content, name, key in cases:
            with self.subTest(key=key):
                path = self.run / 'logs' / name
                path.write_text(content)
                result = d.launch_failures(self.run, [])
                self.assertEqual([f['key'] for f in result], [key])
                path.unlink()

    def test_successful_or_incomplete_preflight_is_not_a_failure(self):
        (self.run / 'logs/map_manifest.log').write_text('/maps/map_manifest.yaml\n')
        (self.run / 'logs/input_check.log').write_text('Sensor delivery check passed for 4 required topics\n')
        self.assertEqual(d.launch_failures(self.run, []), [])
        partial = self.root / 'sany_lidar_20261030_103000'
        (partial / 'logs').mkdir(parents=True)
        self.assertEqual(d.choose_path(partial, preserve_run=True), partial)
        with patch.object(d, 'bag_rows', return_value=[]):
            report = d.analyze(partial, None, '指定路径', [])
        self.assertFalse(report['has_launch_evidence'])
        self.assertIn('可能仍在启动或提前退出', '\n'.join(report['warnings']))

    def test_explicit_event_is_not_replaced_by_live(self):
        event = self.event(100)
        self.event(200)
        self.assertEqual(d.select(self.args(path=str(event / 'summary.json')), [self.p])[0], event)

    def test_automatic_run_path_only_selects_its_own_events(self):
        own = self.event(100)
        self.event(200, 99, self.root / 'other')
        self.assertEqual(d.choose_path(self.run), own)

    def test_explicit_run_cli_keeps_early_faults_when_newer_own_event_exists(self):
        event_ns = int(dt.datetime(2026, 10, 28, 12, 0, 10, tzinfo=dt.timezone.utc).timestamp() * 10**9)
        event = self.event(event_ns)
        (event / 'timeline.csv').write_text('wall_ns,raw,published\n' + f'{event_ns},0,0\n')
        self.log.write_text(
            'I1028 12:00:01.000000 42] SAFETY_DIAGNOSTIC sequence=1 raw_runstate=4 published_runstate=4 localization_uncertainty=0\n'
            'I1028 12:00:02.000000 42] SAFETY_DIAGNOSTIC sequence=2 raw_runstate=-9 published_runstate=-9 localization_lost=1 output_age_sec=0.01\n'
            'I1028 12:00:10.000000 42] SAFETY_DIAGNOSTIC sequence=3 raw_runstate=0 published_runstate=0\n')
        selected, process, _, _ = d.select(self.args(path=str(self.run)), [self.p])
        self.assertEqual(selected, self.run)
        self.assertIsNone(process)
        self.assertEqual(d.choose_path(self.run), event)
        with patch.object(d, 'low_impact', return_value={}), patch('sys.stdout', new_callable=io.StringIO) as output:
            status = d.main([str(self.run), '--log-mb', '64', '--json'])
            report = json.loads(output.getvalue())
        self.assertEqual(status, 0)
        self.assertEqual(Path(report['selected']), self.run)
        self.assertEqual(report['log_samples'], 3)
        self.assertEqual([f['code'] for f in report['findings']], [4, -9])
        self.assertTrue(any('首次已观测门控关闭' in item['text'] for item in report['priority_summary']))
        with patch.object(d, 'low_impact', return_value={}), patch('sys.stdout', new_callable=io.StringIO) as output:
            status = d.main([str(event), '--log-mb', '64', '--json'])
            event_report = json.loads(output.getvalue())
        self.assertEqual(status, 0)
        self.assertEqual(Path(event_report['selected']), event)
        self.assertEqual(event_report['log_samples'], 1)
        self.assertEqual(event_report['findings'], [])

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

    def driver_report(self, active=True, has_evidence=True):
        return dict(schema=1, kind='livox_driver_diagnosis', selected='/runs/lidar_test',
                    active=active, has_evidence=has_evidence,
                    roles=[dict(role='front', device=dict(evidence_scope='stored_snapshot',
                        collected_at='2026-10-30T10:00:00+08:00', age_sec=120))],
                    findings=[dict(code='device_not_normal', severity='error', role='front',
                        message='历史快照 work_state=4 (Error)',
                        evidence=[dict(path='/runs/lidar_test/device.json', field='parameters.work_state',
                                       time='2026-10-30T10:00:00+08:00')])],
                    warnings=['当前未刷新；启动流量证据不能证明当前持续出数。'])

    def test_driver_discovery_prefers_configured_install_over_sibling_source(self):
        script = self.root / 'source/lightning-lm/scripts/diagnose_localization.py'
        sibling = self.root / 'source/sdk/livox_sdk/src/livox_ros_driver2/scripts/diagnose_lidar.py'
        installed = self.root / 'vehicle/lightning_lm_ws/install/livox_ros_driver2_core/lib/livox_ros_driver2_core/diagnose_lidar.py'
        for path in (sibling, installed):
            path.parent.mkdir(parents=True)
            path.write_text('')
        config = self.root / '.config/sany'
        config.mkdir(parents=True)
        (config / 'workspace.bash').write_text('export SANY_WS="$HOME/vehicle"\n')
        self.assertEqual(d.find_driver_diagnostic({}, self.root, script), installed)
        self.assertEqual(d.find_driver_diagnostic({'SANY_WS': ''}, self.root, script), sibling)
        self.assertEqual(d.find_driver_diagnostic({'SANY_WS': '', 'LIGHTNING_LM_INSTALL_SETUP':
                         str(installed.parents[3] / 'setup.bash')}, self.root, script), installed)

    def test_driver_subprocess_uses_cache_only_and_keeps_report(self):
        data = self.driver_report()
        completed = subprocess.CompletedProcess([], 0, json.dumps(data), '')
        with patch.object(d, 'find_driver_diagnostic', return_value=self.root / 'diagnose_lidar.py'), \
                patch.object(d.subprocess, 'run', return_value=completed) as run:
            upstream = d.driver_diagnosis(str(self.root))
        self.assertEqual(upstream['status'], 'available')
        self.assertEqual(upstream['report'], data)
        self.assertEqual(run.call_args.args[0], [d.sys.executable, '-B', str(self.root / 'diagnose_lidar.py'),
                                               '--json', '--no-refresh', '--root', str(self.root)])
        self.assertEqual(run.call_args.kwargs['timeout'], 10)
        text = '\n'.join(d.render_driver_diagnosis(upstream))
        self.assertIn('当前驱动运行实例（设备证据未刷新）', text)
        self.assertIn('stored_snapshot', text)
        self.assertIn('距采集=120 秒', text)
        self.assertIn('当前未刷新', text)
        self.assertIn('parameters.work_state', text)

    def test_driver_history_survives_many_startup_findings(self):
        report = self.driver_report()
        report['findings'] = [dict(code='startup', severity='error', message=f'启动异常 {i}', evidence=[])
                              for i in range(12)]
        report['findings'] += [dict(code='device_history', stage='device_history', severity='info',
            role='right', message='右雷达历史 Normal → Error → Normal；不代表当前状态',
            evidence=[dict(path='/runs/lidar/diagnostics/right_state.jsonl', line=22,
                           time='2026-10-30T14:00:00+08:00', seq=10, steady_ns=100)])]
        text = '\n'.join(d.render_driver_diagnosis(dict(message='读取完成', report=report)))
        self.assertIn('启动异常 9', text)
        self.assertNotIn('启动异常 10', text)
        self.assertIn('右雷达历史 Normal → Error → Normal；不代表当前状态', text)
        self.assertIn('/runs/lidar/diagnostics/right_state.jsonl:line=22 @ 2026-10-30T14:00:00+08:00', text)
        self.assertNotIn('steady_ns', text)
        self.assertIn('设备历史摘要单独保留', text)

    def test_driver_missing_old_failed_or_malformed_does_not_break_localization(self):
        self.log.write_text('I1028 21:10:58.502960 42] SAFETY_DIAGNOSTIC raw_runstate=0 published_runstate=0\n')
        cases = [subprocess.CompletedProcess([], 2, '', 'unrecognized arguments: --no-refresh'),
                 subprocess.CompletedProcess([], 0, json.dumps(dict(schema=0, kind='livox_driver_diagnosis')), ''),
                 subprocess.CompletedProcess([], 2, json.dumps(dict(schema=1, kind='livox_driver_diagnosis',
                                                    has_evidence=False, error='没有雷达包')), ''),
                 subprocess.CompletedProcess([], 0, json.dumps(self.driver_report() | {'roles': [None]}), ''),
                 subprocess.CompletedProcess([], 0, json.dumps(self.driver_report() | {'value': float('nan')}), ''),
                 subprocess.TimeoutExpired('driver', 10), OSError('unavailable')]
        for response in cases:
            with self.subTest(response=response), \
                    patch.object(d, 'low_impact', return_value={}), \
                    patch.object(d, 'select', return_value=(self.run, None, '历史包', [])), \
                    patch.object(d, 'find_driver_diagnostic', return_value=self.root / 'driver.py'), \
                    patch.object(d.subprocess, 'run', **({'side_effect': response} if isinstance(response, Exception)
                                                       else {'return_value': response})), \
                    patch('sys.stdout', new_callable=io.StringIO) as output:
                self.assertEqual(d.main(['--json']), 0)
                report = json.loads(output.getvalue())
                self.assertTrue(report['has_safety_evidence'])
                self.assertIn(report['upstream_lidar']['status'], ('unsupported', 'unavailable', 'failed'))
        with patch.object(d, 'find_driver_diagnostic', return_value=None), patch.object(d.subprocess, 'run') as run:
            self.assertEqual(d.driver_diagnosis()['status'], 'unavailable')
            run.assert_not_called()

    def test_explicit_localization_path_never_queries_current_driver(self):
        self.log.write_text('I1028 21:10:58.502960 42] SAFETY_DIAGNOSTIC raw_runstate=0 published_runstate=0\n')
        with patch.object(d, 'low_impact', return_value={}), patch.object(d, 'driver_diagnosis') as driver, \
                patch('sys.stdout', new_callable=io.StringIO) as output:
            self.assertEqual(d.main([str(self.run), '--json']), 0)
            report = json.loads(output.getvalue())
        driver.assert_not_called()
        self.assertEqual(report['evidence_scope'], 'explicit')
        self.assertEqual(report['upstream_lidar']['status'], 'not_requested')

    def test_auto_history_zero_is_not_current_health_with_active_driver_error(self):
        self.log.write_text('I1028 21:10:58.502960 42] SAFETY_DIAGNOSTIC raw_runstate=0 published_runstate=0\n')
        upstream = dict(status='available', message='缓存证据', report=self.driver_report())
        with patch.object(d, 'low_impact', return_value={}), \
                patch.object(d, 'select', return_value=(self.run, None, '历史包', [])), \
                patch.object(d, 'driver_diagnosis', return_value=upstream), \
                patch('sys.stdout', new_callable=io.StringIO) as output:
            self.assertEqual(d.main(['--color', 'always']), 0)
            text = output.getvalue()
        self.assertNotIn('\x1b[1;32m', text)
        self.assertIn('当前未发现定位进程，包内状态码不代表当前状态', text)
        self.assertIn('当前驱动运行实例', text)
        self.assertIn('历史快照 work_state=4 (Error)', text)

    def test_auto_no_localization_package_can_report_driver_only(self):
        upstream = dict(status='available', message='缓存证据', report=self.driver_report())
        with patch.object(d, 'low_impact', return_value={}), \
                patch.object(d, 'select', side_effect=ValueError('无定位包')), \
                patch.object(d, 'driver_diagnosis', return_value=upstream), \
                patch('sys.stdout', new_callable=io.StringIO) as output:
            self.assertEqual(d.main(['--json']), 0)
            report = json.loads(output.getvalue())
        self.assertFalse(report['has_safety_evidence'])
        self.assertEqual(report['evidence_scope'], 'unavailable')
        self.assertIn('定位证据不可用: 无定位包', d.render(report))
        self.assertIn('历史快照 work_state=4', d.render(report))

    def test_active_localization_and_historical_driver_keep_separate_scopes(self):
        self.log.write_text('I1028 21:10:58.502960 42] SAFETY_DIAGNOSTIC raw_runstate=4 published_runstate=4 yaw_condition=2\n')
        upstream = dict(status='available', message='缓存证据', report=self.driver_report(active=False))
        with patch.object(d, 'low_impact', return_value={}), \
                patch.object(d, 'select', return_value=(self.run, self.p, '运行中的定位实例', [])), \
                patch.object(d, 'driver_diagnosis', return_value=upstream), \
                patch('sys.stdout', new_callable=io.StringIO) as output:
            self.assertEqual(d.main(['--json']), 0)
            report = json.loads(output.getvalue())
        self.assertEqual(report['evidence_scope'], 'live_run')
        self.assertEqual(report['process']['pid'], 42)
        self.assertEqual(report['latest']['raw_runstate'], 4)
        self.assertIn('历史驱动包（当前未发现雷达 worker）', d.render(report))

    def test_auto_without_evidence_remains_exit_two(self):
        upstream = dict(status='no_evidence', message='无证据',
                        report=self.driver_report(active=False, has_evidence=False) | {'findings': []})
        with patch.object(d, 'low_impact', return_value={}), \
                patch.object(d, 'select', side_effect=ValueError('无定位包')), \
                patch.object(d, 'driver_diagnosis', return_value=upstream), \
                patch('sys.stdout', new_callable=io.StringIO) as output:
            self.assertEqual(d.main([]), 2)
        self.assertIn('历史驱动包（当前未发现雷达 worker）', output.getvalue())
        self.assertIn('不能据此确认当前设备或数据流正常', output.getvalue())

    def test_explicit_selection_or_analysis_errors_are_not_hidden_by_driver(self):
        for argv in ([str(self.root / 'missing'), '--json'], ['--pid', '99', '--json']):
            with self.subTest(argv=argv), patch.object(d, 'low_impact', return_value={}), \
                    patch.object(d, 'select', side_effect=ValueError('指定路径或 PID 无效')), \
                    patch.object(d, 'driver_diagnosis') as driver, \
                    patch('sys.stdout', new_callable=io.StringIO) as output:
                self.assertEqual(d.main(argv), 2)
                self.assertIn('error', json.loads(output.getvalue()))
                driver.assert_not_called()
        with patch.object(d, 'low_impact', return_value={}), \
                patch.object(d, 'select', return_value=(self.run, self.p, '当前定位', [])), \
                patch.object(d, 'analyze', side_effect=ValueError('损坏的定位包')), \
                patch.object(d, 'driver_diagnosis') as driver, \
                patch('sys.stdout', new_callable=io.StringIO) as output:
            self.assertEqual(d.main(['--json']), 2)
            self.assertEqual(json.loads(output.getvalue())['error'], '损坏的定位包')
            driver.assert_not_called()

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
        self.assertIn('位姿航向增量与 IMU 机体系 z 角速度积分', out[0]['advice'])
        self.assertIn('CAN 用于平移检查', out[0]['advice'])
        self.assertIn('缺少累计转角 turn 及窗口/位姿分支', out[0]['advice'])
        self.assertIn('无法重建确切动态阈值', out[0]['advice'])
        self.assertNotIn('IMU/CAN 积分', out[0]['advice'])

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

    def test_context_does_not_cross_observed_normal_between_short_losses(self):
        lines = [
            'I1028 12:00:00.100000 42] localization failed! map_checked=1 map_overlap=0.71',
            'I1028 12:00:00.200000 42] SAFETY_DIAGNOSTIC sequence=1 raw_runstate=-9 localization_lost=1',
            'I1028 12:00:00.250000 42] localization pose and map outputs disabled first',
            'I1028 12:00:00.300000 42] SAFETY_DIAGNOSTIC sequence=2 raw_runstate=0',
            'I1028 12:00:00.400000 42] localization failed! map_checked=1 map_overlap=0.72',
            'I1028 12:00:00.500000 42] SAFETY_DIAGNOSTIC sequence=3 raw_runstate=-9 localization_lost=1',
            'I1028 12:00:00.550000 42] localization pose and map outputs disabled second',
            'I1028 12:00:00.600000 42] SAFETY_DIAGNOSTIC sequence=4 raw_runstate=0',
            'I1028 12:00:00.700000 42] localization failed! map_checked=1 map_overlap=0.73',
        ]
        self.log.write_text('\n'.join(lines) + '\n')
        report = d.analyze(self.run, None, '指定路径', [])
        first, second = report['findings']
        self.assertIn('map_overlap=0.71', first['onset_context'][0]['text'])
        self.assertIn('map_overlap=0.72', second['onset_context'][0]['text'])
        self.assertTrue(first['onset_context'][1]['text'].endswith('first'))
        self.assertTrue(second['onset_context'][1]['text'].endswith('second'))
        self.assertEqual(first['context'], [])
        self.assertEqual(second['context'], [])

    def test_context_does_not_bridge_a_safety_sequence_gap(self):
        lines = [
            'I1028 12:00:00.100000 42] localization failed! map_checked=1 map_overlap=0.71',
            'I1028 12:00:00.200000 42] SAFETY_DIAGNOSTIC sequence=1 raw_runstate=-9 localization_lost=1',
            'I1028 12:00:00.250000 42] localization failed! map_checked=1 map_overlap=0.72',
            'I1028 12:00:00.300000 42] SAFETY_DIAGNOSTIC sequence=4 raw_runstate=-9 localization_lost=1',
            'I1028 12:00:00.400000 42] SAFETY_DIAGNOSTIC sequence=5 raw_runstate=0',
        ]
        self.log.write_text('\n'.join(lines) + '\n')
        first, second = d.analyze(self.run, None, '指定路径', [])['findings']
        self.assertEqual(first['state'], 'observation_gap')
        self.assertIn('map_overlap=0.71', first['onset_context'][0]['text'])
        self.assertEqual(first['context'], [])
        self.assertEqual(second['onset_context'], [])
        self.assertEqual(second['context'], [])

    def test_priority_summary_precedes_protections_and_distinguishes_output_phases(self):
        with (self.run / 'config.yaml').open('a') as config:
            config.write('relocalization:\n  backend: solid_kiss\n  map_min_overlap_ratio: 0.92\n'
                         '  kiss_matcher:\n    map_min_overlap_ratio: 0.99\n'
                         '  solid_kiss:\n    map_min_overlap_ratio: 0.75\n')
        lines = [
            'I1028 12:00:00.000000 42] SAFETY_DIAGNOSTIC sequence=1 raw_runstate=4 localization_uncertainty=0 output_age_sec=inf',
            'I1028 12:00:01.000000 42] SAFETY_DIAGNOSTIC sequence=2 raw_runstate=0',
            'I1028 12:00:02.000000 42] localization failed! score: 2.28838, map_checked=1 map_overlap=0.728835 raw_odom_valid=1',
            'I1028 12:00:02.010000 42] SAFETY_DIAGNOSTIC sequence=3 raw_runstate=-9 localization_lost=1 output_age_sec=0.005',
            'I1028 12:00:02.020000 42] localization pose and map outputs disabled',
            'I1028 12:00:02.300000 42] SOLID_KISS_REJECT candidate=91 inliers=9 reason=too_few_robust_inliers',
            'I1028 12:00:02.400000 42] GLOBAL_RELOCALIZATION[solid_kiss] accepted: candidate=14 overlap=0.98',
            'I1028 12:00:02.500000 42] SAFETY_DIAGNOSTIC sequence=4 raw_runstate=4 localization_consistency=0 localization_uncertainty=1 covariance_valid=1 unknown_tolerated=0 optional_can_unavailable=0 yaw_pass_age_sec=inf output_age_sec=0.008',
            'I1028 12:00:03.000000 42] SAFETY_DIAGNOSTIC sequence=5 raw_runstate=0',
            'I1028 12:00:03.100000 42] drop stale lidar correction before matching: age=0.31 hard_deadline=0.3',
        ]
        self.log.write_text('\n'.join(lines) + '\n')
        report = d.analyze(self.run, None, '指定路径', [])
        summary = '\n'.join(item['text'] for item in report['priority_summary'])
        self.assertIn('启动无有效定位输出', summary)
        self.assertIn('首次已观测门控关闭: 10-28 12:00:02.010000', summary)
        self.assertIn('map_overlap=0.728835 < relocalization.solid_kiss.map_min_overlap_ratio=0.75', summary)
        self.assertIn('运行配置判据', summary)
        self.assertIn('已证实未满足项: map_overlap', summary)
        self.assertNotIn('0.92', summary)
        self.assertNotIn('0.99', summary)
        self.assertIn('reason=too_few_robust_inliers', summary)
        self.assertIn('重定位接受', summary)
        self.assertIn('有定位输出后的运动证据 UNKNOWN', summary)
        self.assertIn('12:00:03.000000 已观测到原始码 0', summary)
        rendered = d.render(report)
        self.assertLess(rendered.index('优先结论 / 问题链'), rendered.index('运行保护线索'))

    def test_matching_summary_prefers_logged_thresholds_and_lists_all_failures(self):
        clue = dict(time='t=1', text='localization failed! map_checked=1 raw_odom_valid=1 map_points=100 '
                    'map_inside_xy=0.8 map_overlap=0.7 map_gravity_cos=0.9 '
                    'map_min_inside_xy=0.9 map_min_overlap=0.8 map_min_gravity_cos=0.95')
        text, config_used = d.matching_summary(clue, {'map_min_overlap_ratio': .65, 'overlap_key': 'config.overlap'})
        self.assertFalse(config_used)
        self.assertIn('运行日志判据: map_overlap=0.7 < map_min_overlap=0.8', text)
        self.assertIn('已证实未满足项: map_inside_xy、map_overlap、map_gravity_cos', text)
        self.assertNotIn('运行配置判据', text)

    def test_recent_yaw_fail_summary_keeps_four_real_shaped_segments_ahead_of_old_loss(self):
        def row(clock, code, **extra):
            stamp = dt.datetime.fromisoformat('2026-10-29T' + clock).replace(tzinfo=dt.timezone.utc)
            result = self.log_row(stamp.timestamp(), code, **extra)
            result['_time'] = '10-29 ' + clock
            return result
        rows = [row('17:29:21.998689', -9, localization_lost=1, output_age_sec=.0052619),
                row('17:44:03.471648', 0)]
        for start, end, error in (
                ('18:45:42.800266', '18:45:42.820221', .0489574),
                ('18:45:42.900221', '18:45:42.920220', .0455395),
                ('19:24:28.811980', '19:24:28.911988', .0583416),
                ('19:52:05.714780', '19:52:05.814821', .0423858)):
            rows += [row(start, 4, localization_consistency=2, localization_uncertainty=1,
                         yaw_condition=2, error_yaw_rad=error), row(end, 0)]
        findings = d.fault_episodes(rows, {})
        summary = d.priority_summary(findings, {})
        self.assertIn('航向运动一致性 FAIL: 已读取 4 个日志判据段', summary[0]['text'])
        self.assertIn('最近段 10-29 19:52:05.714780', summary[0]['text'])
        self.assertIn('残差 2.42853°', summary[0]['text'])
        self.assertIn('已观测到原始码 0: 10-29 19:52:05.814821', summary[0]['text'])
        self.assertEqual(summary[0]['source'], findings[-1]['source'])
        self.assertTrue(any('首次已观测门控关闭' in item['text'] for item in summary[1:]))

    def test_priority_separates_first_severe_gate_clear_and_later_raw_recovery(self):
        rows = [self.log_row(1, -9, output_stale=1, radar_all_failed=1, localization_lost=0),
                self.log_row(1.02, -9, output_stale=1, localization_lost=1),
                self.log_row(1.2, -9, output_stale=1, localization_lost=0), self.log_row(1.5, 0),
                self.log_row(2, -9, radar_all_failed=1, localization_lost=0),
                self.log_row(2.02, -9, localization_lost=1),
                self.log_row(2.4, -9, output_stale=1, localization_lost=0),
                self.log_row(2.5, 4, localization_consistency=0, localization_uncertainty=1,
                             unknown_tolerated=0, optional_can_unavailable=0), self.log_row(3, 0)]
        summary = '\n'.join(f['text'] for f in d.priority_summary(d.fault_episodes(rows, {}), {}))
        self.assertIn('首次已观测 -9: t=1；当时判据: 定位输出超时、所有必需雷达均不可用', summary)
        self.assertIn('首次已观测门控关闭: t=1.02', summary)
        self.assertIn('最长已观测门控关闭段: t=2.02 → t=2.4；门控关闭判据已解除', summary)
        self.assertIn('该时刻原始码仍为 -9；后续已观测到原始码 0: t=3', summary)
        self.assertNotIn('未在本段结束时确认恢复为 0', summary)

    def test_priority_includes_severe_without_gate_and_sync_role_recovery(self):
        rows = [self.log_row(1, -9, radar_all_failed=1, localization_lost=0), self.log_row(2, 0),
                self.log_row(3, 2, time_sync_fault_mask=1, bounded_policy_enabled=1,
                             bounded_sync_warning_mask=0, non_sync_fault=0, non_sync_all_pass=1),
                self.log_row(27.1, 0)]
        summary = '\n'.join(f['text'] for f in d.priority_summary(d.fault_episodes(rows, {}), {}))
        self.assertIn('首次已观测 -9: t=1', summary)
        self.assertNotIn('门控关闭', summary)
        self.assertIn('雷达同步故障 2: front，已读取 1 个判据段；最近段 t=3', summary)
        self.assertIn('已观测到原始码 0: t=27.1', summary)
        for ending in (self.log_row(27.1, 4, localization_uncertainty=0), None):
            changed = rows[:-1] + ([ending] if ending else [])
            text = '\n'.join(f['text'] for f in d.priority_summary(d.fault_episodes(changed, {}), {}))
            sync = next(line for line in text.splitlines() if '雷达同步故障 2:' in line)
            self.assertNotIn('已观测到原始码 0:', sync)
            self.assertIn('当时未恢复为原始码 0' if ending else '尚未观测到该段恢复', sync)

    def test_queue_lag_near_severe_is_timed_corroboration_not_root_cause(self):
        lines = [
            'I1028 12:00:00.000000 42] SAFETY_DIAGNOSTIC sequence=1 raw_runstate=0',
            'I1028 12:00:01.000000 42] SAFETY_DIAGNOSTIC sequence=2 raw_runstate=-9 output_stale=1 radar_all_failed=1 localization_lost=0',
            'W1028 12:00:01.001000 42] sensor queue lag reached 0.709954 sec; keeping fresh lidar admissible while the bounded queue evicts stale input',
            'I1028 12:00:01.002000 42] sensor queue lag recovered to 0 sec; resuming lidar input',
            'W1028 12:00:01.004000 42] drop stale lidar correction before matching: age=0.682182 sec, hard_deadline=0.3',
            'I1028 12:00:01.020000 42] SAFETY_DIAGNOSTIC sequence=3 raw_runstate=-9 output_stale=1 localization_lost=1',
            'I1028 12:00:01.800000 42] SAFETY_DIAGNOSTIC sequence=4 raw_runstate=0 localization_lost=0',
            'W1028 12:00:01.900000 42] sensor queue lag reached 9 sec; after healthy recovery',
            'W1028 12:00:01.901000 42] drop stale lidar correction before matching: age=9 sec, hard_deadline=0.3',
        ]
        self.log.write_text('\n'.join(lines) + '\n')
        report = d.analyze(self.run, None, '指定路径', [])
        summaries = [s for s in report['priority_summary'] if '处理链线索' in s['text']]
        self.assertEqual(len(summaries), 1)  # one queue record, despite several -9 predicates
        text = summaries[0]['text']
        self.assertIn('与 -9 起点 10-28 12:00:01.000000 相近', text)
        self.assertIn('10-28 12:00:01.001000 处理进度落后 0.709954s', text)
        self.assertIn('最新入队与最新处理的传感器时间差', text)
        self.assertIn('10-28 12:00:01.004000 过期 LIO 校正被丢弃，age=0.682182s', text)
        self.assertIn('不证明最早物理根因', text)
        self.assertNotIn('处理进度落后 9s', text)
        self.assertEqual(summaries[0]['source'].count('run_loc_online.stderr.log:byte='), 2)
        for finding in report['findings']:
            self.assertFalse(any('after healthy recovery' in c['text'] for c in finding['onset_context']))
        self.assertEqual(report['runtime_protections'][0]['count'], 2)  # full-run protection statistics remain separate

    def test_queue_recovery_alone_does_not_invent_processing_lag(self):
        self.log.write_text(
            'I1028 12:00:01.000000 42] SAFETY_DIAGNOSTIC sequence=1 raw_runstate=-9 output_stale=1\n'
            'I1028 12:00:01.001000 42] sensor queue lag recovered to 0 sec; resuming lidar input\n'
            'W1028 12:00:01.002000 42] drop stale lidar correction before matching: age=0.4 sec, hard_deadline=0.3\n'
            'I1028 12:00:02.000000 42] SAFETY_DIAGNOSTIC sequence=2 raw_runstate=0\n')
        report = d.analyze(self.run, None, '指定路径', [])
        self.assertFalse(any('处理进度落后' in s['text'] for s in report['priority_summary']))
        self.assertTrue(any('sensor queue lag recovered' in c['text'] for c in report['findings'][0]['context']))

    def test_priority_summary_does_not_invent_yaw_fail_for_unknown_or_uncertainty(self):
        for extra in (
                dict(localization_consistency=0, yaw_condition=0, unknown_tolerated=0,
                     optional_can_unavailable=0, yaw_pass_age_sec='inf', output_age_sec=.01),
                dict(localization_uncertainty=2, covariance_valid=1, major_sigma=.7)):
            rows = [self.log_row(1, 4, **extra), self.log_row(2, 0)]
            summary = d.priority_summary(d.fault_episodes(rows, {'major_sigma_limit_m': .5}), {})
            self.assertFalse(any('航向运动一致性 FAIL' in item['text'] for item in summary))

    def test_matching_summary_does_not_attribute_zero_points_or_unchecked_defaults(self):
        for extra, expected in [('map_checked=0 map_points=0', '未确认执行'),
                                ('map_checked=1 map_points=0', '有效评估点数为 0'),
                                ('map_checked=1', '缺少有效评估点证据')]:
            with self.subTest(extra=extra):
                clue = dict(time='t=1', text='localization failed! map_overlap=0 map_min_overlap=0.75 ' + extra)
                text, config_used = d.matching_summary(clue, {'map_min_overlap_ratio': .9})
                self.assertIn(expected, text)
                self.assertNotIn('已证实未满足项', text)
                self.assertFalse(config_used)

    def test_missing_matching_configuration_does_not_invent_overlap_threshold(self):
        config = d.matching_thresholds(self.run / 'config.yaml')
        clue = dict(time='t=1', text='localization failed! map_checked=1 map_overlap=0.728835')
        text, used = d.matching_summary(clue, config)
        self.assertFalse(used)
        self.assertNotIn('运行配置判据', text)
        self.assertNotIn('已证实未满足项', text)

    def test_explicit_nonfinite_map_points_cannot_use_legacy_positive_overlap(self):
        for points in ('nan', 'inf', '-inf'):
            with self.subTest(points=points):
                clue = dict(time='t=1', text=f'localization failed! map_checked=1 map_points={points} '
                            'map_overlap=0.7 map_min_overlap=0.75')
                text, used = d.matching_summary(clue, {'map_min_overlap_ratio': .9})
                self.assertIn('缺少有效评估点证据', text)
                self.assertNotIn('已证实未满足项', text)
                self.assertFalse(used)

    def test_matching_threshold_domains_distinguish_ratios_from_cosine(self):
        for metric, limit, valid, invalid in (
                ('map_inside_xy', 'map_min_inside_xy', (0, 1), (-.1, 1.1)),
                ('map_overlap', 'map_min_overlap', (0, 1), (-.1, 1.1)),
                ('map_gravity_cos', 'map_min_gravity_cos', (-1, 1), (-1.1, 1.1))):
            for threshold in valid + invalid:
                with self.subTest(metric=metric, threshold=threshold):
                    clue = dict(time='t=1', text='localization failed! map_checked=1 map_points=100 '
                                f'{metric}=0 {limit}={threshold}')
                    text, used = d.matching_summary(clue, {})
                    self.assertEqual('运行日志判据: ' + metric in text, threshold in valid)
                    self.assertFalse(used)

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


class ConfirmedLossPolicyTest(unittest.TestCase):
    def row(self, code, **values):
        return dict(raw_runstate=code, published_runstate=code,
                    policy_version=d.CONFIRMED_LOSS_POLICY, _time='now', _source='test', _scope='test', **values)

    def test_startup_grace_and_expiration_are_distinct_from_runtime_stale(self):
        pending = d.classify(self.row(4, startup_pending=1, startup_expired=0, startup_age_sec=7.9, output_age_sec='inf'), {})
        expired = d.classify(self.row(-9, startup_pending=1, startup_expired=1, startup_age_sec=8.1, output_stale=1), {})
        self.assertEqual([x['key'] for x in pending], ['startup_pending'])
        self.assertEqual([x['key'] for x in expired], ['startup_expired'])
        self.assertTrue(all(x['phase'] == 'startup' for x in pending + expired))

    def test_runtime_stale_is_degraded_and_does_not_claim_confirmed_loss(self):
        values = d.classify(self.row(4, output_stale=1, output_age_sec=.6, match_stale=1,
                           valid_match_age_sec=.7, confirmed_tracking_lost=0, localization_lost=0), {})
        self.assertEqual([x['key'] for x in values], ['output_stale_degraded', 'match_stale_degraded'])
        self.assertTrue(all(x['code'] == 4 for x in values))
        confirmed = d.classify(self.row(-9, output_stale=1, confirmed_tracking_lost=1, localization_lost=1), {})
        self.assertEqual([x['key'] for x in confirmed], ['confirmed_tracking_lost'])

    def test_one_rejection_closes_gate_without_establishing_loss(self):
        values = d.classify(self.row(4, startup_pending=0, output_gate_open=0,
                                    latest_match_valid=0, confirmed_tracking_lost=0, lost_frames=1), {})
        self.assertEqual([x['key'] for x in values], ['match_rejected_unconfirmed'])
        self.assertIn('尚未达到', values[0]['title'])

    def test_legacy_minus9_stale_keeps_legacy_meaning_but_new_version_does_not(self):
        legacy = d.classify(dict(raw_runstate=-9, output_stale=1, output_age_sec=.4), {})
        modern = d.classify(self.row(-9, output_stale=1, output_age_sec=.4), {})
        self.assertEqual(legacy[0]['key'], 'output_stale')
        self.assertEqual(modern[0]['key'], 'policy_state_incomplete')

    def test_unknown_future_policy_is_not_silently_treated_as_legacy(self):
        future = d.classify(dict(raw_runstate=-9, output_stale=1,
                                policy_version='localization-dtc-v99'), {})
        self.assertEqual(future[0]['key'], 'unsupported_policy')

    def test_suppressed_partial_radar_preserves_original_fault_and_strict_code(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            (run / 'logs').mkdir()
            (run / 'config.yaml').write_text('diagnostics: {}\n')
            (run / 'logs/run_loc_online.stderr.log').write_text(
                'I1031 12:00:00.000000 42] SAFETY_DIAGNOSTIC sequence=1 policy_version=' + d.CONFIRMED_LOSS_POLICY +
                ' raw_runstate=2 strict_runstate=2 selected_runstate=0 published_runstate=0 time_sync_fault_mask=4'
                ' partial_radar_suppressed=1 localization_healthy_for_degraded=1 output_continuous=1'
                ' output_age_sec=0.01 valid_match_age_sec=0.1 latest_match_valid=1 localization_consistency=1 localization_uncertainty=1\n')
            report = d.analyze(run, None, '指定路径', [])
        self.assertEqual(report['findings'][0]['code'], 2)
        self.assertEqual(report['findings'][0]['published'], 0)
        self.assertEqual(report['latest']['selected_runstate'], 0)
        self.assertEqual(report['partial_radar_suppressions'][0]['evidence']['strict_runstate'], 2)
        self.assertIn('不等于雷达故障已修复', d.render(report))

    def test_no_candidate_output_has_no_fake_overlap_or_three_pose_claim(self):
        report = dict(selection='指定路径', selected='/evidence', process=None, has_match_evidence=True,
                      match_failure_image=dict(status='generated', phase='startup', stage='search',
                      candidate_available=False, candidate_id=-1, retrieval_score=None,
                      reason='no_descriptor_candidates', visualization_mode='sensor_only', comparison={},
                      image='/evidence/comparison.png', html='/evidence/index.html', wall_ns=1000000000,
                      sensor_stamp_ns=1000000000, original_sample_replay_available=False))
        text = d.render(report)
        self.assertIn('初始定位 / search', text)
        self.assertIn('不叠加地图、不计算地图重合率', text)
        self.assertNotIn('命中=None', text)
        self.assertNotIn('初值 / NDT 结果', text)

    def test_incident_csv_policy_marker_preserves_new_startup_semantics(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'timeline.csv'
            path.write_text('wall_ns,raw,published,policy_version,startup_pending,startup_expired,startup_age_sec,output_stale\n'
                            '1793400000000000000,-9,-9,' + d.CONFIRMED_LOSS_POLICY + ',1,1,8.1,1\n')
            rows = d.timeline_rows(path, [])
        self.assertEqual(rows[0]['policy_version'], d.CONFIRMED_LOSS_POLICY)
        self.assertEqual(d.classify(rows[0], {})[0]['key'], 'startup_expired')


class MatchFailureSelectionTest(unittest.TestCase):
    def setUp(self):
        from test_match_failure_birdseye import fixture
        self.fixture = fixture
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.run = self.root / 'run'
        self.run.mkdir()
        self.frames = self.root / 'match_failures'
        self.report = dict(selected=str(self.run), resolved_run=str(self.run),
                           process=dict(run=str(self.run), pid=42, start_wall_ns=1000))

    def tearDown(self):
        self.tmp.cleanup()

    def rendered(self):
        return subprocess.CompletedProcess([], 0, json.dumps(dict(status='generated', image='test.png')))

    def test_only_latest_complete_same_run_and_pid_lifetime_is_selected(self):
        old, _ = self.fixture(self.frames, self.run, 900)
        expected, _ = self.fixture(self.frames, self.run, 2000)
        self.fixture(self.frames, self.root / 'other', 3000)
        self.fixture(self.frames, self.run, 4000, pid=43)
        partial, value = self.fixture(self.frames, self.run, 5000)
        value['status'] = 'capturing'
        (partial / 'summary.json').write_text(json.dumps(value))
        with patch.object(d.subprocess, 'run', return_value=self.rendered()) as run:
            result = d.match_failure_image(self.report)
        self.assertEqual(result['selection'], 'latest_same_run_frame')
        self.assertEqual(run.call_args.args[0][-1], str(expected))
        self.assertIn('-s', run.call_args.args[0])
        self.assertEqual(run.call_args.kwargs['timeout'], 30)
        self.assertEqual(run.call_args.kwargs['env']['OPENBLAS_NUM_THREADS'], '1')
        self.assertEqual(run.call_args.kwargs['env']['PYTHONNOUSERSITE'], '1')

    def test_explicit_incident_does_not_borrow_later_frame(self):
        event = self.root / 'incidents/event_1500_42'
        event.mkdir(parents=True)
        self.report.update(selected=str(event), process=None,
                           capture_windows=[dict(path=str(event), start_ns=1000, end_ns=1800)])
        expected, _ = self.fixture(self.frames, self.run, 1500)
        self.fixture(self.frames, self.run, 2000)
        with patch.object(d.subprocess, 'run', return_value=self.rendered()) as run:
            result = d.match_failure_image(self.report)
        self.assertEqual(result['selection'], 'latest_same_run_frame_in_event')
        self.assertEqual(run.call_args.args[0][-1], str(expected))
        self.report['capture_windows'] = []
        with patch.object(d.subprocess, 'run', side_effect=AssertionError('no timeframe')):
            self.assertEqual(d.match_failure_image(self.report)['status'], 'unavailable')

    def test_copied_run_uses_recorded_full_identity_not_basename(self):
        self.report['process'] = None
        original = self.root / 'original/run'
        (self.run / 'run_metadata.txt').write_text('run_dir=' + str(original) + '\n')
        expected, _ = self.fixture(self.frames, original, 2000)
        self.fixture(self.frames, self.root / 'other/run', 3000)
        with patch.object(d.subprocess, 'run', return_value=self.rendered()) as run:
            d.match_failure_image(self.report)
        self.assertEqual(run.call_args.args[0][-1], str(expected))

    def test_timeout_or_dependencies_do_not_abort_diagnosis(self):
        self.fixture(self.frames, self.run, 2000)
        with patch.object(d.subprocess, 'run', side_effect=subprocess.TimeoutExpired('figure', 30)):
            self.assertEqual(d.match_failure_image(self.report)['status'], 'timeout')
        response = subprocess.CompletedProcess([], 2, json.dumps(dict(status='dependency_missing', message='scipy')))
        with patch.object(d.subprocess, 'run', return_value=response):
            self.assertEqual(d.match_failure_image(self.report)['status'], 'dependency_missing')
        with patch.object(d.subprocess, 'run', side_effect=AssertionError('disabled')):
            self.assertEqual(d.match_failure_image(self.report, enabled=False)['status'], 'disabled')

    def test_direct_failure_package_preserves_scope(self):
        evidence, _ = self.fixture(self.frames, self.run, 2000)
        self.assertEqual(d.choose_path(evidence), evidence.resolve())
        self.report.update(selected=str(evidence), process=None)
        with patch.object(d.subprocess, 'run', return_value=self.rendered()) as run:
            result = d.match_failure_image(self.report)
        self.assertEqual(result['selection'], 'selected_frame')
        self.assertEqual(run.call_args.args[0][-1], str(evidence))

    def test_protected_quota_status_explains_missing_image(self):
        self.frames.mkdir()
        value = dict(schema=1, kind='lightning_match_failure_collector', pid=42, run_dir=str(self.run),
                     session_start_wall_ns=1100, running=True, received=2, admitted=2, saved=0,
                     rejected_disk=2, rejected_map=0, writer_errors=0, protection_days=30, disk_mb=2048,
                     last_error='protected evidence quota full')
        (self.frames / 'collector_status.json').write_text(json.dumps(value))
        result = d.match_failure_image(self.report)
        self.assertEqual(result['status'], 'not_found')
        self.assertEqual(result['collector']['status'], 'observed')
        self.assertEqual(result['collector']['counters']['rejected_disk'], 2)
        self.assertEqual(result['collector']['last_error'], 'protected evidence quota full')

    def test_collector_rejects_cross_run_reused_pid_and_historical_scope(self):
        self.frames.mkdir()
        path = self.frames / 'collector_status.json'
        value = dict(schema=1, kind='lightning_match_failure_collector', pid=42,
                     run_dir=str(self.root / 'other'), session_start_wall_ns=1100, running=True, saved=99)
        path.write_text(json.dumps(value))
        self.assertEqual(d.match_failure_collector(self.report)['status'], 'unavailable')
        value.update(run_dir=str(self.run), session_start_wall_ns=900)
        path.write_text(json.dumps(value))
        self.assertEqual(d.match_failure_collector(self.report)['status'], 'unavailable')
        value['session_start_wall_ns'] = 1100
        path.write_text(json.dumps(value))
        self.report['evidence_scope'] = 'historical'
        self.assertEqual(d.match_failure_collector(self.report)['status'], 'not_requested')

    def test_disabled_collector_is_reported_as_configuration_not_runtime_fact(self):
        (self.run / 'config.yaml').write_text('diagnostics:\n  match_failure_capture:\n    enabled: false\n')
        result = d.match_failure_collector(self.report)
        self.assertEqual(result['status'], 'configured_disabled')

    def test_custom_capture_root_from_run_metadata_keeps_identity_checks(self):
        custom = self.root / 'custom_capture'
        expected, _ = self.fixture(custom, self.run, 2000)
        self.fixture(custom, self.root / 'foreign', 3000)
        (self.run / 'run_metadata.txt').write_text('match_failure_dir=' + str(custom) + '\n')
        (custom / 'collector_status.json').write_text(json.dumps(dict(schema=1,
            kind='lightning_match_failure_collector', pid=42, run_dir=str(self.run),
            session_start_wall_ns=1100, running=True, saved=1)))
        with patch.object(d.subprocess, 'run', return_value=self.rendered()) as run:
            result = d.match_failure_image(self.report)
        self.assertEqual(run.call_args.args[0][-1], str(expected))
        self.assertEqual(result['collector']['status'], 'observed')
        self.assertEqual(result['collector']['counters']['saved'], 1)

    def test_latest_observation_remains_distinct_from_retained_snapshot(self):
        self.frames.mkdir()
        value = dict(schema=2, kind='lightning_match_failure_collector', pid=42, run_dir=str(self.run),
                     session_start_wall_ns=1100, running=True, saved=3, dropped_busy=2,
                     episode=dict(id=4, admitted=3, limit=3, suppressed_at_limit=12),
                     latest_observation=dict(wall_ns=9999, sensor_stamp_ns=9900, phase='startup', stage='search',
                     rejection_reason='no_candidates', candidate_available=False, candidate_id=-1,
                     retrieval_score=None, map_evaluated=False, map_overlap=0))
        (self.frames / 'collector_status.json').write_text(json.dumps(value))
        result = d.match_failure_image(self.report)
        self.assertEqual(result['status'], 'not_found')
        self.assertEqual(result['collector']['latest_observation']['wall_ns'], 9999)
        self.assertEqual(result['collector']['episode']['admitted'], 3)


if __name__ == '__main__':
    unittest.main()
