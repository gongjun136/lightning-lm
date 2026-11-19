'''diagnostic_summary 纯计算摘要层的行为测试（只用合成字典）。'''
import ast
import copy
import datetime as dt
import json
import math
import os
import re
import sys
import tempfile
from pathlib import Path
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import diagnostic_summary as ds  # noqa: E402
import importlib  # noqa: E402
from unittest import mock  # noqa: E402


def live(**over):
    report = {
        'schema': 2, 'selection': '自动选择最新运行', 'selected': 'runs/example', 'resolved_run': 'example',
        'process': {'pid': 1}, 'evidence_scope': 'live_run', 'has_safety_evidence': True, 'latest_age_sec': 1.0,
        'latest': {'_time': '00:00:01.000', '_source': 'safety.log:100', '_scope': 'live_run', 'raw_runstate': 0,
                   'strict_runstate': 0, 'selected_runstate': 0, 'published_runstate': 0, 'policy_version': 3},
        'current_findings': [], 'findings': [], 'warnings': []}
    report.update(over)
    return report


def finding(code, key, evidence=None, **over):
    item = {'code': code, 'key': key, 'title': str(key) + ' 标题', 'confidence': '判据已确认', 'evidence': evidence or {},
            'advice': '', 'time': '00:00:01', 'source': 'safety.log:1', 'scope': 'live_run', 'published': code}
    item.update(over)
    return item


def episode(code, key, start, evidence=None, **over):
    item = finding(code, key, evidence)
    item.update({'start_ns': start, 'last_ns': start + 10, 'end_ns': None, 'first_time': 't%d' % start,
                 'last_time': 't%d' % (start + 10), 'end_time': None, 'state': 'open', 'samples': 3,
                 'capture_coverage': 'range_overlap', 'event_paths': [], 'input_timing_evidence': [],
                 'onset_context': [], 'context': []})
    item.update(over)
    return item


def delay(source_gap, arrival_gap, when='t1'):
    return {'category': 'delivery_delay', 'stream': 'lidar', 'topic': '/points', 'time': when, 'source': 'log:9',
            'source_gap_sec': source_gap, 'arrival_gap_sec': arrival_gap, 'threshold_sec': 0.5}


def lines_with(text, needle):
    return [line for line in text.split(ds.NL) if needle in line]


def complex_report():
    return live(
        thresholds={'output_timeout_sec': 0.5},
        current_findings=[finding(4, 'major_sigma', {'major_sigma': float('nan')}),
                          finding('x', None, {'nested': {'a': [1, 2]}}, confidence='证据不足')],
        findings=[episode(2, 'frame_supply', 10, {'lidar_frame_age_sec': '0.2'},
                          input_timing_evidence=[delay(0.1, float('inf'))]),
                  episode(-9, 'startup_expired', 5, {'startup_age_sec': 31}, phase='startup', state='observation_gap')],
        warnings=['截断'], match_snapshots={'records': [{'saved_path': 'snap/a.pcd', 'current_status': 'kept'}]},
        input_timing={'status': 'unavailable', 'limitation': '未启用'})


class FreshnessTest(unittest.TestCase):
    def state(self, **over):
        return ds.build_summary(live(**over))['freshness']

    def test_five_second_boundary(self):
        self.assertTrue(self.state(latest_age_sec=5.0)['is_current'])
        self.assertTrue(self.state(latest_age_sec=0)['is_current'])
        self.assertTrue(self.state(latest_age_sec='4.5')['is_current'])
        stale = self.state(latest_age_sec=5.0001)
        self.assertEqual(stale['state'], 'stale')
        self.assertIn('当前未知', stale['text'])

    def test_future_missing_and_nonfinite_are_not_fresh(self):
        self.assertEqual(self.state(latest_age_sec=-0.2)['state'], 'clock_anomaly')
        for bad in (None, 'nan', 'inf', float('nan'), 'abc'):
            result = self.state(latest_age_sec=bad)
            self.assertFalse(result['is_current'], bad)
            self.assertIn('当前未知', result['text'])

    def test_requires_process_safety_and_live_scope(self):
        self.assertFalse(self.state(process=None)['is_current'])
        self.assertEqual(self.state(has_safety_evidence=False)['state'], 'unavailable')
        self.assertFalse(self.state(evidence_scope='unavailable')['is_current'])
        self.assertEqual(self.state(evidence_scope='historical')['state'], 'offline')

    def test_is_fresh_live_shared_rule(self):
        self.assertTrue(ds.is_fresh_live(live()))
        self.assertFalse(ds.is_fresh_live(live(latest_age_sec=6)))
        self.assertFalse(ds.is_fresh_live(live(evidence_scope='explicit')))


class RunstateTest(unittest.TestCase):
    def test_offline_zero_only_describes_last_record(self):
        text = ds.render_summary(live(evidence_scope='explicit', process=None))
        self.assertIn('包内末条记录', text)
        self.assertIn('所读最新记录（非当前）', text)
        self.assertNotIn('当前记录：', text)
        self.assertIn('只说明该记录为 0', text)

    def test_tolerated_raw4_published0_not_healthy(self):
        report = live(current_findings=[finding(4, 'match_rejected_grace_suppressed',
                                                {'match_failure_grace_age_sec': 0.12, 'raw_runstate': 4})])
        report['latest'] = dict(report['latest'], raw_runstate='4', strict_runstate=4, match_failure_grace_suppressed=True)
        verdict = ds.build_summary(report)['latest_record']['verdict']
        self.assertIn('不是健康', verdict)
        self.assertNotIn('只说明该记录为 0', verdict)
        text = ds.render_summary(report)
        self.assertIn('原始码 4', text)
        self.assertIn('上报码 0', text)
        self.assertIn('匹配失败容忍已用时 实测 0.12 s / 阈值未记录', text)
        self.assertIn('200ms 仅为最大容忍', text)

    def test_unknown_code_and_missing_rule_visible(self):
        report = live(current_findings=[finding(7, 'future_key', {'alpha': 1}, confidence='规则缺失')])
        report['latest'] = dict(report['latest'], raw_runstate=7)
        self.assertIn('未知运行码', ds.build_summary(report)['latest_record']['verdict'])
        text = ds.render_summary(report)
        for needle in ('未知码 7', 'future_key', '规则缺失', 'alpha=1', '不能视为正常'):
            self.assertIn(needle, text)

    def test_zero_with_findings_not_healthy(self):
        report = live(current_findings=[finding(4, 'major_sigma', {'major_sigma': 0.3})])
        self.assertIn('不能判为健康', ds.build_summary(report)['latest_record']['verdict'])

    def test_startup_pending_kept_as_risk(self):
        report = live(current_findings=[finding(4, 'startup_pending', {'startup_age_sec': 3, 'startup_grace_sec': 30,
                                                                       'startup_pending': True}, phase='startup')])
        report['latest'] = dict(report['latest'], raw_runstate=4, strict_runstate=4, selected_runstate=4, published_runstate=4)
        text = ds.render_summary(report)
        self.assertIn('启动已用时 实测 3 s / 阈值 30 s（证据）', text)
        self.assertIn('不等于无害', text)
        self.assertIn('该记录为故障状态', text)


class FindingFormatTest(unittest.TestCase):
    def test_multiple_independent_code4_criteria(self):
        report = live(thresholds={'output_timeout_sec': 0.5, 'major_sigma_limit_m': 0.2}, current_findings=[
            finding(4, 'match_stale_degraded', {'valid_match_age_sec': 1.25, 'localization_output_max_lidar_age_sec': 1.0,
                                                'match_stale': 1}),
            finding(4, 'output_stale_degraded', {'output_age_sec': 'inf', 'first_valid_output_missing': True}),
            finding(4, 'major_sigma', {'major_sigma': 0.35, 'covariance_valid': True})])
        text = ds.render_summary(report)
        self.assertIn('有效匹配年龄 实测 1.25 s / 阈值 1 s（证据）', text)
        self.assertIn('输出年龄 实测 +inf / 阈值 0.5 s（有效配置）', text)
        self.assertIn('可能尚无首次有效输出', text)
        self.assertIn('主轴 sigma 实测 0.35 m / 阈值 0.2 m（有效配置）', text)
        self.assertEqual(text.count('[C'), 3)

    def test_threshold_precedence_missing_and_no_relation(self):
        report = live(thresholds={'radar_frame_timeout_sec': 1.0}, current_findings=[finding(2, 'frame_supply', {
            'lidar_frame_age_sec': 0.05, 'radar_frame_timeout_sec': 0.5, 'lidar_frame_recovery_count': 3})])
        text = ds.render_summary(report)
        self.assertIn('雷达帧年龄 实测 0.05 s / 阈值 0.5 s（证据）', text)
        self.assertIn('恢复确认帧 实测 3 帧 / 阈值未记录', text)
        self.assertIn('不推断大小关系', text)
        line = lines_with(text, '雷达帧年龄')[0]
        self.assertNotIn('>', line)
        self.assertNotIn('超过', line)

    def test_yaw_units_and_dynamic_threshold(self):
        report = live(current_findings=[
            finding(4, 'yaw_sigma_rad', {'yaw_sigma_rad': 0.1, 'yaw_sigma_limit_deg': 1.5}),
            finding(4, 'yaw_condition', {'yaw_evidence_available': True, 'yaw_published_error_rad': 0.1,
                                         'yaw_tolerance_deg': 3})])
        text = ds.render_summary(report)
        self.assertIn('航向 sigma 实测 5.729577951 deg（原值 0.1 rad） / 阈值 1.5 deg（证据）', text)
        self.assertIn('发布航向误差 实测 0.1 rad / 阈值未记录', text)
        self.assertIn('yaw_tolerance_deg 不是最终阈值', text)
        self.assertNotIn('阈值 3', text)

    def test_matching_rejection_score_and_iterations(self):
        rejection = {'reason_code': 'low_score', 'title': '匹配分数低', 'time': 't', 'source': 'log:2', 'confidence': 0.42,
                     'evidence': {'iterations': 30, 'ndt_converged': False, 'min_confidence': 0.6,
                                  'consecutive_rejected_matches': 2, 'lost_frame_threshold': 5}}
        report = live(current_findings=[finding(4, 'match_rejected_unconfirmed', {
            'latest_match_valid': False, 'lost_frames': 2, 'lost_frame_threshold': 5}, matching_rejection=rejection)])
        text = ds.render_summary(report)
        for needle in ('拒绝项 low_score 匹配分数低', '迭代数=30', 'NDT收敛=false', '匹配分数=0.42 / 门槛=0.6',
                       '连续拒绝=2 / 丢失阈值=5', '匹配分数不是可信等级', '判据已确认',
                       '连续丢失帧 实测 2 帧 / 阈值 5 帧（证据）'):
            self.assertIn(needle, text)

    def test_missing_key_keeps_insufficient_evidence(self):
        text = ds.render_summary(complex_report())
        self.assertIn('未知码(x)', text)
        self.assertIn('未知判据', text)
        self.assertIn('证据不足', text)


class HistoryTest(unittest.TestCase):
    def test_history_does_not_become_current(self):
        report = live(findings=[episode(-9, 'confirmed_tracking_lost', 100, {'lost_frames': 10, 'lost_frame_threshold': 5})])
        summary = ds.build_summary(report)
        self.assertIn('只说明该记录为 0', summary['latest_record']['verdict'])
        self.assertEqual(summary['current_findings'], [])
        text = ds.render_summary(report)
        self.assertIn('当前记录的判据（0 项）', text)
        self.assertIn('未见结束（不代表当前仍在故障）', text)

    def test_episode_states_are_distinct(self):
        states = ['open', 'recovered', 'predicate_changed', 'observation_gap']
        eps = [episode(4, 'k%d' % i, i * 100, state=st) for i, st in enumerate(states)]
        text = ds.render_summary(live(findings=eps))
        for value in ds.STATE_TEXT.values():
            self.assertIn(value, text)
        gap_line = [line for line in lines_with(text, '[H4]') if '码' in line][0]
        self.assertIn('不能判为恢复', gap_line)
        self.assertNotIn('回 0', gap_line)

    def test_first_code4_onset_visible_despite_later_minus9(self):
        eps = [episode(4, 'match_stale_degraded', 100, {'valid_match_age_sec': 2.0},
                       onset_evidence={'valid_match_age_sec': 1.5, 'localization_output_max_lidar_age_sec': 1.0}),
               episode(-9, 'confirmed_tracking_lost', 200, {'lost_frames': 9, 'lost_frame_threshold': 5}),
               episode(4, 'major_sigma', 300, {'major_sigma': 1.0})]
        text = ds.render_summary(live(findings=eps, event_status='triggered', trigger_time='t100', trigger_reasons=['4']))
        onset = [line for line in text.split(ds.NL) if line.startswith('  码 4 [H1]')]
        self.assertEqual(len(onset), 1)
        self.assertIn('有效匹配年龄 实测 1.5 s / 阈值 1 s（证据）', onset[0])
        self.assertIn('首次标题未记录', onset[0])
        self.assertIn('触发原因：4', text)

    def test_grouping_limit_counts_and_index(self):
        eps = [episode(4, 'key%d' % i, i * 100) for i in range(8)]
        eps.append(episode(4, 'key7', 50, source='timeline:1'))
        report = live(findings=eps)
        summary = ds.build_summary(report)
        self.assertEqual(len(summary['history']['groups']), 8)
        self.assertEqual(summary['history']['omitted_groups'], 2)
        text = ds.render_summary(report)
        self.assertIn('另有 2 种未展开', text)
        self.assertIn('不代表物理故障次数', text)
        self.assertIn('同类 2 段', text)
        self.assertIn('H8、H9', ds.render_index(summary))
        index = ds.render_index(summary)
        for i in range(1, 10):
            self.assertIn('[H%d]' % i, index)


class TimingTest(unittest.TestCase):
    def test_delivery_delay_same_pair_without_dds_claim(self):
        ep = episode(4, 'match_stale_degraded', 100, input_timing_evidence=[
            delay(0.1, 0.9), delay(0.1, 2.3, 't2'), {'category': 'callback_silence', 'stream': 'lidar'}])
        text = ds.render_summary(live(findings=[ep]))
        self.assertIn('雷达延迟交付 2 条', text)
        self.assertIn('source_gap=0.1 s / arrival_gap(回调到达间隔)=2.3 s / 取证阈值=0.5 s，时间 t2', text)
        self.assertIn('[H1.t2]', text)
        self.assertIn('已观察到延迟交付，环节未定', text)
        for candidate in ds.DELAY_CANDIDATES:
            self.assertIn(candidate, text)
        self.assertIn('不是绝对网络单向延迟', text)
        self.assertIn('不能确认 DDS', text)
        self.assertNotIn('DDS 根因', text)
        self.assertNotIn('已确认 DDS', text)

    def test_source_gap_only_is_not_dds_or_hardware_stop(self):
        ep = episode(2, 'frame_supply', 100, input_timing_evidence=[
            {'category': 'source_gap', 'stream': 'lidar', 'source_gap_sec': 1.2, 'arrival_gap_sec': 1.2, 'threshold_sec': 0.5}])
        text = ds.render_summary(live(findings=[ep]))
        self.assertIn('源间隔超限 1 条', text)
        self.assertIn('不等于硬件停发', text)
        self.assertNotIn('环节未定', text)
        self.assertNotIn('DDS', text)

    def test_numeric_and_startup_not_attributed_to_global_timing(self):
        timing = {'status': 'available', 'counts': {'delivery_delay': 4}, 'events': [delay(0.1, 3.0)], 'observed_records': 10}
        eps = [episode(4, 'major_sigma', 100, {'major_sigma': 0.5}),
               episode(-9, 'startup_expired', 200, {'startup_age_sec': 31, 'startup_grace_sec': 30}, phase='startup')]
        text = ds.render_summary(live(findings=eps, input_timing=timing))
        self.assertNotIn('环节未定', text)
        self.assertIn('所示段无关联输入时序证据', text)
        self.assertIn('不参与归因', text)
        self.assertIn('启动阶段', text)
        self.assertIn('startup_expired', text)


class OtherEvidenceTest(unittest.TestCase):
    def test_upstream_driver_independent_even_without_safety(self):
        report = live(has_safety_evidence=False, latest={}, localization_error='未找到 Safety 日志', upstream_lidar={
            'status': 'ok', 'message': '驱动报告', 'report': {
                'active': True, 'has_evidence': True, 'warnings': ['驱动报告过旧'],
                'findings': [{'severity': 'error', 'role': 'front', 'stage': 'udp', 'message': '历史组帧超时', 'evidence': {'n': 3}}],
                'roles': [{'role': 'front', 'device': {'evidence_scope': 'historical', 'collected_at': 'c1', 'age_sec': 900}}]}})
        text = ds.render_summary(report)
        for needle in ('未找到 Safety 日志', '[U1] error front/udp：历史组帧超时', '[UW1] 驱动报告过旧',
                       '不拼接为因果', 'active 或退出码 0 不等于健康', '运行码未记录'):
            self.assertIn(needle, text)
        self.assertNotIn('环节未定', text)

    def test_all_risk_categories_visible_and_indexed(self):
        report = live(
            warnings=['合同版本不一致', '配置失效：参数未加载', '日志截断'],
            launch_findings=[{'title': '启动失败', 'stage': 'configure', 'evidence': {'rc': 1}, 'source': 'launch.log:3', 'advice': '检查'}],
            partial_radar_suppressions=[{'time': 't5', 'scope': 'live_run', 'evidence': {'device': 'rear'}, 'source': 's:1'}],
            coexisting_radar_conditions=[{'key': 'radar_time_sync', 'device': 'left', 'last_time': 't6', 'scope': 'live_run',
                                          'evidence': {}, 'source': 's:2'}],
            runtime_protections=[{'key': 'grace', 'title': '匹配失败容忍', 'count': 4, 'first_time': 't1', 'last_time': 't2',
                                  'evidence': {}, 'source': 's:3', 'text': '容忍 4 次'}],
            runtime_build={'status': 'mismatch', 'compiled': {'source_revision': 'abc', 'source_dirty': True, 'built_at': 'b'},
                           'expected_contract': 'c2'},
            motion_evidence={'status': 'UNKNOWN', 'unknown_advisory': True},
            match_failure_image={'status': 'disabled', 'message': '已用 --no-match-image 禁用', 'collector': {
                'status': 'error', 'running': False, 'fresh': False, 'last_error': '配额已满', 'counters': {'rejected_disk': 7}}},
            input_timing={'status': 'available', 'counts': {}, 'events': [], 'truncated_events': 12, 'lost_by_topic': {'/points': 3}},
            collector_totals={'dropped': 2})
        text = ds.render_summary(report)
        for needle in ('[W1] 合同版本不一致', '[W2] 配置失效', '[W3] 日志截断', '启动失败×1 [L1]', '[S1]',
                       'radar_time_sync×1 [X1]', 'grace×4 [P1]', '运行时保护 1', '不匹配', '未知不等于通过', 'rejected_disk=7',
                       '配额已满', '时序事件截断 12', '按话题丢失 /points=3', 'dropped=2', '--no-match-image'):
            self.assertIn(needle, text, needle)
        index = ds.render_index(ds.build_summary(report))
        for ref in ('[W1]', '[W2]', '[W3]', '[L1]', '[S1]', '[X1]', '[P1]'):
            self.assertIn(ref, index)


class PurityTest(unittest.TestCase):
    def test_input_not_mutated_deterministic_and_json_safe(self):
        report = complex_report()
        before = copy.deepcopy(report)
        first = ds.render_summary(report)
        second = ds.render_summary(report)
        summary = ds.build_summary(report)
        index = ds.render_index(summary)
        self.assertEqual(report, before)
        self.assertEqual(first, second)
        json.dumps(summary, allow_nan=False)
        self.assertEqual(json.dumps(summary, sort_keys=True), json.dumps(ds.build_summary(report), sort_keys=True))
        self.assertIn('snap/a.pcd', first)
        self.assertIn('[M1]', index)
        self.assertIn('[H1.t1]', index)
        self.assertTrue(all(len(line) < 800 for line in index.split(ds.NL)))

    def test_default_output_is_concise(self):
        lines = ds.render_summary(live()).split(ds.NL)
        self.assertLessEqual(len(lines), 15)
        self.assertNotIn(chr(27), ds.render_summary(live()))
        self.assertIn(chr(27), ds.render_summary(live(), color=True))

    def test_module_is_pure_computation(self):
        import translation_diagnostics as td
        for module, allowed in ((ds, {'json', 'math', 'translation_diagnostics'}),
                                (td, {'__future__', 'math', 'collections', 'typing'})):
            with self.subTest(module=module.__name__):
                with open(module.__file__, encoding='utf-8') as handle:
                    tree = ast.parse(handle.read())
                imported = set()
                for node in ast.walk(tree):
                    if isinstance(node, ast.Import):
                        imported.update(alias.name.split('.')[0] for alias in node.names)
                    elif isinstance(node, ast.ImportFrom):
                        imported.add((node.module or '').split('.')[0])
                    elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                        self.assertNotIn(node.func.id, ('open', 'exec', 'eval', 'compile', '__import__'))
                self.assertTrue(imported <= allowed, imported)


class ColorTest(unittest.TestCase):
    def latest_line(self, report):
        return next(line for line in ds.render_summary(report, color=True).splitlines() if '最新记录：' in line)

    def test_current_fault_normal_and_suppression_have_distinct_colors(self):
        self.assertTrue(self.latest_line(live()).startswith('\x1b[1;32m'))
        fault = live()
        fault['latest'].update(raw_runstate=4, strict_runstate=4, selected_runstate=4, published_runstate=4)
        self.assertTrue(self.latest_line(fault).startswith('\x1b[1;31m'))
        fault['latest'].update(selected_runstate=0, published_runstate=0, match_failure_grace_suppressed=1)
        self.assertTrue(self.latest_line(fault).startswith('\x1b[1;33m'))
        zero_suppressed = live()
        zero_suppressed['latest']['match_failure_grace_suppressed'] = 1
        self.assertTrue(self.latest_line(zero_suppressed).startswith('\x1b[1;33m'))

    def test_offline_stale_and_incomplete_records_are_never_green(self):
        for report in (live(evidence_scope='explicit'), live(latest_age_sec=6),
                       live(latest={'raw_runstate': 0}), live(latest_age_sec=-1)):
            self.assertTrue(self.latest_line(report).startswith('\x1b[1;33m'))
        offline_fault = live(evidence_scope='explicit')
        offline_fault['latest'].update(raw_runstate=-9, published_runstate=-9)
        self.assertTrue(self.latest_line(offline_fault).startswith('\x1b[1;33m'))

    def test_colors_preserve_text_and_highlight_observation_and_threshold(self):
        report = rich_report()
        colored = ds.render_summary(report, color=True)
        plain = ds.render_summary(report)
        self.assertEqual(re.sub(r'\x1b\[[0-9;]*m', '', colored), plain)
        self.assertNotIn('\x1b', plain)
        self.assertIn('\x1b[1m实测 0.61 s / 阈值 0.5 s', colored)
        self.assertIn('\x1b[1;36m历史判据', colored)
        self.assertIn('\x1b[33m  [W1]', colored)
        self.assertNotIn('\x1b', json.dumps(ds.build_summary(report), ensure_ascii=False))


class WiringTest(unittest.TestCase):
    def read(self, name):
        path = os.path.join(HERE, name)
        if not os.path.exists(path):
            self.skipTest(name + ' 不在本目录')
        with open(path, encoding='utf-8') as handle:
            return handle.read()

    def test_cli_wiring(self):
        src = self.read('diagnose_localization.py')
        for needle in ('import diagnostic_summary', 'diagnostic_summary.build_summary(report)',
                       'diagnostic_summary.render_summary(report', 'diagnostic_summary.render_index(diag_summary)',
                       'diagnostic_summary.is_fresh_live(report)', 'refresh_live_status(report)',
                       'fresh=0 <= age <= 10', 'onset_title='):
            self.assertIn(needle, src)

    def test_analyze_keeps_negative_age(self):
        src = self.read('diagnose_localization.py')
        func = next((node for node in ast.walk(ast.parse(src))
                     if isinstance(node, ast.FunctionDef) and node.name == 'analyze'), None)
        self.assertIsNotNone(func)
        body = ast.get_source_segment(src, func)
        self.assertIn('age < 0', body)
        self.assertNotIn('max(0, (time.time_ns() - row_ns', body)

    def test_refresh_is_read_only(self):
        src = self.read('diagnose_localization.py')
        func = next((node for node in ast.walk(ast.parse(src))
                     if isinstance(node, ast.FunctionDef) and node.name == 'refresh_live_status'), None)
        self.assertIsNotNone(func)
        body = ast.get_source_segment(src, func)
        for forbidden in ('rclpy', 'subprocess', 'open(', '.write', 'sleep', 'unlink'):
            self.assertNotIn(forbidden, body)

    def test_contract_lists_module(self):
        self.assertIn('scripts/diagnostic_summary.py', self.read('check_diagnostic_contract.py'))


NOW_NS = 2000000000 * 10 ** 9
PROC_INFO = {'pid': 42, 'start_ticks': 1000, 'exe': '/opt/loc/bin/node', 'run': 'runs/r1',
             'log': 'runs/r1/safety.log', 'config': 'cfg.yaml', 'start_wall_ns': 1}


def rich_report():
    shared = dict(delay(0.1, 2.0, 't101'), observed_wall_ns=101)
    other = dict(delay(0.1, 1.1, 't102'), observed_wall_ns=102)
    gap = {'category': 'source_gap', 'stream': 'lidar', 'topic': '/points', 'source': 'log:9', 'time': 't103',
           'observed_wall_ns': 103, 'source_gap_sec': 1.5, 'threshold_sec': 0.5}
    eps = [
        episode(2, 'frame_supply', 100, {'lidar_frame_age_sec': 0.2, 'radar_frame_timeout_sec': 0.5},
                title='雷达帧恢复确认中', onset_title='雷达帧供给中断', state='recovered',
                onset_evidence={'lidar_frame_age_sec': 0.61, 'radar_frame_timeout_sec': 0.5},
                input_timing_evidence=[shared, gap]),
        episode(4, 'match_stale_degraded', 105, {'valid_match_age_sec': 1.4}, state='predicate_changed',
                onset_evidence={'valid_match_age_sec': 1.2, 'localization_output_max_lidar_age_sec': 1.0},
                input_timing_evidence=[shared, other]),
        episode(-9, 'confirmed_tracking_lost', 120, {'lost_frames': 12, 'lost_frame_threshold': 10},
                state='observation_gap', capture_coverage='no_retained_overlap'),
        episode(4, 'major_sigma', 130, {'major_sigma': 0.4, 'major_sigma_limit_m': 0.3}),
        episode(-9, 'startup_expired', 10, {'startup_age_sec': 31, 'startup_grace_sec': 30}, phase='startup',
                state='recovered'),
        episode(2, 'frame_supply', 102, {'lidar_frame_age_sec': 0.2}, last_ns=105, source='timeline:5',
                input_timing_evidence=[shared]),
    ]
    report = live(
        findings=eps, event_status='triggered', trigger_time='t100', trigger_reasons=['2', '4'],
        current_findings=[finding(4, 'match_stale_degraded', {'valid_match_age_sec': 1.3,
                                                              'localization_output_max_lidar_age_sec': 1.0}),
                          finding(4, 'major_sigma', {'major_sigma': 0.4, 'major_sigma_limit_m': 0.3})],
        warnings=['合同版本不一致', '配置失效：参数未加载', '日志截断'],
        input_timing={'status': 'available', 'counts': {'delivery_delay': 2, 'source_gap': 1},
                      'events': [shared, other, gap], 'observed_records': 500},
        upstream_lidar={'status': 'ok', 'message': '驱动报告', 'report': {
            'active': False, 'findings': [{'severity': 'warn', 'role': 'front', 'stage': 'frame', 'message': '组帧超时'}]}},
        match_failure_image={'status': 'ok', 'image': 'img/a.png', 'html': 'img/a.html', 'synthetic_fixture': False,
                             'check_score_matches_recorded': False, 'original_sample_replay_available': False},
        launch_findings=[{'title': '启动失败', 'stage': 'configure'}],
        coexisting_radar_conditions=[{'key': 'radar_time_sync', 'device': 'left'},
                                     {'key': 'radar_time_sync', 'device': 'rear'}],
        runtime_build={'status': 'unverified'})
    report['latest'] = dict(report['latest'], raw_runstate=4, strict_runstate=4, selected_runstate=4, published_runstate=4)
    return report


class CompactAccuracyTest(unittest.TestCase):
    def test_rich_report_compact_and_deduplicated(self):
        report = rich_report()
        text = ds.render_summary(report)
        lines = text.split(ds.NL)
        self.assertLessEqual(len(lines), 45, text)
        for needle in ('去重 3 条', '雷达延迟交付 2 条', '[H2.t1]', '源间隔超限 1 条', '[H1.t2]', '[W1] 合同版本不一致',
                       '[W2] 配置失效', '[W3] 日志截断', 'radar_time_sync×2', '启动失败×1 [L1]',
                       '历史驱动包（当前未发现雷达 worker）', '运行构建未验证', '未验证前仅作可视化参考',
                       '起始时刻未找到保留采集覆盖（不代表整段未采集）', '不代表物理故障次数'):
            self.assertIn(needle, text)
        self.assertEqual(text.count('环节未定'), 1)
        self.assertEqual(text.count('DDS 传输/接收队列'), 1)
        onset = lines_with(text, '码 2 [H1] 首个触发')
        self.assertEqual(len(onset), 1)
        self.assertIn('雷达帧供给中断', onset[0])
        self.assertNotIn('恢复确认中', onset[0])
        self.assertIn('雷达帧年龄 实测 0.61 s / 阈值 0.5 s（证据）', onset[0])
        index = ds.render_index(ds.build_summary(report))
        for needle in ('[H6]', '[H2.t1] 同 [H1.t1]', '[W3]', '[X2]', '[L1]'):
            self.assertIn(needle, index)

    def test_coverage_wording_is_onset_only(self):
        report = live(findings=[episode(4, 'a', 100, capture_coverage='range_overlap'),
                                episode(4, 'b', 200, capture_coverage='no_retained_overlap')])
        text = ds.render_summary(report)
        index = ds.render_index(ds.build_summary(report))
        self.assertIn('起始时刻未找到保留采集覆盖（不代表整段未采集）', text)
        self.assertIn('起始时刻落在保留采集范围内（不保证整段完整）', index)
        self.assertNotIn('保留采集覆盖该段', text + index)

    def test_numeric_precision_keeps_trigger_values(self):
        report = live(current_findings=[
            finding(4, 'match_stale_degraded', {'valid_match_age_sec': 0.5123456, 'localization_output_max_lidar_age_sec': 0.5}),
            finding(4, 'output_stale', {'output_age_sec': 0.5000001, 'output_timeout_sec': 0.5})])
        text = ds.render_summary(report)
        self.assertIn('有效匹配年龄 实测 0.5123456 s / 阈值 0.5 s（证据）', text)
        self.assertIn('输出年龄 实测 0.5000001 s / 阈值 0.5 s（证据）', text)
        self.assertEqual(ds.fmt(1700000000123456789), '1700000000123456789')
        self.assertEqual(ds.fmt('1700000000123456789'), '1700000000123456789')

    def test_invalid_evidence_threshold_not_replaced_by_config(self):
        report = live(thresholds={'localization_output_max_lidar_age_sec': 1.0}, current_findings=[
            finding(4, 'match_stale_degraded', {'valid_match_age_sec': 0.7, 'localization_output_max_lidar_age_sec': 'nan'})])
        text = ds.render_summary(report)
        self.assertIn('证据中阈值为 非有限值(nan)，不可验证', text)
        self.assertNotIn('阈值 1 s', text)

    def test_missing_runstates_not_called_all_zero(self):
        report = live(latest={'_time': 't', '_source': 's:1', 'raw_runstate': 0})
        verdict = ds.build_summary(report)['latest_record']['verdict']
        self.assertIn('已记录的运行码为 0', verdict)
        self.assertIn('上报码', verdict)
        self.assertNotIn('各运行码为 0', verdict)

    def test_scope_wording_for_history_and_upstream(self):
        text = ds.render_summary(live(evidence_scope='historical', process=None, upstream_lidar={
            'status': 'ok', 'report': {'active': False, 'findings': []}}))
        self.assertIn('当前未发现定位进程，包内状态码不代表当前状态', text)
        self.assertIn('历史驱动包（当前未发现雷达 worker）', text)
        self.assertIn('不能据此确认当前设备或数据流正常', text)
        active = ds.render_summary(live(upstream_lidar={'status': 'ok', 'report': {'active': True, 'findings': []}}))
        self.assertIn('当前驱动运行实例', active)
        self.assertIn('不能据此确认当前设备或数据流正常', active)

    def test_image_flags_kept(self):
        text = ds.render_summary(live(match_failure_image={
            'status': 'ok', 'image': 'img/a.png', 'synthetic_fixture': True,
            'check_score_matches_recorded': False, 'original_sample_replay_available': True}))
        for needle in ('img/a.png', 'synthetic_fixture=true', '未验证前仅作可视化参考'):
            self.assertIn(needle, text)


def refresh_report():
    old = {'_wall_ns': NOW_NS - 30 * 10 ** 9, '_time': 'old', '_source': 'safety.log:1', 'raw_runstate': 4,
           'strict_runstate': 4, 'selected_runstate': 4, 'published_runstate': 4}
    return live(process=dict(PROC_INFO), latest=old, latest_age_sec=30.0,
                current_findings=[finding(4, 'output_stale', {'output_age_sec': 30})],
                findings=[episode(4, 'output_stale', 100)],
                input_timing={'status': 'available', 'counts': {'delivery_delay': 1}, 'events': [delay(0.1, 2.0)]},
                warnings=['最新诊断记录已过去 30.0 秒，不能作为当前实时状态。', '其他提示'],
                thresholds={'output_timeout_sec': 0.5})


class LiveRefreshTest(unittest.TestCase):
    def setUp(self):
        try:
            self.dl = importlib.import_module('diagnose_localization')
        except Exception as exc:  # 缺少运行依赖时跳过，不伪造通过
            self.skipTest('无法导入 diagnose_localization：%s' % exc)

    def _refresh(self, report, proc_lists, rows):
        reads = []
        procs = iter(proc_lists)

        def fake_rows(path, budget, warnings, *args, **kwargs):
            reads.append((str(path), budget, warnings))
            warnings.append('补读窗口限制：尾部 0MiB，序号缺口 1')
            return rows, {'limit_from_log': 2}

        def fake_classify(row, limits):
            return [{'code': 4, 'key': 'refreshed'}] if row.get('raw_runstate') else []

        dl = self.dl
        with mock.patch.object(dl, 'find_processes', side_effect=lambda *a, **k: next(procs)) as find, mock.patch.object(dl, 'log_rows', side_effect=fake_rows), mock.patch.object(dl, 'classify', side_effect=fake_classify), mock.patch.object(dl, 'row_ns', side_effect=lambda row: row.get('_wall_ns')):
            dl.refresh_live_status(report, now_ns=NOW_NS)
        return reads, find.call_count

    def assert_unknown(self, report, reason):
        meta = report['live_status_refresh']
        self.assertEqual(meta['status'], 'failed')
        self.assertIn(reason, meta['reason'])
        self.assertEqual(report['latest']['_time'], 'old')
        fresh = ds.build_summary(report)['freshness']
        self.assertFalse(fresh['is_current'])
        self.assertIn('补读失败', fresh['text'])

    def test_new_record_after_slow_analysis(self):
        report = refresh_report()
        history = copy.deepcopy((report['findings'], report['input_timing']))
        new = {'_wall_ns': NOW_NS - 10 ** 9, '_time': 'new', '_source': 'safety.log:900', 'raw_runstate': 0,
               'strict_runstate': 0, 'selected_runstate': 0, 'published_runstate': 0}
        reads, finds = self._refresh(report, [[dict(PROC_INFO)], [dict(PROC_INFO)]], [report['latest'], new])
        self.assertEqual(report['live_status_refresh']['status'], 'refreshed')
        self.assertEqual(report['latest']['_time'], 'new')
        self.assertEqual(report['latest_age_sec'], 1.0)
        self.assertEqual(report['current_findings'], [])
        self.assertEqual(report['analyzed_latest']['latest']['_time'], 'old')
        self.assertEqual((report['findings'], report['input_timing']), history)
        self.assertEqual(finds, 2)
        self.assertEqual(len(reads), 1)
        self.assertEqual(reads[0][0], PROC_INFO['log'])
        self.assertTrue(128 * 1024 <= reads[0][1] <= 256 * 1024)
        self.assertIsNot(reads[0][2], report['warnings'])
        self.assertIn('补读窗口限制', report['live_status_refresh']['limitations'][0])
        self.assertFalse(any('补读窗口限制' in w for w in report['warnings']))
        self.assertTrue(report['warnings'][0].startswith('历史分析快照：'))
        fresh = ds.build_summary(report)['freshness']
        self.assertTrue(fresh['is_current'])
        self.assertIn('有界补读', fresh['text'])

    def test_pid_reuse_is_unknown(self):
        report = refresh_report()
        reads, _ = self._refresh(report, [[dict(PROC_INFO, start_ticks=2)]], [])
        self.assertEqual(reads, [])
        self.assert_unknown(report, 'PID')

    def test_rejection_evidence_does_not_move_to_a_new_record(self):
        for advance in (0, 1):
            with self.subTest(record_advanced=bool(advance)):
                report = refresh_report()
                rejection = {'reason_code': 'ndt_not_converged', 'ndt_iterations': 20}
                report['current_findings'] = [finding(4, 'refreshed', matching_rejection=rejection)]
                row = dict(report['latest'], _wall_ns=report['latest']['_wall_ns'] + advance)
                self._refresh(report, [[dict(PROC_INFO)], [dict(PROC_INFO)]], [row])
                current = report['current_findings'][0]
                if advance:
                    self.assertNotIn('matching_rejection', current)
                else:
                    self.assertEqual(current['matching_rejection'], rejection)

    def test_exited_process_is_unknown(self):
        report = refresh_report()
        self._refresh(report, [[]], [])
        self.assert_unknown(report, '已退出')

    def test_log_switch_after_read_is_unknown(self):
        report = refresh_report()
        new = {'_wall_ns': NOW_NS - 10 ** 9, '_time': 'new', 'raw_runstate': 0}
        self._refresh(report, [[dict(PROC_INFO)], [dict(PROC_INFO, log='runs/r2/safety.log')]], [new])
        self.assert_unknown(report, '切换')

    def test_missing_or_future_time_is_unknown(self):
        report = refresh_report()
        self._refresh(report, [[dict(PROC_INFO)], [dict(PROC_INFO)]], [{'_time': 'x', 'raw_runstate': 0}])
        self.assert_unknown(report, '缺少有效时间')
        report = refresh_report()
        future = {'_wall_ns': NOW_NS + 3 * 10 ** 9, '_time': 'future', 'raw_runstate': 0}
        self._refresh(report, [[dict(PROC_INFO)], [dict(PROC_INFO)]], [future])
        self.assert_unknown(report, '未来')

    def test_explicit_package_never_refreshes(self):
        report = dict(refresh_report(), evidence_scope='explicit')
        reads, finds = self._refresh(report, [], [])
        self.assertEqual((reads, finds), ([], 0))
        self.assertEqual(report['live_status_refresh']['status'], 'skipped')
        self.assertEqual(report['latest']['_time'], 'old')

    def test_real_log_timestamp_and_runtime_threshold_are_read(self):
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / 'safety.log'
            when = dt.datetime.fromtimestamp(NOW_NS / 1e9 - .25).astimezone()
            text = ('LOCALIZATION_DTC_CONFIG output_timeout_sec=0.7\n' +
                    when.strftime('I%m%d %H:%M:%S.%f') +
                    ' 42] SAFETY_DIAGNOSTIC sequence=1 raw_runstate=0 published_runstate=0 localization_consistency=1\n')
            log.write_text(text, encoding='utf-8')
            original = log.read_bytes()
            proc = dict(PROC_INFO, log=str(log), run=directory, start_wall_ns=NOW_NS-60*10**9)
            report = dict(refresh_report(), process=proc, has_safety_evidence=False,
                          motion_evidence={'status': 'FAIL'})
            with mock.patch.object(self.dl, 'find_processes', return_value=[proc]), \
                    mock.patch.object(self.dl, 'classify', wraps=self.dl.classify) as classify:
                self.dl.refresh_live_status(report, now_ns=NOW_NS)
            self.assertEqual(report['live_status_refresh']['status'], 'refreshed')
            self.assertAlmostEqual(report['latest_age_sec'], .25, places=5)
            self.assertTrue(report['has_safety_evidence'])
            self.assertTrue(ds.is_fresh_live(report))
            self.assertEqual(classify.call_args.args[1]['output_timeout_sec'], .7)
            self.assertEqual(report['thresholds']['output_timeout_sec'], .7)
            self.assertEqual(report['motion_evidence']['status'], 'PASS')
            self.assertEqual(report['motion_evidence']['evidence']['_time'], report['latest']['_time'])
            self.assertEqual(report['analyzed_latest']['motion_evidence']['status'], 'FAIL')
            self.assertEqual(log.read_bytes(), original)

    def test_last_record_without_time_does_not_reuse_previous_healthy_record(self):
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / 'safety.log'
            when = dt.datetime.fromtimestamp(NOW_NS / 1e9 - .25).astimezone()
            log.write_text(when.strftime('I%m%d %H:%M:%S.%f') +
                           ' 42] SAFETY_DIAGNOSTIC raw_runstate=0 published_runstate=0\n' +
                           'bad_time SAFETY_DIAGNOSTIC raw_runstate=4 published_runstate=4\n', encoding='utf-8')
            proc = dict(PROC_INFO, log=str(log), run=directory, start_wall_ns=NOW_NS-60*10**9)
            report = dict(refresh_report(), process=proc)
            with mock.patch.object(self.dl, 'find_processes', return_value=[proc]):
                self.dl.refresh_live_status(report, now_ns=NOW_NS)
            self.assert_unknown(report, '缺少有效时间')


if __name__ == '__main__':
    unittest.main()
