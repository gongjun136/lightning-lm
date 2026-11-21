"""只用合成输入，验证逐路归因边界、历史证据与一帧恢复展示。"""
import copy
import unittest
import diagnose_localization as diagnosis
import diagnostic_summary as summary
from radar_input_diagnostics import describe_radar_inputs
from test_diagnostic_summary import live, episode


def sample():
    row = dict(policy_version='localization-dtc-v9-match-grace', raw_runstate=2,
               lidar_frame_started=1, lidar_frame_seen=1, lidar_frame_fault=1,
               lidar_frame_age_sec=.4, lidar_frame_recovery_count=0,
               radar_rx_schema=1, radar_rx_available=1, radar_rx_consistent=1,
               radar_rx_after_safety_sec=.001, radar_rx_ids='0|1|2', radar_rx_primary=0,
               radar_rx_safety_accepted=3, radar_rx_snapshot_accepted=3,
               radar_rx_timeout_sec=.3, radar_rx_min_lidars=1, radar_rx_tolerance_sec=.02,
               radar_rx_reorder_sec=.1, radar_rx_pending=1, radar_rx_pending_present='1|2',
               radar_rx_pending_missing=0, radar_rx_pending_progress_sec=0,
               _time='00:00:01', _source='example.log:1', _scope='live_run')
    for i in range(3):
        row.update({f'radar_rx_{i}_{k}': v for k, v in dict(
            topic=f'/lidar/{i}', samples=5, window_samples=2, age_sec=.4 if i == 0 else .01,
            lag_sec=.05, max_lag_sec=.06, source_gap_sec=.1, arrival_gap_sec=.1,
            stamp_state='advancing', assembly='accepted').items()})
    return row


class RadarInputDiagnosticsTest(unittest.TestCase):
    def test_primary_silence_is_shown_without_dds_root_cause(self):
        row = sample()
        findings = diagnosis.classify(row, dict(radar_frame_timeout_sec=.3, radar_frame_recovery_frames=1))
        text = summary.render_summary(live(current_findings=findings))
        self.assertIn('雷达 [0] 未到达或回调已超时', text)
        self.assertIn('主雷达是合格组帧的必需来源', text)
        self.assertIn('首个新鲜合格组帧清除码 2', text)
        self.assertIn('时钟同步未验证', text)
        self.assertIn('50.00 ms', text)
        self.assertIn('缺少 [0]', text)

    def test_delay_and_tolerance_are_distinct(self):
        row = sample()
        row['radar_rx_0_age_sec'] = .01
        row['radar_rx_2_arrival_gap_sec'] = .51
        row['radar_rx_2_lag_sec'] = .45
        row['radar_rx_2_assembly'] = 'late_bucket'
        text = '\n'.join(describe_radar_inputs(row))
        self.assertIn('雷达 [2] 有相邻消息延迟交付线索', text)
        self.assertIn('所属时间桶已排出', text)
        row['radar_rx_2_arrival_gap_sec'] = .1
        row['radar_rx_2_assembly'] = 'tolerance'
        text = '\n'.join(describe_radar_inputs(row))
        self.assertIn('同桶时间跨度超容差', text)
        self.assertNotIn('有相邻消息延迟交付线索', text)

    def test_all_silent_and_no_data_do_not_blame_one_lidar(self):
        row = sample()
        for i in range(3): row[f'radar_rx_{i}_samples'] = 0
        self.assertIn('不能锁定单一路', '\n'.join(describe_radar_inputs(row)))
        del row['radar_rx_0_samples']
        row['radar_rx_0_age_sec'] = float('nan')
        self.assertIn('雷达 [0] 接收证据不完整', '\n'.join(describe_radar_inputs(row)))

    def test_old_busy_unknown_and_concurrent_snapshots(self):
        self.assertIn('旧日志无法补算', '\n'.join(describe_radar_inputs({})))
        row = sample()
        row['radar_rx_available'] = 0
        self.assertIn('不沿用上一次状态', '\n'.join(describe_radar_inputs(row)))
        row['radar_rx_available'] = 1
        for key, value in [('radar_rx_consistent', 0), ('radar_rx_after_safety_sec', .06)]:
            modified = {**row, key: value}
            text = '\n'.join(describe_radar_inputs(modified))
            self.assertIn('不认定为触发瞬间原因', text)
            self.assertIn('逐路结论：不归因', text)
        row['radar_rx_schema'] = 99
        self.assertIn('版本不支持', '\n'.join(describe_radar_inputs(row)))
        row = {**sample(), 'radar_rx_snapshot_accepted':4}
        self.assertIn('逐路结论：不归因', '\n'.join(describe_radar_inputs(row)))

    def test_negative_lag_and_frozen_stamps_are_not_healthy(self):
        row = sample()
        row.update(radar_rx_1_lag_sec=-1, radar_rx_1_stamp_state='rollback')
        text = '\n'.join(describe_radar_inputs(row))
        self.assertIn('接收时间差为负', text)
        self.assertIn('源时间戳回退', text)
        self.assertIn('-1000.00 ms', text)

    def test_causal_gate_missing_primary_and_auxiliary_threshold(self):
        row = sample()
        row.update(radar_rx_0_age_sec=.01, radar_rx_2_arrival_gap_sec=.5)
        text = '\n'.join(describe_radar_inputs(row))
        self.assertIn('待组帧桶缺主雷达 0', text)
        self.assertIn('因果未闭合', text)
        self.assertIn('辅助雷达缺失单独不构成码 2', text)
        row.update(radar_rx_pending_present=0, radar_rx_pending_missing='1|2', radar_rx_min_lidars=3)
        self.assertIn('来源数 1/3，未满足数量门槛', '\n'.join(describe_radar_inputs(row)))

    def test_empty_window_has_no_fake_peak(self):
        row = sample()
        row['radar_rx_0_window_samples'] = 0
        text = '\n'.join(describe_radar_inputs(row))
        self.assertIn('不适用（窗口内无回调）', text)
        self.assertIn('此后无回调', text)

    def test_recovery_count_in_history_is_explicitly_onset(self):
        row = sample()
        f = episode(2, 'frame_supply', 1, row, onset_evidence=row, state='recovered',
                    end_time='00:00:02', end_source='sample.log:2')
        text = summary.render_summary(live(findings=[f], thresholds={'radar_frame_recovery_frames':1}))
        self.assertIn('触发时 恢复确认帧', text)
        self.assertIn('原始码回 0 首见 00:00:02', text)
        self.assertIn('回 0 依据：sample.log:2', text)

    def test_history_uses_onset_not_recovered_route_state(self):
        onset = sample()
        latest = {**onset, 'radar_rx_0_age_sec': .01}
        f = episode(2, 'frame_supply', 1, latest, onset_evidence=onset)
        report = live(findings=[f])
        before = copy.deepcopy(report)
        text = summary.render_summary(report)
        self.assertIn('以下为触发时输入快照', text)
        self.assertIn('雷达 [0] 未到达或回调已超时', text)
        self.assertEqual(report, before)

    def test_flat_log_round_trip_and_single_id(self):
        row = diagnosis.fields('radar_rx_schema=1 radar_rx_available=1 radar_rx_consistent=1 '
            'radar_rx_after_safety_sec=0 radar_rx_ids=0 radar_rx_primary=0 radar_rx_timeout_sec=0.3 '
            'radar_rx_0_samples=0 radar_rx_0_age_sec=nan radar_rx_0_lag_sec=nan')
        self.assertIn('雷达 0（主', '\n'.join(describe_radar_inputs(row)))
        row['radar_rx_ids'] = '0|bad'
        self.assertIn('列表缺失或无效', '\n'.join(describe_radar_inputs(row)))

    def test_timeline_link_requires_unique_matching_sample(self):
        log = {**sample(), '_wall_ns':1_020_000_000}
        row = {k:v for k,v in log.items() if not k.startswith('radar_rx_')}
        row.update(_wall_ns=1_000_000_000, _source='timeline.csv:2')
        before = copy.deepcopy(row)
        linked = diagnosis.link_radar_input_evidence([row], [log])[0]
        self.assertEqual(linked['radar_rx_link_source'], 'example.log:1')
        self.assertEqual(linked['_source'], 'timeline.csv:2')
        self.assertEqual(row, before)
        ambiguous = diagnosis.link_radar_input_evidence([row], [log, dict(log)])[0]
        self.assertIn('候选 2 条', '\n'.join(describe_radar_inputs(ambiguous)))
        self.assertIn('非事件时间线原记录', '\n'.join(describe_radar_inputs(linked)))
        self.assertIn('独立的取证上限', '\n'.join(describe_radar_inputs(linked)))
        for logs in ([{**log, '_wall_ns':1_100_000_000}],
                     [{**log, 'lidar_frame_age_sec':.41}],
                     [{**log, 'lidar_frame_recovery_count':1}], [log, dict(log)]):
            self.assertNotIn('radar_rx_schema', diagnosis.link_radar_input_evidence([row], logs)[0])
        rounded = {**log, 'lidar_frame_age_sec':.31925}
        self.assertIn('radar_rx_schema', diagnosis.link_radar_input_evidence(
            [{**row, 'lidar_frame_age_sec':.319250073}], [rounded])[0])


if __name__ == '__main__':
    unittest.main()
