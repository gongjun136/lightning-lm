"""观察流与真实诊断入口的隔离回归。"""
import datetime as dt
import json
from pathlib import Path
import tempfile
import unittest
import diagnose_localization as diagnose
import diagnostic_summary as summary
from localization_log_events import any_event, has_event, CONTEXT_EVENTS


class ObservationIntegrationTest(unittest.TestCase):
    def test_raw_log_flags_and_large_counters_survive(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);log=root/'run.log'
            log.write_text('I1119 12:00:00.000000 1 core.cc:1] event=LOC_MAP_CONSISTENCY_CONFIG tracking_map_overlap_check_enabled=0 map_consistency_enabled=1 map_min_overlap=0.6\n'
                           'W1119 12:00:01.000000 1 core.cc:2] event=LOC_TRACKING_LOW_OVERLAP accepted=1 map_checked=1 map_overlap_check_enabled=0 map_overlap=0.59 map_min_overlap=0.6 observation_streak_frames=3 observation_streak_sec=0.2 observation_total_low_frames=9007199254740993 observation_capture_suppressed=2\n',encoding='utf-8')
            context=[];warnings=[]
            rows,config=diagnose.log_rows(log,1024*1024,warnings,dt.datetime(2026,11,19),context)
            self.assertEqual(rows,[])
            evidence=diagnose.overlap_observation_report(config['_map_overlap_policy'],context,root,None)
            self.assertIs(evidence['policy']['map_overlap_check_enabled'],False)
            self.assertTrue(evidence['latest']['evidence']['is_observation'])
            self.assertEqual(evidence['latest']['evidence']['total_low_frames'],9007199254740993)
            json.dumps(evidence,allow_nan=False)

    def test_configuration_alone_cannot_claim_running_gate(self):
        report={'schema':2,'selected':'run','latest':{},'warnings':[],
                'matching_thresholds':{'map_overlap_check_enabled':False,'map_min_overlap_ratio':.6}}
        rendered=summary.render_summary(report)
        self.assertIn('重合率检查：未知',rendered)
        self.assertNotIn('启动日志确认',rendered)
        report['overlap_observations']={'policy':{'map_overlap_check_enabled':False,'source':'startup_log','map_min_overlap_ratio':.6}}
        self.assertIn('启动日志确认',summary.render_summary(report))

    def test_new_events_have_no_accidental_legacy_match(self):
        self.assertFalse(has_event('unrelated legacy text','LOC_TRACKING_LOW_OVERLAP'))
        self.assertFalse(any_event('unrelated legacy text',CONTEXT_EVENTS))
        self.assertTrue(has_event('event=LOC_TRACKING_LOW_OVERLAP accepted=1','LOC_TRACKING_LOW_OVERLAP'))

    def test_collector_pid_and_kind_must_match(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);folder=root/'overlap_observations';folder.mkdir()
            path=folder/'collector_status.json'
            data={'schema':1,'kind':'tracking_overlap_observation_collector','pid':12,'saved':2,'suppressed':8}
            path.write_text(json.dumps(data))
            self.assertEqual(diagnose.overlap_observation_report({},[],root,{'pid':13})['collector']['status'],'unavailable')
            actual=diagnose.overlap_observation_report({},[],root,{'pid':12})
            self.assertEqual(actual['collector']['counters']['saved'],2)
            self.assertIsNone(actual['policy']['map_overlap_check_enabled'])
            data['kind']='lightning_match_failure_collector';path.write_text(json.dumps(data))
            self.assertEqual(diagnose.overlap_observation_report({},[],root,None)['collector']['status'],'unavailable')


if __name__=='__main__':unittest.main()
