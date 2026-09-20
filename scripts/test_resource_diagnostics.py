#!/usr/bin/env python3
"""Unit and live Linux checks for diagnostic-only resource accounting."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import yaml

import monitor_process_resources as monitor
import summarize_run_config as config_summary


class ResourceTests(unittest.TestCase):
    def test_stat_with_spaces_and_parentheses(self):
        fields=['0']*50
        fields[0]='R';fields[2]='321';fields[11]='200';fields[12]='100'
        fields[17]='6';fields[19]='999';fields[36]='3'
        row=monitor.parse_stat('42 (worker ) extra) '+' '.join(fields))
        self.assertEqual(row['name'],'worker ) extra')
        self.assertEqual((row['cpu_ticks'],row['start_ticks'],row['last_cpu']),(300,999,3))
        self.assertEqual(row['pgrp'],321)

    def test_cpu_equivalents_and_reused_tid(self):
        old=dict(pid=1,start_ticks=10,name='worker',state='R',last_cpu=2,
                 allowed_cpus='0-3',cpu_ticks=100,sched_wait_ns=500)
        new=dict(old,cpu_ticks=250,sched_wait_ns=1500)
        result=monitor.thread_interval(new,old,2,100)
        self.assertAlmostEqual(result['cpu_core_equivalents'],.75)
        self.assertEqual(result['sched_wait_ns_delta'],1000)
        result=monitor.thread_interval(dict(new,start_ticks=20),old,2,100)
        self.assertIsNone(result['cpu_core_equivalents'])
        self.assertIsNone(result['sched_wait_ns_delta'])

    def test_missing_and_reset_counters_are_null(self):
        self.assertIsNone(monitor.counter_delta({'a':2},{'a':3},'a'))
        self.assertIsNone(monitor.counter_delta({'a':2},{},'a'))
        self.assertEqual(monitor.cpu_list_count('0-3,6,8-9'),7)

    def test_host_cpu_and_iowait(self):
        row=monitor.host_intervals({'cpu0':(200,60,10)},{'cpu0':(100,20,0)})[0]
        self.assertEqual(row['busy_percent'],50)
        self.assertEqual(row['iowait_percent'],10)
        self.assertIsNone(monitor.host_intervals({'cpu0':(200,60,10)},{})[0]['busy_percent'])

    def test_formal_config_matches_conservative(self):
        repo=Path(__file__).resolve().parent.parent
        path=repo/'config/reproduction/multi_lidar/sany_3livox/sany_3lidar_localization_solid.yaml'
        declared=config_summary.summarize(path,{'LIGHTNING_LM_NDT_MAX_POINTS':'3600'})
        self.assertTrue(declared['lio_point_budget_enabled'])
        load=declared['adaptive_load']
        self.assertEqual(load['lio_point_budgets'],[2200,1804,1500])
        self.assertEqual(load['tracking_lidar_count'],3)
        self.assertFalse(load['rotate_secondary_lidars'])
        self.assertEqual(load['hard_latency_sec'],.3)
        self.assertEqual(declared['requested_effective_budget']['ndt_max_points'],3600)
        self.assertEqual(declared['compute_budget']['ndt_max_points'],3500)
        with tempfile.TemporaryDirectory() as temp:
            generated=Path(temp)/'candidate.yaml'
            subprocess.run([sys.executable,str(repo/'scripts/prepare_compute_pruning_config.py'),
                            '--base',str(path),'--output',str(generated)],check=True,capture_output=True)
            result=config_summary.summarize(generated,{})
            self.assertEqual(result['adaptive_load'],declared['adaptive_load'])

    def test_three_lidar_uses_current_map_without_fixed_transform(self):
        repo=Path(__file__).resolve().parent.parent
        path=repo/'config/reproduction/multi_lidar/sany_3livox/sany_3lidar_localization_solid.yaml'
        config=yaml.safe_load(path.read_text(encoding='utf-8'))
        self.assertIs(config['output']['fixed_map_transform']['enabled'],False)

    def test_system_cpu_ranking_handles_pid_reuse(self):
        old={1:dict(start_ticks=10,cpu_ticks=100),2:dict(start_ticks=20,cpu_ticks=80),9:{}}
        current={1:dict(pid=1,start_ticks=10,cpu_ticks=250,name='worker',pgrp=1,state='R',threads=3),
                 2:dict(pid=2,start_ticks=21,cpu_ticks=900,name='reused',pgrp=2,state='S',threads=1)}
        result=monitor.system_process_intervals(current,old,2,100)
        self.assertEqual(len(result['top']),1)
        self.assertAlmostEqual(result['top'][0]['cpu_core_equivalents'],.75)
        self.assertEqual(result['unpaired_processes'],1)
        self.assertEqual(result['vanished_processes'],1)

    @unittest.skipUnless(sys.platform=='linux','requires real Linux /proc')
    def test_find_executable_below_launcher(self):
        # ros2 run is a Python parent; the sampler must identify its executable child.
        wrapper=subprocess.Popen([sys.executable,'-c',
            'import subprocess\np=subprocess.Popen(["sleep","2"])\nprint(p.pid,flush=True)\np.wait()'],
            start_new_session=True,stdout=subprocess.PIPE,text=True)
        try:
            child_pid=int(wrapper.stdout.readline())
            target=monitor.discover_target(Path('/proc'),wrapper.pid,'sleep')
            self.assertIsNotNone(target)
            self.assertEqual(target['pid'],child_pid)
            self.assertNotEqual(target['pid'],wrapper.pid)
            wrapper.wait(timeout=5)
        finally:
            wrapper.stdout.close()
            if wrapper.poll() is None:
                os.killpg(wrapper.pid,15)
                wrapper.wait(timeout=5)

    @unittest.skipUnless(sys.platform=='linux','requires real Linux /proc')
    def test_live_process_group_and_interval(self):
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp)/'resources.jsonl'
            child=subprocess.Popen([sys.executable,'-c',
                'import time\nend=time.monotonic()+4\nwhile time.monotonic()<end: sum(range(5000))'],start_new_session=True)
            sampler=None
            try:
                exe=Path(os.readlink(f'/proc/{child.pid}/exe')).name
                self.assertIsNone(monitor.discover_target(Path('/proc'),child.pid+1000000,exe))
                sampler=subprocess.Popen([sys.executable,monitor.__file__,'--process-group',str(child.pid),
                    '--executable',exe,'--interval-sec','.5','--output',str(out)],stderr=subprocess.PIPE,text=True)
                child.wait(timeout=10)
                _,error=sampler.communicate(timeout=10)
                self.assertEqual(sampler.returncode,0,error)
                records=[json.loads(line) for line in out.read_text().splitlines()]
                self.assertEqual(records[0]['pid'],child.pid)
                self.assertEqual(records[0]['schema_version'],2)
                self.assertEqual(len(records[0]['executable_sha256']),64)
                samples=[r for r in records if r['type']=='sample']
                self.assertGreaterEqual(len(samples),3)
                self.assertIsNone(samples[0]['cpu_core_equivalents'])
                self.assertGreater(max(r['cpu_core_equivalents'] or 0 for r in samples),.1)
                self.assertTrue(all(r['pid']==child.pid for r in samples))
                self.assertTrue(all(r['memory']['VmRSS_kb'] is not None for r in samples))
                self.assertTrue(all('wchan' in thread for r in samples for thread in r['threads']))
                self.assertTrue(all(len(r['system_process_cpu']['top'])<=20 for r in samples))
                self.assertTrue(any(row['pid']==child.pid and row['cpu_core_equivalents']>.1
                                    for r in samples for row in r['system_process_cpu']['top']))
                self.assertEqual(records[-1]['type'],'end')
            finally:
                if child.poll() is None:child.terminate();child.wait(timeout=5)
                if sampler and sampler.poll() is None:sampler.terminate();sampler.wait(timeout=5)


if __name__=='__main__':unittest.main()
