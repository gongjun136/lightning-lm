#!/usr/bin/env python3
"""Standalone Linux compile/smoke for the bounded asynchronous diagnostic sink."""
import csv,os,re,subprocess,tempfile,unittest
from pathlib import Path

@unittest.skipUnless(os.name=='posix','requires Linux compiler and syscalls')
class CausalTraceTests(unittest.TestCase):
    def test_multithread_footer_and_no_overwrite(self):
        repo=Path(__file__).resolve().parent.parent
        with tempfile.TemporaryDirectory() as temp:
            work=Path(temp);binary=work/'trace_test'
            subprocess.run(['g++','-std=c++17','-pthread','-I'+str(repo/'src'),
                str(repo/'src/utils/causal_trace.cc'),str(repo/'src/test/causal_trace_driver.cc'),
                '-o',str(binary)],check=True)
            path=work/'trace.csv';env=dict(os.environ,LIGHTNING_LM_CAUSAL_TRACE_PATH=str(path))
            subprocess.run([str(binary),'on'],env=env,check=True)
            text=path.read_text();rows=list(csv.DictReader(v for v in text.splitlines() if not v.startswith('#')))
            footer=re.search(r'# received=(\d+) written=(\d+) dropped=(\d+) io_errors=(\d+)',text)
            self.assertIsNotNone(footer)
            received,written,dropped,errors=map(int,footer.groups())
            self.assertEqual(received,4000);self.assertEqual(received,written+dropped)
            self.assertEqual(errors,0);self.assertEqual(len(rows),written);self.assertGreater(written,0)
            self.assertEqual(len({r['sequence'] for r in rows}),written)
            self.assertTrue(all(r['source']=='test' and r['after_v']=='2' for r in rows))
            subprocess.run([str(binary),'off'],env=env,check=True)  # Existing file disables trace.
            self.assertEqual(path.read_text(),text)
            env.pop('LIGHTNING_LM_CAUSAL_TRACE_PATH')
            subprocess.run([str(binary),'off'],env=env,check=True)

if __name__=='__main__':unittest.main()
