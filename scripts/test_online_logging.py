#!/usr/bin/env python3
"""Real glog integration; block its output fd, not a mock queue callback."""
import os
from pathlib import Path
import re
import selectors
import subprocess
import tempfile
import time
import unittest


@unittest.skipUnless(os.name == 'posix', 'requires Linux and glog')
class OnlineLoggingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.temp.name) / 'logging_test'
        repo = Path(__file__).resolve().parent.parent
        subprocess.run(['g++', '-O2', '-std=c++17', '-pthread', '-I'+str(repo/'src'),
                        str(repo/'src/test/online_logging_driver.cc'),
                        str(repo/'src/utils/online_logging.cc'),
                        str(repo/'src/utils/causal_trace.cc'),
                        '-lglog', '-lgflags', '-o', str(cls.binary)], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def check_records(self, text):
        rows = [line for line in text.splitlines() if ' async_seq=' in line]
        seq = [int(re.search(r'async_seq=(\d+)', line)[1]) for line in rows]
        self.assertEqual(seq, list(range(len(seq))))
        self.assertEqual(text.count('TEST_WARNING'), 1)
        self.assertEqual(text.count('TEST_ERROR'), 1)
        # No per-producer reversal, and no torn/interleaved records.
        for t in range(4):
            items = [int(m[1]) for line in rows
                     if (m := re.search(rf'producer={t} item=(\d+)', line))]
            self.assertEqual(items, sorted(set(items)))
        summary = re.search(r'ASYNC_LOG_SUMMARY received=(\d+) written=(\d+) dropped=(\d+) critical_dropped=(\d+) io_errors=(\d+) truncated=(\d+)', text)
        self.assertIsNotNone(summary)
        received, written, dropped, critical, errors, truncated = map(int, summary.groups())
        self.assertEqual(received, written + dropped)
        self.assertEqual(written, len(rows))
        self.assertEqual((critical, errors, truncated), (0, 0, 0))
        return dropped

    def test_normal_order_and_shutdown_drain(self):
        result = subprocess.run([str(self.binary)], capture_output=True, text=True, timeout=15)
        self.assertEqual(result.returncode, 0)
        self.check_records(result.stderr)

    def test_blocked_output_does_not_block_producers_or_flush(self):
        proc = subprocess.Popen([str(self.binary), 'flood'], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        try:
            # Do not read stderr: the actual writer blocks once the pipe fills.
            with selectors.DefaultSelector() as selector:
                selector.register(proc.stdout, selectors.EVENT_READ)
                self.assertTrue(selector.select(3), 'producers blocked on logging output')
            done = proc.stdout.readline().decode()
            self.assertIn('PRODUCERS_DONE', done)
            time.sleep(0.3)
            self.assertIsNone(proc.poll(), 'shutdown should still be draining the blocked writer')
            _, text = proc.communicate(timeout=15)
            self.assertEqual(proc.returncode, 0)
            dropped = self.check_records(text.decode())
            self.assertGreater(dropped, 0)
            print(done.strip(), 'dropped=', dropped)
        finally:
            if proc.poll() is None:
                proc.kill()
                proc.communicate()

    def test_synchronous_control_reproduces_blocking(self):
        proc = subprocess.Popen([str(self.binary), 'sync'], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        try:
            with selectors.DefaultSelector() as selector:
                selector.register(proc.stdout, selectors.EVENT_READ)
                self.assertFalse(selector.select(0.5), 'control did not fill stderr pipe')
            done, _ = proc.communicate(timeout=15)
            self.assertIn(b'PRODUCERS_DONE', done)
            self.assertEqual(proc.returncode, 0)
        finally:
            if proc.poll() is None:
                proc.kill()
                proc.communicate()

    def test_write_failure_does_not_retry_on_producer(self):
        with open('/dev/full', 'wb') as target:
            result = subprocess.run([str(self.binary), 'flood'], stdout=subprocess.PIPE,
                                    stderr=target, timeout=5)
        self.assertEqual(result.returncode, 0)
        self.assertIn(b'PRODUCERS_DONE', result.stdout)


if __name__ == '__main__':
    unittest.main()
