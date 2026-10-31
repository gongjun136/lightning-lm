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

import compute_profile_report


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
                        str(repo/'src/utils/compute_profiling.cc'),
                        '-lglog', '-lgflags', '-o', str(cls.binary)], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def check_records(self, text):
        rows = [line for line in text.splitlines() if re.match(r'^[IWEF]\d{4} ', line)]
        self.assertTrue(rows, 'no glog records written')
        self.assertNotIn('async_seq=', text)
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

    def profile_run(self, root, mode='profile', **options):
        (root / 'results').mkdir(exist_ok=True)
        (root / 'run_metadata.txt').write_text('compute_profile=1\ncompute_profile_storage=dedicated-v1\n')
        env = dict(os.environ, LIGHTNING_LM_COMPUTE_PROFILE='1',
                   LIGHTNING_LM_COMPUTE_PROFILE_PATH=str(root / compute_profile_report.PROFILE))
        env.update(options.pop('env', {}))
        return subprocess.run([str(self.binary), mode], cwd=root, env=env,
                              capture_output=True, text=True, timeout=15, **options)

    def test_full_frame_profile_is_complete_and_not_duplicated_to_operational_log(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result = self.profile_run(root)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.check_records(result.stderr)
            self.assertIn('COMPUTE_BENCH_SUMMARY module=lidar_end_to_end', result.stderr)
            for tag in ('COMPUTE_BENCH_FRAME', 'LIO_BENCH_FRAME', 'COMPUTE_BENCH_EVENT', 'module=imu_dr_pipeline'):
                self.assertNotIn(tag, result.stderr)
            report = compute_profile_report.validate(root, True)
            self.assertEqual(report['status'], 'complete', report)
            self.assertEqual(report['records'], 200)
            self.assertEqual(set(report['modules'].values()), {25})
            text = (root / compute_profile_report.PROFILE).read_text()
            self.assertIn('frame_id=24 outer_wall_ms=1.25 outer_thread_cpu_ms=0.75 outer_process_cpu_ms_concurrent=2.5', text)
            self.assertNotIn('lidar_end_to_end', text)
            self.assertNotIn('TEST_WARNING', text)

    def test_existing_profile_fails_startup_without_overwriting_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'results').mkdir()
            path = root / compute_profile_report.PROFILE
            path.write_text('retained evidence\n')
            result = self.profile_run(root)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('cannot create compute profile', result.stderr)
            self.assertEqual(path.read_text(), 'retained evidence\n')

    def test_profile_truncation_and_abrupt_exit_are_not_reported_as_complete(self):
        for mode, problem in [('profile_long', 'truncated'), ('profile_crash', 'missing_or_invalid_footer')]:
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                result = self.profile_run(root, mode)
                self.assertEqual(result.returncode, 77 if mode == 'profile_crash' else 0, result.stderr)
                report = compute_profile_report.validate(root, True)
                self.assertEqual(report['status'], 'incomplete')
                self.assertIn(problem, report['problems'])

    def test_profile_disk_write_failure_is_visible_without_blocking_producers(self):
        import resource
        import signal
        def limit_file_size():
            signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
            resource.setrlimit(resource.RLIMIT_FSIZE, (512, 512))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result = self.profile_run(root, preexec_fn=limit_file_size)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('PRODUCERS_DONE', result.stdout)
            self.assertIn('COMPUTE_PROFILE_INCOMPLETE', result.stderr)
            self.assertEqual(compute_profile_report.validate(root, True)['status'], 'incomplete')

    def test_stalled_profile_writer_does_not_block_operational_logging(self):
        # Interpose only fflush on the new regular profile file. Its header is
        # flushed normally; stall the first writer batch without production hooks.
        shim_source = r'''
#define _GNU_SOURCE
#include <dlfcn.h>
#include <fcntl.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
int fflush(FILE *stream) {
    static int (*real_flush)(FILE *);
    static atomic_int matches;
    if (!real_flush) real_flush = dlsym(RTLD_NEXT, "fflush");
    const char *path = getenv("LIGHTNING_LM_COMPUTE_PROFILE_PATH");
    if (stream && path) {
        char fd_path[64], target[4096];
        snprintf(fd_path, sizeof(fd_path), "/proc/self/fd/%d", fileno(stream));
        ssize_t size = readlink(fd_path, target, sizeof(target)-1);
        if (size >= 0) {
            target[size] = 0;
            if (!strcmp(path, target) && atomic_fetch_add(&matches, 1) == 1) {
                int marker = open(getenv("PROFILE_PAUSE_MARKER"), O_CREAT|O_WRONLY, 0600);
                if (marker >= 0) close(marker);
                usleep(3000000);
            }
        }
    }
    return real_flush(stream);
}
'''
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, shim = root/'slow_profile.c', root/'slow_profile.so'
            source.write_text(shim_source)
            subprocess.run(['gcc', '-shared', '-fPIC', str(source), '-ldl', '-o', str(shim)], check=True)
            marker = root/'writer_paused'
            (root/'results').mkdir()
            env = dict(os.environ, LD_PRELOAD=str(shim), PROFILE_PAUSE_MARKER=str(marker),
                       LIGHTNING_LM_COMPUTE_PROFILE='1',
                       LIGHTNING_LM_COMPUTE_PROFILE_PATH=str(root / compute_profile_report.PROFILE))
            proc = subprocess.Popen([str(self.binary), 'profile_flood'], cwd=root, env=env,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            try:
                with selectors.DefaultSelector() as selector:
                    selector.register(proc.stdout, selectors.EVENT_READ)
                    self.assertTrue(selector.select(2), 'producers waited for the profile writer')
                self.assertIn(b'PRODUCERS_DONE', proc.stdout.readline())
                self.assertTrue(marker.is_file(), 'profile writer was not actually stalled')
                self.assertIsNone(proc.poll(), 'profile writer should still be paused')
                operational = b''
                deadline = time.monotonic() + 0.5
                with selectors.DefaultSelector() as selector:
                    selector.register(proc.stderr, selectors.EVENT_READ)
                    while b'TEST_ERROR' not in operational and time.monotonic() < deadline:
                        if selector.select(max(0, deadline-time.monotonic())):
                            operational += os.read(proc.stderr.fileno(), 4096)
                self.assertIn(b'TEST_ERROR', operational, 'operational logs waited for the profile writer')
                _, remaining = proc.communicate(timeout=15)
                text = (operational + remaining).decode()
                self.assertEqual(proc.returncode, 0, text)
            finally:
                if proc.poll() is None:
                    proc.kill()
                    proc.communicate()
            self.check_records(text)
            report = compute_profile_report.validate(root, True)
            self.assertGreater(report['footer']['dropped'], 0)
            self.assertIn('COMPUTE_PROFILE_INCOMPLETE', text)
            self.assertEqual(report['status'], 'incomplete')


if __name__ == '__main__':
    unittest.main()
