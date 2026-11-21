"""Exercise affinity readback against synthetic /proc snapshots; no live mutation."""
from pathlib import Path
import tempfile
import unittest

from verify_localization_affinity import verify


class AffinityTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.process = self.root / '123'
        (self.process / 'task').mkdir(parents=True)
        (self.process / 'environ').write_bytes(b'')
        self.isolated = self.root / 'isolated'
        self.isolated.write_text('10-11')
        for tid, name, mask in [(123, 'loc_executor', '11'), (124, 'loc_sensor', '10'),
                                (125, 'loc_output', '11'), (126, 'dds_receiver', '0-9')]:
            self.thread(tid, name, mask)

    def tearDown(self):
        self.tmp.cleanup()

    def thread(self, tid, name, mask):
        task = self.process / 'task' / str(tid)
        task.mkdir(exist_ok=True)
        (task / 'status').write_text(f'Name:\t{name}\nCpus_allowed_list:\t{mask}\n')

    def check(self):
        return verify(123, True, proc_root=self.root, isolated_path=self.isolated)

    def candidate(self):
        (self.process / 'environ').write_bytes(
            b'LIGHTNING_LM_MATCH_CPUS=11\0LIGHTNING_LM_NDT_THREADS=1\0')
        self.thread(127, 'loc_match', '11')

    def test_single_match_shares_output_core(self):
        self.candidate()
        self.assertTrue(self.check()['passed'])
        self.assertEqual(self.check()['match_cpus'], [11])

    def test_legacy_four_thread_background_remains_supported(self):
        for tid in range(127, 131):
            self.thread(tid, 'loc_match', '0-9')
        self.assertTrue(self.check()['passed'])

    def test_explicit_shared_core_comparison(self):
        (self.process / 'environ').write_bytes(
            b'LIGHTNING_LM_MATCH_CPUS=0-9\0LIGHTNING_LM_NDT_THREADS=4\0')
        for tid in range(127, 131):
            self.thread(tid, 'loc_match', '0-9')
        self.assertTrue(self.check()['passed'])

    def test_wrong_match_core_or_unexpected_worker_fails(self):
        self.candidate()
        self.thread(127, 'loc_match', '0-9')
        self.assertFalse(self.check()['passed'])
        self.thread(127, 'loc_match', '11')
        self.thread(128, 'loc_match', '11')
        self.assertFalse(self.check()['passed'])

    def test_background_cannot_leak_onto_output_core(self):
        self.candidate()
        self.thread(126, 'dds_receiver', '11')
        self.assertFalse(self.check()['passed'])

    def test_missing_configuration_or_match_is_not_success(self):
        self.assertFalse(self.check()['passed'])
        self.candidate()
        (self.process / 'environ').unlink()
        self.assertFalse(self.check()['passed'])

    def dual_lio(self, count=2):
        (self.process / 'environ').write_bytes(
            b'LIGHTNING_LM_MATCH_CPUS=0-9\0LIGHTNING_LM_NDT_THREADS=4\0'
            b'LIGHTNING_LM_LIO_WORKER_CPUS=11\0' + f'LIGHTNING_LM_LIO_THREADS={count}\0'.encode())
        for tid in range(127, 131):
            self.thread(tid, 'loc_match', '0-9')

    def test_dual_lio_requires_one_worker_on_cpu11(self):
        self.dual_lio()
        self.assertFalse(self.check()['passed'])
        self.thread(131, 'loc_lio_worker', '11')
        self.assertTrue(self.check()['passed'])
        self.thread(131, 'loc_lio_worker', '0-9')
        self.assertFalse(self.check()['passed'])
        self.thread(131, 'loc_lio_worker', '10-11')
        self.assertFalse(self.check()['passed'])

    def test_worker_count_and_ndt_overlap_are_checked(self):
        self.dual_lio(1)
        self.assertTrue(self.check()['passed'])
        self.thread(131, 'loc_lio_worker', '11')
        self.assertFalse(self.check()['passed'])
        self.dual_lio()
        self.assertTrue(self.check()['passed'])
        self.thread(132, 'loc_lio_worker', '11')
        self.assertFalse(self.check()['passed'])
        (self.process / 'environ').write_bytes(
            b'LIGHTNING_LM_MATCH_CPUS=11\0LIGHTNING_LM_LIO_WORKER_CPUS=11\0')
        self.assertFalse(self.check()['passed'])


if __name__ == '__main__':
    unittest.main()
