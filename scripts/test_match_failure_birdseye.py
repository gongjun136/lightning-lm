#!/usr/bin/env python3
"""Evidence integrity, exact 3D scoring and optional visualization regressions."""
import json
import hashlib
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch

import match_failure_birdseye as b


def fixture(root, run=None, wall_ns=1000, pid=42):
    evidence = root / ('failure_' + str(wall_ns))
    evidence.mkdir(parents=True)
    identity = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    shifted = identity.copy()
    shifted[3] = 1
    # Third point shares XY with the map but is 2 m above it. Fourth point is
    # display-only and must never dilute the original validation denominator.
    points = [(0, 0, 0, 0, 0, 1), (1, 0, 0, 1, 1, 1),
              (0, 0, 2, 2, 2, 1), (3, 0, 0, 0, 3, 0)]
    (evidence / 'scan.bin').write_bytes(b''.join(struct.pack('<fffIII', *p) for p in points))
    (evidence / 'map_local.bin').write_bytes(struct.pack('<ffffff', 0, 0, 0, 1, 0, 0))
    summary = dict(schema=1, kind=b.KIND, status='complete', run_dir=str(run or root), wall_ns=wall_ns,
                   sensor_stamp_ns=wall_ns - 1, pid=pid, session='test',
                   build=dict(commit='test'), map=dict(local_points=2, bounds_min=[0, 0, 0], bounds_max=[1, 0, 0],
                   coverage='all_possible_neighbors_for_saved_scan_and_all_three_poses', identity_sha256='test'),
                   poses=dict(guess=shifted, ndt=identity, check=identity),
                   criterion=dict(evaluated=True, overlap=2 / 3, points=3, nearest_neighbor_distance_m=.5,
                                  min_overlap=.75, sampling_stride=1),
                   scan=dict(source_points=4, saved_points=4, validation_points=3, display_downsampled=False,
                             record_bytes=24, encoding='little_endian_xyz_f32_lidar_u32_source_index_u32_validation_u32'),
                   files=dict(scan='scan.bin', map='map_local.bin'),
                   files_sha256=dict(scan=hashlib.sha256((evidence / 'scan.bin').read_bytes()).hexdigest(),
                                     map=hashlib.sha256((evidence / 'map_local.bin').read_bytes()).hexdigest()))
    (evidence / 'summary.json').write_text(json.dumps(summary), encoding='utf-8')
    return evidence, summary


class BirdseyeTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.evidence, self.summary = fixture(self.root)

    def tearDown(self):
        self.tmp.cleanup()

    def save(self):
        (self.evidence / 'summary.json').write_text(json.dumps(self.summary), encoding='utf-8')

    def dependencies(self):
        try:
            return b._dependencies()
        except ImportError:
            self.skipTest('optional plotting dependencies unavailable')

    def test_3d_score_uses_original_validation_sample_not_xy_or_display(self):
        np, tree, _ = self.dependencies()
        scan, cloud, mask, poses, threshold = b.load_points(self.evidence, self.summary, np)
        result = b.compare(scan, cloud, mask, poses, threshold, np, tree)
        self.assertEqual(result['check']['validation_points'], 3)
        self.assertEqual(result['check']['hit_points'], 2)
        self.assertEqual(result['check']['overlap'], 2 / 3)
        self.assertEqual(result['guess']['hit_points'], 1)

    def test_incomplete_manifest_is_not_evidence(self):
        self.summary['status'] = 'capturing'
        self.save()
        self.assertFalse(b.is_evidence(self.evidence))
        self.assertEqual(b.render_evidence(self.evidence)['status'], 'unavailable')
        self.assertFalse((self.evidence / 'birdseye').exists())

    def test_path_escape_rejected_before_import_or_write(self):
        self.summary['files']['scan'] = '../scan.bin'
        self.save()
        with patch.object(b, '_dependencies', side_effect=AssertionError('must not import')):
            result = b.render_evidence(self.evidence)
        self.assertEqual(result['status'], 'unavailable')
        self.assertFalse((self.evidence / 'birdseye').exists())

    def test_invalid_count_and_partial_binary_are_rejected(self):
        np, _, _ = self.dependencies()
        self.summary['criterion']['points'] = 4
        with self.assertRaisesRegex(ValueError, 'sample'):
            b.load_points(self.evidence, self.summary, np)
        (self.evidence / 'map_local.bin').write_bytes(b'x')
        with self.assertRaisesRegex(ValueError, 'partial'):
            b.load_points(self.evidence, self.summary, np)

    def test_missing_dependency_is_explicit_nonfatal_status(self):
        with patch.object(b, '_dependencies', side_effect=ImportError('scipy')):
            result = b.render_evidence(self.evidence)
        self.assertEqual(result['status'], 'dependency_missing')
        self.assertIn('scipy', result['message'])

    def test_height_filter_does_not_change_score_and_cache_hashes_inputs(self):
        self.dependencies()
        result = b.render_evidence(self.evidence, z_max=.1)
        self.assertEqual(result['status'], 'generated', result)
        self.assertEqual(result['comparison']['check']['overlap'], 2 / 3)
        self.assertTrue(result['check_score_matches_recorded'])
        self.assertFalse(result['display']['height_filter_affects_score'])
        self.assertGreater(Path(result['image']).stat().st_size, 1000)
        self.assertIn('data:image/png;base64,', Path(result['html']).read_text())
        with patch.object(b, '_dependencies', side_effect=AssertionError('cache must not import')):
            self.assertEqual(b.render_evidence(self.evidence, z_max=.1)['status'], 'cached')
        self.summary['criterion']['overlap'] = .25
        self.save()
        rerendered = b.render_evidence(self.evidence, z_max=.1)
        self.assertEqual(rerendered['status'], 'generated')
        self.assertFalse(rerendered['check_score_matches_recorded'])
        self.assertNotEqual(result['fingerprint'], rerendered['fingerprint'])

    def test_zero_local_neighbours_is_valid_zero_overlap(self):
        np, tree, _ = self.dependencies()
        (self.evidence / 'map_local.bin').write_bytes(b'')
        self.summary['map']['local_points'] = 0
        scan, cloud, mask, poses, threshold = b.load_points(self.evidence, self.summary, np)
        result = b.compare(scan, cloud, mask, poses, threshold, np, tree)
        self.assertEqual(result['check']['overlap'], 0)

    def test_oversized_binary_and_nonfinite_json_are_rejected(self):
        with patch.object(b, 'MAX_MAP_BYTES', 12):
            self.assertEqual(b.render_evidence(self.evidence)['status'], 'unavailable')
        (self.evidence / 'summary.json').write_text('{"x":NaN}')
        self.assertFalse(b.is_evidence(self.evidence))

    def test_hash_mismatch_is_rejected_before_dependency_import(self):
        data = bytearray((self.evidence / 'scan.bin').read_bytes())
        data[0] ^= 1
        (self.evidence / 'scan.bin').write_bytes(data)
        with patch.object(b, '_dependencies', side_effect=AssertionError('bad hash')):
            result = b.render_evidence(self.evidence)
        self.assertEqual(result['status'], 'unavailable')
        self.assertIn('SHA256', result['message'])

    def test_unevaluated_guard_labels_exploratory_score(self):
        self.dependencies()
        content = bytearray((self.evidence / 'scan.bin').read_bytes())
        for offset in range(20, len(content), 24):
            content[offset:offset + 4] = struct.pack('<I', 0)
        (self.evidence / 'scan.bin').write_bytes(content)
        self.summary['files_sha256']['scan'] = hashlib.sha256(content).hexdigest()
        self.summary['criterion'].update(evaluated=False, points=0)
        self.summary['scan']['validation_points'] = 0
        self.save()
        result = b.render_evidence(self.evidence)
        self.assertEqual(result['status'], 'generated', result)
        self.assertIsNone(result['check_score_matches_recorded'])
        self.assertEqual(result['comparison']['check']['score_basis'], 'saved_display_sample_only')

    def test_evaluated_but_zero_samples_is_not_claimed_as_reproduction(self):
        self.dependencies()
        content = bytearray((self.evidence / 'scan.bin').read_bytes())
        for offset in range(20, len(content), 24):
            content[offset:offset + 4] = struct.pack('<I', 0)
        (self.evidence / 'scan.bin').write_bytes(content)
        self.summary['files_sha256']['scan'] = hashlib.sha256(content).hexdigest()
        self.summary['criterion'].update(evaluated=True, points=0, sample_replay_available=False)
        self.summary['scan']['validation_points'] = 0
        self.summary['synthetic_fixture'] = True
        self.save()
        result = b.render_evidence(self.evidence)
        self.assertEqual(result['status'], 'generated', result)
        self.assertIsNone(result['check_score_matches_recorded'])
        self.assertTrue(result['synthetic_fixture'])
        self.assertIn('SYNTHETIC TEST FIXTURE', Path(result['html']).read_text())

    def test_registration_rejection_does_not_claim_optimized_ndt_pose(self):
        self.dependencies()
        self.summary['reason'] = 'registration_rejected'
        self.summary['poses']['ndt'] = self.summary['poses']['guess'].copy()
        self.save()
        result = b.render_evidence(self.evidence)
        self.assertEqual(result['status'], 'generated', result)
        self.assertEqual(result['comparison']['ndt']['candidate_status'], 'unconfirmed')
        self.assertIn('may remain initial guess', result['comparison']['ndt']['label'])
        self.assertNotEqual(result['comparison']['ndt']['label'], 'NDT result')
        self.assertEqual(result['recorded_reason'], 'registration_rejected')

    def test_zero_sample_map_bounds_reason_does_not_prove_outside_map(self):
        self.dependencies()
        content = bytearray((self.evidence / 'scan.bin').read_bytes())
        for offset in range(20, len(content), 24):
            content[offset:offset + 4] = struct.pack('<I', 0)
        (self.evidence / 'scan.bin').write_bytes(content)
        self.summary['files_sha256']['scan'] = hashlib.sha256(content).hexdigest()
        self.summary['reason'] = 'map_bounds'
        self.summary['criterion'].update(evaluated=True, points=0, sample_replay_available=False)
        self.summary['scan']['validation_points'] = 0
        self.save()
        result = b.render_evidence(self.evidence)
        self.assertEqual(result['status'], 'generated', result)
        self.assertEqual(result['recorded_reason'], 'map_bounds')
        self.assertIn('rejection predicate not established', result['display_reason'])
        self.assertNotEqual(result['reason'], 'map_bounds')
        self.assertIn('Acceptance candidate', result['comparison']['check']['label'])
        self.assertIn('no usable overlap samples', result['comparison']['check']['label'])
        self.assertEqual(result['comparison']['check']['score_basis'], 'saved_display_sample_only')
        self.assertIsNone(result['check_score_matches_recorded'])


if __name__ == '__main__':
    unittest.main()
