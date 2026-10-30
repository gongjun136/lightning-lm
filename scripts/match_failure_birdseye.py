#!/usr/bin/env python3
"""Render saved match-failure evidence; no ROS, playback, or live subscriptions."""
import argparse
import base64
import hashlib
import html
import json
import math
import os
from pathlib import Path
import sys
import tempfile

SCHEMA = 1
KIND = 'lightning_match_failure'
RENDER_VERSION = 2
MAX_SUMMARY_BYTES = 1024 * 1024
MAX_SCAN_POINTS = 52000
MAX_MAP_BYTES = 64 * 1024 * 1024
MAX_CANDIDATES = 256
POSES = ('guess', 'ndt', 'check')
NOTE = 'Red means no map neighbour within the 3D distance threshold; it does not identify a dynamic object.'


def read_summary(directory):
    """The writer atomically publishes a complete directory, never a partial JSON."""
    directory = Path(directory)
    path = directory / 'summary.json'
    if directory.name.startswith('.') or directory.is_symlink() or path.is_symlink():
        raise ValueError('evidence directory/summary must not be a symbolic link')
    if path.stat().st_size > MAX_SUMMARY_BYTES:
        raise ValueError('summary exceeds 1 MiB')
    def invalid(value):
        raise ValueError('nonfinite JSON value: ' + value)
    value = json.loads(path.read_text(encoding='utf-8'), parse_constant=invalid)
    if not isinstance(value, dict) or value.get('kind') != KIND or value.get('schema') != SCHEMA:
        raise ValueError('unsupported match-failure evidence schema')
    if value.get('status') != 'complete':
        raise ValueError('match-failure evidence is not complete')
    return value


def is_evidence(directory):
    try:
        return read_summary(directory).get('kind') == KIND
    except (OSError, ValueError, TypeError):
        return False


def _finite(value, name, positive=False):
    if isinstance(value, bool):
        raise ValueError(name + ' must be numeric')
    number = float(value)
    if not math.isfinite(number) or (positive and number <= 0):
        raise ValueError(name + ' is invalid')
    return number


def _file(directory, name, expected, limit):
    if name != expected:
        raise ValueError('unsupported evidence filename: ' + str(name))
    path = directory / expected
    if path.is_symlink() or path.resolve().parent != directory.resolve():
        raise ValueError('evidence file escapes its directory')
    size = path.stat().st_size
    if size < 0 or (size == 0 and expected != 'map_local.bin') or size > limit:
        raise ValueError(expected + ' is empty or exceeds the visualization limit')
    return path, size


def _fingerprint(paths, options):
    digest = hashlib.sha256(json.dumps(options, sort_keys=True).encode())
    hashes = {}
    for path in paths:
        digest.update(path.name.encode())
        own_digest = hashlib.sha256()
        with path.open('rb') as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(chunk)
                own_digest.update(chunk)
        hashes[path.name] = own_digest.hexdigest()
    return digest.hexdigest(), hashes


def _dependencies():
    # Set defaults before numpy imports: the checker must not create a large
    # BLAS pool or use the localization process's isolated CPUs.
    for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[name] = '1'
    import numpy as np
    from scipy.spatial import cKDTree
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return np, cKDTree, plt


def load_points(directory, summary, np):
    scan_path, scan_size = _file(directory, summary.get('files', {}).get('scan'), 'scan.bin', MAX_SCAN_POINTS * 24)
    map_path, map_size = _file(directory, summary.get('files', {}).get('map'), 'map_local.bin', MAX_MAP_BYTES)
    if scan_size % 24 or map_size % 12:
        raise ValueError('point file has a partial record')
    scan_meta = summary.get('scan', {})
    if scan_meta.get('record_bytes') != 24 or scan_meta.get('encoding') != 'little_endian_xyz_f32_lidar_u32_source_index_u32_validation_u32':
        raise ValueError('unsupported scan encoding')
    dtype = np.dtype([('xyz', '<f4', (3,)), ('lidar', '<u4'), ('source_index', '<u4'), ('validation', '<u4')])
    scan = np.fromfile(scan_path, dtype=dtype)
    cloud = np.fromfile(map_path, dtype='<f4').reshape(-1, 3)
    if not np.isfinite(scan['xyz']).all() or not np.isfinite(cloud).all():
        raise ValueError('point evidence contains nonfinite coordinates')
    validation = scan['validation'] == 1
    if np.any(scan['validation'] > 1):
        raise ValueError('invalid validation flag')
    if len(scan) != scan_meta.get('saved_points') or int(validation.sum()) != scan_meta.get('validation_points'):
        raise ValueError('scan count does not match summary')
    if len(cloud) != summary.get('map', {}).get('local_points'):
        raise ValueError('map count does not match summary')
    criterion = summary.get('criterion', {})
    if int(validation.sum()) != criterion.get('points'):
        raise ValueError('original overlap validation sample is missing/incomplete')
    if not criterion.get('evaluated') and validation.any():
        raise ValueError('unevaluated criterion unexpectedly has validation samples')
    if summary.get('map', {}).get('coverage') != 'all_possible_neighbors_for_saved_scan_and_all_three_poses':
        raise ValueError('local map does not guarantee all 3D neighbours for the three saved poses')
    distance = _finite(criterion.get('nearest_neighbor_distance_m'), 'distance threshold', positive=True)
    if distance > 10:
        raise ValueError('distance threshold exceeds visualization limit')
    poses = {}
    for name in POSES:
        matrix = np.asarray(summary.get('poses', {}).get(name), dtype=np.float64)
        if matrix.size != 16 or not np.isfinite(matrix).all():
            raise ValueError('missing/invalid pose: ' + name)
        matrix = matrix.reshape(4, 4)
        if not np.allclose(matrix[3], [0, 0, 0, 1], atol=1e-7) or not np.allclose(matrix[:3, :3].T @ matrix[:3, :3], np.eye(3), atol=1e-4) or abs(np.linalg.det(matrix[:3, :3]) - 1) > 1e-4:
            raise ValueError('pose is not a rigid transform: ' + name)
        poses[name] = matrix
    return scan, cloud, validation, poses, distance


def compare(scan, cloud, validation, poses, distance, np, tree_type, bounds=None):
    tree = tree_type(cloud)
    result = {}
    for name, matrix in poses.items():
        transformed = scan['xyz'].astype(np.float64) @ matrix[:3, :3].T + matrix[:3, 3]
        # The live PCL guard transforms in double, then casts the query to
        # float32 and compares a float32 squared distance with the double limit.
        queries = transformed.astype(np.float32)
        _, indices = tree.query(queries, k=1, workers=1)
        hit = np.zeros(len(scan), dtype=bool)
        present = indices < len(cloud)
        deltas = queries[present] - cloud[indices[present]]
        hit[present] = np.sum(deltas * deltas, axis=1, dtype=np.float32) <= distance * distance
        if bounds is not None:
            low, high = bounds
            hit &= np.all((transformed >= low - distance) & (transformed <= high + distance), axis=1)
        score_mask = validation if validation.any() else np.ones(len(scan), dtype=bool)
        count = int(score_mask.sum())
        result[name] = dict(xyz=transformed, hit=hit, validation_points=count,
                            hit_points=int(hit[score_mask].sum()), overlap=float(hit[score_mask].mean()),
                            score_basis='original_validation_sample' if validation.any() else 'saved_display_sample_only',
                            display_points=len(scan))
    return result


def _output_directory(evidence, output):
    target = Path(output).expanduser() if output else evidence / 'birdseye'
    # No output paths are taken from the untrusted evidence manifest.
    if target.is_symlink():
        raise ValueError('output directory must not be a symbolic link')
    if not output and target.resolve().parent != evidence.resolve():
        raise ValueError('output directory escapes evidence')
    target.mkdir(parents=True, exist_ok=True)
    for name in ('comparison.png', 'index.html', 'analysis.json'):
        if (target / name).is_symlink():
            raise ValueError('output file must not be a symbolic link')
    return target.resolve()


def _atomic(path, data):
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name + '.', dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as stream:
            stream.write(data)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def render_evidence(directory, output_dir=None, z_min=None, z_max=None):
    """Return a machine-readable status even when optional plotting is unavailable."""
    directory = Path(directory).expanduser()
    try:
        summary = read_summary(directory)
        directory = directory.resolve()
        if z_min is not None:
            z_min = _finite(z_min, 'z_min')
        if z_max is not None:
            z_max = _finite(z_max, 'z_max')
        if z_min is not None and z_max is not None and z_min >= z_max:
            raise ValueError('z_min must be below z_max')
        scan_path, _ = _file(directory, summary.get('files', {}).get('scan'), 'scan.bin', MAX_SCAN_POINTS * 24)
        map_path, _ = _file(directory, summary.get('files', {}).get('map'), 'map_local.bin', MAX_MAP_BYTES)
        fingerprint, hashes = _fingerprint([directory / 'summary.json', scan_path, map_path],
                                           dict(renderer=RENDER_VERSION, z_min=z_min, z_max=z_max))
        for key, name in (('scan', 'scan.bin'), ('map', 'map_local.bin')):
            if hashes[name] != summary.get('files_sha256', {}).get(key):
                raise ValueError('evidence SHA256 mismatch or missing: ' + name)
        output = _output_directory(directory, output_dir)
        cached_path = output / 'analysis.json'
        if cached_path.is_file() and cached_path.stat().st_size <= MAX_SUMMARY_BYTES:
            try:
                cached = json.loads(cached_path.read_text(encoding='utf-8'))
                if cached.get('fingerprint') == fingerprint and all((output / name).is_file() for name in ('comparison.png', 'index.html')):
                    cached.update(status='cached', cached=True, evidence=str(directory), image=str(output / 'comparison.png'), html=str(output / 'index.html'), analysis=str(cached_path))
                    cached['provenance'].update(summary=str(directory / 'summary.json'), scan=str(scan_path), map=str(map_path))
                    return cached
            except (OSError, ValueError):
                pass
        try:
            np, tree_type, plt = _dependencies()
        except ImportError as exc:
            return dict(status='dependency_missing', evidence=str(directory), message='缺少离线绘图依赖 numpy/scipy/matplotlib；诊断继续，未生成图。' + str(exc))
        scan, cloud, validation, poses, distance = load_points(directory, summary, np)
        bounds = (np.asarray(summary['map'].get('bounds_min'), dtype=np.float64),
                  np.asarray(summary['map'].get('bounds_max'), dtype=np.float64))
        if any(value.shape != (3,) or not np.isfinite(value).all() for value in bounds) or np.any(bounds[0] > bounds[1]):
            raise ValueError('invalid original map bounds')
        comparison = compare(scan, cloud, validation, poses, distance, np, tree_type, bounds)
        evaluated = summary['criterion'].get('evaluated') is True
        replay = evaluated and bool(validation.any())
        original = _finite(summary['criterion'].get('overlap'), 'recorded overlap') if evaluated else None
        exact = (abs(comparison['check']['overlap'] - original) <= max(1e-6, 0.01 / int(validation.sum()))) if replay else None
        result = dict(schema=1, kind='match_failure_birdseye', status='generated', cached=False,
                      fingerprint=fingerprint, evidence=str(directory), wall_ns=summary.get('wall_ns'),
                      sensor_stamp_ns=summary.get('sensor_stamp_ns'), pid=summary.get('pid'), run_dir=summary.get('run_dir'),
                      reason=summary.get('reason'), synthetic_fixture=summary.get('synthetic_fixture') is True,
                      original_sample_replay_available=replay, poses=summary['poses'],
                      provenance=dict(summary=str(directory / 'summary.json'), scan=str(scan_path), map=str(map_path),
                                      map_identity_sha256=summary.get('map', {}).get('identity_sha256'), build=summary.get('build'),
                                      renderer_version=RENDER_VERSION, distance_dimension=3, map_coverage=summary['map']['coverage'],
                                      verified_file_sha256=hashes, score_arithmetic='double_transform_float32_query_and_squared_distance'),
                      criterion=summary['criterion'], original_check_overlap=original,
                      check_score_matches_recorded=exact, comparison={name: {k: v for k, v in values.items() if k not in ('xyz', 'hit')}
                                                                  for name, values in comparison.items()},
                      display=dict(z_min=z_min, z_max=z_max, height_filter_affects_score=False,
                                   source_scan_points=summary['scan'].get('source_points'), saved_scan_points=len(scan),
                                   validation_points=int(validation.sum()), display_downsampled=summary['scan'].get('display_downsampled')),
                      warning=NOTE, image=str(output / 'comparison.png'), html=str(output / 'index.html'), analysis=str(cached_path))
        def visible(points):
            mask = np.ones(len(points), dtype=bool)
            if z_min is not None:
                mask &= points[:, 2] >= z_min
            if z_max is not None:
                mask &= points[:, 2] <= z_max
            return mask
        shown_map = cloud[visible(cloud)]
        # Only the map display is thinned; the NN test always uses every saved
        # map point and every original validation sample.
        shown_map = shown_map[::max(1, math.ceil(len(shown_map) / 100000))]
        all_xy = np.concatenate([v['xyz'][:, :2] for v in comparison.values()])
        lower, upper = all_xy.min(axis=0) - 2, all_xy.max(axis=0) + 2
        fig, axes = plt.subplots(1, 3, figsize=(18, 7.3))
        # Reserve explicit header/legend space across matplotlib 3.5+; equal
        # axes can otherwise push their titles into a multiline suptitle.
        fig.subplots_adjust(left=.04, right=.99, top=.79, bottom=.17, wspace=.20)
        titles = dict(guess='Initial guess', ndt='NDT result', check='Actual overlap-check pose')
        try:
            for axis, name in zip(axes, POSES):
                value = comparison[name]
                shown = visible(value['xyz'])
                axis.scatter(shown_map[:, 0], shown_map[:, 1], s=.35, c='#aeb4bc', label='Map', rasterized=True)
                for match, color, label in [(True, '#168559', 'Map neighbour <= threshold'), (False, '#d63d42', 'No map neighbour')]:
                    points = value['xyz'][shown & (value['hit'] == match)]
                    axis.scatter(points[:, 0], points[:, 1], s=2, c=color, label=label, rasterized=True)
                origin = poses[name][:3, 3]
                axis.scatter([origin[0]], [origin[1]], marker='^', s=35, color='#203852', label='Pose origin')
                axis.set(xlim=(lower[0], upper[0]), ylim=(lower[1], upper[1]), xlabel='Map X [m]', ylabel='Map Y [m]',
                         title=f'{titles[name]}\n3D overlap {value["overlap"]:.4%} ({value["hit_points"]}/{value["validation_points"]})')
                axis.set_aspect('equal', adjustable='box')
                axis.grid(alpha=.15)
            handles, labels = axes[0].get_legend_handles_labels()
            fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(.5, .025), ncol=4, fontsize=9)
            score_note = (f'recorded check {original:.4%}, reproduced' if exact else
                          f'WARNING: recorded check {original:.4%}, recomputation differs' if replay else
                          'No original validation sample: scores use saved display sample only')
            fixture_note = 'SYNTHETIC TEST FIXTURE | ' if summary.get('synthetic_fixture') else ''
            fig.suptitle(f'{fixture_note}Match-failure comparison | sensor stamp {summary.get("sensor_stamp_ns")}\n'
                         f'3D NN threshold {distance:g} m | {score_note}\n{NOTE}', fontsize=11, y=.98)
            import io
            buffer = io.BytesIO()
            fig.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
        finally:
            plt.close(fig)
        png = buffer.getvalue()
        payload = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)
        document = ('<!doctype html><meta charset="utf-8"><title>Match failure comparison</title>'
                    '<style>body{font:16px system-ui;margin:24px;color:#203040}img{width:100%;height:auto}pre{white-space:pre-wrap}</style>'
                    '<h1>' + fixture_note + 'Match failure: three saved poses</h1><p>' + html.escape(NOTE) + '</p><p>' + html.escape(score_note) + '</p>'
                    '<p>Height filtering affects display only. Reported scores use all original validation points and 3D distances.</p>'
                    '<img alt="Map and scan comparison at initial, NDT and actual check poses" src="data:image/png;base64,' + base64.b64encode(png).decode() + '">'
                    '<details><summary>Provenance and scores</summary><pre>' + html.escape(payload) + '</pre></details>')
        _atomic(output / 'comparison.png', png)
        _atomic(output / 'index.html', document.encode('utf-8'))
        _atomic(cached_path, payload.encode('utf-8'))
        return result
    except (OSError, ValueError, TypeError, KeyError, RuntimeError) as exc:
        return dict(status='unavailable', evidence=str(directory), message='失败匹配证据无法绘图: ' + str(exc))


def main(argv=None):
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('evidence', help='complete match-failure directory')
    parser.add_argument('--output-dir')
    parser.add_argument('--z-min', type=float, help='map-frame display filter only')
    parser.add_argument('--z-max', type=float, help='map-frame display filter only')
    args = parser.parse_args(argv)
    # Standalone use observes the same isolated-core policy as --diagnose.
    from diagnose_localization import low_impact
    try:
        scheduling = low_impact()
    except RuntimeError as exc:
        print(json.dumps(dict(status='unavailable', message=str(exc)), ensure_ascii=False))
        return 2
    result = render_evidence(args.evidence, args.output_dir, args.z_min, args.z_max)
    result['scheduling'] = scheduling
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
    return 0 if result['status'] in ('generated', 'cached') else 2


if __name__ == '__main__':
    sys.exit(main())
