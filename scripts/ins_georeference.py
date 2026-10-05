#!/usr/bin/env python3
"""Rigid rear-axle trajectory alignment and immutable task/map asset migration.

Requires numpy, scipy, PyYAML. TUM: seconds x y z qx qy qz qw.
Never estimate scale or silently sort/deduplicate trajectories.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from scipy.spatial.transform import Rotation, Slerp
import yaml


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_tum(path):
    a = np.loadtxt(path, ndmin=2)
    if a.shape[1] != 8 or len(a) < 2 or not np.isfinite(a).all():
        raise ValueError('TUM must contain finite rows: t x y z qx qy qz qw')
    if (np.diff(a[:, 0]) <= 0).any():
        raise ValueError('TUM timestamps must strictly increase; repair a separate copy explicitly')
    if not np.allclose(np.linalg.norm(a[:, 4:], axis=1), 1, atol=1e-5):
        raise ValueError('TUM quaternions must be unit length')
    return a


def fit_rigid(source, target, minimum_span):
    x, y = source - source.mean(0), target - target.mean(0)
    # A straight path does not constrain roll around the path.
    if min(np.linalg.svd(x, compute_uv=False)[1], np.linalg.svd(y, compute_uv=False)[1]) / np.sqrt(len(x)) < minimum_span:
        raise ValueError('trajectory is stationary/near-collinear; collect a route with turns')
    u, _, vt = np.linalg.svd(x.T @ y)
    d = np.eye(3)
    d[2, 2] = np.linalg.det(vt.T @ u.T)
    r = vt.T @ d @ u.T
    return Rotation.from_matrix(r), target.mean(0) - r @ source.mean(0)


def residuals(source, target, rotation, translation, source_q, target_q):
    distance = np.linalg.norm(rotation.apply(source) + translation - target, axis=1)
    angles = (Rotation.from_quat(target_q).inv() * rotation * source_q).magnitude() * 180 / np.pi
    return dict(count=len(distance), position_rmse_m=float(np.sqrt(np.mean(distance**2))),
                position_p95_m=float(np.percentile(distance, 95)), position_max_m=float(max(distance)),
                attitude_rmse_deg=float(np.sqrt(np.mean(angles**2))), attitude_max_deg=float(max(angles)))


def fit(args):
    a, b = read_tum(args.source), read_tum(args.target)
    a[:, 0] += args.time_offset_sec
    idx = np.searchsorted(a[:, 0], b[:, 0], side='right')
    valid = (idx > 0) & (idx < len(a))
    idx = np.clip(idx, 1, len(a)-1)
    valid &= a[idx, 0] - a[idx-1, 0] <= args.max_gap_sec
    b, idx = b[valid], idx[valid]
    if len(b) < 30:
        raise ValueError('need at least 30 synchronized pairs with valid interpolation coverage')
    ratio = (b[:, 0] - a[idx-1, 0]) / (a[idx, 0] - a[idx-1, 0])
    x = (1-ratio[:, None])*a[idx-1, 1:4] + ratio[:, None]*a[idx, 1:4]
    q = Slerp(a[:, 0], Rotation.from_quat(a[:, 4:]))(b[:, 0])
    split = int(len(b)*0.7)
    r, t = fit_rigid(x[:split], b[:split, 1:4], args.min_span_m)
    train = residuals(x[:split], b[:split, 1:4], r, t, q[:split], b[:split, 4:])
    hold = residuals(x[split:], b[split:, 1:4], r, t, q[split:], b[split:, 4:])
    accepted = all(v['position_rmse_m'] <= args.max_rmse_m and v['position_max_m'] <= args.max_error_m
                   and v['attitude_max_deg'] <= args.max_angle_deg for v in [train, hold])
    with open(args.site) as f:
        geo = yaml.safe_load(f)['ins_only']['georeference']
    if not geo['confirmed'] or geo['datum'] != 'WGS84' or not np.isfinite(geo['origin_llh']).all():
        raise ValueError('use the confirmed site configuration that produced the target ENU trajectory')
    result = dict(version=1, accepted=bool(accepted), convention='target_from_source',
                  reference='rear_axle', target_frame='fixed_enu', georeference=geo,
                  target_from_source=dict(translation_m=t.tolist(), quaternion_xyzw=r.as_quat().tolist()),
                  source_sha256=digest(args.source), target_sha256=digest(args.target),
                  site_sha256=digest(args.site), time_offset_sec=args.time_offset_sec,
                  matching=dict(max_gap_sec=args.max_gap_sec, excluded_target_rows=int(len(read_tum(args.target))-len(b))),
                  validation=dict(split='first 70% fit, final 30% held out; no scale adjustment',
                                  training=train, held_out=hold,
                                  limits=dict(rmse_m=args.max_rmse_m, max_error_m=args.max_error_m, max_angle_deg=args.max_angle_deg)))
    with open(args.output, 'x') as f:
        yaml.safe_dump(result, f, sort_keys=False)
    print(json.dumps(result['validation'], indent=2))
    if not accepted:
        raise SystemExit('Alignment rejected; diagnostic transform saved with accepted: false')


def load_transform(path):
    with open(path) as f:
        c = yaml.safe_load(f)
    if c.get('accepted') is not True or c.get('convention') != 'target_from_source':
        raise ValueError('transform must be accepted and have target_from_source convention')
    q = np.asarray(c['target_from_source']['quaternion_xyzw'], dtype=float)
    t = np.asarray(c['target_from_source']['translation_m'], dtype=float)
    if q.shape != (4,) or t.shape != (3,) or not np.isfinite(q).all() or not np.isfinite(t).all() or abs(np.linalg.norm(q)-1) > 1e-6:
        raise ValueError('invalid rigid transform')
    return Rotation.from_quat(q), t


def transform_tum(args):
    r, t = load_transform(args.transform)
    a = read_tum(args.source)
    a[:, 1:4] = r.apply(a[:, 1:4]) + t
    a[:, 4:] = (r * Rotation.from_quat(a[:, 4:])).as_quat()
    with open(args.output, 'x') as f:
        np.savetxt(f, a, fmt='%.12f')


def transform_tasks(args):
    r, t = load_transform(args.transform)
    with open(args.source, newline='') as f:
        reader = csv.DictReader(f)
        names, rows = reader.fieldnames, list(reader)
    required = ['x', 'y', 'z', 'roll', 'pitch', 'yaw']
    if not names or not set(required).issubset(names) or not rows:
        raise ValueError('task CSV requires x,y,z,roll,pitch,yaw columns; preserves other columns')
    for row in rows:
        values = np.array([float(row[k]) for k in required])
        if not np.isfinite(values).all():
            raise ValueError('nonfinite task coordinate')
        p = r.apply(values[:3]) + t
        # ROS fixed-axis XYZ: roll, pitch, yaw; no CHC heading interpretation here.
        e = (r * Rotation.from_euler('xyz', values[3:], degrees=args.angles == 'degrees')).as_euler('xyz', degrees=args.angles == 'degrees')
        row.update({k: format(v, '.12g') for k, v in zip(required, np.r_[p, e])})
    with open(args.output, 'x', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)


def transform_grid(args):
    r, t = load_transform(args.transform)
    if np.linalg.norm(r.apply([0, 0, 1]) - [0, 0, 1]) > 1e-6:
        raise ValueError('tilted 3D transform cannot preserve a 2D occupancy grid; regenerate from transformed scans and ray origins')
    source = Path(args.source)
    with source.open() as f:
        c = yaml.safe_load(f)
    x, y, yaw = c['origin']
    p = r.apply([x, y, args.source_ground_z]) + t
    c['origin'] = [float(p[0]), float(p[1]), float(yaw + r.as_euler('xyz')[2])]
    image = (source.parent / c['image']).resolve()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(image, output / image.name)
    c['image'] = image.name
    with (output / 'map.yaml').open('x') as f:
        yaml.safe_dump(c, f, sort_keys=False)
    shutil.copyfile(args.transform, output / 'georeference.yaml')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('fit')
    p.add_argument('--source', required=True, help='local rear-axle TUM')
    p.add_argument('--target', required=True, help='accepted fixed-ENU rear-axle TUM')
    p.add_argument('--site', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--time-offset-sec', type=float, default=0)
    p.add_argument('--max-gap-sec', type=float, default=.1)
    p.add_argument('--min-span-m', type=float, default=1)
    p.add_argument('--max-rmse-m', type=float, required=True)
    p.add_argument('--max-error-m', type=float, required=True)
    p.add_argument('--max-angle-deg', type=float, required=True)
    p.set_defaults(func=fit)
    for name, fn in [('tum', transform_tum), ('tasks', transform_tasks), ('grid', transform_grid)]:
        p = commands.add_parser(name)
        p.add_argument('--source', required=True)
        p.add_argument('--transform', required=True)
        p.add_argument('--output', required=True)
        p.set_defaults(func=fn)
        if name == 'tasks':
            p.add_argument('--angles', choices=['degrees', 'radians'], required=True)
        if name == 'grid':
            p.add_argument('--source-ground-z', type=float, required=True)
    args = parser.parse_args()
    for k, v in vars(args).items():
        if isinstance(v, float) and (not np.isfinite(v) or (k not in ['time_offset_sec', 'source_ground_z'] and v <= 0)):
            parser.error(f'invalid {k}')
    args.func(args)


if __name__ == '__main__':
    main()
