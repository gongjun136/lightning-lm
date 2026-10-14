#!/usr/bin/env python3
"""Convert a trajectory's rigid sensor reference point to rear-axle FLU.

Preserves world frame, timestamps, and map normalization. Never rebases a
trajectory or estimates extrinsics. T_AB maps coordinates in B into A.
"""
import argparse
from decimal import Decimal, ROUND_HALF_EVEN
import hashlib
import json
from pathlib import Path

import numpy as np
import yaml


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def quaternion(value):
    q = np.asarray(value, dtype=float)
    if q.shape != (4,) or not np.isfinite(q).all() or abs(np.linalg.norm(q) - 1) > 1e-5:
        raise ValueError('quaternion must be a finite unit xyzw quaternion')
    return q / np.linalg.norm(q)


def multiply(a, b):
    return np.r_[a[3] * b[:3] + b[3] * a[:3] + np.cross(a[:3], b[:3]),
                 a[3] * b[3] - np.dot(a[:3], b[:3])]


def rotate(q, v):
    return v + 2 * np.cross(q[:3], np.cross(q[:3], v) + q[3] * v)


def transform(config, name):
    value = config.get(name, {})
    if value.get('confirmed') is not True or not value.get('evidence'):
        raise ValueError(f'{name}: confirmed extrinsic and evidence are required')
    t = np.asarray(value.get('translation_m'), dtype=float)
    if t.shape != (3,) or not np.isfinite(t).all():
        raise ValueError(f'{name}: three finite translation values are required')
    return quaternion(value.get('quaternion_xyzw')), t


def source_from_body(config):
    if config.get('schema_version') != 1:
        raise ValueError('schema_version must equal 1')
    source = config.get('source_frame')
    if source == 'rear_axle_flu':
        if config.get('published_reference_point_confirmed') is not True or not config.get('published_reference_point_evidence'):
            raise ValueError('published rear-axle point/body convention requires confirmation and evidence')
        return np.array([0., 0., 0., 1.]), np.zeros(3)
    if source not in ('imu', 'primary_lidar'):
        raise ValueError('source_frame must be imu, primary_lidar, or rear_axle_flu')
    q_bl, t_bl = transform(config, 'T_rear_axle_lidar')
    q_lb = np.r_[-q_bl[:3], q_bl[3]]
    t_lb = -rotate(q_lb, t_bl)
    if source == 'primary_lidar':
        return q_lb, t_lb
    q_il, t_il = transform(config, 'T_imu_lidar')
    return multiply(q_il, q_lb), t_il + rotate(q_il, t_lb)


def normalize(source, config_path, recording_metadata, output):
    with Path(config_path).open(encoding='utf-8') as f:
        c = yaml.safe_load(f)
    with Path(recording_metadata).open(encoding='utf-8') as f:
        metadata = json.load(f) if str(recording_metadata).endswith('.json') else yaml.safe_load(f)
    q_sb, t_sb = source_from_body(c)
    for key in ('source_recording_id', 'time_basis', 'recorded_intervals'):
        if not metadata.get(key):
            raise ValueError(f'recording metadata missing {key}; do not infer recording duration from successful poses')
    for key in ('world_frame_id', 'trajectory_kind'):
        if not c.get(key):
            raise ValueError(f'config missing {key}')
    if c['trajectory_kind'] not in ('frontend', 'optimized_keyframes', 'optimized_dense_reconstructed', 'localization'):
        raise ValueError('invalid trajectory_kind')
    intervals = metadata['recorded_intervals']
    for interval in intervals:
        if not all(type(interval.get(k)) is int for k in ('start_ns', 'end_ns')) or interval['start_ns'] >= interval['end_ns'] or not isinstance(interval.get('epoch'), str) or not interval['epoch'].strip():
            raise ValueError('recorded_intervals require integer start_ns < end_ns and a named string epoch')
    if len({x['epoch'] for x in intervals}) != 1:
        raise ValueError('TUM has no epoch column: export one clock epoch per physical file')
    source_world = c.get('source_world_frame_id', c['world_frame_id'])
    if metadata.get('world_frame_id') and metadata['world_frame_id'] != source_world:
        raise ValueError('input metadata world_frame_id differs from declared source world; numeric poses cannot be relabeled')
    if metadata.get('body_frame_id') and metadata['body_frame_id'] != c['source_frame']:
        raise ValueError('input metadata body_frame_id differs from source_frame; use frame-free recording metadata for raw mapping IMU output')
    offset_ns = c.get('timestamp_offset_ns')
    if type(offset_ns) is not int:
        raise ValueError('timestamp_offset_ns must be an explicit integer')
    world_transform = c.get('world_transform', {})
    if world_transform.get('enabled'):
        q_tw, t_tw = transform(c, 'world_transform')
        if not source_world or source_world == c['world_frame_id']:
            raise ValueError('world_transform requires distinct, explicit source and target world_frame_id')
    else:
        if source_world != c['world_frame_id']:
            raise ValueError('different source/target worlds require an explicit confirmed world_transform')
        q_tw, t_tw = np.array([0., 0., 0., 1.]), np.zeros(3)
    rows, timestamps, time_tokens = [], [], []
    with Path(source).open(encoding='utf-8') as f:
        for number, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) != 8:
                raise ValueError(f'TUM line {number}: expected 8 columns')
            v = np.asarray(fields, dtype=float)
            if not np.isfinite(v).all():
                raise ValueError(f'TUM line {number}: nonfinite value')
            q = quaternion(v[4:])
            stamp = int((Decimal(fields[0]) * 1_000_000_000).to_integral_value(rounding=ROUND_HALF_EVEN)) + offset_ns
            if timestamps and stamp <= timestamps[-1]:
                raise ValueError('TUM timestamps must strictly increase; never silently sort/deduplicate')
            if not any(x['start_ns'] <= stamp <= x['end_ns'] for x in intervals):
                raise ValueError('TUM timestamp outside recorded interval: check receive/header time bases')
            v[1:4] += rotate(q, t_sb)
            v[4:] = multiply(q, q_sb)
            v[1:4] = rotate(q_tw, v[1:4]) + t_tw
            v[4:] = multiply(q_tw, v[4:])
            rows.append(v)
            timestamps.append(stamp)
            time_tokens.append(fields[0] if offset_ns == 0 else format(Decimal(stamp) / 1_000_000_000, '.9f'))
    out = Path(output)
    sidecar = out.with_suffix('.metadata.json')
    if out.exists() or sidecar.exists():
        raise ValueError('output or metadata exists; choose a new artifact path')
    out.parent.mkdir(parents=True, exist_ok=True)
    metadata.update(schema_version=1, world_frame_id=c['world_frame_id'], body_frame_id='rear_axle_flu',
                    trajectory_kind=c['trajectory_kind'], clock_epoch=intervals[0]['epoch'],
                    source_sha256=sha256(source), normalization_config_sha256=sha256(config_path),
                    source_frame=c['source_frame'], source_world_frame_id=source_world,
                    applied_timestamp_offset_ns=offset_ns,
                    world_transform_applied=bool(world_transform.get('enabled')),
                    convention='T_world_body=T_world_source*T_source_body; no rebasing',
                    source_from_body={'quaternion_xyzw': q_sb.tolist(), 'translation_m': t_sb.tolist()},
                    extrinsic_evidence={k: c[k] for k in ('T_imu_lidar', 'T_rear_axle_lidar',
                        'published_reference_point_evidence', 'world_transform') if k in c})
    with out.open('x', encoding='utf-8', newline='\n') as f:
        for token, row in zip(time_tokens, rows):
            f.write(token + ' ' + ' '.join(f'{v:.12f}' for v in row[1:]) + '\n')
    metadata['trajectory_sha256'] = sha256(out)
    with sidecar.open('x', encoding='utf-8', newline='\n') as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write('\n')
    return {'trajectory': str(out), 'metadata': str(sidecar), 'rows': len(rows)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', required=True)
    p.add_argument('--config', required=True)
    p.add_argument('--recording-metadata', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    try:
        print(json.dumps(normalize(a.source, a.config, a.recording_metadata, a.output), ensure_ascii=False))
    except (ValueError, KeyError, OSError, TypeError) as error:
        p.exit(2, f'accuracy_normalize_tum: {error}\n')


if __name__ == '__main__':
    main()
