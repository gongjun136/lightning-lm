#!/usr/bin/env python3
"""Compare fixed extrinsics on held-out stationary ROS 2 bags without ICP fitting.

Run after sourcing ROS Humble. Samples cloud frames 30..39 for each LiDAR.
Reports conditional reciprocal local-plane residuals, overlap and ground diagnostics;
these do not certify full six-DOF calibration or absolute accuracy.
"""
import argparse
import hashlib
import importlib.util
from itertools import combinations
import json
from pathlib import Path
import sqlite3

import numpy as np
from rclpy.serialization import deserialize_message
from rosidl_runtime_py.utilities import get_message
import yaml


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_clouds(bag, config):
    files = list(bag.glob('*.db3'))
    if len(files) != 1:
        raise ValueError('Expected one SQLite split in this assessment: ' + str(bag))
    message_type = get_message('sensor_msgs/msg/PointCloud2')
    imu_type = get_message('sensor_msgs/msg/Imu')
    formats = {1:'i1', 2:'u1', 3:'i2', 4:'u2', 5:'i4', 6:'u4', 7:'f4', 8:'f8'}
    clouds = {}
    with sqlite3.connect(files[0].resolve().as_uri() + '?mode=ro&immutable=1', uri=True) as data:
        topics = dict(data.execute('SELECT name,id FROM topics'))
        imu_topic = config['multi_lidar']['topics']['imu_0']
        accelerations = []
        for payload, in data.execute('SELECT data FROM messages WHERE topic_id=? ORDER BY timestamp LIMIT 200 OFFSET 600', (topics[imu_topic],)):
            msg = deserialize_message(payload, imu_type)
            accelerations.append([msg.linear_acceleration.x, msg.linear_acceleration.y, msg.linear_acceleration.z])
        mean_acceleration = np.mean(accelerations, axis=0)
        if not np.isfinite(mean_acceleration).all() or np.linalg.norm(mean_acceleration) < 0.1:
            raise ValueError('Insufficient primary IMU acceleration for ground search')
        # Config rotation is T_imu_from_lidar. Only a broad ground-search direction;
        # this diagnostic does not validate the inherited LiDAR/IMU calibration.
        rotation_imu_lidar = np.asarray(config['fasterlio']['extrinsic_R']).reshape(3,3)
        upward = rotation_imu_lidar.T @ mean_acceleration
        upward /= np.linalg.norm(upward)
        for i in range(4):
            topic = config['multi_lidar']['topics']['lidar_' + str(i)]
            points = []
            for payload, in data.execute('SELECT data FROM messages WHERE topic_id=? ORDER BY timestamp LIMIT 10 OFFSET 30', (topics[topic],)):
                msg = deserialize_message(payload, message_type)
                if msg.row_step != msg.width * msg.point_step or msg.is_bigendian:
                    raise ValueError('Unsupported padded/big-endian cloud')
                dtype = np.dtype({'names':[f.name for f in msg.fields], 'formats':[formats[f.datatype] for f in msg.fields],
                                  'offsets':[f.offset for f in msg.fields], 'itemsize':msg.point_step})
                cloud = np.frombuffer(msg.data, dtype=dtype)
                xyz = np.column_stack([cloud[key] for key in ['x','y','z']]).astype(float)
                radius = np.linalg.norm(xyz, axis=1)
                points.append(xyz[np.isfinite(xyz).all(axis=1) & (radius >= 4) & (radius <= 45)])
            if len(points) != 10:
                raise ValueError('Insufficient stationary sample frames: ' + topic)
            clouds[i] = np.concatenate(points)
    return clouds, upward


def ground_plane(points, upward, qc, seed):
    points = points[np.linalg.norm(points, axis=1) <= 25]
    points = qc.voxel_sample(points, 0.15, 12000, seed)
    rng = np.random.default_rng(seed)
    best = np.zeros(len(points), dtype=bool)
    for _ in range(350):
        sample = points[rng.choice(len(points), 3, replace=False)]
        normal = np.cross(sample[1]-sample[0], sample[2]-sample[0])
        norm = np.linalg.norm(normal)
        if norm < 1e-8: continue
        normal /= norm
        if np.dot(normal, upward) < 0: normal *= -1
        distance = -np.dot(normal, sample[0])
        if not 1 < abs(distance) < 5 or np.dot(normal, upward) < np.cos(np.deg2rad(15)):
            continue
        inliers = np.abs(points @ normal + distance) < 0.10
        if inliers.sum() > best.sum(): best = inliers
    if best.sum() < 100:
        return {'status':'insufficient_ground_plane_points'}
    center = points[best].mean(axis=0)
    _, _, basis = np.linalg.svd(points[best]-center, full_matrices=False)
    normal = basis[-1]
    if np.dot(normal, upward) < 0: normal *= -1
    return {'status':'measured', 'normal':normal.tolist(), 'distance_m':float(-np.dot(normal, center)),
            'inliers':int(best.sum()), 'sample_points':len(points)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bag', action='append', required=True, type=Path)
    parser.add_argument('--baseline-config', required=True, type=Path)
    parser.add_argument('--new-config', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    module_path = Path(__file__).parent / 'reproduction/multi_lidar/sany_4livox/analyze_sany_map_cross_source.py'
    spec = importlib.util.spec_from_file_location('map_qc', module_path)
    qc = importlib.util.module_from_spec(spec); spec.loader.exec_module(qc)
    configs = {key:yaml.safe_load(path.read_text()) for key, path in [('baseline',args.baseline_config), ('new',args.new_config)]}
    if configs['baseline']['multi_lidar']['topics'] != configs['new']['multi_lidar']['topics']:
        raise ValueError('Different sensor topic mappings')
    if configs['new']['multi_lidar']['primary_lidar_id'] != 0:
        raise ValueError('This assessment expects lidar0 as primary')
    results, aggregate = [], {'baseline':[], 'new':[]}
    for bag in args.bag:
        clouds, upward = read_clouds(bag, configs['new'])
        result = {'bag':bag.name, 'read_only_processing_path':str(bag), 'configurations':{}}
        result['ground_search_upward_in_primary_frame'] = upward.tolist()
        for key, config in configs.items():
            transformed = {}
            for i, cloud in clouds.items():
                transform = np.array(config['multi_lidar']['extrinsics']['lidar' + str(i)]['T']).reshape(4,4)
                transformed[i] = qc.voxel_sample(cloud @ transform[:3,:3].T + transform[:3,3], 0.15, 30000, 20261010+i)
            pairs, residuals = {}, []
            for left, right in combinations(range(4), 2):
                lr, _, lo = qc.directional_planar_residual(transformed[left], transformed[right], 0.5, 6000, 20260711+left*10+right)
                rl, _, ro = qc.directional_planar_residual(transformed[right], transformed[left], 0.5, 6000, 20260711+right*10+left)
                values = np.concatenate([lr,rl]); residuals.append(values)
                pairs[f'{left}-{right}'] = {'point_to_plane_m':qc.distribution(values), 'overlap_left_to_right':lo, 'overlap_right_to_left':ro}
            values = np.concatenate(residuals); aggregate[key].append(values)
            planes = {str(i):ground_plane(cloud, upward, qc, 20261010+i) for i,cloud in transformed.items()}
            for i,plane in planes.items():
                if plane['status'] == 'measured' and planes['0']['status'] == 'measured':
                    plane['normal_difference_from_primary_deg'] = float(np.degrees(np.arccos(np.clip(np.dot(plane['normal'], planes['0']['normal']),-1,1))))
                    plane['plane_distance_difference_from_primary_m'] = float(plane['distance_m'] - planes['0']['distance_m'])
            result['configurations'][key] = {'aggregate_point_to_plane_m':qc.distribution(values), 'pairs':pairs, 'ground_planes_in_primary_frame':planes}
        results.append(result)
        print(json.dumps({'bag':bag.name, **{k:v['aggregate_point_to_plane_m'] for k,v in result['configurations'].items()}}), flush=True)
    output = {'method':'held_out_stationary_clouds_fixed_extrinsics_no_ICP_refinement',
              'config_sha256':{'baseline':sha256(args.baseline_config), 'new':sha256(args.new_config)},
              'script_sha256':sha256(Path(__file__)), 'qc_script_sha256':sha256(module_path),
              'parameters':{'frames_per_sensor':10, 'frame_offset':30, 'raw_range_m':[4,45], 'voxel_m':0.15,
                            'max_correspondence_distance_m':0.5, 'max_pairs_per_direction':6000, 'max_points_per_sensor':30000,
                            'ground_search':{'range_max_m':25, 'max_sample_points':12000, 'ransac_iterations':350,
                                             'inlier_threshold_m':0.10, 'max_angle_from_imu_deg':15,
                                             'plane_distance_m':[1,5], 'seed_base':20261010,
                                             'upward_source':'primary mean IMU acceleration using inherited extrinsic_R; not an IMU calibration verdict'}},
              'aggregate':{k:qc.distribution(np.concatenate(v)) for k,v in aggregate.items()}, 'bags':results,
              'interpretation_limit':'Conditional reciprocal planar correspondences within 0.5 m; scene selection and vegetation affect statistics. Ground planes can cover different sloped surfaces. No absolute truth, no six-DOF observability certification; no extrinsic fitting performed.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
