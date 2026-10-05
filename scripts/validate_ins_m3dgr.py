#!/usr/bin/env python3
"""Read-only M3DGR GNSS audit and GeographicLib vs independent ECEF check.

This validates numerical conversion, not GNSS/CGI absolute positioning accuracy.
Dependencies: rosbags, numpy. No ROS installation is needed by the reader.
"""
import argparse
from collections import Counter
import csv
import io
import json
from pathlib import Path
import subprocess

import numpy as np
from rosbags.highlevel import AnyReader


def ecef(llh):
    lat, lon, h = np.asarray(llh, dtype=float).T
    lat, lon = np.deg2rad(lat), np.deg2rad(lon)
    a, f = 6378137., 1/298.257223563
    e2 = f*(2-f)
    n = a/np.sqrt(1-e2*np.sin(lat)**2)
    return np.column_stack(((n+h)*np.cos(lat)*np.cos(lon), (n+h)*np.cos(lat)*np.sin(lon),
                            (n*(1-e2)+h)*np.sin(lat)))


def independent_enu(llh, origin):
    lat, lon = np.deg2rad(origin[:2])
    r = np.array([[-np.sin(lon), np.cos(lon), 0],
                  [-np.sin(lat)*np.cos(lon), -np.sin(lat)*np.sin(lon), np.cos(lat)],
                  [np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)]])
    return (ecef(llh)-ecef([origin])) @ r.T


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bag', required=True)
    p.add_argument('--projector', required=True, help='built ins_project_llh executable')
    p.add_argument('--ground-truth', help='optional M3DGR TUM; only audit, do not align or score')
    p.add_argument('--output', required=True, help='new directory')
    a = p.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    rows, connections = [], []
    with AnyReader([Path(a.bag)]) as reader:
        connections = [dict(topic=c.topic, type=c.msgtype, count=c.msgcount) for c in reader.connections]
        selected = [c for c in reader.connections if c.topic == '/ublox_driver/receiver_pvt']
        if len(selected) != 1:
            raise ValueError('requires one M3DGR /ublox_driver/receiver_pvt connection')
        for c, stamp, raw in reader.messages(connections=selected):
            m = reader.deserialize(raw, c.msgtype)
            rows.append([stamp*1e-9, m.time.week, m.time.tow, m.latitude, m.longitude, m.altitude,
                         m.height_msl, m.fix_type, int(m.valid_fix), int(m.diff_soln), m.carr_soln,
                         m.h_acc, m.v_acc, m.vel_e, m.vel_n, -m.vel_d])
    with (out/'gnss.csv').open('x', newline='') as f:
        w = csv.writer(f)
        w.writerow(['bag_stamp','gps_week','gps_tow','latitude','longitude','ellipsoid_height','height_msl',
                    'fix_type','valid_fix','diff_soln','carr_soln','h_acc','v_acc','ve','vn','vu'])
        w.writerows(rows)
    data = np.asarray(rows)
    valid = np.isfinite(data[:, 3:6]).all(1) & (data[:, 8] == 1)
    llh = data[valid, 3:6]
    if len(llh) < 2:
        raise ValueError('not enough finite valid GNSS LLH values')
    origin = llh[0]
    stream = io.StringIO()
    np.savetxt(stream, llh, fmt='%.16g')
    run = subprocess.run([a.projector, *map(str, origin)], input=stream.getvalue(), text=True,
                         capture_output=True, check=True)
    actual = np.loadtxt(io.StringIO(run.stdout), ndmin=2)
    expected = independent_enu(llh, origin)
    errors = np.linalg.norm(actual-expected, axis=1)
    np.savetxt(out/'enu_comparison.csv', np.c_[data[valid, 0], actual, expected, errors], delimiter=',',
               header='stamp,e,n,u,independent_e,independent_n,independent_u,error_m', comments='')
    report = dict(scope='numerical coordinate conversion only; not CGI or GNSS accuracy',
        bag=str(Path(a.bag).resolve()), bag_size_bytes=Path(a.bag).stat().st_size,
        source_topic='/ublox_driver/receiver_pvt', height='ellipsoidal altitude, not height_msl',
        connections=connections, rows=len(rows), finite_valid_rows=int(valid.sum()),
        carrier_solution_counts=dict(Counter(map(str, data[:, 10].astype(int)))),
        fixed_origin_llh=origin.tolist(), max_conversion_difference_m=float(max(errors)),
        rms_conversion_difference_m=float(np.sqrt(np.mean(errors**2))),
        acceptance_tolerance_m=1e-7, passed=bool(max(errors)<1e-7))
    if a.ground_truth:
        gt = np.loadtxt(a.ground_truth, ndmin=2)
        report['ground_truth_audit'] = dict(path=str(Path(a.ground_truth).resolve()), rows=len(gt),
            adjacent_duplicate_stamps=int(np.count_nonzero(np.diff(gt[:, 0]) == 0)),
            backward_stamps=int(np.count_nonzero(np.diff(gt[:, 0]) < 0)),
            all_identity_quaternions=bool(np.allclose(gt[:, 4:], [0,0,0,1])),
            note='No full-attitude accuracy claim and no preprocessing of original data')
    (out/'report.json').write_text(json.dumps(report, indent=2))
    print(json.dumps({k:v for k,v in report.items() if k!='connections'}, indent=2))
    if not report['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
