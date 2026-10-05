#!/usr/bin/env python3
"""Black-box synthetic tests for frame migration and held-out rejection."""
import argparse
import csv
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.spatial.transform import Rotation
import yaml


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--map-binary')
    args = parser.parse_args()
    tool = str(Path(__file__).with_name('ins_georeference.py'))
    with tempfile.TemporaryDirectory(prefix='ins_georeference_test_') as temporary:
        root = Path(temporary)
        def run(*arguments, good=True):
            result = subprocess.run([sys.executable, tool, *map(str, arguments)], capture_output=True, text=True)
            assert (result.returncode == 0) == good, result.stdout + result.stderr
            return result
        times = np.arange(401)*.01+100
        theta = np.linspace(0, 2*np.pi, len(times))
        xyz = np.c_[20*np.cos(theta), 10*np.sin(theta), np.sin(3*theta)]
        r = Rotation.from_euler('xyz', [4, -6, 35], degrees=True)
        translation = np.array([12., 30., -3.])
        source = np.c_[times-.013, xyz, np.tile([0,0,0,1], (len(times),1))]
        target = np.c_[times, r.apply(xyz)+translation, np.tile(r.as_quat(), (len(times),1))]
        np.savetxt(root/'source.tum', source)
        np.savetxt(root/'target.tum', target)
        site = dict(ins_only=dict(georeference=dict(confirmed=True, datum='WGS84', origin_llh=[30,120,0])))
        (root/'site.yaml').write_text(yaml.safe_dump(site))
        common = ['fit','--source',root/'source.tum','--target',root/'target.tum','--site',root/'site.yaml',
                  '--time-offset-sec','.013','--max-rmse-m','.02','--max-error-m','.05','--max-angle-deg','.1']
        run(*common, '--output', root/'transform.yaml')
        fit = yaml.safe_load((root/'transform.yaml').read_text())
        assert fit['accepted'] and fit['validation']['held_out']['position_max_m'] < 1e-10
        assert np.linalg.norm(np.array(fit['target_from_source']['translation_m'])-translation)<1e-10
        run(*common, '--output', root/'transform.yaml', good=False)
        run('tum','--source',root/'source.tum','--transform',root/'transform.yaml','--output',root/'converted.tum')
        actual = np.loadtxt(root/'converted.tum')
        assert np.allclose(actual[:,1:], target[:,1:], atol=1e-10)
        assert np.allclose(actual[:,0], source[:,0]), 'frame conversion must not alter time'
        with (root/'tasks.csv').open('w') as f:
            f.write('id,x,y,z,roll,pitch,yaw\npark,1,2,3,0,0,0\n')
        run('tasks','--source',root/'tasks.csv','--transform',root/'transform.yaml', '--output',root/'tasks_enu.csv','--angles','degrees')
        with (root/'tasks_enu.csv').open() as f:
            row = next(csv.DictReader(f))
        assert row['id'] == 'park'
        assert np.allclose([float(row[k]) for k in ['x','y','z']], r.apply([1,2,3])+translation)
        assert np.allclose([float(row[k]) for k in ['roll','pitch','yaw']], [4,-6,35])
        # A fit to the early segment must not conceal drift in the held-out tail.
        target[300:,1] += 1
        np.savetxt(root/'target.tum', target)
        run(*common, '--output', root/'bad.yaml', good=False)
        assert not yaml.safe_load((root/'bad.yaml').read_text())['accepted']
        run('tum','--source',root/'source.tum','--transform',root/'bad.yaml','--output',root/'bad.tum',good=False)
        straight = source.copy()
        straight[:,1:4] = np.c_[np.linspace(0,50,len(times)), np.zeros((len(times),2))]
        np.savetxt(root/'source.tum', straight)
        run(*common, '--output', root/'line.yaml', good=False)
        np.savetxt(root/'source.tum', source)
        # PGM cell data is preserved for planar transforms; tilted transforms reject.
        (root/'map.pgm').write_bytes(b'P5\n2 2\n255\n'+bytes([0,128,255,0]))
        (root/'map.yaml').write_text(yaml.safe_dump(dict(image='map.pgm', resolution=.1, origin=[1,2,.2])))
        run('grid','--source',root/'map.yaml','--transform',root/'transform.yaml','--output',root/'tilted','--source-ground-z','0',good=False)
        planar = dict(fit)
        planar['target_from_source'] = dict(translation_m=[10,20,3],quaternion_xyzw=Rotation.from_euler('z',90,degrees=True).as_quat().tolist())
        (root/'planar.yaml').write_text(yaml.safe_dump(planar))
        run('grid','--source',root/'map.yaml','--transform',root/'planar.yaml','--output',root/'planar','--source-ground-z','0')
        grid = yaml.safe_load((root/'planar/map.yaml').read_text())
        assert np.allclose(grid['origin'], [8,21,.2+np.pi/2])
        assert (root/'planar/map.pgm').read_bytes() == (root/'map.pgm').read_bytes()
        if args.map_binary:
            (root/'map.pcd').write_text('VERSION .7\nFIELDS x y z intensity\nSIZE 4 4 4 4\nTYPE F F F F\nCOUNT 1 1 1 1\nWIDTH 3\nHEIGHT 1\nVIEWPOINT 0 0 0 1 0 0 0\nPOINTS 3\nDATA ascii\n1 2 3 10\n4 5 6 20\n7 8 9 30\n')
            command = [args.map_binary, str(root/'map.pcd'),str(root/'transform.yaml'), str(root/'source.tum'),str(root/'newmap')]
            result = subprocess.run(command, capture_output=True, text=True)
            assert result.returncode == 0, result.stderr
            assert (root/'newmap/tiled/index.txt').exists() and (root/'newmap/map.pcd').exists()
            assert subprocess.run(command, capture_output=True).returncode != 0, 'must reject overwriting map'
        print('PASS: rigid transform, time offset, held-out drift rejection, degeneracy, task/TUM/grid migration, overwrite protection' + (', PCD/tiled export' if args.map_binary else ''))


if __name__ == '__main__':
    main()
