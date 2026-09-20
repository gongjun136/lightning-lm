#!/usr/bin/env python3
"""Record declared compute settings; frame telemetry remains runtime authority."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import yaml


def summarize(path, environ):
    raw = path.read_bytes()
    config = yaml.safe_load(raw)
    budget = config.get('compute_budget', {})
    load = config.get('multi_lidar', {}).get('adaptive_load', {})
    overrides = {key: environ.get('LIGHTNING_LM_' + key.upper()) for key in
                 ('lio_threads', 'ndt_threads', 'ndt_max_points', 'solid_icp_workers', 'solid_worker_nice')}
    return dict(config_sha256=hashlib.sha256(raw).hexdigest(), compute_budget=budget,
                environment_overrides=overrides,
                requested_effective_budget={key: int(value) if value else budget.get(key) for key, value in overrides.items()},
                adaptive_load=load,
                lio_point_budget_enabled=bool(load.get('enabled') and load.get('lio_point_budgets')),
                note='Declared settings, not observed CPU occupancy; LIO_BENCH_FRAME records per-frame budgets. Null means unspecified.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = summarize(args.config, os.environ)
    with args.output.open('x', encoding='utf-8') as out:
        json.dump(result, out, ensure_ascii=False, indent=2)
    print('LIO point budgets:', result['adaptive_load'].get('lio_point_budgets', 'DISABLED'))
    print('Requested compute budget:', json.dumps(result['requested_effective_budget']))
    print('Configuration SHA256:', result['config_sha256'])
    if not result['lio_point_budget_enabled']:
        print('WARNING: LIO hard point budgets are NOT enabled; NDT_MAX_POINTS does not limit LIO.')


if __name__ == '__main__':
    main()
