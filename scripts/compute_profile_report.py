#!/usr/bin/env python3
"""Validate dedicated timing evidence after shutdown; never copy timing records."""
import argparse
from collections import Counter
import json
from pathlib import Path
import re
import sys

PROFILE = 'results/compute_profile.log'
STATUS = 'results/compute_profile_status.json'
HEADER = b'# COMPUTE_PROFILE schema=1 sampling=all sequence=profile_seq\n'
TIMING = re.compile(rb'\] ((COMPUTE_BENCH_[A-Z_]+|LIO_BENCH_FRAME)\b.*)')
COUNTERS = ('received', 'written', 'dropped', 'critical_dropped', 'io_errors', 'truncated')


def metadata(run):
    try:
        return dict(line.split('=', 1) for line in (Path(run) / 'run_metadata.txt').read_text().splitlines() if '=' in line)
    except OSError:
        return {}


def validate(run, enabled):
    """Full sequential scan only after the algorithm/writers have stopped."""
    path = Path(run) / PROFILE
    report = dict(schema=1, enabled=enabled, path=PROFILE, records=0,
                  status='disabled' if not enabled else 'missing', problems=[], modules={})
    if not enabled:
        if path.exists():
            report.update(status='incomplete', problems=['unexpected_profile_when_disabled'])
        return report
    if not path.is_file():
        report['problems'].append('profile_file_missing')
        return report
    footer = None
    sequence_errors = 0
    unknown = 0
    modules = Counter()
    with path.open('rb') as stream:
        if stream.readline() != HEADER:
            report['problems'].append('missing_schema_header')
        for line in stream:
            if line.startswith(b'COMPUTE_PROFILE_SUMMARY '):
                if footer is not None:
                    report['problems'].append('duplicate_footer')
                footer = dict((k.decode(), int(v)) for k, v in re.findall(rb'(\w+)=(\d+)', line))
                if not line.endswith(b'\n'):
                    report['problems'].append('partial_footer')
                continue
            match = TIMING.search(line)
            if not match:
                unknown += 1
                continue
            sequence = re.search(rb'\bprofile_seq=(\d+)\b', line)
            if not sequence or int(sequence[1]) != report['records'] or not line.endswith(b'\n'):
                sequence_errors += 1
            if footer is not None:
                report['problems'].append('records_after_footer')
            module = re.search(rb'\bmodule=([^\s]+)', match[1])
            modules[(module[1] if module else match[2]).decode(errors='replace')] += 1
            report['records'] += 1
    report.update(footer=footer, modules=dict(modules), sequence_errors=sequence_errors,
                  file_bytes=path.stat().st_size, file_mtime_ns=path.stat().st_mtime_ns)
    if unknown:
        report['problems'].append('unrecognized_or_partial_records')
    if sequence_errors:
        report['problems'].append('sequence_gaps_or_partial_records')
    if footer is None or any(key not in footer for key in COUNTERS):
        report['problems'].append('missing_or_invalid_footer')
    else:
        if footer['written'] != report['records'] or footer['received'] != footer['written'] + footer['dropped']:
            report['problems'].append('counter_mismatch')
        report['problems'].extend(key for key in ('dropped', 'critical_dropped', 'io_errors', 'truncated') if footer[key])
    if report['records'] == 0:
        report['problems'].append('no_timing_records')
    operational = Path(run) / 'logs/run_loc_online.stderr.log'
    if operational.is_file():
        with operational.open('rb') as stream:
            stream.seek(max(0, operational.stat().st_size - 64 * 1024))
            statuses = [line for line in stream if b'COMPUTE_PROFILE_STATUS path=' in line and b'COMPUTE_PROFILE_SUMMARY ' in line]
        if not statuses:
            report['problems'].append('operational_status_missing')
        else:
            counters = dict((k.decode(), int(v)) for k, v in re.findall(rb'(\w+)=(\d+)', statuses[-1]))
            report['operational_footer'] = {key: counters.get(key) for key in COUNTERS}
            # The operational status is emitted after fclose; the file's own
            # footer cannot describe an error encountered while closing it.
            if any(counters.get(key) != (footer or {}).get(key) for key in COUNTERS):
                report['problems'].append('post_close_counter_mismatch')
    report['status'] = 'incomplete' if report['problems'] else 'complete'
    return report


def inspect(run, live=False):
    """Use the final validation cache only while it matches the evidence file."""
    run = Path(run)
    meta = metadata(run)
    enabled = meta.get('compute_profile') == '1'
    if meta.get('compute_profile_storage') != 'dedicated-v1':
        return dict(status='legacy_unverified', enabled=enabled, path=PROFILE)
    path = run / PROFILE
    if not enabled:
        return dict(status='incomplete' if path.exists() else 'disabled', enabled=False, path=PROFILE)
    try:
        stat = path.stat()
        if (run / STATUS).stat().st_size <= 1024 * 1024:
            report = json.loads((run / STATUS).read_text())
            if (isinstance(report, dict) and report.get('schema') == 1 and report.get('enabled') is True and
                    report.get('path') == PROFILE and report.get('status') in ('complete', 'incomplete') and
                    report.get('file_bytes') == stat.st_size and report.get('file_mtime_ns') == stat.st_mtime_ns):
                return report
    except (OSError, ValueError, TypeError):
        pass
    return dict(status=('recording' if live else 'unverified') if path.is_file() else 'missing',
                enabled=True, path=PROFILE)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    args = parser.parse_args()
    report = validate(args.run, metadata(args.run).get('compute_profile') == '1')
    (args.run / STATUS).write_text(json.dumps(report, indent=2) + '\n')
    print('compute_profile_record_count=' + str(report['records']))
    print('compute_profile_integrity=' + report['status'])
    if report['status'] not in ('complete', 'disabled'):
        print('WARNING: compute profile evidence is incomplete; see ' + str(args.run / STATUS), file=sys.stderr)
    return 0  # Evidence failure must not replace the algorithm's exit status.


if __name__ == '__main__':
    sys.exit(main())
