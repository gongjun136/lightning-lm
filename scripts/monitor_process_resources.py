#!/usr/bin/env python3
"""Low-rate Linux /proc sampling. CPU equivalents are CPU seconds / wall seconds.

Only observes the run_loc_online executable in the launcher's process group;
never attaches to another deployment and never changes affinity or scheduling.
No third-party Python dependency. Missing counters are null, not zero.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import sys
import threading
import time

ENV_KEYS = ('OMP_NUM_THREADS', 'OMP_WAIT_POLICY', 'OMP_PROC_BIND', 'OMP_PLACES',
            'GOMP_SPINCOUNT', 'LIGHTNING_LM_LIO_THREADS', 'LIGHTNING_LM_NDT_THREADS',
            'LIGHTNING_LM_NDT_MAX_POINTS', 'LIGHTNING_LM_CPU_AFFINITY',
            'LIGHTNING_LM_SOLID_ICP_WORKERS', 'LIGHTNING_LM_SOLID_CPU_AFFINITY')


def parse_stat(text):
    left, right = text.index('('), text.rindex(')')
    fields = text[right + 2:].split()
    return dict(pid=int(text[:left].strip()), name=text[left + 1:right], state=fields[0],
                pgrp=int(fields[2]), cpu_ticks=int(fields[11]) + int(fields[12]),
                threads=int(fields[17]), start_ticks=int(fields[19]),
                last_cpu=int(fields[36]))


def read_status(path):
    result = {}
    for line in path.read_text().splitlines():
        if ':' in line:
            key, value = line.split(':', 1)
            result[key] = value.strip()
    return result


def optional_int(text):
    return int(text.split()[0]) if text else None


def cpu_list_count(text):
    cpus = set()
    for item in text.split(','):
        if not item: continue
        bounds = item.split('-')
        cpus.update(range(int(bounds[0]), int(bounds[-1]) + 1))
    return len(cpus)


def read_thread(path):
    result = parse_stat((path / 'stat').read_text())
    status = read_status(path / 'status')
    result.update(allowed_cpus=status.get('Cpus_allowed_list'),
                  voluntary_switches=optional_int(status.get('voluntary_ctxt_switches')),
                  involuntary_switches=optional_int(status.get('nonvoluntary_ctxt_switches')))
    try:
        run_ns, wait_ns, slices = map(int, (path / 'schedstat').read_text().split()[:3])
        result.update(sched_run_ns=run_ns, sched_wait_ns=wait_ns, sched_slices=slices)
    except (OSError, ValueError):
        result.update(sched_run_ns=None, sched_wait_ns=None, sched_slices=None)
    try: result['wchan'] = (path / 'wchan').read_text().strip() or None
    except OSError: result['wchan'] = None
    return result


def counter_delta(current, previous, key):
    a, b = current.get(key), previous.get(key)
    return a - b if a is not None and b is not None and a >= b else None


def thread_interval(current, previous, elapsed, ticks_per_sec):
    # TIDs can be reused; an unpaired new thread has no interval estimate.
    previous = previous if previous and previous['start_ticks'] == current['start_ticks'] else {}
    ticks = counter_delta(current, previous, 'cpu_ticks')
    result = dict(tid=current['pid'], start_ticks=current['start_ticks'], name=current['name'],
                  state=current['state'], last_cpu=current['last_cpu'], allowed_cpus=current['allowed_cpus'],
                  cpu_core_equivalents=ticks / ticks_per_sec / elapsed if ticks is not None and elapsed > 0 else None)
    result['wchan'] = current.get('wchan')  # Read-time snapshot, not interval attribution.
    for key in ('sched_run_ns', 'sched_wait_ns', 'sched_slices', 'voluntary_switches', 'involuntary_switches'):
        result[key + '_delta'] = counter_delta(current, previous, key)
    return result


def read_host_cpus(proc):
    result = {}
    for line in (proc / 'stat').read_text().splitlines():
        fields = line.split()
        if fields and fields[0].startswith('cpu') and fields[0][3:].isdigit():
            values = list(map(int, fields[1:9]))  # guest time is already in user/nice
            result[fields[0]] = (sum(values), values[3], values[4])
    return result


def host_intervals(current, previous):
    rows = []
    for cpu, values in current.items():
        old = previous.get(cpu)
        delta = [a-b for a,b in zip(values, old)] if old else None
        valid = delta is not None and delta[0] > 0 and min(delta) >= 0
        rows.append(dict(cpu=cpu, busy_percent=100*(delta[0]-delta[1]-delta[2])/delta[0] if valid else None,
                         iowait_percent=100*delta[2]/delta[0] if valid else None))
    return rows


def discover_target(proc, pgrp, executable):
    matches = []
    for path in proc.iterdir():
        if not path.name.isdigit(): continue
        try:
            stat = parse_stat((path / 'stat').read_text())
            if stat['pgrp'] == pgrp and Path(os.readlink(path / 'exe')).name == executable:
                matches.append(stat)
        except (OSError, ValueError, IndexError):
            continue  # Process exited during discovery, or /proc is restricted.
    if len(matches) > 1:
        raise RuntimeError('multiple target executables in the launch process group')
    return matches[0] if matches else None


def read_system_processes(proc):
    """Only /proc stat counters and comm; no command lines or environments."""
    rows = {}; missing = 0
    for path in proc.iterdir():
        if not path.name.isdigit(): continue
        try:
            row = parse_stat((path / 'stat').read_text())
            rows[row['pid']] = row
        except (OSError, ValueError, IndexError): missing += 1
    return rows, missing


def system_process_intervals(current, previous, elapsed, ticks_per_sec, top_count=20):
    rows = []; unpaired = 0
    for pid, row in current.items():
        old = previous.get(pid)
        if old is None or old['start_ticks'] != row['start_ticks']:
            unpaired += 1
            continue
        delta = counter_delta(row, old, 'cpu_ticks')
        if delta is None or elapsed <= 0: continue
        rows.append(dict(pid=pid, start_ticks=row['start_ticks'], name=row['name'],
                         pgrp=row['pgrp'], state=row['state'], threads=row['threads'],
                         cpu_core_equivalents=delta/ticks_per_sec/elapsed))
    rows.sort(key=lambda r:(-r['cpu_core_equivalents'],r['pid']))
    return dict(top=rows[:top_count], ranked_processes=len(rows), unpaired_processes=unpaired,
                vanished_processes=len(set(previous)-set(current)), top_count=top_count)


def process_snapshot(proc, pid, include_pss):
    path = proc / str(pid)
    stamp = time.monotonic()
    stat = parse_stat((path / 'stat').read_text())
    if stat['state'] in ('Z', 'X'):
        raise ProcessLookupError('target has exited')
    status = read_status(path / 'status')
    threads = {}; missed = 0
    for task in (path / 'task').iterdir():
        try: threads[int(task.name)] = read_thread(task)
        except (OSError, ValueError, IndexError): missed += 1
    memory = {key + '_kb': optional_int(status.get(key)) for key in ('VmRSS', 'VmHWM', 'VmSize', 'VmSwap', 'RssAnon', 'RssFile')}
    pss_error = None
    if include_pss:
        try:
            rollup = read_status(path / 'smaps_rollup')
            memory.update({key + '_kb': optional_int(rollup.get(key)) for key in ('Pss', 'Private_Clean', 'Private_Dirty', 'SwapPss')})
        except OSError as exc: pss_error = type(exc).__name__
    try: io = {k:int(v) for k,v in read_status(path / 'io').items()}
    except OSError: io = {}
    system_processes, system_missing = read_system_processes(proc)
    # Protect against PID reuse across the multi-file read.
    if parse_stat((path / 'stat').read_text())['start_ticks'] != stat['start_ticks']:
        raise RuntimeError('target PID reused during sample')
    return dict(monotonic_s=stamp, stat=stat, memory=memory, threads=threads, missing_threads=missed,
                allowed_cpus=status.get('Cpus_allowed_list', ''), io=io, pss_error=pss_error,
                host_cpus=read_host_cpus(proc), pss_sampled=include_pss,
                system_processes=system_processes, system_missing=system_missing)


def sample_interval(current, previous, ticks_per_sec):
    elapsed = current['monotonic_s'] - previous['monotonic_s'] if previous else 0
    ticks = counter_delta(current['stat'], previous['stat'], 'cpu_ticks') if previous else None
    cores = ticks / ticks_per_sec / elapsed if ticks is not None and elapsed > 0 else None
    allowed = cpu_list_count(current['allowed_cpus'])
    threads = [thread_interval(row, previous['threads'].get(tid) if previous else None, elapsed, ticks_per_sec)
               for tid,row in current['threads'].items()]
    return dict(type='sample', wall_time_ns=time.time_ns(), monotonic_s=current['monotonic_s'], elapsed_s=elapsed,
                pid=current['stat']['pid'], process_start_ticks=current['stat']['start_ticks'],
                cpu_core_equivalents=cores, cpu_percent_one_core=cores*100 if cores is not None else None,
                cpu_percent_allowed_capacity=cores/allowed*100 if cores is not None and allowed else None,
                allowed_cpus=current['allowed_cpus'], allowed_cpu_count=allowed,
                thread_count=current['stat']['threads'], sampled_thread_count=len(threads),
                active_thread_count=sum((r['cpu_core_equivalents'] or 0)>0 for r in threads),
                runnable_thread_snapshot_count=sum(r['state']=='R' for r in threads),
                missing_thread_samples=current['missing_threads'], memory=current['memory'],
                pss_sampled=current['pss_sampled'], pss_error=current['pss_error'],
                io_delta={k:counter_delta(current['io'],previous['io'],k) if previous else None for k in current['io']},
                host_cpus=host_intervals(current['host_cpus'],previous['host_cpus'] if previous else {}),threads=threads,
                system_process_cpu=system_process_intervals(current['system_processes'],
                    previous['system_processes'] if previous else {},elapsed,ticks_per_sec),
                system_process_read_misses=current['system_missing'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--process-group',type=int,required=True)
    parser.add_argument('--executable',default='run_loc_online')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--interval-sec',type=float,default=1.0)
    parser.add_argument('--pss-interval-sec',type=float,default=30.0)
    parser.add_argument('--startup-timeout-sec',type=float,default=30.0)
    args=parser.parse_args()
    if (not 0.5<=args.interval_sec<=60 or not 5<=args.pss_interval_sec<=600 or
            not 0<args.startup_timeout_sec<=600 or args.process_group<=0):
        parser.error('positive process group, interval in [0.5,60], PSS interval in [5,600] required')
    stop=threading.Event()
    for sig in (signal.SIGINT,signal.SIGTERM):signal.signal(sig,lambda *_:stop.set())
    proc=Path('/proc');ticks=os.sysconf('SC_CLK_TCK');target=None
    deadline=time.monotonic()+args.startup_timeout_sec
    # Exclusive file creation avoids overwriting a previous diagnostic run.
    with args.output.open('x',encoding='utf-8') as output:
        def emit(row):
            output.write(json.dumps(row,ensure_ascii=False,separators=(',',':'),allow_nan=False)+'\n');output.flush()
        while not stop.is_set() and time.monotonic()<deadline:
            target=discover_target(proc,args.process_group,args.executable)
            if target:break
            stop.wait(.2)
        if target is None:
            emit(dict(type='error',reason='target_not_found',process_group=args.process_group))
            return 2
        path=proc/str(target['pid']);exe=os.readlink(path/'exe')
        with (path/'exe').open('rb') as source:
            digest=hashlib.file_digest(source,'sha256').hexdigest() if hasattr(hashlib,'file_digest') else None
            if digest is None:
                sha=hashlib.sha256()
                for chunk in iter(lambda:source.read(1024*1024),b''):sha.update(chunk)
                digest=sha.hexdigest()
        try:
            env=dict(item.split('=',1) for item in (path/'environ').read_bytes().decode(errors='replace').split('\0') if '=' in item)
            selected_env={k:env.get(k) for k in ENV_KEYS}
        except OSError: selected_env=None
        try: schedstats=(proc/'sys/kernel/sched_schedstats').read_text().strip()
        except OSError:schedstats=None
        emit(dict(type='metadata',schema_version=2,pid=target['pid'],process_start_ticks=target['start_ticks'],
                  process_group=args.process_group,executable=exe,executable_sha256=digest,
                  ticks_per_sec=ticks,interval_sec=args.interval_sec,pss_interval_sec=args.pss_interval_sec,
                  kernel_sched_schedstats=schedstats,selected_environment=selected_env,
                  system_process_scope='top20 by paired CPU delta; comm only, no cmdline; new/exited processes not attributed; wchan is a read-time snapshot',
                  semantics='CPU equivalents=process CPU delta/wall delta; active threads are not occupied cores; last_cpu is a snapshot; sched_wait=runqueue wait, not mutex/sleep or pipeline queue wait; zero schedstat may mean disabled'))
        previous=None;next_pss=0;samples=0
        while not stop.is_set():
            before=time.monotonic();cpu_before=time.process_time()
            try:current=process_snapshot(proc,target['pid'],before>=next_pss)
            except (OSError,ValueError,IndexError) as exc:
                emit(dict(type='end',reason='target_unavailable',error=type(exc).__name__,samples=samples));return 0
            if current['stat']['start_ticks']!=target['start_ticks']:
                emit(dict(type='end',reason='target_pid_reused',samples=samples));return 0
            if current['pss_sampled']:next_pss=before+args.pss_interval_sec
            row=sample_interval(current,previous,ticks)
            row['sampler_read_wall_ms']=(time.monotonic()-before)*1000
            row['sampler_read_cpu_ms']=(time.process_time()-cpu_before)*1000
            row['sampler_process_cpu_s']=time.process_time()
            emit(row);previous=current;samples+=1
            stop.wait(max(0,args.interval_sec-(time.monotonic()-before)))
        emit(dict(type='end',reason='stopped',samples=samples))
    return 0


if __name__=='__main__':
    try:sys.exit(main())
    except (OSError,RuntimeError) as exc:
        print(f'resource monitor failed: {exc}',file=sys.stderr);sys.exit(2)
