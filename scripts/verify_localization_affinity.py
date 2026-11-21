#!/usr/bin/env python3
"""Verify actual Linux thread masks and boot isolation; never changes affinity."""
import argparse
import json
from pathlib import Path


def cpus(text):
    result = set()
    for item in text.strip().split(","):
        if not item:
            continue
        ends = item.split("-")
        first, last = int(ends[0]), int(ends[-1])
        result.update(range(first, last + 1))
    return result


def verify(pid, require_isolation=False, *, proc_root=Path('/proc'),
           isolated_path=Path('/sys/devices/system/cpu/isolated')):
    expected = {"loc_sensor": {10}, "loc_output": {11}, "loc_executor": {11}}
    background = set(range(10))
    errors = []
    threads = []
    found = set()
    process = proc_root / str(pid)
    match_cpus = background
    single_match = False
    try:
        environment = dict(item.split('=', 1) for item in
                           (process / 'environ').read_bytes().decode(errors='surrogateescape').split('\0')
                           if '=' in item)
        if 'LIGHTNING_LM_MATCH_CPUS' in environment:
            match_cpus = cpus(environment['LIGHTNING_LM_MATCH_CPUS'])
            if not match_cpus or 10 in match_cpus:
                errors.append('匹配 CPU 配置为空或与传感器 CPU10 重叠')
        single_match = environment.get('LIGHTNING_LM_NDT_THREADS') == '1'
    except (OSError, ValueError) as exc:
        errors.append(f'无法核对目标进程的匹配 CPU 配置: {exc}')
    expected['loc_match'] = match_cpus
    for task in sorted((process / "task").iterdir()):
        try:
            status = dict(line.split(":", 1) for line in (task / "status").read_text().splitlines())
        except FileNotFoundError:
            continue
        name = status["Name"].strip()
        allowed = cpus(status["Cpus_allowed_list"])
        desired = expected.get(name, background)
        valid = allowed == desired if name in expected else bool(allowed) and allowed <= background
        if not valid:
            errors.append(f"TID {task.name} {name}: unexpected mask {sorted(allowed)}")
        if name in expected:
            if name in found and (name != 'loc_match' or single_match):
                errors.append(f"duplicate critical role: {name}")
            found.add(name)
        threads.append(dict(tid=int(task.name), name=name, cpus=sorted(allowed), valid=valid))
    for name in expected.keys() - found:
        errors.append(f"missing critical role: {name}")
    isolated = cpus(isolated_path.read_text())
    if require_isolation and not {10, 11} <= isolated:
        errors.append("boot scheduler isolation for CPUs 10-11 is absent")
    return dict(pid=pid, isolated_cpus=sorted(isolated), match_cpus=sorted(match_cpus),
                single_match_requested=single_match, threads=threads, errors=errors,
                passed=not errors,
                limitation="Affinity/isolation readback does not establish deadlines; collect loaded timing and IRQ deltas.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pid", type=int)
    parser.add_argument("--require-isolation", action="store_true")
    args = parser.parse_args()
    result = verify(args.pid, args.require_isolation)
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["passed"] else 1)
