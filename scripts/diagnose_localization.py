#!/usr/bin/env python3
"""Field diagnosis with optional cached figures; no ROS node or localization changes."""
import argparse
from collections import Counter, deque
import csv
import datetime as dt
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time

MIB = 1024 * 1024
PROC = Path('/proc')
ALIASES = dict(raw='raw_runstate', published='published_runstate',
               selected='selected_runstate', strict='strict_runstate',
               consistency='localization_consistency', uncertainty='localization_uncertainty',
               unobservable='pose_unobservable', yaw_sigma='yaw_sigma_rad')
LIMITS = ('major_sigma_limit_m', 'lateral_sigma_limit_m', 'yaw_sigma_limit_deg',
          'longitudinal_tolerance_m', 'lateral_tolerance_m', 'yaw_tolerance_deg',
          'distance_tolerance_ratio', 'angle_tolerance_ratio', 'max_input_gap_sec',
          'consistency_unknown_grace_sec', 'output_timeout_sec', 'startup_grace_sec',
          'motion_window_sec', 'evaluation_timeout_sec', 'config_version')
FAULT_CODES = (2, 4, -9)
CONTRACT_PATH = Path(__file__).with_name('diagnose_contract.json')
LAUNCH_LOGS = ('map_manifest.log', 'input_check.log', 'config_resolution.log',
               'launcher.log', 'launch.log', 'launch_compute_config.log')


def color_enabled(mode='auto', stream=None, env=None):
    stream = sys.stdout if stream is None else stream
    env = os.environ if env is None else env
    if mode != 'auto':
        return mode == 'always'
    return stream.isatty() and 'NO_COLOR' not in env and env.get('TERM') != 'dumb'


def row_ns(row):
    # Keep integer nanoseconds, especially for coverage boundaries.
    for key in ('_wall_ns', 'wall_ns'):
        try:
            return int(row[key])
        except (KeyError, TypeError, ValueError, OverflowError):
            pass
    return None


def read_text(path, limit=MIB):
    with Path(path).open('rb') as f:
        return f.read(limit).decode('utf-8', errors='replace')


def key_values(path):
    try:
        return dict(line.split('=', 1) for line in read_text(path).splitlines() if '=' in line)
    except OSError:
        return {}


def number(value):
    try:
        n = float(value)
        return n if math.isfinite(n) else None
    except (ValueError, TypeError):
        return None


def fields(text):
    result = {}
    for key, value in re.findall(r'([A-Za-z_]\w*)=([^\s,]+)', text):
        n = number(value)
        result[key] = n if n is not None else value
    return result


def cpu_set(text):
    values = set()
    for item in text.strip().split(','):
        if item:
            ends = [int(x) for x in item.split('-')]
            values.update(range(ends[0], ends[-1] + 1))
    return values


def low_impact():
    result = {}
    try:
        os.nice(max(0, 10 - os.getpriority(os.PRIO_PROCESS, 0)))
        result['nice'] = os.getpriority(os.PRIO_PROCESS, 0)
        isolated = cpu_set(read_text('/sys/devices/system/cpu/isolated'))
        allowed = os.sched_getaffinity(0)
        shared = allowed - isolated
        if not shared:
            raise RuntimeError('诊断脚本仅允许使用隔离核；请从普通终端运行，以免干扰定位。')
        os.sched_setaffinity(0, shared)
        result.update(cpus=sorted(shared), isolated=sorted(isolated))
    except (AttributeError, OSError) as exc:
        result['note'] = '调度限制未完全应用: ' + str(exc)
    return result


def argument(args, name):
    for i, arg in enumerate(args):
        if arg.startswith(name + '='):
            return arg.split('=', 1)[1]
        if arg == name and i + 1 < len(args):
            return args[i + 1]
    return ''


def find_processes(proc=PROC):
    """Linux comm is loc_executor after pthread_setname_np; never match comm."""
    found = []
    if not proc.is_dir():
        return found
    try:
        boot = int(re.search(r'^btime (\d+)$', read_text(proc / 'stat'), re.M)[1])
    except (OSError, TypeError, ValueError):
        boot = 0
    for entry in proc.iterdir():
        if not entry.name.isdigit():
            continue
        try:
            exe = os.readlink(entry / 'exe').removesuffix(' (deleted)')
            if Path(exe).name != 'run_loc_online':
                continue  # excludes ros2 run, taskset, shells and this analyzer
            args = (entry / 'cmdline').read_bytes().decode(errors='replace').split('\0')
            cwd = Path(os.readlink(entry / 'cwd'))
            def resolved(flag):
                value = argument(args, flag)
                return str((cwd / value).resolve()) if value else ''
            config, tum = resolved('--config'), resolved('--output_tum')
            archive = resolved('--incident_dir') or (str(Path(tum).parent / 'incidents') if tum else '')
            log = os.readlink(entry / 'fd/2')
            log = log if Path(log).is_file() else ''
            run = str(Path(config).parent) if config else (str(Path(log).parent.parent) if log else '')
            # Start ticks distinguish restarted processes when selecting an instance.
            start_ticks = int(read_text(entry / 'stat').rsplit(')', 1)[1].split()[19])
            found.append(dict(pid=int(entry.name), comm=read_text(entry / 'comm').strip(),
                              exe=exe, config=config, log=log, run=run, archive=archive,
                              match_failure_root=resolved('--match_failure_dir'),
                              start_ticks=start_ticks,
                              start_wall_ns=int((boot + start_ticks / os.sysconf('SC_CLK_TCK')) * 1e9) if boot else 0))
        except (OSError, ValueError, IndexError):
            continue  # process exited or inaccessible; no root required
    return sorted(found, key=lambda p: p['start_ticks'], reverse=True)


def workspace_config(env=None, home=None):
    env = dict(os.environ if env is None else env)
    home = Path.home() if home is None else Path(home)
    # Parse only literal exports/variable references, never execute a config file.
    config = dict(HOME=str(home))
    try:
        for line in read_text(home / '.config/sany/workspace.bash').splitlines():
            match = re.match(r'^export\s+(\w+)=(.+?)\s*$', line)
            if match:
                value = match[2].strip('"\'')
                value = re.sub(r'\$\{(\w+)\}|\$(\w+)',
                               lambda m: config.get(m[1] or m[2], m[0]), value)
                if '$' not in value and '`' not in value:
                    config[match[1]] = value
    except OSError:
        pass
    config.update(env)
    return config


def default_roots(env=None, home=None):
    home = Path.home() if home is None else Path(home)
    config = workspace_config(env, home)
    roots = [config.get('LIGHTNING_LM_OUT_ROOT')]
    for key in ('SANY_WS', 'LIGHTNING_LM_WS'):
        if config.get(key):
            roots.append(str(Path(config[key]) / 'runs'))
    roots += [str(home / 'Sany/runs'), '/home/nvidia/project/gj_ws/runs']
    return list(dict.fromkeys(Path(p).expanduser().resolve() for p in roots if p))


def find_driver_diagnostic(env=None, home=None, script=None):
    """Check known install/source locations only; no ROS discovery or recursive scan."""
    config = workspace_config(env, home)
    script = Path(__file__ if script is None else script).resolve()
    prefixes = [Path(p) for p in config.get('AMENT_PREFIX_PATH', '').split(os.pathsep) if p]
    if config.get('LIGHTNING_LM_INSTALL_SETUP'):
        setup = Path(config['LIGHTNING_LM_INSTALL_SETUP']).expanduser()
        prefixes += [setup.parent, setup.parent.parent]
    for key in ('SANY_WS', 'LIGHTNING_LM_WS'):
        if config.get(key):
            workspace = Path(config[key]).expanduser()
            prefixes += [workspace / 'lightning_lm_ws/install', workspace / 'install']
    prefixes += [script.parents[2], script.parents[2].parent, script.parents[1] / 'install']
    relative = Path('lib/livox_ros_driver2_core/diagnose_lidar.py')
    candidates = [p / suffix for p in prefixes for suffix in (relative, Path('livox_ros_driver2_core') / relative)]
    candidates.append(script.parents[2] / 'sdk/livox_sdk/src/livox_ros_driver2/scripts/diagnose_lidar.py')
    return next((p.resolve() for p in candidates if p.is_file()), None)


def driver_diagnosis(root=None):
    """Delegate driver rules; exit 0/has_evidence means evidence, never health."""
    script = find_driver_diagnostic()
    if script is None:
        return dict(status='unavailable', message='未找到独立雷达诊断脚本；上游状态未检查。')
    result = dict(script=str(script))
    command = [sys.executable, '-B', str(script), '--json', '--no-refresh']
    if root:
        command += ['--root', str(root)]
    try:
        def reject_nonfinite(value):
            raise ValueError('雷达报告包含非有限数值: ' + value)
        completed = subprocess.run(command, capture_output=True, text=True, encoding='utf-8',
                                   errors='replace', timeout=10, check=False)
        data = json.loads(completed.stdout, parse_constant=reject_nonfinite)
        if not isinstance(data, dict) or data.get('schema') != 1 or data.get('kind') != 'livox_driver_diagnosis':
            return dict(result, status='unsupported', message='雷达诊断接口版本不兼容；请更新独立雷达诊断脚本。')
        if data.get('error'):
            return dict(result, status='unavailable', message='雷达诊断未取得证据: ' + str(data['error']))
        if (completed.returncode not in (0, 2) or not isinstance(data.get('active'), bool)
                or not isinstance(data.get('has_evidence'), bool)
                or not isinstance(data.get('findings'), list) or not isinstance(data.get('warnings'), list)
                or not all(isinstance(f, dict) for f in data['findings'])
                or not isinstance(data.get('roles'), list)
                or not all(isinstance(r, dict) and isinstance(r.get('device'), dict) for r in data['roles'])
                or (completed.returncode == 2 and data['has_evidence'])):
            return dict(result, status='failed', message='雷达诊断返回了不完整或不一致的结果；上游状态未知。')
        return dict(result, status='available' if data['has_evidence'] else 'no_evidence', report=data,
                    message='独立雷达诊断仅读取本机已有证据（--no-refresh），未刷新设备状态；与定位包未建立同次采集关联。')
    except (OSError, ValueError, subprocess.TimeoutExpired) as exc:
        return dict(result, status='failed', message='雷达诊断调用失败或接口不可读，上游状态未知: ' + str(exc))


def render_driver_diagnosis(upstream):
    lines = ['上游雷达独立诊断: ' + upstream['message']]
    data = upstream.get('report')
    if not data:
        return lines
    scope = '当前驱动运行实例（设备证据未刷新）' if data['active'] else '历史驱动包（当前未发现雷达 worker）'
    lines += ['  范围: ' + scope, '  雷达数据包: ' + str(data.get('selected'))]
    for role in data.get('roles', []):
        device = role.get('device', {})
        lines.append(f'  {role.get("role")}: 设备证据={device.get("evidence_scope", "未知")}，'
                     f'采集={device.get("collected_at")}，距采集={device.get("age_sec")} 秒')
    findings = data['findings']
    history = [f for f in findings if f.get('stage') == 'device_history']
    other = [f for f in findings if f.get('stage') != 'device_history']
    # Retain each device's history/coverage summary even when lifecycle and
    # startup failures fill the ordinary ten-item summary.
    for finding in other[:10] + history:
        lines.append(f'  [{finding.get("severity", "未知")}] {finding.get("role") or "生命周期"}: {finding.get("message", "未知判据")}')
        evidence = finding.get('evidence', [])
        if finding.get('stage') == 'device_history':
            locations = []
            for item in evidence:
                if isinstance(item, dict) and item.get('path'):
                    position = (f':line={item["line"]}' if item.get('line') is not None
                                else f':byte={item["byte_offset"]}' if item.get('byte_offset') is not None else '')
                    locations.append(str(item['path']) + position + (' @ ' + str(item['time']) if item.get('time') else ''))
            lines.append('    证据: ' + ' | '.join(dict.fromkeys(locations)))
        else:
            lines.append('    证据: ' + json.dumps(evidence, ensure_ascii=False))
    if len(other) > 10:
        lines.append(f'  共 {len(findings)} 条雷达记录；普通判据展示前 10 条，设备历史摘要单独保留；--json 保留完整独立报告。')
    elif not findings:
        lines.append('  已读取证据中没有列出雷达异常；不能据此确认当前设备或数据流正常。')
    lines += ['  雷达提示: ' + str(w) for w in data['warnings']]
    return lines


def summary_at(path):
    try:
        value = json.loads(read_text(path / 'summary.json'))
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError):
        return {}


def stamp(path):
    summary = summary_at(path)
    if number(summary.get('trigger_wall_ns')) is not None:
        return int(summary['trigger_wall_ns'])
    meta = key_values(path / 'run_metadata.txt')
    try:
        return int(dt.datetime.fromisoformat(meta['started_at']).timestamp() * 1e9)
    except (KeyError, ValueError):
        # Before preflight succeeds, the launcher has not written metadata yet.
        # Its standard directory name records launch time, unlike copy/mtime.
        match = re.fullmatch(r'sany_lidar_(\d{8}_\d{6})', path.name)
        if match:
            try:
                return int(dt.datetime.strptime(match[1], '%Y%m%d_%H%M%S').astimezone().timestamp() * 1e9)
            except ValueError:
                pass
        return path.stat().st_mtime_ns


def is_run_directory(path):
    return path.is_dir() and ((path / 'logs/run_loc_online.stderr.log').is_file() or
        any((path / 'logs' / name).is_file() for name in LAUNCH_LOGS) or
        (re.fullmatch(r'sany_lidar_\d{8}_\d{6}', path.name) is not None and (path / 'logs').is_dir()))


def launch_failures(run, warnings):
    """Explain recorded preflight errors without assuming every tool error is a bad map."""
    findings = []
    for name in LAUNCH_LOGS:
        path = run / 'logs' / name
        if not path.is_file():
            continue
        try:
            with path.open('rb') as stream:
                size = path.stat().st_size
                stream.seek(max(0, size - 256 * 1024))
                if size > 256 * 1024:
                    stream.readline()
                    warnings.append(f'启动日志 {name} 仅读取末尾 256 KiB。')
                while stream.tell() < size:
                    offset = stream.tell()
                    line = stream.readline(min(65536, size - offset)).decode(errors='replace').strip()
                    key = title = advice = None
                    layer = 'map_preflight' if name == 'map_manifest.log' else 'input_preflight' if name == 'input_check.log' else 'launcher'
                    missing = re.search(r'error while loading shared libraries:\s*([^:\s]+):', line)
                    if missing:
                        key, layer = 'shared_library_missing', 'runtime_loader'
                        title = f'启动工具动态库加载失败: {missing[1]}'
                        advice = '检查安装空间、动态库搜索路径和工具链接配置；此错误发生在工具执行前，不能据此判断地图损坏或传感器故障。'
                    elif 'Map manifest rejected:' in line:
                        key, title = 'map_preflight_rejected', '地图预检工具校验未通过'
                        advice = '按原始拒绝原因检查地图清单或运行配置；工具拒绝不等于地图文件损坏。'
                    elif 'ERROR: sensor delivery check timed out' in line:
                        key, title = 'input_check_timeout', '定位输入检查超时，算法尚未通过启动输入门槛'
                        advice = '核对 input_check.log 的各 topic 接收计数及本次 DDS 环境；接收缺失不能单独区分设备停流与订阅环境问题。'
                    elif 'ERROR: sensor delivery check:' in line:
                        key, title = 'input_check_error', '定位输入检查执行失败'
                        advice = '按原始错误检查检查器参数或运行环境，不将执行失败直接解释为传感器故障。'
                    elif 'SANY vehicle configuration rejected:' in line or 'vehicle configuration resolution failed' in line:
                        key, title = 'configuration_rejected', '车辆/定位配置解析未通过'
                        advice = '按原始错误检查所选车辆配置与必填字段；缺少详细原因时不猜测具体配置项。'
                    elif '[ros2run]: Process exited with failure' in line or 'map manifest validation failed;' in line:
                        key, title = 'preflight_tool_failed', '启动预检工具失败，现有记录未确定具体原因'
                        advice = '核对同一日志中更早的错误和工具运行环境；不能仅凭日志名称认定地图有问题。'
                    if key:
                        findings.append(dict(key=key, stage=layer, title=title, advice=advice,
                                             source=f'{path}:byte={offset}', evidence=line))
                        break  # The first explicit error is more specific than a trailing ros2run status.
        except OSError as error:
            warnings.append(f'启动日志读取失败 {path}: {error}')
    return findings


def same_run(event, run, pid=None):
    meta = key_values(event / 'run_metadata.txt')
    if meta.get('run_dir'):
        return Path(meta['run_dir']).resolve() == Path(run).resolve()
    # Older/manual captures lack metadata. Require both PID and config identity.
    if pid and event.name.endswith('_' + str(pid)):
        try:
            return (event / 'config.yaml').read_bytes() == (Path(run) / 'config.yaml').read_bytes()
        except OSError:
            pass
    return False


def events_in(root):
    if not root.is_dir():
        return []
    return [p for p in root.glob('event_*') if p.is_dir() and
            ((p / 'summary.json').is_file() or (p / 'timeline.csv').is_file())]


def choose_path(path, run=None, pid=None, start_wall_ns=0, preserve_run=False):
    path = Path(path).expanduser().resolve()
    if not path.exists():
        raise ValueError('路径不存在: ' + str(path))
    if path.is_file() and path.name in ('summary.json', 'summary.html', 'timeline.csv', 'metadata.yaml'):
        path = path.parent
    if path.name.startswith('segment_') and (path.parent / 'summary.json').is_file():
        path = path.parent
    if path.is_file() or (path / 'summary.json').is_file() or (path / 'timeline.csv').is_file():
        return path
    candidates = events_in(path) + events_in(path / 'incidents')
    if is_run_directory(path):
        if preserve_run:
            return path
        candidates += [e for e in events_in(path.parent / 'incidents') if same_run(e, path, pid)]
        candidates += events_in(path / 'results/incidents')
        candidates = [e for e in candidates if same_run(e, path, pid)]
        return max(candidates, key=stamp) if candidates else path
    if run:
        candidates = [p for p in candidates if same_run(p, run, pid) and stamp(p) >= start_wall_ns]
    else:
        runs = {p.parent for p in path.glob('*/run_metadata.txt') if is_run_directory(p.parent)}
        runs.update(p.parent.parent for name in LAUNCH_LOGS for p in path.glob('*/logs/' + name))
        runs.update(p for p in path.glob('sany_lidar_*') if is_run_directory(p))
        candidates += list(runs)
    if candidates:
        return max(candidates, key=stamp)
    if (path / 'metadata.yaml').is_file() or list(path.glob('*.mcap')) or list(path.glob('*.db3')):
        return path
    raise ValueError('没有找到事件包或定位运行包: ' + str(path))


def select(args, processes=None, roots=None):
    if args.path:
        return choose_path(args.path, preserve_run=True), None, '指定路径', []
    processes = find_processes() if processes is None else processes
    if args.pid:
        processes = [p for p in processes if p['pid'] == args.pid]
        if not processes:
            raise ValueError('指定 PID 不是可读取的定位实例；通过 exe 识别，不依赖线程名。')
    if processes:
        p = processes[0]
        warnings = [f'发现 {len(processes)} 个定位实例，选择启动最晚的 PID {p["pid"]}；可用 --pid 指定。'] if len(processes) > 1 else []
        try:
            chosen = choose_path(p['archive'], p['run'], p['pid'], p.get('start_wall_ns', 0)) if p['archive'] else None
        except ValueError:
            chosen = None
        if chosen is None:
            chosen = Path(p['run'])
            warnings.append('当前定位实例尚无匹配事件包，分析本次运行日志；未使用历史实例的事件。')
        return chosen, p, '运行中的定位实例', warnings
    roots = [Path(args.root).expanduser().resolve()] if args.root else (default_roots() if roots is None else roots)
    for root in roots:
        try:
            return choose_path(root), None, '无定位进程，默认目录最新包', []
        except ValueError:
            pass
    raise ValueError('未找到定位进程或数据包。检查过: ' + ', '.join(map(str, roots)) + '；也可直接传入包路径。')


def thresholds(config):
    """Only diagnostics scalar values are needed; no PyYAML/ROS dependency."""
    result = {}
    try:
        active = False
        for line in read_text(config).splitlines():
            if line.startswith('diagnostics:'):
                active = True
            elif line and not line[0].isspace() and not line.startswith('#'):
                active = False
            if active:
                match = re.match(r'^  (\w+):\s*([^#]+)', line)
                if match and match[1] in LIMITS:
                    value = match[2].strip().strip('"\'')
                    result[match[1]] = number(value) if number(value) is not None else value
    except OSError:
        pass
    return result


def matching_thresholds(config):
    """Read explicit scalar paths only; do not invent missing YAML defaults."""
    scalars, parents = {}, []
    try:
        for line in read_text(config).splitlines():
            match = re.match(r'^( *)(\w+):\s*([^#]*)', line)
            if not match:
                continue
            indent, key, value = len(match[1]), match[2], match[3].strip().strip('"\'')
            while parents and parents[-1][0] >= indent:
                parents.pop()
            path = tuple(p[1] for p in parents) + (key,)
            if not value:
                parents.append((indent, key))
            else:
                scalars[path] = value
    except OSError:
        return {}
    backend = scalars.get(('relocalization', 'backend'))
    overlap_path = ('relocalization', backend, 'map_min_overlap_ratio')
    if overlap_path not in scalars:
        overlap_path = ('relocalization', 'map_min_overlap_ratio')
    return dict(source=str(config), backend=backend,
                map_min_overlap_ratio=number(scalars.get(overlap_path)) if backend else None,
                overlap_key='.'.join(overlap_path),
                min_tracking_confidence=number(scalars.get(('lidar_loc', 'min_tracking_confidence'),
                                                           scalars.get(('lidar_loc', 'min_init_confidence')))))


def record_runtime_protection(protections, key, line, row):
    """Aggregate two protection markers without retaining an unbounded row list."""
    if key == 'stale_lidar_correction':
        title = '过期 LIO 校正被丢弃'
        evidence = {k: row[k] for k in ('age', 'hard_deadline') if k in row}
    else:
        title = '自适应雷达负载级别变化'
        evidence = {k: row[k] for k in ('processing_ms', 'latency_ms', 'tracking_healthy') if k in row}
        transition = re.search(r'adaptive lidar load step (\d+)\s*->\s*(\d+)', line)
        if transition:
            evidence.update(from_step=int(transition[1]), to_step=int(transition[2]))
    entry = protections.setdefault(key, dict(key=key, title=title, count=0,
        first_time=row['_time'], first_ns=row_ns(row), first_source=row['_source']))
    entry['count'] += 1
    entry.update(last_time=row['_time'], last_ns=row_ns(row), evidence=evidence,
                 text=line.strip(), source=row['_source'])


def is_boundary_context(text):
    return any(marker in text for marker in (
        'localization failed!', 'localization pose and map outputs disabled',
        'localization pose and map outputs recovered', 'initial pose is not configured',
        'sensor queue lag reached')) or (
        'GLOBAL_RELOCALIZATION[' in text and '] accepted' in text)


def log_rows(path, budget, warnings, reference=None, context=None,
             protections=None, protection_window=None):
    """Bounded tail, original byte offsets as evidence; does not slurp the log."""
    result, config = [], {}
    if not path or not Path(path).is_file():
        return result, config
    recent_context = deque(context, maxlen=1000) if context is not None else None
    boundary_context = [c for c in context or [] if is_boundary_context(c['text'])]
    boundary_truncated = len(boundary_context) > 128
    boundary_context = boundary_context[:128]
    context_truncated = False
    path = Path(path)
    try:
        with path.open('rb') as f:
            # Runtime config is near startup; read at most 1 MiB independently.
            head = f.read(MIB).decode(errors='replace')
            for line in head.splitlines():
                if 'LOCALIZATION_DTC_CONFIG' in line:
                    config.update(fields(line))
                if 'LOCALIZATION_DIAGNOSTIC_BUILD' in line:
                    config['_build'] = fields(line)
            size = os.fstat(f.fileno()).st_size
            start = max(0, size - budget)
            f.seek(start)
            if start:
                f.readline()
                warnings.append(f'日志仅检查末尾 {budget // MIB} MiB；更早的故障未覆盖，可增大 --log-mb。')
            while f.tell() < size:
                offset = f.tell()
                line = f.readline(min(128 * 1024, size - offset)).decode(errors='replace')
                modern = 'SAFETY_DIAGNOSTIC' in line
                protection_key = None
                if protections is not None:
                    if 'drop stale lidar correction before matching:' in line:
                        protection_key = 'stale_lidar_correction'
                    elif 'adaptive lidar load step ' in line:
                        protection_key = 'adaptive_lidar_load_step'
                contextual = context is not None and any(marker in line for marker in (
                    'initial pose is not configured', 'GLOBAL_RELOCALIZATION[',
                    'SOLID_KISS_REJECT', 'failed to load SOLiD ICP cloud',
                    'failed to build SOLiD+KISS target', 'failed to prepare SOLiD+KISS query',
                    'ignore CAN wheel speed:', 'propagate over long imu interval:',
                    'sensor queue lag ', 'drop stale lidar correction before matching:',
                    'localization failed!', 'localization pose and map outputs disabled',
                    'localization pose and map outputs stopped',
                    'localization pose and map outputs recovered'))
                if not modern and 'SAFETY_TIME_SYNC_SHADOW' not in line and not contextual and not protection_key:
                    continue
                row = fields(line)
                if not modern:
                    row['raw_runstate'] = row.get('published_runstate')
                match = re.match(r'\w(\d{4})\s+([\d:.]+)', line)
                row['_time'] = f'{match[1][:2]}-{match[1][2:]} {match[2]}' if match else '时间未知'
                if match and reference:
                    try:
                        dates = [dt.datetime.strptime(f'{y}{match[1]} {match[2]}', '%Y%m%d %H:%M:%S.%f').replace(tzinfo=reference.tzinfo)
                                 for y in (reference.year - 1, reference.year, reference.year + 1)]
                        closest = min(dates, key=lambda d: abs((d - reference).total_seconds()))
                        row['_wall_ns'] = int(closest.timestamp() * 1e9)
                    except ValueError:
                        pass
                row['_source'] = f'{path}:byte={offset}'
                row['_scope'] = '运行日志（有界尾部）' if modern else '旧日志上报码（原始码未知）'
                if protection_key and not modern:
                    if protection_window is not None:
                        stamp = row_ns(row)
                        if stamp is None:
                            warning = '保护日志缺少可解析时间戳，未计入指定事件窗口。'
                            if warning not in warnings:
                                warnings.append(warning)
                            continue
                        if not protection_window or not protection_window[0] <= stamp <= protection_window[1]:
                            continue
                    record_runtime_protection(protections, protection_key, line, row)
                    if not contextual:
                        continue
                if contextual and not modern:
                    # Retain the latest bounded set of actual log clues, never infer a
                    # hardware cause solely from a matching keyword.
                    if len(recent_context) == recent_context.maxlen:
                        context_truncated = True
                    clue = dict(time=row['_time'], wall_ns=row.get('_wall_ns'),
                                source=row['_source'], byte_offset=offset, text=line.strip()[:1200])
                    recent_context.append(clue)
                    if is_boundary_context(clue['text']):
                        if len(boundary_context) < 128:
                            boundary_context.append(clue)
                        else:
                            boundary_truncated = True
                    continue
                if number(row.get('raw_runstate')) is None:
                    warning = '忽略缺少有效故障码的日志记录；可能存在格式变化或未写完的行。'
                    if warning not in warnings:
                        warnings.append(warning)
                    continue
                result.append(row)
    except OSError as exc:
        warnings.append('日志读取失败: ' + str(exc))
    if context is not None:
        merged = {c['source']: c for c in boundary_context + list(recent_context)}
        context[:] = sorted(merged.values(), key=lambda c: (c.get('wall_ns') or 0, c.get('byte_offset', 0)))
    if context_truncated:
        warnings.append('关联日志超过 1000 条，仅保留已检查日志范围内最新 1000 条，并额外保留最多 128 条最早关键边界线索；其他较早线索可能未覆盖。')
    if boundary_truncated:
        warnings.append('关键边界线索超过 128 条，额外缓存仅保留最早 128 条；其余仅在最近 1000 条关联日志内保留，部分故障起点可能未覆盖。')
    # Avoid double-counting shadow and unconditional diagnostics from the same run.
    if any(r['_scope'] == '运行日志（有界尾部）' for r in result):
        result = [r for r in result if r['_scope'] == '运行日志（有界尾部）']
    gaps = sum(1 for a, b in zip(result, result[1:])
               if number(a.get('sequence')) is not None and number(b.get('sequence')) is not None
               and b['sequence'] != a['sequence'] + 1)
    if gaps:
        warnings.append(f'Safety 日志序号有 {gaps} 处不连续；缺口期间的短暂故障和恢复过程无法保证覆盖。')
    return result, config


def timeline_rows(path, warnings):
    result = []
    if not path.is_file():
        return result
    try:
        with path.open(newline='') as f:
            for line, row in enumerate(csv.DictReader(f), 2):
                if line > 30000:
                    warnings.append('时间线超过 30000 行，仅检查前 30000 行。')
                    break
                if number(row.get('raw')) is None or not row.get('wall_ns') or None in row or None in row.values():
                    warnings.append('忽略时间线尚未写完/格式损坏的一行。')
                    continue
                clean = {ALIASES.get(k, k): (number(v) if number(v) is not None else v) for k, v in row.items()}
                clean['_wall_ns'] = int(row['wall_ns'])
                clean['_time'] = wall_time(row.get('wall_ns'))
                clean['_source'] = f'{path}:{line}'
                clean['_scope'] = '事件时间线'
                result.append(clean)
    except (OSError, csv.Error, ValueError) as exc:
        warnings.append('时间线读取失败: ' + str(exc))
    return result


def wall_time(ns):
    try:
        return dt.datetime.fromtimestamp(int(ns) / 1e9).astimezone().isoformat(timespec='milliseconds')
    except (ValueError, TypeError, OSError, OverflowError):
        return '时间未知'


def recent_context_clues(clues, limit):
    selected = []
    for marker in ('SOLID_KISS_REJECT', 'GLOBAL_RELOCALIZATION['):
        clue = next((c for c in reversed(clues) if marker in c['text']), None)
        if clue is not None and clue not in selected:
            selected.append(clue)
    selected = selected[-limit:]
    selected += [c for c in reversed(clues) if c not in selected][:max(0, limit - len(selected))]
    return [c for c in clues if c in selected]


def classify(row, limits):
    """Mirror Safety decision precedence; never equate UNKNOWN with FAIL."""
    code = number(row.get('raw_runstate'))
    causes = []
    def add(key, title, evidence, advice, confidence='判据已确认'):
        values = {k: row[k] for k in evidence if k in row}
        causes.append(dict(code=int(code), key=key, title=title, evidence=values,
                           advice=advice, confidence=confidence, time=row.get('_time'),
                           source=row.get('_source'), scope=row.get('_scope'),
                           published=row.get('published_runstate')))
    if code == 4:
        uncertainty = number(row.get('localization_uncertainty'))
        if row.get('pose_unobservable') == 1:
            add('unobservable', '定位被判为不可观测', ['pose_unobservable', 'covariance_valid'], '检查匹配约束与重定位前后状态。')
        else:
            for metric, limit, factor in [('major_sigma', 'major_sigma_limit_m', 1),
                                           ('lateral_sigma', 'lateral_sigma_limit_m', 1),
                                           ('yaw_sigma_rad', 'yaw_sigma_limit_deg', 180 / math.pi)]:
                value, boundary = number(row.get(metric)), number(limits.get(limit))
                if row.get('covariance_valid') == 1 and value is not None and boundary is not None and value * factor > boundary:
                    add(metric, f'不确定度超限: {metric}={value * factor:.6g} > {boundary:g}' + ('°' if factor != 1 else ' m'),
                        [metric, 'covariance_valid', 'localization_uncertainty'], '检查协方差与实际误差证据，不能仅凭轨迹连续判误报。')
            if uncertainty == 0:
                add('uncertainty_unknown', '不确定度检查为 UNKNOWN（无有效结论）', ['localization_uncertainty', 'covariance_valid'], '检查协方差是否有效，以及启动/重定位状态。')
            elif uncertainty == 2 and not causes:
                add('uncertainty_fail', '不确定度检查失败，现有字段不足以还原超限项', ['localization_uncertainty', 'covariance_valid', 'major_sigma', 'lateral_sigma', 'yaw_sigma_rad'], '需要本次实际配置和更完整判据。')
        for key, title, evidence in [('yaw_condition', '航向运动一致性检查失败', ['yaw_condition', 'error_yaw_rad']),
                                     ('translation_condition', '平移运动一致性检查失败', ['translation_condition', 'error_x', 'error_y'])]:
            if row.get(key) == 2:
                add(key, title, evidence, '比较位姿增量与 IMU/CAN 积分；运动阈值包含增量比例项，不是不确定度阈值。')
                if key == 'yaw_condition':
                    causes[-1]['advice'] = '比较位姿航向增量与 IMU 机体系 z 角速度积分；CAN 用于平移检查。旧字段缺少累计转角 turn 及窗口/位姿分支，无法重建确切动态阈值；这不是协方差不确定度阈值。'
                    if number(row.get('error_yaw_rad')) is not None:
                        causes[-1]['evidence']['error_yaw_deg'] = row['error_yaw_rad'] * 180 / math.pi
        if row.get('localization_consistency') == 2 and row.get('yaw_condition') != 2 and row.get('translation_condition') != 2:
            add('consistency_fail', '运动一致性检查失败，缺少平移/航向细分状态', ['localization_consistency'], '需要同次运行的细分判据或完整事件。')
        latched = row.get('confirmed_failure') == 1 and row.get('localization_consistency') != 2 and uncertainty != 2
        if latched:
            add('latched_failure', '先前确认的失败尚未重新通过', ['confirmed_failure', 'localization_consistency', 'localization_uncertainty'], '向前检查首次 FAIL；UNKNOWN 不会清除已确认失败。')
        if not latched and row.get('localization_consistency') == 0 and row.get('unknown_tolerated') == 0 and row.get('optional_can_unavailable') == 0:
            add('consistency_unknown', '运动一致性证据不足，未满足 UNKNOWN 宽容条件',
                ['localization_consistency', 'yaw_condition', 'translation_condition', 'yaw_pass_age_sec',
                 'output_age_sec', 'imu_arrival_age_sec', 'imu_stamp_lag_sec', 'unknown_tolerated'],
                '检查 IMU 时间区间覆盖、时间戳跳变与输出间隔；当前 IMU 到达新鲜也不证明历史积分区间完整。')
            window = number(limits.get('motion_window_sec'))
            if window is not None:
                causes[-1]['evidence']['motion_window_sec'] = window
            if (str(row.get('yaw_pass_age_sec')) == 'inf' and uncertainty == 1 and
                    row.get('covariance_valid') == 1 and number(row.get('output_age_sec')) is not None):
                causes[-1]['advice'] = '尚无可用航向 PASS，核查是否处于运动评估窗口建立/重建阶段；也可能是 IMU 历史区间覆盖不足，不能仅凭当前新鲜判定根因。'
            grace = number(limits.get('consistency_unknown_grace_sec'))
            age = number(row.get('yaw_pass_age_sec'))
            if grace is not None and age is not None and age > grace:
                causes[-1]['title'] += f'；最近航向 PASS 已过去 {age:.6g}s > {grace:g}s'
            gap, imu_age = number(limits.get('max_input_gap_sec')), number(row.get('imu_arrival_age_sec'))
            if gap is not None and imu_age is not None and imu_age > gap:
                causes[-1]['title'] += f'；主 IMU 有效证据年龄 {imu_age:.6g}s > {gap:g}s'
        if str(row.get('output_age_sec')) == 'inf':
            for cause in causes:
                if cause['key'] in ('uncertainty_unknown', 'consistency_unknown'):
                    cause['phase'] = 'startup'
                    cause['title'] = '启动后尚无有效定位输出；' + cause['title']
                    cause['evidence']['output_age_sec'] = row['output_age_sec']
                    cause['advice'] = '先核对初始位姿、全局重定位候选拒绝原因及接受时间，再核对首个有效定位输出。'
    elif code == -9:
        for key, title in [('output_stale', '定位输出超时'), ('localization_lost', '定位输出门控关闭/定位丢失'), ('radar_all_failed', '所有必需雷达均不可用')]:
            if row.get(key) == 1:
                startup = key == 'output_stale' and str(row.get('output_age_sec')) == 'inf'
                if startup:
                    title = '启动后尚无有效定位输出，超过启动宽限期'
                add(key, title, [key, 'output_age_sec', 'imu_arrival_age_sec', 'can_arrival_age_sec',
                                'imu_stamp_lag_sec', 'can_stamp_lag_sec', 'match_stale', 'match_age_sec', 'lost_frames'],
                    '检查初始位姿及全局重定位接受时间。' if startup else '沿输入、处理队列、匹配、发布时间线检查最早异常。')
                causes[-1]['phase'] = 'startup' if startup else 'runtime'
                ages = [number(row.get(k)) for k in ('imu_arrival_age_sec', 'can_arrival_age_sec')]
                gap = number(limits.get('max_input_gap_sec'))
                if not startup and gap is not None and all(a is not None and a > gap for a in ages):
                    causes[-1]['advice'] = 'IMU 与 CAN 证据同时变旧；优先检查共享回调调度和输入链路，现有证据不能确定 DDS 或硬件根因。'
        if not causes and number(row.get('publish_end_age')) is not None and row['publish_end_age'] > 0.5:
            add('publish_gap', '同时观测到发布完成间隔增大', ['publish_end_age', 'publishing', 'imu_callback_age', 'lidar_callback_age', 'sensor_pending', 'localization_pending'], '需 Safety 详细判据确认 -9 分支；发布返回不等于 DDS 下游收到。', '相关线索')
    elif code == 2:
        mask = number(row.get('time_sync_fault_mask'))
        if mask:
            blocking = int(mask)
            if row.get('bounded_policy_enabled') == 1:
                blocking &= ~int(number(row.get('bounded_sync_warning_mask')) or 0)
            roles = [name for bit, name in enumerate(('front', 'rear', 'left', 'right')) if blocking & (1 << bit)]
            for role in roles:
                add('radar_time_sync', '雷达时间同步诊断未通过: ' + role,
                    ['time_sync_fault_mask', 'bounded_sync_warning_mask', 'bounded_policy_enabled', 'non_sync_fault', 'non_sync_all_pass'],
                    '查看对应雷达 diagnostic_state 的原生 DTC 和 PTP 证据；同步判据失败不等于偏差数值一定超限。')
                causes[-1]['device'] = role
        if row.get('non_sync_fault') == 1 or row.get('non_sync_all_pass') == 0:
            add('radar_other', '雷达非时间同步诊断失败或证据未齐', ['non_sync_fault', 'non_sync_all_pass'], '检查各必需雷达的状态与 HMS 证据。')
    if code in (2, 4, -9) and not causes:
        add('missing_details', '记录了故障码，但缺少判据字段，无法确定具体触发分支',
            ['raw_runstate', 'published_runstate', 'major_sigma', 'lateral_sigma', 'yaw_sigma_rad', 'covariance_valid'],
            '补充同次运行日志或升级后生成的完整事件；不能直接认定误报。', '证据不足')
    elif code is not None and code not in (0, 2, 4, -9):
        add('unknown_code', '发现未支持的 Safety 码，诊断规则需要同步', ['raw_runstate', 'published_runstate'],
            '检查运行库版本与故障码定义，并更新诊断规则；不能按正常状态处理。', '规则缺失')
    return causes


def bag_rows(path, byte_limit, warnings):
    """Fallback for bare bags. Filter diagnostic topics, never deserialize clouds."""
    result = []
    try:
        import rosbag2_py
        from rclpy.serialization import deserialize_message
        from rosidl_runtime_py.utilities import get_message
    except ImportError as exc:
        warnings.append('裸 rosbag 需要 ROS 2 及消息环境；请使用 diagnose_localization.sh。' + str(exc))
        return result
    files = [path] if path.is_file() else list(path.glob('*.mcap')) + list(path.glob('*.db3'))
    if sum(p.stat().st_size for p in files) > byte_limit:
        warnings.append(f'裸 bag 大于 {byte_limit // MIB} MiB 扫描预算；未扫描。可用 --bag-mb 增大，或提供事件/运行目录。')
        return result
    try:
        storage = 'sqlite3' if any(p.suffix == '.db3' for p in files) else 'mcap'
        reader = rosbag2_py.SequentialReader()
        reader.open(rosbag2_py.StorageOptions(uri=str(path), storage_id=storage), rosbag2_py.ConverterOptions('', ''))
        topics = {t.name: t.type for t in reader.get_all_topics_and_types()}
        selected = [t for t, typ in topics.items() if typ.endswith('/Safety')]
        if not selected:
            warnings.append('裸 bag 不包含 Safety 消息；无法重建原始故障判据。')
            return result
        reader.set_filter(rosbag2_py.StorageFilter(topics=selected))
        classes = {t: get_message(topics[t]) for t in selected}
        deadline = time.monotonic() + 15
        last = {}
        while reader.has_next():
            if time.monotonic() > deadline or len(result) >= 30000:
                warnings.append('裸 bag 达到 15 秒/30000 状态变化预算，结果不完整。')
                break
            topic, data, ns = reader.read_next()
            state = int(deserialize_message(data, classes[topic]).runstate)
            if last.get(topic) == state:
                continue
            last[topic] = state
            # The bag contains the reported code only; masked raw faults are unknowable.
            result.append(dict(raw_runstate=state, published_runstate=state, _time=wall_time(ns),
                               _wall_ns=int(ns), _source=f'{path}:{topic}@{ns}', _scope='裸 bag 上报码（原始码未知）'))
        warnings.append('裸 bag 的 Safety 只有上报码，不包含细分判据；关闭上报时不能据此推断原始故障。')
    except Exception as exc:
        warnings.append('rosbag 离线读取失败: ' + str(exc))
    return result


def cause_key(finding):
    return (finding['code'], finding['key'], finding.get('device'), finding.get('phase'))


def fault_episodes(rows, limits):
    """Keep every selected-predicate interval and device, including brief faults."""
    episodes, scopes = [], {}
    for row in rows:
        scopes.setdefault(row.get('_scope'), []).append(row)
    metrics = ('error_x', 'error_y', 'error_yaw_rad', 'output_age_sec',
               'imu_arrival_age_sec', 'major_sigma', 'lateral_sigma', 'yaw_sigma_rad')
    def score(finding):
        return max((abs(number(finding['evidence'].get(k)) or 0) for k in metrics), default=0)
    for scope_rows in scopes.values():
        active = {}
        previous = None
        for row in scope_rows:
            if previous and number(previous.get('sequence')) is not None and number(row.get('sequence')) is not None and row['sequence'] != previous['sequence'] + 1:
                for episode in active.values():
                    episode['state'] = 'observation_gap'
                active = {}
            previous = row
            current = {cause_key(f): f for f in classify(row, limits)}
            for key in list(active):
                if key not in current:
                    episode = active.pop(key)
                    episode['end_time'] = row.get('_time')
                    episode['end_ns'] = row_ns(row)
                    episode['next_raw'] = row.get('raw_runstate')
                    episode['end_localization_lost'] = row.get('localization_lost')
                    episode['state'] = 'recovered' if row.get('raw_runstate') == 0 else 'predicate_changed'
                    if episode['start_ns'] is not None and episode['end_ns'] is not None:
                        episode['observed_duration_sec'] = max(0, (episode['end_ns'] - episode['start_ns']) / 1e9)
            for key, finding in current.items():
                if key not in active:
                    episode = dict(finding, samples=0, first_time=finding['time'],
                                   start_ns=row_ns(row), last_ns=row_ns(row), end_ns=None,
                                   end_time=None, state='open', onset_evidence=finding['evidence'],
                                   onset_source=finding['source'], onset_advice=finding['advice'])
                    active[key] = episode
                    episodes.append(episode)
                episode = active[key]
                episode['samples'] += 1
                episode['last_time'] = finding['time']
                episode['last_ns'] = row_ns(row)
                # Preserve a representative worst sample, not a healthy-looking
                # last sample after input recovered but output is still stale.
                if score(finding) >= score(episode):
                    episode.update(finding)
    return episodes


def context_bounds(finding, rows):
    start = finding['start_ns']
    if start is None:
        return None, None
    end = finding.get('end_ns') or finding.get('last_ns') or start
    lower, upper = start - 2e9, end + 2e9
    previous = None
    for row in rows:
        stamp = row_ns(row)
        if stamp is None:
            continue
        gap = (previous is not None and number(previous.get('sequence')) is not None and
               number(row.get('sequence')) is not None and row['sequence'] != previous['sequence'] + 1)
        if number(row.get('raw_runstate')) == 0 or gap:
            if stamp < start or (gap and stamp == start):
                lower = max(lower, stamp)
            elif stamp >= start:
                upper = min(upper, (row_ns(previous) or start) if gap else stamp)
        previous = row
    return lower, upper


def matching_summary(clue, matching):
    values = fields(clue['text'])
    text = f'首次失锁前匹配记录: {clue["time"]}，map_checked={values.get("map_checked", "未知")}，map_overlap={values.get("map_overlap", "未知")}，raw_odom_valid={values.get("raw_odom_valid", "未知")}。'
    points, overlap = number(values.get('map_points')), number(values.get('map_overlap'))
    if values.get('map_checked') != 1:
        return text + ' 地图一致性未确认执行，不判断其失败项。', False
    if points == 0:
        return text + ' 地图一致性有效评估点数为 0，未获得可用比例证据。', False
    if not (points is not None and points > 0) and not ('map_points' not in values and overlap is not None and overlap > 0):
        return text + ' 缺少有效评估点证据，不从默认比例归因。', False
    comparisons, failed, config_used = [], [], False
    for metric, limit in (('map_inside_xy', 'map_min_inside_xy'),
                          ('map_overlap', 'map_min_overlap'),
                          ('map_gravity_cos', 'map_min_gravity_cos')):
        value, threshold = number(values.get(metric)), number(values.get(limit))
        basis = '运行日志判据'
        if threshold is None and metric == 'map_overlap' and limit not in values:
            threshold = matching.get('map_min_overlap_ratio')
            limit = matching.get('overlap_key', limit)
            basis = '运行配置判据'
        minimum = -1 if metric == 'map_gravity_cos' else 0
        if value is None or threshold is None or not minimum <= threshold <= 1:
            continue
        config_used |= basis == '运行配置判据'
        relation = '<' if value < threshold else '>='
        comparisons.append(f'{basis}: {metric}={value:g} {relation} {limit}={threshold:g}')
        if value < threshold:
            failed.append(metric)
    if comparisons:
        text += ' ' + '；'.join(comparisons) + '。'
        text += ' 已证实未满足项: ' + '、'.join(failed) + '。' if failed else ' 已记录的这些判据未见超限；其他未记录条件不作判断。'
    text += ' 不据此推断车辆出图或其他物理原因。'
    return text, config_used


def priority_summary(findings, matching):
    """A short evidence chain; proximity never establishes a physical cause."""
    findings = [f for f in findings if f['scope'] == '运行日志（有界尾部）'] or findings
    summary = []
    def add(text, source):
        summary.append(dict(text=text, source=source))
    yaw_failures = [f for f in findings if f['key'] == 'yaw_condition' and f['evidence'].get('yaw_condition') == 2]
    if yaw_failures:
        last = max(yaw_failures, key=lambda f: f['start_ns'] or 0)
        error = number(last['evidence'].get('error_yaw_deg'))
        residual = f'{error:.6g}°' if error is not None else '未记录'
        recovery = ('已观测到原始码 0: ' + last['end_time']) if last.get('state') == 'recovered' else (
            '后续存在观测缺口，恢复时刻未知' if last.get('state') == 'observation_gap' else '未观测到本段恢复为 0')
        kind = '日志判据段' if last['scope'] == '运行日志（有界尾部）' else '判据段'
        add(f'航向运动一致性 FAIL: 已读取 {len(yaw_failures)} 个{kind}（不等于物理故障次数）；'
            f'最近段 {last["first_time"]}，残差 {residual}（样例 {last["time"]}）；{recovery}。', last['source'])
    startup = next((f for f in findings if f.get('phase') == 'startup'), None)
    severe = [f for f in findings if f['code'] == -9 and f['start_ns'] is not None]
    if severe:
        first_severe = min(severe, key=lambda f: f['start_ns'])
        initial = [f['title'] for f in severe if f['start_ns'] == first_severe['start_ns']]
        add(f'首次已观测 -9: {first_severe["first_time"]}；当时判据: {"、".join(initial)}。', first_severe['onset_source'])
        seen_queue = set()
        for finding in sorted(severe, key=lambda f: f['start_ns']):
            queue = next((c for c in finding.get('onset_context', []) if 'sensor queue lag reached' in c['text']), None)
            if queue is None or queue['source'] in seen_queue:
                continue
            lag = re.search(r'sensor queue lag reached\s+([\d.eE+-]+)\s+sec', queue['text'])
            lag_sec = number(lag[1]) if lag else None
            if lag_sec is None:
                continue
            seen_queue.add(queue['source'])
            text = (f'与 -9 起点 {finding["first_time"]} 相近的处理链线索: {queue["time"]} '
                    f'处理进度落后 {lag_sec:g}s（最新入队与最新处理的传感器时间差）。')
            source = queue['source']
            stale = next((c for c in finding.get('onset_context', []) if 'drop stale lidar correction before matching:' in c['text']), None)
            if stale:
                values = fields(stale['text'])
                text += (f' {stale["time"]} 过期 LIO 校正被丢弃，age={values.get("age", "未知")}s，'
                         f'hard_deadline={values.get("hard_deadline", "未知")}s。')
                source += '；' + stale['source']
            text += ' 这是时间相近的相关证据，不证明最早物理根因；不能据此断定 CPU 阻塞、DDS 故障或全部雷达硬件故障。'
            add(text, source)
    losses = [f for f in findings if f['key'] == 'localization_lost' and f['start_ns'] is not None]
    if losses:
        first = min(losses, key=lambda f: f['start_ns'])
        add(f'首次已观测门控关闭: {first["first_time"]}，localization_lost=1；'
            f'output_age_sec={first["onset_evidence"].get("output_age_sec", "未知")}。', first['onset_source'])
        clue = next((c for c in first.get('onset_context', []) if 'localization failed!' in c['text']), None)
        if clue:
            text, config_used = matching_summary(clue, matching)
            add(text, clue['source'] + ('；配置: ' + matching['source'] if config_used else ''))
        longest = max(losses, key=lambda f: (f.get('end_ns') or f.get('last_ns') or f['start_ns']) - f['start_ns'])
        if longest is not first:
            next_raw = number(longest.get('next_raw'))
            if longest.get('state') == 'recovered':
                state = '门控关闭判据已解除，且同条记录原始码为 0'
            elif longest.get('state') == 'predicate_changed':
                cleared = (longest.get('end_localization_lost') == 0 or next_raw in (0, 2, 4))
                state = ('门控关闭判据已解除' if cleared else '门控关闭判据不再出现，解除状态缺字段确认')
                state += f'；该时刻原始码仍为 {next_raw:g}' if next_raw is not None else '；该时刻原始码未知'
            elif longest.get('state') == 'observation_gap':
                state = '后续存在观测缺口，门控判据解除时刻未知'
            else:
                state = '截至该段末保留记录，尚未观测到门控关闭判据解除'
            recovery = [f for f in findings if f.get('state') == 'recovered' and f.get('end_ns') is not None and
                        f['end_ns'] >= (longest.get('end_ns') or longest.get('last_ns') or longest['start_ns'])]
            if longest.get('state') != 'recovered' and recovery:
                state += '；后续已观测到原始码 0: ' + min(recovery, key=lambda f: f['end_ns'])['end_time']
            add(f'最长已观测门控关闭段: {longest["first_time"]} → {longest.get("end_time") or longest.get("last_time")}；{state}。', longest['onset_source'])
        for marker, title in (('SOLID_KISS_REJECT', '该段最近保留的重定位拒绝'), ('] accepted', '该段重定位接受')):
            clue = next((c for c in reversed(longest.get('context', [])) if marker in c['text']), None)
            if clue:
                values = fields(clue['text'])
                detail = ', '.join(f'{k}={values[k]}' for k in ('candidate', 'inliers', 'reason', 'overlap') if k in values)
                add(f'{title}: {clue["time"]}，{detail}。', clue['source'])
    sync = [f for f in findings if f['code'] == 2 and f['key'] == 'radar_time_sync']
    for role in sorted({f.get('device', '未知') for f in sync}):
        episodes = [f for f in sync if f.get('device', '未知') == role]
        last = max(episodes, key=lambda f: f['start_ns'] or 0)
        if last.get('state') == 'recovered':
            ending = '已观测到原始码 0: ' + last['end_time']
        elif last.get('state') == 'predicate_changed':
            ending = '同步判据段结束: ' + last['end_time'] + '；当时未恢复为原始码 0'
        elif last.get('state') == 'observation_gap':
            ending = '后续存在观测缺口，恢复时刻未知'
        else:
            ending = '尚未观测到该段恢复'
        add(f'雷达同步故障 2: {role}，已读取 {len(episodes)} 个判据段；最近段 {last["first_time"]}；{ending}。', last['onset_source'])
    if startup:
        add(f'启动无有效定位输出: {startup["first_time"]} → {startup.get("end_time") or "尚未见结束"}。', startup['onset_source'])
    motion = [f for f in findings if f['key'] == 'consistency_unknown' and
              number(f['onset_evidence'].get('output_age_sec')) is not None and
              str(f['onset_evidence'].get('yaw_pass_age_sec')) == 'inf']
    if motion:
        last = max(motion, key=lambda f: f['start_ns'] or 0)
        ending = (' → ' + last['end_time'] + ' 已观测到原始码 0') if last.get('state') == 'recovered' else '；未据此确认恢复'
        add(f'有定位输出后的运动证据 UNKNOWN: {last["first_time"]}{ending}；尚无可用航向 PASS，窗口建立/重建与历史 IMU 覆盖仍需区分。', last['onset_source'])
    return summary


def capture_windows(path, run, process, summary, rows, warnings):
    """Summary envelopes show possible coverage, never message completeness."""
    candidates = [path] if summary else []
    if process:
        root = Path(process['archive']) if process.get('archive') else run / 'results/incidents'
        candidates += [p for p in events_in(root) if same_run(p, run, process['pid'])]
    elif not summary and run.is_dir():
        candidates += [p for root in (run / 'results/incidents', run.parent / 'incidents')
                       for p in events_in(root) if same_run(p, run)]
    candidates = list(dict.fromkeys(candidates))
    if len(candidates) > 256:
        candidates = sorted(candidates, key=stamp, reverse=True)[:256]
        warnings.append('保留事件超过 256 个，仅核查最近 256 个摘要的覆盖范围。')
    windows = []
    for event in candidates:
        s = summary if event == path else summary_at(event)
        trigger = number(s.get('trigger_wall_ns'))
        if trigger is None:
            continue
        # Diagnose the overlap using the recorded Safety stream, not a cloud
        # stream whose range may be longer while diagnostics were omitted.
        ranges = [(stream.get('first_relative_sec'), stream.get('last_relative_sec'))
                  for stream in s.get('streams', [])
                  if str(stream.get('type', '')).endswith('/Safety') and stream.get('written', 0)]
        for first, last in ranges:
            if number(first) is not None and number(last) is not None:
                windows.append(dict(path=str(event), start_ns=int(trigger + float(first) * 1e9),
                                    end_ns=int(trigger + float(last) * 1e9), status=s.get('status'),
                                    evidence='Safety 首末记录范围，中间完整性未验证'))
        if event == path and rows:
            times = [row_ns(r) for r in rows if row_ns(r) is not None]
            if times:
                windows.append(dict(path=str(event), start_ns=min(times), end_ns=max(times),
                                    status=s.get('status'), evidence='时间线首末记录范围，中间完整性未验证'))
    return windows


def match_failure_collector(report):
    """Only a live process can authorize reading its current collector status."""
    process = report.get('process')
    if not process or report.get('evidence_scope') in ('explicit', 'historical'):
        return dict(status='not_requested', message='历史/指定包未拼接当前采集器状态。')
    run = Path(process['run'])
    roots = [run / 'match_failures', run / 'results/match_failures', run.parent / 'match_failures']
    declared_root = key_values(run / 'run_metadata.txt').get('match_failure_dir')
    if declared_root:
        roots.insert(0, Path(declared_root))
    if process.get('match_failure_root'):
        roots.insert(0, Path(process['match_failure_root']))
    for root in dict.fromkeys(roots):
        path = root / 'collector_status.json'
        if root.is_symlink() or path.is_symlink() or not path.is_file():
            continue
        try:
            if path.stat().st_size > MIB:
                continue
            value = json.loads(read_text(path))
            if not isinstance(value, dict) or value.get('schema') != 1 or value.get('kind') != 'lightning_match_failure_collector':
                continue
            if value.get('pid') != process['pid'] or not value.get('run_dir') or Path(value['run_dir']).resolve() != run.resolve():
                continue
            if int(value.get('session_start_wall_ns', 0)) < process.get('start_wall_ns', 0):
                continue
            counters = {}
            for key in ('received', 'admitted', 'saved', 'suppressed', 'dropped_busy', 'dropped_queue',
                        'dropped_oversize', 'rejected_disk', 'rejected_map', 'writer_errors'):
                count = number(value.get(key))
                if count is not None and count >= 0:
                    counters[key] = int(count)
            age = max(0, (time.time_ns() - path.stat().st_mtime_ns) / 1e9)
            return dict(status='observed', source=str(path), running=value.get('running') is True,
                        status_age_sec=round(age, 3), fresh=age <= 10, counters=counters,
                        last_error=str(value.get('last_error', ''))[:1024],
                        protection_days=number(value.get('protection_days')), disk_mb=number(value.get('disk_mb')),
                        limitation='采集器进程累计计数，不代表所选快照损坏；rejected_disk 含配额/剩余空间不足。')
        except (OSError, ValueError, TypeError, OverflowError):
            continue
    # This is a declared setting, not a claim that a missing writer is healthy.
    config = Path(process.get('config') or run / 'config.yaml')
    try:
        section = False
        level = None
        for line in read_text(config).splitlines():
            stripped = line.strip()
            indent = len(line) - len(line.lstrip())
            if stripped == 'match_failure_capture:':
                section, level = True, indent
                continue
            if section and stripped and not stripped.startswith('#') and indent <= level:
                section = False
            if section and re.fullmatch(r'enabled:\s*(false|False|FALSE|0)(\s*#.*)?', stripped):
                return dict(status='configured_disabled', source=str(config), message='运行配置快照关闭了失败匹配采集；未找到同次采集器状态。')
    except OSError:
        pass
    return dict(status='unavailable', message='未找到与当前运行/PID/启动时间匹配的采集器状态，无法确认是否启用或写入成功。')


def match_failure_image(report, explicit=None, output_dir=None, enabled=True):
    result = _match_failure_image(report, explicit, output_dir, enabled)
    result['collector'] = match_failure_collector(report)
    return result


def _match_failure_image(report, explicit=None, output_dir=None, enabled=True):
    """Select one complete, same-run frame; rendering is bounded and optional."""
    if not enabled:
        return dict(status='disabled', message='已跳过失败匹配俯视图。')
    import match_failure_birdseye as bird
    chosen = Path(report['selected']) if report.get('selected') else None
    process = report.get('process')
    run = Path(report.get('resolved_run') or process.get('run')) if (report.get('resolved_run') or process) else None
    direct = Path(explicit).expanduser() if explicit else chosen
    if direct and direct.name == 'summary.json':
        direct = direct.parent
    evidence = None
    selection = 'explicit_frame' if explicit else 'selected_frame'
    if direct and bird.is_evidence(direct):
        evidence = direct
    elif explicit:
        return dict(status='unavailable', message='指定路径不是完整的失败匹配证据包: ' + str(direct))
    elif run:
        # A copied run can retain its original absolute identity in metadata.
        # Never associate merely by directory basename or reused PID.
        run_meta = key_values(run / 'run_metadata.txt')
        identities = {str(run.resolve())}
        if run_meta.get('run_dir'):
            identities.add(str(Path(run_meta['run_dir']).resolve()))
        window = None
        if chosen and chosen != run and not process:
            windows = [w for w in report.get('capture_windows', []) if Path(w['path']).resolve() == chosen.resolve()]
            if not windows:
                return dict(status='unavailable', message='指定事件缺少时间范围，未使用同次运行其他时刻的匹配图。')
            window = (min(w['start_ns'] for w in windows), max(w['end_ns'] for w in windows))
        roots = [run / 'match_failures', run / 'results/match_failures', run.parent / 'match_failures']
        if run_meta.get('match_failure_dir'):
            roots.insert(0, Path(run_meta['match_failure_dir']))
        if chosen and chosen != run:
            roots.insert(0, chosen / 'match_failures')
        if process and process.get('match_failure_root'):
            roots.insert(0, Path(process['match_failure_root']))
        candidates, checked, invalid = [], 0, 0
        for root in dict.fromkeys(roots):
            if not root.is_dir() or root.is_symlink():
                continue
            # Names contain creation timestamps; cap discovery before reading
            # JSON, and report truncation instead of silently claiming coverage.
            directories = sorted((p for p in root.iterdir() if p.is_dir() and not p.is_symlink()),
                                 key=lambda p: p.name, reverse=True)
            for item in directories[:bird.MAX_CANDIDATES - checked]:
                checked += 1
                try:
                    summary = bird.read_summary(item)
                    if summary.get('synthetic_fixture') is True:
                        continue
                    if not summary.get('run_dir') or str(Path(summary['run_dir']).resolve()) not in identities:
                        continue
                    wall = int(summary['wall_ns'])
                    if process and (summary.get('pid') != process['pid'] or wall < process.get('start_wall_ns', 0)):
                        continue
                    if window and not window[0] <= wall <= window[1]:
                        continue
                    candidates.append((wall, item))
                except (OSError, ValueError, TypeError, KeyError):
                    invalid += 1
            if checked >= bird.MAX_CANDIDATES:
                break
        if candidates:
            _, evidence = max(candidates, key=lambda pair: pair[0])
            selection = 'latest_same_run_frame_in_event' if window else 'latest_same_run_frame'
        else:
            return dict(status='not_found', checked=checked, unreadable_or_incomplete=invalid,
                        discovery_limited=checked >= bird.MAX_CANDIDATES,
                        message='未找到所选运行/事件范围内完整的失败匹配点云证据；日志和轨迹不能重建当时动态物体。')
    if evidence is None:
        return dict(status='not_found', message='没有可关联的失败匹配证据。')
    # Avoid user-site NumPy shadowing the distro SciPy/matplotlib ABI on Orin.
    command = [sys.executable, '-s', str(Path(__file__).with_name('match_failure_birdseye.py')), str(evidence)]
    if output_dir:
        command += ['--output-dir', str(output_dir)]
    environment = dict(os.environ)
    environment.update(PYTHONNOUSERSITE='1', PYTHONIOENCODING='utf-8', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
    try:
        completed = subprocess.run(command, capture_output=True, text=True, encoding='utf-8', errors='replace',
                                   timeout=30, env=environment)
        if len(completed.stdout) > MIB:
            raise ValueError('绘图返回内容超过限制')
        result = json.loads(completed.stdout)
        if not isinstance(result, dict) or result.get('status') not in ('generated', 'cached', 'unavailable', 'dependency_missing'):
            raise ValueError('绘图返回格式无效')
        result['selection'] = selection
        return result
    except subprocess.TimeoutExpired:
        return dict(status='timeout', evidence=str(evidence), message='俯视图计算超过 30 秒，已终止绘图；故障诊断结果保留。')
    except (OSError, ValueError, TypeError) as exc:
        return dict(status='unavailable', evidence=str(evidence), message='俯视图生成失败，诊断继续: ' + str(exc))


def runtime_contract(logged_limits, warnings):
    built = logged_limits.get('_build', {})
    try:
        expected = json.loads(read_text(CONTRACT_PATH))['contract_id']
    except (OSError, ValueError, KeyError):
        expected = None
        warnings.append('诊断合同文件缺失或无效，无法检查脚本与运行库的兼容性。')
    if not built:
        status = 'unverified'
        warnings.append('运行日志没有编译版本标记，无法确认新采集修复已生效；run_metadata 的 git_commit 只代表启动时源码。')
    elif expected is None or built.get('contract_id') != expected or number(built.get('schema')) != 2:
        status = 'mismatch'
        warnings.append('运行库与脚本的诊断合同不匹配；结果按现有字段解释，不能视为最新版完整验收。')
    else:
        status = 'matched'
    return dict(status=status, expected_contract=expected, compiled=built)


def analyze(path, process, selection, warnings, log_mb=32, bag_mb=128):
    warnings = list(warnings)
    summary = summary_at(path) if path.is_dir() else {}
    meta = key_values(path / 'run_metadata.txt') if path.is_dir() else {}
    run = path
    if process:
        run = Path(process['run'])
    elif not is_run_directory(path) and (meta.get('run_dir') or summary.get('kind') == 'lightning_match_failure' and summary.get('run_dir')):
        original = Path(meta.get('run_dir') or summary['run_dir'])
        relocated = path.parent.parent / original.name
        run = relocated if (relocated / 'logs/run_loc_online.stderr.log').is_file() else original
    config = path / 'config.yaml' if path.is_dir() else path.parent / 'config.yaml'
    launch_findings = launch_failures(run, warnings) if run == path and not summary else []
    limits = thresholds(config)
    matching = matching_thresholds(config)
    log = process.get('log') if process else str(run / 'logs/run_loc_online.stderr.log')
    try:
        reference = dt.datetime.fromisoformat(meta['started_at'])
    except (KeyError, ValueError):
        reference = dt.datetime.fromtimestamp(int(summary.get('trigger_wall_ns', time.time_ns())) / 1e9).astimezone()
    rows = timeline_rows(path / 'timeline.csv', warnings) if path.is_dir() else []
    protection_window = None
    protection_scope = '已扫描运行日志（有界尾部，不代表整个运行）'
    if not process and path.is_dir() and ((path / 'summary.json').is_file() or (path / 'timeline.csv').is_file()):
        times = [row_ns(r) for r in rows if row_ns(r) is not None]
        protection_window = (min(times), max(times)) if times else ()
        protection_scope = '指定事件时间线首末范围内已扫描日志（不扩展窗口）'
        if not times:
            warnings.append('缺少指定事件时间线范围，未统计运行保护线索；不使用窗口外日志补充次数。')
    context, protections = [], {}
    logdata, logged_limits = log_rows(log, int(log_mb * MIB), warnings, reference, context,
                                      protections, protection_window)
    if not limits:
        limits = {k: v for k, v in logged_limits.items() if k in LIMITS}
        warnings.append('缺少配置快照，使用能找到的运行配置日志；未知阈值不作超限判断。')
    build = runtime_contract(logged_limits, warnings)
    if rows and not all(k in rows[-1] for k in ('translation_condition', 'yaw_condition', 'output_stale', 'time_sync_fault_mask')):
        warnings.append('事件时间线是旧格式或缺少细分字段；需同次运行日志补充，不能仅凭包头源码版本认定已升级。')
    if summary and not process:
        times = [number(r.get('wall_ns')) for r in rows if number(r.get('wall_ns')) is not None]
        if times:
            logdata = [r for r in logdata if min(times) - 1e9 <= r.get('_wall_ns', -1) <= max(times) + 1e9]
        else:
            logdata = []
        if not logdata:
            warnings.append('未找到与指定事件时间窗口对应的详细运行日志；只使用包内证据。')
    if not rows and not logdata and not launch_findings:
        rows = bag_rows(path, int(bag_mb * MIB), warnings)
    if summary:
        if summary.get('status') != 'complete':
            warnings.append('事件状态 ' + str(summary.get('status')) + '；采集中/采集缺失不等于算法故障。')
        if summary.get('error'):
            warnings.append('采集错误: ' + summary['error'])
        if summary.get('capture_loss_seen'):
            warnings.append('capture_loss_seen 是进程累计采集丢弃标记，不等于本事件全部丢失。')
        if 'startup_without_output' in summary.get('reasons', []):
            ages = [number(r.get('publish_end_age')) for r in rows]
            if any(a is not None and a >= 0 for a in ages):
                warnings.append('startup_without_output 与窗口内已有发布完成记录并存；旧版并发时间采样可误触发此标签，不能据此认定启动无输出。')
        if process:
            warnings.append('同时显示本事件和当前运行日志尾部的故障；日志中的较早故障不一定属于此事件，见每项时间与来源。')
    findings = fault_episodes(rows + logdata, limits)
    scope_rows = {}
    for row in rows + logdata:
        scope_rows.setdefault(row.get('_scope'), []).append(row)
    windows = capture_windows(path, run, process, summary, rows, warnings)
    for finding in findings:
        start = finding['start_ns']
        finding['event_paths'] = sorted(set(w['path'] for w in windows
            if start is not None and w['start_ns'] <= start <= w['end_ns']))
        finding['capture_coverage'] = ('range_overlap' if finding['event_paths'] else
                                     'no_retained_overlap' if start is not None else 'unknown_time')
        finding['capture_coverage_basis'] = 'episode_start'
        # Nearby facts are corroboration, not a claim about the physical root cause.
        lower, upper = context_bounds(finding, scope_rows.get(finding['scope'], []))
        nearby_context = [c for c in context if start is not None and c['wall_ns'] is not None
                          and lower <= c['wall_ns'] <= upper]
        onset_context = []
        for marker in ('localization failed!', 'localization pose and map outputs disabled'):
            candidates = [c for c in nearby_context if c['wall_ns'] <= start + 2e9 and marker in c['text']]
            before = [c for c in candidates if c['wall_ns'] <= start]
            after = [c for c in candidates if c['wall_ns'] > start]
            clue = (before[-1] if before else None) if marker == 'localization failed!' else (after[0] if after else before[-1] if before else None)
            if clue is not None:
                onset_context.append(clue)
        if finding['code'] == -9:
            for marker in ('sensor queue lag reached', 'drop stale lidar correction before matching:'):
                candidates = [c for c in nearby_context if c['wall_ns'] <= start + 2e9 and marker in c['text']]
                if candidates:
                    onset_context.append(min(candidates, key=lambda c: c['wall_ns']))
        finding['onset_context'] = [c for c in nearby_context if c in onset_context]
        finding['context'] = recent_context_clues([c for c in nearby_context if c not in onset_context], 6)
    absent = sum(f['capture_coverage'] == 'no_retained_overlap' for f in findings if f['scope'] == '运行日志（有界尾部）')
    if absent:
        warnings.append(f'{absent} 个日志判据段的起始时刻未找到保留事件范围覆盖；可能尚在写入、触发被抑制或已轮转，不代表从未采集。')
    if not rows and not logdata:
        warnings.append('未产生算法 Safety 判据；已单独解释启动前失败日志。' if launch_findings else
                        '没有可解析的 Safety 判据记录；此次不能给出故障原因。')
        if is_run_directory(path) and not launch_findings:
            warnings.append('所选运行目录尚无算法判据和明确入口错误；可能仍在启动或提前退出，不能只凭文件缺失认定故障。')
    all_rows = rows + logdata
    coexisting_radar = {}
    for row in all_rows:
        if number(row.get('raw_runstate')) != -9:
            continue
        mask = int(number(row.get('time_sync_fault_mask')) or 0)
        if row.get('bounded_policy_enabled') == 1:
            mask &= ~int(number(row.get('bounded_sync_warning_mask')) or 0)
        conditions = [('radar_time_sync', role) for bit, role in enumerate(('front', 'rear', 'left', 'right'))
                      if mask & (1 << bit)]
        if row.get('non_sync_fault') == 1:
            conditions.append(('radar_other', None))
        for key, device in conditions:
            previous = coexisting_radar.get((key, device))
            if previous and (previous['last_ns'] or 0) > (row_ns(row) or 0):
                continue
            coexisting_radar[key, device] = dict(
                key=key, device=device, raw_runstate=-9, last_ns=row_ns(row),
                last_time=row.get('_time'), source=row.get('_source'), scope=row.get('_scope'),
                evidence={k: row[k] for k in ('time_sync_fault_mask', 'bounded_policy_enabled',
                          'bounded_sync_warning_mask', 'non_sync_fault') if k in row})
    latest = max(enumerate(all_rows), key=lambda pair: (row_ns(pair[1]) or 0, pair[0]))[1] if all_rows else {}
    age = max(0, (time.time_ns() - row_ns(latest)) / 1e9) if process and row_ns(latest) else None
    if age is not None and age > 5:
        warnings.append(f'最新诊断记录已过去 {age:.1f} 秒，不能作为当前实时状态。')
    dropped = [number(r.get('collector_dropped')) for r in rows]
    dropped = [v for v in dropped if v is not None]
    return dict(schema=2, selected=str(path), resolved_run=str(run), selection=selection, process=process,
                has_match_evidence=summary.get('kind') == 'lightning_match_failure' and summary.get('status') == 'complete',
                event_status=summary.get('status'), trigger_time=wall_time(summary.get('trigger_wall_ns')),
                trigger_reasons=summary.get('reasons', []), thresholds=limits,
                latest={k: latest.get(k) for k in ('_time', 'raw_runstate', 'published_runstate', '_source', '_scope')},
                latest_age_sec=age, current_findings=classify(latest, limits), runtime_build=build,
                capture_windows=windows, collector_totals=summary.get('collector_totals', {}),
                window_collector_drop_delta=max(dropped) - min(dropped) if dropped else None,
                event_samples=dict(Counter(str(int(r['raw_runstate'])) for r in rows if number(r.get('raw_runstate')) is not None)),
                log_samples=len(logdata), findings=findings, warnings=warnings,
                priority_summary=priority_summary(findings, matching), matching_thresholds=matching,
                runtime_protections=list(protections.values()), runtime_protection_scope=protection_scope,
                coexisting_radar_conditions=list(coexisting_radar.values()),
                launch_findings=launch_findings, has_launch_evidence=bool(launch_findings),
                has_safety_evidence=bool(rows or logdata),
                limitation='输出的是故障判据及相关线索，不是已证明的物理根因；日志/时间线样本数不是故障次数或时间占比。')


def render(report, color=False, all_events=False):
    def paint(text, style):
        return f'\033[{style}m{text}\033[0m' if color else text
    lines = [paint('定位现场自动诊断', '1;36'), f'选择: {report["selection"]}', f'数据包: {report["selected"]}']
    if report['process']:
        p = report['process']
        lines.append(f'定位: PID {p["pid"]}, 线程名 {p["comm"]}, exe={p["exe"]}')
    lines.append(f'检查器调度: {report.get("scheduling", {})}')
    if report.get('evidence_scope') == 'historical':
        lines.append('定位证据范围: 历史包；当前未发现定位进程，包内状态码不代表当前状态。')
    elif report.get('evidence_scope') == 'explicit':
        lines.append('定位证据范围: 指定包（离线）。')
    for finding in report.get('launch_findings', []):
        prefix = ('入口阻断证据' if not report.get('has_safety_evidence') and not report.get('process') else
                  '启动预检历史失败记录（也存在算法实例或判据，不据此认定当前仍阻断）')
        lines += [paint(f'{prefix}: {finding["title"]}；层级={finding["stage"]}', '1;31'),
                  '  原始错误: ' + finding['evidence'], '  证据: ' + finding['source'],
                  '  定位方向: ' + finding['advice']]
    if report.get('upstream_lidar'):
        lines += render_driver_diagnosis(report['upstream_lidar'])
    match_image = report.get('match_failure_image', {})
    if match_image:
        lines.append(paint('失败匹配三姿态俯视图: ' + str(match_image.get('status')), '1;36'))
        collector = match_image.get('collector', {})
        if collector.get('status') == 'observed':
            lines.append('  失败证据采集器: ' + ('运行标记为 true' if collector.get('running') else '运行标记为 false') +
                         f'；状态文件距今 {collector.get("status_age_sec")} 秒；累计计数=' + json.dumps(collector.get('counters', {}), ensure_ascii=False))
            lines.append('  ' + collector['limitation'])
            lines.append('  采集状态证据: ' + collector['source'])
            if not collector.get('fresh'):
                lines.append(paint('  采集器状态文件已陈旧，不能视为当前实时状态。', '33'))
            if collector.get('counters', {}).get('rejected_disk', 0):
                lines.append(paint('  曾因配额/剩余空间不足拒绝保存；保护期内旧证据不会为新事件让位，缺图不能解释为没有匹配失败。', '33'))
            if collector.get('last_error'):
                lines.append(paint('  采集器最近错误: ' + collector['last_error'], '33'))
        elif collector.get('status') in ('configured_disabled', 'unavailable'):
            lines.append('  ' + collector.get('message', ''))
        if match_image.get('image'):
            scope = {'explicit_frame': '明确指定的失败帧', 'selected_frame': '所选失败帧',
                     'latest_same_run_frame_in_event': '指定事件时间范围内最新失败帧',
                     'latest_same_run_frame': '同次运行最新保留失败帧'}.get(match_image.get('selection'), '未知')
            lines += ['  采集时间: ' + wall_time(match_image.get('wall_ns')) + '；扫描时间: ' + wall_time(match_image.get('sensor_stamp_ns')),
                      '  失败原因: ' + str(match_image.get('reason')),
                      '  范围: ' + scope + '；这是保留帧，不保证对应最新 Safety 故障。',
                      '  原始证据: ' + str(match_image.get('evidence')),
                      '  图片: ' + match_image['image'], '  可离线打开: ' + match_image['html'],
                      '  初值 / NDT 结果 / 实际检查位姿；按三维最近邻计分后俯视投影，红点不等于动态物体。']
            if match_image.get('synthetic_fixture'):
                lines.append(paint('  合成测试数据，不是现场事故证据。', '33'))
            for name, title in (('guess', '初值'), ('ndt', 'NDT 结果'), ('check', '实际检查位姿')):
                value = match_image.get('comparison', {}).get(name, {})
                lines.append(f'  {title}: 命中={value.get("hit_points")}/{value.get("validation_points")}，重合率={value.get("overlap")}')
            if match_image.get('check_score_matches_recorded') is False:
                lines.append(paint('  注意: 重算检查分数与记录不一致，不能视为精确复现；请核查证据/版本。', '33'))
            elif match_image.get('original_sample_replay_available') is False:
                lines.append(paint('  原判据未执行或无有效样本；图中比分仅针对保存的显示样本，不能视为原判据复现。', '33'))
            else:
                lines.append('  原检查重合率: ' + str(match_image.get('original_check_overlap')) + '；复算对照: ' + str(match_image.get('check_score_matches_recorded')))
        elif match_image.get('message'):
            lines.append('  ' + match_image['message'])
    if report.get('localization_error'):
        lines.append('定位证据不可用: ' + report['localization_error'])
        return '\n'.join(lines)
    if report.get('has_match_evidence'):
        lines.append('本次指定失败匹配快照；图中三姿态不代表独立真值，动态物体需结合连续帧判断。')
        return '\n'.join(lines)
    if report.get('launch_findings') and not report.get('has_safety_evidence'):
        lines.append('本包没有可解析的算法 Safety 记录；以上为启动入口失败，不是运行中 2/4/-9 判据。')
        lines += [paint('注意: ' + w, '33') for w in report['warnings']]
        return '\n'.join(lines)
    if report['event_status']:
        lines.append(f'事件: {report["trigger_time"]}, {report["event_status"]}, 触发={report["trigger_reasons"]}')
    latest = report['latest']
    age = report.get('latest_age_sec')
    fresh = report.get('evidence_scope') not in ('historical', 'explicit') and (age is None or age <= 5)
    raw = number(latest['raw_runstate'])
    style = '1;31' if fresh and raw is not None and raw != 0 else '1;32' if fresh and raw == 0 else '1;33'
    lines.append(paint(f'最新已读取状态: {latest["_time"]}, 原始={latest["raw_runstate"]}, 上报={latest["published_runstate"]}', style))
    if raw == 0:
        lines.append('最新记录为 0，不代表整个运行无故障；下方保留历史证据。')
    if report.get('priority_summary'):
        lines.append(paint('优先结论 / 问题链（已读取证据，不是物理根因）', '1;36'))
        for item in report['priority_summary']:
            lines.append('  ' + item['text'])
            lines.append('    证据: ' + item['source'])
    build = report.get('runtime_build', {})
    compiled = build.get('compiled', {})
    status = {'matched': '诊断合同匹配', 'unverified': '缺少编译标记，尚未核实',
              'mismatch': '诊断合同不匹配'}.get(build.get('status'), '未知')
    lines.append(paint(f'运行库版本核对: {status}；提交={compiled.get("source_revision", "未知")}；'
                       f'源码含改动={compiled.get("source_dirty", "未知")}；编译={compiled.get("built_at", "未知")}',
                       '36' if build.get('status') == 'matched' else '33'))
    for finding in report.get('current_findings', []):
        lines.append(paint(f'最近记录判据 [{finding["code"]}] {finding["title"]}', '1;31' if fresh else '33'))
        lines.append('  定位方向: ' + finding['advice'])
    if report.get('runtime_protections'):
        lines.append(paint('运行保护线索（独立于 Safety 故障判据，不计入故障段，不等于定位失败）', '33'))
        lines.append('  统计范围: ' + report['runtime_protection_scope'])
        for clue in report['runtime_protections']:
            lines += [f'  {clue["title"]}: {clue["count"]} 条；首末记录: {clue["first_time"]} → {clue["last_time"]}',
                      '  最近读取数值: ' + json.dumps(clue['evidence'], ensure_ascii=False),
                      '  首条位置: ' + clue['first_source'],
                      f'  最近原始日志: {clue["text"]}\n    位置: {clue["source"]}']
    if report.get('coexisting_radar_conditions'):
        lines.append(paint('检查范围内曾与 -9 并存的较低优先级雷达诊断条件（对应记录原始码为 -9；不代表当前仍异常，不计入故障段，不能据此认定启动失败原因）', '33'))
        for clue in report['coexisting_radar_conditions']:
            title = '时间同步诊断未通过: ' + clue['device'] if clue['key'] == 'radar_time_sync' else '非时间同步诊断失败'
            lines += [f'  {title}；最近并存观测: {clue["last_time"]}；来源: {clue["scope"]}',
                      '  数值: ' + json.dumps(clue['evidence'], ensure_ascii=False),
                      '  证据: ' + str(clue['source'])]
    lines.append('生效阈值: ' + json.dumps(report['thresholds'], ensure_ascii=False))
    findings = sorted(report['findings'], key=lambda f: f.get('start_ns') or 0, reverse=True)
    lines.append(paint(f'故障判据段: 共 {len(findings)} 段（不同证据来源分别统计，不等于物理故障次数）', '1;36'))
    groups = Counter((f['code'], f['key'], f.get('device', ''), f.get('phase', '')) for f in findings)
    for (code, key, device, phase), count in groups.items():
        example = next(f for f in findings if (f['code'], f['key'], f.get('device', ''), f.get('phase', '')) == (code, key, device, phase))
        lines.append(f'  [{code}] {example["title"]}: {count} 段')
    shown = findings if all_events else findings[:12]
    if len(shown) < len(findings):
        lines.append(paint(f'以下显示最近 {len(shown)} 段；--all 显示全部，--json 始终包含全部已分析判据段。', '33'))
    for f in shown:
        ended = f.get('end_time') or f.get('last_time') or '时间未知'
        state = '检查范围结束，未见后续状态（不等于已恢复）'
        if f.get('state') == 'observation_gap':
            state = '后续存在观测缺口，结束或恢复时刻未知'
        elif f.get('state') == 'recovered':
            state = '已观测到原始码恢复为 0'
        elif f.get('state') == 'predicate_changed':
            next_code = number(f.get('next_raw'))
            next_label = f'{next_code:g}' if next_code is not None else '未知'
            state = ('判据已变化，原始码仍为 ' if next_code == f['code'] else '原始码切换为 ') + next_label + '（未恢复为 0）'
        lines += [paint(f'历史记录 [{f["code"]}] {f["title"]}（{f["confidence"]}；{f["scope"]}）', '33'),
                  f'  时间段: {f["first_time"]} → {ended}；判据样本={f["samples"]}；样例上报码={f["published"]}',
                  '  时段状态: ' + state,
                  f'  样例时间: {f["time"]}；数值: ' + json.dumps(f['evidence'], ensure_ascii=False),
                  f'  证据: {f["source"]}']
        if f['onset_evidence'] != f['evidence']:
            lines.append('  首次异常数值: ' + json.dumps(f['onset_evidence'], ensure_ascii=False))
        if f['onset_advice'] != f['advice']:
            lines.append('  首次异常定位方向: ' + f['onset_advice'])
        coverage = ('起始时刻落在保留事件范围内（不保证整段覆盖或中间完整）' if f['event_paths'] else
                    '起始时刻未找到保留事件范围覆盖（不代表整段未采集），需核对写入/抑制/轮转' if f['capture_coverage'] == 'no_retained_overlap' else '起始时间证据不足')
        lines.append(paint('  故障起点数据覆盖: ' + coverage, '36' if f['event_paths'] else '33'))
        for clue in f.get('onset_context', []):
            lines.append(f'  起点附近日志（相关线索）: {clue["text"]}\n    位置: {clue["source"]}')
        for clue in recent_context_clues(f['context'], 2):
            lines.append(f'  关联日志: {clue["text"]}\n    位置: {clue["source"]}')
        lines.append('  定位方向: ' + f['advice'])
    if report.get('collector_totals'):
        lines.append('采集器累计统计: ' + json.dumps(report['collector_totals'], ensure_ascii=False))
        lines.append(f'所选时间线窗口内采集丢弃增量: {report.get("window_collector_drop_delta")}（采集丢弃不等于输入丢失）')
    lines += [paint('注意: ' + w, '33') for w in report['warnings']]
    lines.append(report['limitation'])
    return '\n'.join(lines)


def main(argv=None):
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', nargs='?', help='指定事件、运行目录或 rosbag 路径；省略则自动查找')
    parser.add_argument('--pid', type=int, help='选择指定的运行中定位实例')
    parser.add_argument('--root', help='无定位进程时使用的导出根目录')
    parser.add_argument('--json', action='store_true', help='输出结构化 JSON')
    parser.add_argument('--color', choices=('auto', 'always', 'never'), default='auto', help='终端颜色；auto 尊重 NO_COLOR，重定向不着色')
    parser.add_argument('--all', action='store_true', help='文本显示全部判据段；默认显示最近 12 段及全部分类汇总')
    parser.add_argument('--log-mb', type=int, default=32, help='最多检查日志尾部 MiB，默认 32')
    parser.add_argument('--bag-mb', type=int, default=128, help='裸 bag 最大扫描大小 MiB，默认 128')
    parser.add_argument('--match-evidence', help='指定失败匹配证据包；默认仅绘制所选运行/事件范围内最新一帧')
    parser.add_argument('--match-output-dir', help='俯视图输出目录；默认失败证据包内 birdseye')
    parser.add_argument('--no-match-image', action='store_true', help='跳过失败匹配俯视图，保留文字诊断')
    args = parser.parse_args(argv)
    if args.match_evidence:
        if args.path:
            parser.error('path 与 --match-evidence 二选一，避免混合两个证据范围')
        args.path = args.match_evidence
    if not 1 <= args.log_mb <= 1024 or not 1 <= args.bag_mb <= 8192:
        parser.error('log-mb 必须为 1..1024，bag-mb 必须为 1..8192')
    try:
        scheduling = low_impact()
        try:
            path, process, selection, warnings = select(args)
        except (OSError, ValueError) as exc:
            if args.path or args.pid:
                raise
            report = dict(schema=2, selected=None, process=None, selection='未找到可读取的定位包',
                          evidence_scope='unavailable', localization_error=str(exc), has_safety_evidence=False)
        else:
            report = analyze(path, process, selection, warnings, args.log_mb, args.bag_mb)
            report['evidence_scope'] = 'explicit' if args.path else 'live_run' if process else 'historical'
        report['match_failure_image'] = match_failure_image(report, args.match_evidence, args.match_output_dir, not args.no_match_image)
        report['upstream_lidar'] = (dict(status='not_requested', message='指定定位历史包，未拼接本机当前雷达结果；没有已关联的同次驱动证据。')
                                    if args.path else driver_diagnosis(args.root))
        report['scheduling'] = scheduling
        print(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) if args.json else render(report, color_enabled(args.color), args.all))
        has_driver = report['upstream_lidar'].get('report', {}).get('has_evidence') is True
        return 0 if report['has_safety_evidence'] or report.get('has_launch_evidence') or report.get('has_match_evidence') or has_driver else 2
    except (OSError, ValueError, RuntimeError) as exc:
        print(json.dumps(dict(error=str(exc)), ensure_ascii=False) if args.json else '无法完成诊断: ' + str(exc))
        return 2


if __name__ == '__main__':
    sys.exit(main())
