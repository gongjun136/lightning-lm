#!/usr/bin/env python3
"""Master-Orin acquisition: original LiDAR/CGI/CAN ROS streams, never lever arms.

The bag keeps source timestamps and raw typed CGI CAN payloads. Whole-bus CAN
errors/unknown IDs require a separate receive-only capture on the CGI host.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

import yaml

from accuracy_export_bag import FIELDS, PRIMARY_FIELDS, RAW_SUFFIXES, audit_clock, write_json
from resolve_sany_config import load, resolve

REPO = Path(__file__).resolve().parents[1]
ENV_KEYS = ('ROS_DOMAIN_ID', 'ROS_LOCALHOST_ONLY', 'RMW_IMPLEMENTATION',
            'FASTRTPS_DEFAULT_PROFILES_FILE', 'FASTRMW_DEFAULT_PROFILES_FILE',
            'CYCLONEDDS_URI', 'ROS_DISCOVERY_SERVER')
EVIDENCE = ['/slamPoseRaw_topic', '/slamPose_topic', '/PosRes', '/slamState_topic',
            '/SystemState', '/tf', '/tf_static', '/clock', '/rosout', '/parameter_events',
            '/Angle_sensor_topic', '/ImuRealAngle_topic', '/ImuRealPitch_V_topic',
            '/Imu_BodyAcce_topic', '/phi_topic', '/domain_vcu/status', '/safety/throstle']


def configure_transport_environment(environment, site_profile):
    """Reuse the Orin DDS profile even when only sensor nodes are running."""
    if environment.get('RMW_IMPLEMENTATION', 'rmw_fastrtps_cpp') != 'rmw_fastrtps_cpp':
        return
    profile = environment.get('FASTRTPS_DEFAULT_PROFILES_FILE') or environment.get('FASTRMW_DEFAULT_PROFILES_FILE')
    if not profile and site_profile.is_file():
        profile = str(site_profile)
    if profile:
        environment.setdefault('FASTRTPS_DEFAULT_PROFILES_FILE', profile)
        environment.setdefault('FASTRMW_DEFAULT_PROFILES_FILE', profile)


def active_localizers():
    """Read only this user's selected process settings, never arbitrary secrets."""
    result = []
    for p in Path('/proc').glob('[0-9]*/cmdline'):
        try:
            args = [s.decode() for s in p.read_bytes().split(b'\0') if s]
            if not args or Path(args[0]).name != 'run_loc_online':
                continue
            env = dict(s.split(b'=', 1) for s in p.with_name('environ').read_bytes().split(b'\0') if b'=' in s)
            config = next((a.split('=', 1)[1] for a in args if a.startswith('--config=')), None)
            if '--config' in args:
                config = args[args.index('--config') + 1]
            modules = sorted({line.split()[-1] for line in p.with_name('maps').read_text().splitlines()
                              if 'lightning' in line and '.so' in line and line.split()[-1].startswith('/')})
            result.append({'pid': int(p.parent.name), 'argv': args, 'config': config, 'loaded_modules': modules,
                           'environment': {k: env[k.encode()].decode() for k in ENV_KEYS if k.encode() in env}})
        except (OSError, UnicodeError, IndexError):
            continue
    return result


def make_plan(config, prefix, extras):
    multi = config.get('multi_lidar', {})
    inputs = list(multi.get('topics', {}).values()) if multi.get('enabled') else [
        config.get('common', {}).get('livox_lidar_topic') or config.get('common', {}).get('lidar_topic'),
        config.get('common', {}).get('imu_topic')]
    inputs = [t for t in inputs if t]
    if not inputs:
        raise ValueError('配置中没有雷达与 IMU 输入，不能创建空录制任务')
    required_inputs = list(inputs)
    reference = [prefix + '/' + field[0] for field in FIELDS.values()]
    reference += [prefix + '/time', prefix + '/driver_status']
    optional = [prefix + '/' + suffix for suffix in RAW_SUFFIXES if prefix + '/' + suffix not in reference]
    system = config.get('system', {})
    wheel = system.get('wheel_speed_topic', '/SpeThrCAN4_topic')
    if system.get('enable_wheel_speed_observation', False) or system.get('enable_wheel_speed_dr_observation', False):
        inputs.append(wheel)
        # The localization launcher explicitly permits absent wheel samples.
        # Preserve that contract; --extra-topic can make it mandatory for a run.
        optional.append(wheel)
    required = sorted(set(required_inputs + reference + extras))
    exact = sorted(set(required + EVIDENCE + [wheel]))
    if any(not re.fullmatch(r'/[A-Za-z0-9_/]+', t) for t in exact):
        raise ValueError('话题必须是绝对 ROS 名称，不能含空格或正则表达式')
    # Dynamic discovery also records CGI outputs/Livox streams that start later.
    regex = r'^(?:' + '|'.join([re.escape(t) for t in exact] + [
        re.escape(prefix) + r'/.*', r'/livox/.*', r'/sensors/livox/.*',
        r'/localization/.*', r'/diagnostics(?:/.*)?']) + r')$'
    return {'required_topics': required, 'configured_inputs': sorted(set(inputs)),
            'optional_input_topics': sorted(set(optional) - set(required)),
            'explicit_topics': exact, 'topic_regex': regex}


def quality_report(samples, now_ns, now_mono):
    """An engineering readiness report, not a calibrated accuracy certificate."""
    issues, values, stamps = [], {}, []
    for field, (_, _, attrs) in FIELDS.items():
        item = samples.get(field)
        if not item:
            issues.append(field + ':missing')
            continue
        message, received, stamp = item
        stamps.append(stamp)
        raw_values = {a: getattr(message, a) for a in attrs}
        values[field] = {a: v if math.isfinite(float(v)) else None for a, v in raw_values.items()}
        if not message.frame.valid:
            issues.append(field + ':invalid')
        if now_mono - received > .15 or now_ns - stamp > 150_000_000 or stamp <= 0:
            issues.append(field + ':stale')
        if stamp - now_ns > 20_000_000:
            issues.append(field + ':future')
        if any(v is None for v in values[field].values()):
            issues.append(field + ':nonfinite')
    status = values.get('status', {})
    fixed = status.get('system_state') == 2 and status.get('satellite_status') == 4
    if not fixed:
        issues.append('status:not_2_4')
    age = status.get('differential_age_s', math.inf)
    if age is None or not math.isfinite(age) or not 0 <= age <= 2:
        issues.append('differential_age')
    limits = {'position_sigma': (.1, .1, .2), 'attitude_sigma': (.5, .5, .5),
              'velocity_sigma': (.2, .2, .2)}
    for field, bounds in limits.items():
        if field in values and any(v is None or not math.isfinite(v) or v < 0 or v > b
                                  for v, b in zip(values[field].values(), bounds)):
            issues.append(field + ':limit')
    if stamps and max(stamps) - min(stamps) > 20_000_000:
        issues.append('field_skew')
    lat, lon = values.get('latitude', {}).get('latitude_deg'), values.get('longitude', {}).get('longitude_deg')
    if lat is not None and not -90 <= lat <= 90 or lon is not None and not -180 <= lon <= 180:
        issues.append('llh:range')
    return {'received_at_ns': now_ns, 'device_fixed_flag': fixed, 'instant_quality_ok': not issues,
            'issues': issues, 'fields': values, 'header_stamps_ns': {k: v[2] for k, v in samples.items()}}


class Monitor:
    def __init__(self, node, plan, prefix):
        self.node, self.plan, self.prefix = node, plan, prefix
        self.counts, self.last, self.samples, self.subs, self.errors = Counter(), {}, {}, {}, {}
        self.versions, self.good = {}, 0
        self.maximum_good = 0
        self.matcher = re.compile(plan['topic_regex'])
        self.graph = {}

    def discover(self):
        from rclpy.qos import QoSProfile, ReliabilityPolicy, DurabilityPolicy
        from rosidl_runtime_py.utilities import get_message
        for topic, types in self.node.get_topic_names_and_types():
            if not self.matcher.fullmatch(topic):
                continue
            self.graph[topic] = types
            if topic in self.subs:
                continue
            if len(types) != 1:
                self.errors[topic] = '话题存在多个类型：' + ','.join(types)
                continue
            try:
                msg_type = get_message(types[0])
                qos = QoSProfile(depth=16, reliability=ReliabilityPolicy.BEST_EFFORT,
                                 durability=DurabilityPolicy.VOLATILE)
                if topic == '/tf_static':
                    qos.reliability = ReliabilityPolicy.RELIABLE
                    qos.durability = DurabilityPolicy.TRANSIENT_LOCAL
                self.subs[topic] = self.node.create_subscription(
                    msg_type, topic, lambda data, t=topic, typ=msg_type: self.receive(t, typ, data), qos, raw=True)
            except (ImportError, AttributeError, RuntimeError, ValueError) as exc:
                self.errors[topic] = str(exc)

    def receive(self, topic, msg_type, data):
        from rclpy.serialization import deserialize_message
        at = time.monotonic()
        self.counts[topic] += 1
        self.last[topic] = at
        field = next((k for k, (suffix, _, _) in FIELDS.items() if topic == self.prefix + '/' + suffix), None)
        if not field:
            return
        msg = deserialize_message(data, msg_type)
        stamp = msg.frame.header.stamp.sec * 1_000_000_000 + msg.frame.header.stamp.nanosec
        if field in self.samples and stamp <= self.samples[field][2]:
            self.good = 0
            return
        self.samples[field] = (msg, at, stamp)
        report = quality_report(self.samples, self.node.get_clock().now().nanoseconds, at)
        if any(':invalid' in issue or ':nonfinite' in issue or issue == 'status:not_2_4' for issue in report['issues']):
            self.good = 0
        if all(k in self.samples and self.samples[k][2] > self.versions.get(k, 0) for k in PRIMARY_FIELDS):
            self.good = self.good + 1 if report['instant_quality_ok'] else 0
            self.maximum_good = max(self.maximum_good, self.good)
            self.versions = {k: self.samples[k][2] for k in PRIMARY_FIELDS}

    def report(self):
        now = time.monotonic()
        quality = quality_report(self.samples, self.node.get_clock().now().nanoseconds, now)
        if any(':stale' in issue for issue in quality['issues']):
            self.good = 0
        quality.update(consecutive_good=self.good, maximum_consecutive_good=self.maximum_good,
                       reference_ready=self.good >= 10 and quality['instant_quality_ok'])
        return {'message_counts': dict(self.counts), 'missing_optional_inputs': [t for t in self.plan['optional_input_topics']
                if t not in self.last or now - self.last[t] > 2],
                'missing_required': [t for t in self.plan['required_topics']
                if t not in self.last or now - self.last[t] > 2], 'quality': quality,
                'topic_types': self.graph, 'subscription_errors': self.errors}


def snapshot(out, config_path, config, plan, active):
    files = [config_path, REPO / 'workspace.repos', REPO / 'config/sany/calibrations.yaml',
             Path(os.environ.get('SANY_VEHICLE_CONFIG', '/opt/conf/vehicle_config.yaml'))]
    manifest = config.get('map_package', {}).get('manifest')
    if manifest:
        files.append(Path(manifest))
    for proc in active:
        if proc.get('config'):
            p = Path(proc['config'])
            if p.resolve() != config_path:
                files.append(p)
    for key in ENV_KEYS:
        if key.endswith('PROFILES_FILE') and os.environ.get(key):
            files.append(Path(os.environ[key]))
    copied = []
    (out / 'snapshots').mkdir()
    for i, p in enumerate(dict.fromkeys(files)):
        if p.is_file():
            data = p.read_bytes()
            name = f'{i:02d}_{p.name}'
            (out / 'snapshots' / name).write_bytes(data)
            copied.append({'source': str(p), 'snapshot': 'snapshots/' + name, 'sha256': hashlib.sha256(data).hexdigest()})
    (out / 'config.yaml').write_text(yaml.safe_dump(config, sort_keys=False), encoding='utf-8')
    git = subprocess.run(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], capture_output=True, text=True)
    dirty = subprocess.run(['git', '-C', str(REPO), 'status', '--porcelain'], capture_output=True, text=True)
    return {'schema_version': 1, 'role': 'master_sensors', 'source_recording_id': out.name,
            'epoch': out.name, 'time_basis': None,
            'physical_reference_point': 'unverified_device_output', 'online_lever_arm_applied': False,
            'raw_can_bus_capture': False, 'environment': {k: os.environ.get(k) for k in ENV_KEYS},
            'git_head': git.stdout.strip(), 'git_status': dirty.stdout, 'active_localizers': active,
            'snapshots': copied, 'plan': plan, 'clock_start': audit_clock()}


def main(argv=None):
    parser = argparse.ArgumentParser(description='主 Orin 原始传感器录制；保留设备点，离线补杆臂')
    parser.add_argument('--output', type=Path, help='新的输出目录，不能已存在')
    parser.add_argument('--config', type=Path, help='定位实际使用的 YAML；默认从运行进程读取')
    parser.add_argument('--cgi-prefix', default='/cgi430')
    parser.add_argument('--extra-topic', action='append', default=[])
    parser.add_argument('--duration-seconds', type=float, default=0, help='录制器就绪后的秒数，0 表示 Ctrl-C 结束')
    parser.add_argument('--wait-seconds', type=float, default=20, help='等待输入数据的秒数')
    parser.add_argument('--check-only', action='store_true')
    parser.add_argument('--require-fixed', action='store_true', help='开始前要求连续 10 条通过工程门限的新解')
    args = parser.parse_args(argv)
    if args.duration_seconds < 0 or not 0 < args.wait_seconds <= 300:
        parser.error('录制时长须非负，等待时长须在 (0, 300] 秒内')
    if not re.fullmatch(r'/[A-Za-z0-9_/]+', args.cgi_prefix) or args.cgi_prefix.endswith('/'):
        parser.error('CGI 前缀必须是无尾斜线的绝对命名空间')
    active = active_localizers()
    if args.config:
        # An explicit offline/test configuration must not inherit another
        # locator's domain or transport settings merely because it is running.
        active = [p for p in active if p['config'] and Path(p['config']).resolve() == args.config.resolve()]
    if len(active) > 1:
        raise ValueError('发现多个定位进程，请用 --config 明确被评定位器配置')
    for k, v in (active[0]['environment'].items() if active else []):
        if k in os.environ and os.environ[k] != v:
            raise ValueError(f'{k} 与运行定位器不一致：当前 {os.environ[k]}，定位器 {v}')
        os.environ.setdefault(k, v)
    os.environ.setdefault('ROS_DOMAIN_ID', '42')
    configure_transport_environment(os.environ, Path.home() / 'Documents/ros2/fastdds.xml')
    config_path = (args.config or (Path(active[0]['config']) if active and active[0]['config'] else None)
                   or REPO / 'config/sany/localization.yaml').resolve()
    config = load(config_path)
    if 'sany_vehicle_config' in config:
        config = resolve(config_path, os.environ.get('SANY_VEHICLE_CONFIG', '/opt/conf/vehicle_config.yaml'))
    plan = make_plan(config, args.cgi_prefix, args.extra_topic)
    out = args.output or (Path(os.environ.get('LIGHTNING_LM_OUT_ROOT', REPO.parent.parent / 'runs')) /
                          time.strftime('sensors_%Y%m%d_%H%M%S'))
    out = out.resolve()
    if not args.check_only or args.output:
        out.mkdir(parents=True, exist_ok=False)
        metadata = snapshot(out, config_path, config, plan, active)
        write_json(out / 'capture.metadata.json', metadata)
    else:
        metadata = {}
    import rclpy
    from rclpy.signals import SignalHandlerOptions
    rclpy.init(signal_handler_options=SignalHandlerOptions.NO)
    node = rclpy.create_node('localization_sensor_recorder_' + str(os.getpid()))
    monitor = Monitor(node, plan, args.cgi_prefix)
    stopped = False

    def stop(_sig, _frame):
        nonlocal stopped
        stopped = True

    for sig in (signal.SIGINT, signal.SIGTERM):
        signal.signal(sig, stop)
    proc, log, audit = None, None, None
    exit_code, forced = 2, False
    try:
        print(f'检查实际输入，ROS_DOMAIN_ID={os.environ["ROS_DOMAIN_ID"]}，配置={config_path}', flush=True)
        deadline, discovery = time.monotonic() + args.wait_seconds, 0
        while not stopped and time.monotonic() < deadline:
            now = time.monotonic()
            if now >= discovery:
                monitor.discover()
                discovery = now + 1
            rclpy.spin_once(node, timeout_sec=.05)
            report = monitor.report()
            if not report['missing_required'] and (not args.require_fixed or report['quality']['reference_ready']):
                break
        report = monitor.report()
        print('采集前检查：' + json.dumps(report, ensure_ascii=False), flush=True)
        metadata['preflight'] = report
        if report['missing_required'] or (args.require_fixed and not report['quality']['reference_ready']) or stopped:
            print('输入或固定解门禁未通过；请核对 CGI、雷达、轮速和 DDS，未开始录包。', file=sys.stderr)
            return 2
        if args.check_only:
            exit_code = 0
            return 0
        qos = {topic: {'reliability': 'best_effort', 'durability': 'volatile', 'history': 'keep_last', 'depth': 64}
               for topic in set(plan['explicit_topics']) | set(monitor.graph)}
        qos['/tf_static'].update(reliability='reliable', durability='transient_local')
        (out / 'qos.yaml').write_text(yaml.safe_dump(qos, sort_keys=False), encoding='utf-8')
        log = (out / 'recorder.log').open('w', encoding='utf-8')
        audit = (out / 'quality_and_inputs.jsonl').open('w', encoding='utf-8')
        cmd = ['ros2', 'bag', 'record', '--storage', 'mcap', '--output', str(out / 'bag'),
               '--qos-profile-overrides-path', str(out / 'qos.yaml'), '--regex', plan['topic_regex'],
               '--max-bag-size', str(4 * 1024**3), '--max-cache-size', str(256 * 1024**2)]
        metadata.update(recorder_command=cmd, capture_start_ns=time.time_ns(), capture_start_monotonic_ns=time.monotonic_ns())
        write_json(out / 'capture.metadata.json', metadata)
        proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        start, next_audit, next_clock = time.monotonic(), 0, time.monotonic() + 30
        ready = None
        print(f'正在启动录制器：{out}；按一次 Ctrl-C 后等待写盘。', flush=True)
        exit_code = 0
        while not stopped:
            rclpy.spin_once(node, timeout_sec=.05)
            now = time.monotonic()
            if proc.poll() is not None:
                exit_code = 3
                print('录包进程意外退出，请检查 recorder.log。', file=sys.stderr)
                break
            if ready is None:
                # Humble logs this only after constructing the recorder. Start
                # finite runs here, so discovery/startup do not consume them.
                if 'Recording...' in (out / 'recorder.log').read_text(encoding='utf-8', errors='replace'):
                    ready = now
                    metadata['recorder_ready_ns'] = time.time_ns()
                    write_json(out / 'capture.metadata.json', metadata)
                    print('录制器已就绪；正在保存原始数据。', flush=True)
                elif now - start > 30:
                    print('录制器 30 秒内未就绪，请检查 recorder.log。', file=sys.stderr)
                    exit_code = 3
                    break
            if now >= next_audit:
                monitor.discover()
                audit.write(json.dumps(monitor.report(), ensure_ascii=False, allow_nan=False) + '\n')
                audit.flush()
                next_audit = now + 1
                if os.statvfs(out).f_bavail * os.statvfs(out).f_frsize < 2 * 1024**3:
                    print('剩余空间不足 2 GiB，停止录制并保存已有数据。', file=sys.stderr)
                    exit_code = 3
                    break
            if now >= next_clock:
                with (out / 'clock.jsonl').open('a', encoding='utf-8') as f:
                    f.write(json.dumps(audit_clock(), ensure_ascii=False) + '\n')
                next_clock = now + 30
            if args.duration_seconds and ready is not None and now - ready >= args.duration_seconds:
                break
    finally:
        if proc and proc.poll() is None:
            os.killpg(proc.pid, signal.SIGINT)
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                forced = True
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
        if audit:
            audit.close()
        if log:
            log.close()
        if proc:
            metadata.update(recorder_returncode=proc.returncode, forced_shutdown=forced,
                            capture_end_ns=time.time_ns(), capture_end_monotonic_ns=time.monotonic_ns())
            path = out / 'bag/metadata.yaml'
            bag_meta = yaml.safe_load(path.read_text()) if path.is_file() else {}
            counts = {t['topic_metadata']['name']: t['message_count']
                      for t in bag_meta.get('rosbag2_bagfile_information', {}).get('topics_with_message_count', [])}
            missing = [t for t in plan['required_topics'] if not counts.get(t)]
            metadata.update(bag_message_counts=counts, missing_required_in_bag=missing)
            metadata['discovered_topics_without_bag_messages'] = sorted(t for t in monitor.graph if not counts.get(t))
            if forced or proc.returncode not in (0, 130, -signal.SIGINT) or missing:
                exit_code = 3
            print(f'录制结束，返回码={exit_code}，必需话题缺失={missing}，目录={out}', flush=True)
        metadata.update(exit_code=exit_code, final_observation=monitor.report())
        if out.is_dir() and metadata:
            metadata['clock_end'] = audit_clock()
            write_json(out / 'capture.metadata.json', metadata)
        node.destroy_node()
        rclpy.shutdown()
    return exit_code


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, OSError, yaml.YAMLError) as exc:
        print(f'录制失败：{exc}', file=sys.stderr)
        sys.exit(2)
