#!/usr/bin/env python3
"""Fail closed when reviewed Safety/recorder sources change without diagnosis review."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'scripts/diagnose_contract.json'
SOURCES = (
    'src/core/localization/localization_dtc_monitor.cc',
    'src/core/localization/localization_dtc_monitor.h',
    'src/utils/radar_dtc_state.h',
    'src/utils/sany_radar_topology.h',
    'src/utils/localization_dtc_reporting.h',
    'src/utils/incident_recorder.cc',
    'src/utils/incident_recorder.h',
    'src/core/system/loc_system.cc#PublishSafety',
    'src/core/system/loc_system.cc#DtcLocalizationLost',
    'src/core/system/loc_system.cc#DtcInitialization',
    'config/sany/localization.yaml#diagnostics',
    'workspace.repos',
)


def snapshot(root=ROOT):
    result = {}
    for name in SOURCES:
        path, _, section = name.partition('#')
        text = (root / path).read_text(encoding='utf-8').replace('\r\n', '\n')
        if section in ('PublishSafety', 'DtcLocalizationLost'):
            declaration = 'void LocSystem::PublishSafety()' if section == 'PublishSafety' else 'bool LocSystem::DtcLocalizationLost()'
            start = text.index(declaration)
            text = text[start:text.index('\n}\n', start) + 3]
        elif section == 'DtcInitialization':
            start = text.index('    const YAML::Node dtc = root["diagnostics"];')
            text = text[start:text.index('    const auto dtc_qos = ', start)]
        elif section == 'diagnostics':
            lines = text.splitlines(keepends=True)
            start = next(i for i, line in enumerate(lines) if line.startswith('diagnostics:'))
            end = next((i for i in range(start + 1, len(lines))
                        if lines[i].strip() and not lines[i][0].isspace() and not lines[i].startswith('#')), len(lines))
            text = ''.join(lines[start:end])
        result[name] = hashlib.sha256(text.encode()).hexdigest()
    identity = hashlib.sha256(json.dumps(result, sort_keys=True).encode()).hexdigest()
    return dict(schema=1, contract_id=identity, reviewed_sources=result)


def check(root=ROOT, manifest=MANIFEST):
    actual = snapshot(root)
    expected = json.loads(manifest.read_text(encoding='utf-8'))
    changed = [p for p, digest in actual['reviewed_sources'].items()
               if expected.get('reviewed_sources', {}).get(p) != digest]
    if changed or actual != expected:
        raise ValueError('诊断合同未同步: ' + ', '.join(changed or ['manifest']) +
                         '\n请核对诊断解析、判据和回归用例，再运行 check_diagnostic_contract.py --update；不要只更新哈希跳过审查。')
    return actual['contract_id']


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--update', action='store_true', help='完成判据与测试审查后更新已核对源文件快照')
    group.add_argument('--print-id', action='store_true', help='检查后仅输出合同 ID，供构建标记使用')
    group.add_argument('--print-dependencies', action='store_true', help='输出 CMake 重新配置依赖')
    args = parser.parse_args()
    try:
        if args.print_dependencies:
            print(';'.join(str(ROOT / name.split('#')[0]) for name in SOURCES))
        elif args.update:
            MANIFEST.write_text(json.dumps(snapshot(), indent=2) + '\n', encoding='utf-8')
            print('诊断合同已更新；仍需运行诊断回归测试。')
        else:
            identity = check()
            print(identity if args.print_id else 'diagnostic contract OK: ' + identity)
        return 0
    except (OSError, ValueError, KeyError, StopIteration) as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
