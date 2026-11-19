'''定位诊断证据摘要层（纯计算）。

只格式化与汇总 diagnose_localization 生成的 report 字典：不读写文件、不访问网络、
不启动子进程、不创建 ROS 节点或订阅；判据与可信等级沿用 classify() 的结论，
不重算安全逻辑；不修改输入；同一输入输出确定。
'''
from translation_diagnostics import describe_translation
from diagnostic_history_layout import render_history_sections, history_line_style

import json
import math

SUMMARY_VERSION = 1
FRESH_LIMIT_SEC = 5.0
HISTORY_LIMIT = 6
NL = chr(10)
KNOWN_CODES = (0, 2, 4, -9)
RUN_KEYS = (('raw', '原始码'), ('strict', '严格码'), ('selected', '选择码'), ('published', '上报码'))
DELAY_CANDIDATES = ('驱动/SDK 缓存', 'publish 内部阻塞', 'DDS 传输/接收队列', '执行器调度')
MISSING_STAMPS = ('publish 进入/返回时间', 'DDS source/received timestamp', '回调进入时间', '各时间戳的时钟域与跨节点对齐')
GRACE_KEYS = ('match_rejected_grace_suppressed', 'match_failure_grace_inconsistent')
PHASE_TEXT = {'startup': '启动阶段', 'runtime': '运行阶段'}
STATE_TEXT = {
    'open': '该证据范围内未见结束（不代表当前仍在故障）',
    'recovered': '已观测到原始码回 0',
    'predicate_changed': '判据发生变化（未必恢复）',
    'observation_gap': '观测中断（不能判为恢复）',
}
COVERAGE_TEXT = {
    'range_overlap': '起始时刻落在保留采集范围内（不保证整段完整）',
    'no_retained_overlap': '起始时刻未找到保留采集覆盖（不代表整段未采集）',
    'unknown_time': '起始时间未知，无法判断采集覆盖',
}
# 判据 key -> [(标签, 实测字段, 阈值字段, 单位, 换算)]；只格式化已有字段，不重算判据。
SPECS = {
    'frame_supply': [('雷达帧年龄', 'lidar_frame_age_sec', 'radar_frame_timeout_sec', 's', None),
                     ('恢复确认帧', 'lidar_frame_recovery_count', 'radar_frame_recovery_frames', '帧', None)],
    'match_stale_degraded': [('有效匹配年龄', 'valid_match_age_sec', 'localization_output_max_lidar_age_sec', 's', None)],
    'output_stale_degraded': [('输出年龄', 'output_age_sec', 'output_timeout_sec', 's', None)],
    'startup_pending': [('启动已用时', 'startup_age_sec', 'startup_grace_sec', 's', None)],
    'confirmed_tracking_lost': [('连续丢失帧', 'lost_frames', 'lost_frame_threshold', '帧', None)],
    'match_rejected_grace_suppressed': [('匹配失败容忍已用时', 'match_failure_grace_age_sec', 'match_failure_grace_sec', 's', None)],
    'major_sigma': [('主轴 sigma', 'major_sigma', 'major_sigma_limit_m', 'm', None)],
    'lateral_sigma': [('横向 sigma', 'lateral_sigma', 'lateral_sigma_limit_m', 'm', None)],
    'yaw_sigma_rad': [('航向 sigma', 'yaw_sigma_rad', 'yaw_sigma_limit_deg', 'deg', 'rad2deg')],
    'yaw_condition': [('发布航向误差', 'yaw_published_error_rad', 'yaw_limit_rad', 'rad', None),
                      ('原始航向误差', 'yaw_raw_error_rad', 'yaw_limit_rad', 'rad', None)],
}
SPECS['output_stale'] = SPECS['output_stale_degraded']
SPECS['startup_expired'] = SPECS['startup_pending']
SPECS['match_rejected_unconfirmed'] = SPECS['confirmed_tracking_lost']
SPECS['match_failure_grace_inconsistent'] = SPECS['match_rejected_grace_suppressed']
FLAGS = {
    'frame_supply': ['lidar_frame_fault', 'lidar_frame_seen'],
    'match_stale_degraded': ['match_stale'],
    'output_stale_degraded': ['first_valid_output_missing'],
    'startup_pending': ['startup_pending', 'startup_expired'],
    'confirmed_tracking_lost': ['confirmed_tracking_lost', 'localization_lost'],
    'match_rejected_unconfirmed': ['latest_match_valid', 'output_gate_open'],
    'match_rejected_grace_suppressed': ['raw_runstate', 'selected_runstate', 'published_runstate', 'output_gate_open',
                                        'latest_match_valid', 'match_failure_grace_active', 'match_failure_grace_suppressed'],
    'major_sigma': ['covariance_valid'],
    'yaw_condition': ['yaw_evidence_available', 'error_yaw_deg'],
    'translation_condition': ['translation_condition', 'translation_failure_latched',
                              'translation_evidence_origin', 'translation_selected_window'],
}
for _alias, _base in (('output_stale', 'output_stale_degraded'), ('startup_expired', 'startup_pending'),
                      ('match_failure_grace_inconsistent', 'match_rejected_grace_suppressed'),
                      ('lateral_sigma', 'major_sigma'), ('yaw_sigma_rad', 'major_sigma')):
    FLAGS[_alias] = FLAGS[_base]


def d(value):
    return value if isinstance(value, dict) else {}


def lst(value):
    return list(value) if isinstance(value, (list, tuple)) else []


def num(value):
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        try:
            return float(value)
        except OverflowError:
            return None
    if isinstance(value, str):
        try:
            return float(value.strip())
        except ValueError:
            return None
    return None


def finite(value):
    v = num(value)
    return v if v is not None and math.isfinite(v) else None


def _val(value, default=-1.0):
    v = finite(value)
    return default if v is None else v


def txt(value):
    return '未记录' if value is None or value == '' else str(value)


def fmt(value, unit=''):
    if value is None:
        return '未记录'
    if isinstance(value, bool):
        return 'true' if value else 'false'
    v = num(value)
    if v is None:
        return str(value)
    if math.isnan(v):
        return '非有限值(nan)'
    if math.isinf(v):
        return '+inf' if v > 0 else '-inf'
    if isinstance(value, int) or (isinstance(value, str) and value.strip().lstrip('+-').isdigit()):
        text = str(int(value))
    elif v == int(v) and abs(v) < 1e15:
        text = str(int(v))
    else:
        text = '%.10g' % v
    return text + (' ' + unit if unit else '')


def clean(obj):
    '''生成可 JSON 序列化（allow_nan=False）的新结构，不改动原对象。'''
    if isinstance(obj, dict):
        return {str(k): clean(obj[k]) for k in sorted(obj, key=str)}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    if isinstance(obj, float) and not math.isfinite(obj):
        return fmt(obj)
    if obj is None or isinstance(obj, (bool, int, float, str)):
        return obj
    return str(obj)


def truthy(value):
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.strip().lower() in ('true', 'yes', 'on'):
        return True
    v = finite(value)
    return v is not None and v != 0


def present(value):
    if value is None or value is False or value in ('', [], {}):
        return False
    v = num(value)
    return not (v is not None and v == 0)


def short(value, limit=120):
    if isinstance(value, (dict, list, tuple)):
        text = json.dumps(clean(value), ensure_ascii=False, sort_keys=True)
        return text if len(text) <= limit else text[:limit] + '…'
    return fmt(value)


def brief(obj, limit=4):
    if not isinstance(obj, dict):
        return short(obj)
    keys = sorted(obj, key=str)
    parts = [str(k) + '=' + short(obj[k], 60) for k in keys[:limit]]
    if len(keys) > limit:
        parts.append('另 %d 项见索引' % (len(keys) - limit))
    return '，'.join(parts) if parts else '无'


def code_of(value):
    v = finite(value)
    return int(v) if v is not None and v == int(v) else None


def code_text(value):
    c = code_of(value)
    if c is None:
        return '未知码(' + fmt(value) + ')'
    return str(c) if c in KNOWN_CODES else '未知码 ' + str(c)


def freshness(report):
    report = d(report)
    scope = report.get('evidence_scope')
    latest = d(report.get('latest'))
    raw_age = report.get('latest_age_sec')
    age = num(raw_age)
    err = report.get('localization_error')
    if not report.get('has_safety_evidence'):
        state = 'unavailable'
        text = '没有 Safety 证据（范围 ' + txt(scope) + '），状态未知' + ('：' + str(err) if err else '')
    elif scope in ('explicit', 'historical'):
        state, text = 'offline', '离线/历史证据：以下只描述包内末条记录，不代表当前状态'
        if scope == 'historical' and not report.get('process'):
            text = '当前未发现定位进程，包内状态码不代表当前状态；' + text
    elif scope != 'live_run':
        state, text = 'unknown', '证据范围为 ' + txt(scope) + '，当前未知'
    elif not report.get('process'):
        state, text = 'unknown', '未确认定位进程在运行，当前未知'
    elif d(report.get('live_status_refresh')).get('status') == 'failed':
        state = 'unknown'
        text = '分析后实时补读失败：' + txt(d(report.get('live_status_refresh')).get('reason')) + '，当前未知'
    elif age is None:
        state, text = 'unknown', '最新记录缺少可用时间，当前未知'
    elif not math.isfinite(age):
        state, text = 'unknown', '最新记录年龄为非有限值（' + fmt(raw_age) + '），当前未知'
    elif age < 0:
        state, text = 'clock_anomaly', '最新记录时间比本机晚 %g 秒（未来时间/时钟异常），当前未知' % -age
    elif age > FRESH_LIMIT_SEC:
        state, text = 'stale', '最新记录已过去 %g 秒，超过 %g 秒新鲜度上限，当前未知' % (age, FRESH_LIMIT_SEC)
    else:
        state, text = 'current', '当前记录：距今 %g 秒（上限 %g 秒）' % (age, FRESH_LIMIT_SEC)
    if state == 'current' and d(report.get('live_status_refresh')).get('status') == 'refreshed':
        text += '（分析结束后对同一实例有界补读所得；历史段来自更早的分析快照）'
    return {'state': state, 'is_current': state == 'current', 'text': text, 'limit_sec': FRESH_LIMIT_SEC,
            'age_sec': age if age is not None and math.isfinite(age) else None,
            'age_raw': None if raw_age is None else fmt(raw_age), 'scope': scope,
            'time': latest.get('_time'), 'source': latest.get('_source')}


def is_fresh_live(report):
    return freshness(report)['is_current']


def runstate(latest, finding_count):
    values = {k: latest.get(k + '_runstate') for k, _ in RUN_KEYS}
    codes = {k: code_of(v) for k, v in values.items()}
    parts = [label + ' ' + ('未记录' if values[k] is None else code_text(values[k])) for k, label in RUN_KEYS]
    unknown = [label for k, label in RUN_KEYS if values[k] is not None and codes[k] not in KNOWN_CODES]
    raw, pub = codes['raw'], codes['published']
    suppressed = truthy(latest.get('match_failure_grace_suppressed'))
    missing = [label for k, label in RUN_KEYS if values[k] is None]
    zeros = '各运行码为 0' if not missing else '已记录的运行码为 0（' + '、'.join(missing) + '未记录，不推测）'
    if all(v is None for v in values.values()):
        verdict = '运行码未记录，状态未知'
    elif unknown:
        verdict = '存在未知运行码（' + '、'.join(unknown) + '），规则缺失，不能判为正常'
    elif values['raw'] is None:
        verdict = '缺少原始码，不能判为正常'
    elif raw != 0 and pub == 0:
        verdict = '原始码 %d 而上报码 0：判据被容忍/抑制，不是健康' % raw
    elif raw != 0:
        verdict = '该记录为故障状态：原始码 %d，上报码 %s' % (raw, txt(pub))
    elif any(c not in (0, None) for c in codes.values()):
        verdict = '原始码 0 但其他运行码非 0（选择/锁存差异），不能判为正常'
    elif suppressed:
        verdict = zeros + '，但记录 match_failure_grace_suppressed=true，容忍期间不是健康'
    elif finding_count:
        verdict = zeros + '，但该记录仍有 %d 项判据/线索（见下），不能判为健康' % finding_count
    elif missing:
        verdict = zeros + '；不能确认上报状态，不判为健康'
    else:
        verdict = '该条 Safety 记录各运行码为 0（只说明该记录为 0，不保证定位精度或下游已收到）'
    if suppressed and '容忍' not in verdict:
        verdict += '；记录 match_failure_grace_suppressed=true'
    severity = ('normal' if not missing and not unknown and all(c == 0 for c in codes.values())
                and not suppressed and not finding_count else
                'fault' if pub in (2, 4, -9) or (raw in (2, 4, -9) and pub != 0) else 'attention')
    return {'values': clean(values), 'parts': parts, 'verdict': verdict, 'severity': severity}


def threshold(name, evidence, thresholds):
    if name is None:
        return None, '本判据无阈值字段'
    if name in evidence:
        value, origin = evidence.get(name), '证据'
    elif thresholds.get(name) is not None:
        value, origin = thresholds.get(name), '有效配置'
    else:
        return None, '阈值未记录'
    if finite(value) is None:
        return None, '%s中阈值为 %s，不可验证' % (origin, fmt(value))
    return value, origin


def measures(key, evidence, thresholds):
    '''列出实测与实际阈值；不比较大小（恢复确认阶段等情形下比较会误导）。'''
    if key == 'translation_condition':
        return describe_translation(evidence)['lines']
    out = []
    for label, vkey, tkey, unit, conv in SPECS.get(key, ()):
        value = evidence.get(vkey)
        limit, origin = threshold(tkey, evidence, thresholds)
        if value is None and limit is None:
            continue
        if conv == 'rad2deg' and finite(value) is not None:
            shown = fmt(math.degrees(finite(value)), unit) + '（原值 ' + fmt(value, 'rad') + '）'
        else:
            shown = fmt(value, unit if conv is None else 'rad')
        limit_text = origin if limit is None else '阈值 ' + fmt(limit, unit) + '（' + origin + '）'
        out.append(label + ' 实测 ' + shown + ' / ' + limit_text)
    return out


def flags(key, evidence):
    names = FLAGS[key] if key in FLAGS else sorted(evidence, key=str)[:6]
    return [str(n) + '=' + short(evidence[n], 60) for n in names if n in evidence]


def notes(key, evidence, thresholds, finding):
    out = []
    if key == 'translation_condition':
        out.extend(describe_translation(evidence)['warnings'])
    if key == 'frame_supply':
        out.append('恢复确认阶段年龄可能已低于阈值，按记录展示，不推断大小关系')
    if key in ('output_stale', 'output_stale_degraded'):
        age = num(evidence.get('output_age_sec'))
        if truthy(evidence.get('first_valid_output_missing')) or (age is not None and math.isinf(age)):
            out.append('可能尚无首次有效输出')
    if key == 'yaw_condition' and threshold('yaw_limit_rad', evidence, thresholds)[0] is None:
        out.append('动态航向阈值未记录；yaw_tolerance_deg 不是最终阈值')
    if key in GRACE_KEYS:
        out.append('容忍时长以记录/配置为准（200ms 仅为最大容忍）；容忍期间不是健康')
    if key not in SPECS and key not in FLAGS:
        out.append('该判据无专用格式，按 classify 结论原样展示')
    if finding.get('confidence') == '规则缺失' or code_of(finding.get('code')) not in KNOWN_CODES:
        out.append('规则缺失或未知码，不能视为正常')
    if finding.get('phase') == 'startup':
        out.append('启动阶段判据，不等于无害（可能触发停车）')
    return out


def rejection_view(mr):
    ev = d(mr.get('evidence'))
    items = ['拒绝项 ' + fmt(mr.get('reason_code')) + (' ' + txt(mr.get('title')) if mr.get('title') else '')]
    for label, name in (('迭代数', 'iterations'), ('NDT收敛', 'ndt_converged'),
                        ('迭代预算耗尽', 'ndt_budget_exhausted'), ('数值失败', 'ndt_numerical_failure')):
        if name in ev:
            items.append(label + '=' + fmt(ev.get(name)))
    if mr.get('confidence') is not None or 'min_confidence' in ev:
        items.append('匹配分数=' + fmt(mr.get('confidence')) + ' / 门槛=' + fmt(ev.get('min_confidence')))
    if 'consecutive_rejected_matches' in ev:
        items.append('连续拒绝=' + fmt(ev.get('consecutive_rejected_matches')) + ' / 丢失阈值=' + fmt(ev.get('lost_frame_threshold')))
    if ev.get('map_overlap_check_enabled') in (0, 1):
        items.append('重合率检查=' + ('开启' if ev['map_overlap_check_enabled'] else '关闭（仅统计）'))
    if 'map_overlap' in ev:
        items.append('地图重合率=' + fmt(ev.get('map_overlap')) + ' / 配置阈值=' + fmt(ev.get('map_min_overlap')))
    if ev.get('map_checked') == 1 and 'map_gravity_cos' in ev:
        measured, minimum = finite(ev.get('map_gravity_cos')), finite(ev.get('map_min_gravity_cos'))
        items.append('重力余弦=' + fmt(ev.get('map_gravity_cos')) + ' / 实际下限=' + fmt(ev.get('map_min_gravity_cos')))
        if measured is not None and minimum is not None and -1 <= measured <= 1 and -1 <= minimum <= 1:
            items.append('地图/里程计竖直轴夹角=' + fmt(math.degrees(math.acos(measured)), '°') +
                         ' / 实际上限=' + fmt(math.degrees(math.acos(minimum)), '°'))
    return {'text': '，'.join(items), 'time': mr.get('time'), 'source': mr.get('source'),
            'relation': mr.get('relation'), 'snapshots': len(lst(mr.get('snapshots'))), 'evidence': clean(ev),
            'note': '匹配分数不是可信等级；迭代数/耗时不单独说明失败原因'}


def finding_view(finding, thresholds, ref):
    f = d(finding)
    ev = d(f.get('evidence'))
    key = f.get('key') or '未知判据'
    mr = f.get('matching_rejection')
    return {'ref': ref, 'code': code_text(f.get('code')), 'key': key, 'title': f.get('title') or key,
            'confidence': f.get('confidence') or '未标注', 'phase': f.get('phase'), 'time': f.get('time'),
            'source': f.get('source'), 'scope': f.get('scope'), 'published': clean(f.get('published')),
            'advice': f.get('advice'), 'measures': measures(key, ev, thresholds), 'flags': flags(key, ev),
            'notes': notes(key, ev, thresholds, f), 'evidence': clean(ev),
            'matching_rejection': rejection_view(mr) if isinstance(mr, dict) else None}


def timing_view(events, ref):
    '''只用该段已关联的输入时序证据；时间邻近只表示相关。'''
    tagged = [('%s.t%d' % (ref, i + 1), i, d(e)) for i, e in enumerate(events)]
    counts = {}
    for _, _, e in tagged:
        cat = txt(e.get('category'))
        counts[cat] = counts.get(cat, 0) + 1
    lines = []
    if tagged:
        lines.append('关联事件分类：' + '，'.join('%s %d' % (k, counts[k]) for k in sorted(counts)))
    delays = [t for t in tagged if t[2].get('category') == 'delivery_delay' and t[2].get('stream') == 'lidar']
    if delays:
        r, _, top = max(delays, key=lambda t: (_val(t[2].get('arrival_gap_sec')), -t[1]))
        lines.append('雷达延迟交付 %d 条（时间邻近该段 ±1 秒，仅相关）' % len(delays))
        lines.append('最突出：同一对相邻帧 source_gap=%s / arrival_gap=%s / 取证阈值=%s，时间 %s，话题 %s，来源 %s [%s]' % (
            fmt(top.get('source_gap_sec'), 's'), fmt(top.get('arrival_gap_sec'), 's'), fmt(top.get('threshold_sec'), 's'),
            txt(top.get('time')), txt(top.get('topic')), txt(top.get('source')), r))
        lines.append('结论：已观察到延迟交付；环节未定。候选：' + '、'.join(DELAY_CANDIDATES))
        lines.append('缺少中间时间戳：' + '、'.join(MISSING_STAMPS) + '，不能确认 DDS 或其他单一环节')
        lines.append('注：比较的是同一对相邻帧的源间隔与到达间隔，不是绝对网络单向延迟；取证阈值不是 Safety 阈值')
    gaps = [t for t in tagged if t[2].get('category') == 'source_gap']
    if gaps:
        r, _, top = max(gaps, key=lambda t: (_val(t[2].get('source_gap_sec')), -t[1]))
        lines.append('源间隔超限 %d 条，最大 source_gap=%s（取证阈值 %s，流 %s）[%s]：无法区分停发/中途丢弃/跳时，不等于硬件停发' % (
            len(gaps), fmt(top.get('source_gap_sec'), 's'), fmt(top.get('threshold_sec'), 's'), txt(top.get('stream')), r))
    if 'callback_silence' in counts:
        lines.append('callback_silence 只表示回调暂未到达，不说明原因')
    return {'counts': counts, 'lines': lines, 'events': [dict(clean(e), ref=r) for r, _, e in tagged]}


def episode_view(ep, index, thresholds):
    ref = 'H%d' % (index + 1)
    view = finding_view(ep, thresholds, ref)
    onset = d(ep.get('onset_evidence'))
    state, cov = ep.get('state'), ep.get('capture_coverage')
    view.update({
        'index': index, 'device': ep.get('device') or d(ep.get('evidence')).get('device'),
        'state': state, 'state_text': STATE_TEXT.get(state, '状态未知(' + txt(state) + ')'),
        'coverage': COVERAGE_TEXT.get(cov, '采集覆盖未知(' + txt(cov) + ')'),
        'first_time': ep.get('first_time'), 'last_time': ep.get('last_time'), 'end_time': ep.get('end_time'),
        'start_ns': clean(ep.get('start_ns')), 'last_ns': clean(ep.get('last_ns')), 'end_ns': clean(ep.get('end_ns')),
        'samples': clean(ep.get('samples')), 'next_raw': clean(ep.get('next_raw')),
        'onset_measures': measures(view['key'], onset, thresholds) if onset else [],
        'onset_evidence': clean(onset), 'onset_source': ep.get('onset_source'), 'onset_advice': ep.get('onset_advice'),
        'onset_title': ep.get('onset_title'),
        'event_paths': clean(lst(ep.get('event_paths'))), 'onset_context': clean(lst(ep.get('onset_context'))),
        'context': clean(lst(ep.get('context'))), 'timing': timing_view(lst(ep.get('input_timing_evidence')), ref)})
    return view


def history(report, thresholds):
    raws = [d(ep) for ep in lst(report.get('findings'))]
    views = [episode_view(ep, i, thresholds) for i, ep in enumerate(raws)]

    def when(i, keys):
        for k in keys:
            v = finite(raws[i].get(k))
            if v is not None:
                return (1, v, i)
        return (0, 0.0, i)

    groups = {}
    for v in views:
        gk = (v['code'], v['key'], txt(v['device']), txt(v['phase']))
        g = groups.setdefault(gk, {'code': v['code'], 'key': v['key'], 'device': v['device'], 'phase': v['phase'],
                                   'refs': [], 'latest': None, 'rank': None})
        g['refs'].append(v['ref'])
        r = when(v['index'], ('last_ns', 'end_ns', 'start_ns'))
        if g['rank'] is None or r > g['rank']:
            g['latest'], g['rank'] = v['ref'], r
    ordered = sorted(groups.values(), key=lambda g: g['rank'], reverse=True)
    first = {}
    for v in views:
        r = when(v['index'], ('start_ns', 'last_ns', 'end_ns'))
        rank = (0 if r[0] else 1, r[1], r[2])
        if v['code'] not in first or rank < first[v['code']][0]:
            first[v['code']] = (rank, v)
    onsets = []
    for _, v in sorted(first.values(), key=lambda item: item[0]):
        onsets.append({'code': v['code'], 'ref': v['ref'], 'key': v['key'],
                       'title': v['onset_title'] or '判据 ' + v['key'] + '（首次标题未记录）',
                       'time': v['first_time'], 'source': v['onset_source'],
                       'measures': v['onset_measures'], 'from_onset': bool(v['onset_measures'])})
    out_groups = [{k: g[k] for k in ('code', 'key', 'device', 'phase', 'refs', 'latest')} for g in ordered]
    return {'episodes': views, 'groups': out_groups, 'shown': [g['latest'] for g in ordered[:HISTORY_LIMIT]],
            'omitted_groups': max(0, len(ordered) - HISTORY_LIMIT), 'first_onsets': onsets, 'record_count': len(views)}


def global_timing(report):
    it = d(report.get('input_timing'))
    lines = []
    if it.get('status') == 'available':
        counts = d(it.get('counts'))
        text = '，'.join('%s %s' % (k, fmt(counts[k])) for k in sorted(counts, key=str)) or '无事件'
        lines.append('有界日志全局输入时序计数：%s（未关联到所示段的历史线索不参与归因）' % text)
        lines.append('观察记录 %s 条，全局事件 %d 条（详见 --all/--json 的 input_timing）' % (
            fmt(it.get('observed_records')), len(lst(it.get('events')))))
    else:
        lines.append('输入时序不可用（状态 %s）' % txt(it.get('status')))
    if it.get('limitation'):
        lines.append('局限：' + txt(it.get('limitation')))
    gaps = []
    if present(it.get('truncated_events')):
        gaps.append('时序事件截断 ' + short(it.get('truncated_events')))
    lost = d(it.get('lost_by_topic'))
    if any(present(v) for v in lost.values()):
        gaps.append('按话题丢失 ' + brief(lost, 8))
    if present(it.get('event_sequence_gaps')):
        gaps.append('事件序号缺口 ' + short(it.get('event_sequence_gaps')))
    return {'status': it.get('status'), 'lines': lines, 'gaps': gaps}


def match_failure(report):
    m = d(report.get('match_failure_image'))
    col = d(m.get('collector'))
    lines = []
    if m:
        lines.append('匹配失败图像：状态 %s；%s%s' % (txt(m.get('status')), txt(m.get('message')),
                                              ('；原因 ' + txt(m.get('reason'))) if m.get('reason') else ''))
        if m.get('image') or m.get('html'):
            lines.append('图像 %s；HTML %s' % (txt(m.get('image')), txt(m.get('html'))))
        if col:
            lines.append('图像采集器：状态 %s，运行 %s，新鲜 %s，最近错误 %s；计数 %s' % (
                txt(col.get('status')), fmt(col.get('running')), fmt(col.get('fresh')), txt(col.get('last_error')),
                brief(d(col.get('counters')), 8)))
        lines.append('分数与记录一致 %s；原始样本可重放 %s；可视化模式 %s' % (
            fmt(m.get('check_score_matches_recorded')), fmt(m.get('original_sample_replay_available')),
            txt(m.get('visualization_mode'))))
    else:
        lines.append('匹配失败图像信息缺失')
    lines.append('注：没有图像不代表没有匹配失败；--no-match-image 时不生成/读取匹配图像缓存')
    snaps = d(report.get('match_snapshots'))
    records = [dict(clean(d(r)), ref='M%d' % (i + 1)) for i, r in enumerate(lst(snaps.get('records')))]
    lines.append('匹配快照 %d 条；保留的已存日志 %s%s' % (len(records), short(snaps.get('retained_saved_logs')),
                                                ('；局限：' + txt(snaps.get('limitation'))) if snaps.get('limitation') else ''))
    for r in records[:3]:
        lines.append('[%s] %s，状态 %s' % (r['ref'], txt(r.get('saved_path')), txt(r.get('current_status'))))
    if len(records) > 3:
        lines.append('另 %d 条快照见索引 M4~M%d' % (len(records) - 3, len(records)))
    compact = []
    if m:
        has_file = m.get('image') or m.get('html')
        if has_file:
            compact.append('图像 %s（HTML %s，状态 %s）' % (txt(m.get('image')), txt(m.get('html')), txt(m.get('status'))))
        elif m.get('reason') or m.get('message'):
            compact.append('未提供图像：状态 %s，%s；无图不代表无匹配失败' % (
                txt(m.get('status')), txt(m.get('reason') or m.get('message'))))
        losses = any(present(v) for k, v in d(col.get('counters')).items()
                     if str(k).startswith(('rejected', 'dropped')) or k == 'writer_errors')
        if col and (col.get('last_error') or col.get('fresh') is False or losses
                    or col.get('status') not in (None, 'ok', 'running', 'idle', 'observed', 'not_requested')):
            compact.append('图像采集器：状态 %s，新鲜 %s，最近错误 %s；计数 %s' % (
                txt(col.get('status')), fmt(col.get('fresh')), txt(col.get('last_error')), brief(d(col.get('counters')), 8)))
        if truthy(m.get('synthetic_fixture')):
            compact.append('synthetic_fixture=true：合成测试图，不是事故现场图')
        if has_file and not (m.get('check_score_matches_recorded') is True and m.get('original_sample_replay_available') is True):
            compact.append('图像分数/阈值与记录一致性=%s，原始匹配样本可重放=%s：未验证前仅作可视化参考' % (
                fmt(m.get('check_score_matches_recorded')), fmt(m.get('original_sample_replay_available'))))
    if records:
        compact.append('匹配快照 %d 条（首个 %s）[%s]' % (
            len(records), txt(records[0].get('saved_path')), _refs_short([r['ref'] for r in records])))
    return {'lines': lines, 'snapshots': records, 'compact': compact}


def upstream(report):
    u = d(report.get('upstream_lidar'))
    rep = d(u.get('report'))
    findings = [dict(clean(d(f)), ref='U%d' % (i + 1)) for i, f in enumerate(lst(rep.get('findings')))]
    warns = [{'ref': 'UW%d' % (i + 1), 'text': txt(w)} for i, w in enumerate(lst(rep.get('warnings')))]
    active = rep.get('active')
    if active is True:
        kind = '当前驱动运行实例'
    elif active is False:
        kind = '历史驱动包（当前未发现雷达 worker）'
    else:
        kind = '驱动运行状态未记录'
    lines = [kind + '；状态 %s；%s；active=%s，有证据=%s' % (txt(u.get('status')), txt(u.get('message')),
                                              fmt(rep.get('active')), fmt(rep.get('has_evidence')))]
    for f in findings:
        lines.append('[%s] %s %s/%s：%s' % (f['ref'], txt(f.get('severity')), txt(f.get('role')),
                                         txt(f.get('stage')), txt(f.get('message'))))
    lines.extend('[%s] %s' % (w['ref'], w['text']) for w in warns)
    for role in lst(rep.get('roles')):
        role = d(role)
        dev = d(role.get('device'))
        lines.append('角色 %s：证据范围 %s，采集于 %s，距今 %s' % (txt(role.get('role')), txt(dev.get('evidence_scope')),
                                                    txt(dev.get('collected_at')), fmt(dev.get('age_sec'), 's')))
    lines.append('注：独立的驱动侧证据，与定位包无同次关联，不拼接为因果；active 或退出码 0 不等于健康')
    if not findings:
        lines.append('无驱动异常条目：不能据此确认当前设备或数据流正常')
    return {'lines': lines, 'findings': findings, 'warnings': warns, 'roles': clean(lst(rep.get('roles'))),
            'kind': kind, 'status': u.get('status')}


def _launch(it):
    return '阶段 %s：%s；证据 %s；来源 %s；建议 %s' % (txt(it.get('stage')), txt(it.get('title')), brief(it.get('evidence')),
                                            txt(it.get('source')), txt(it.get('advice')))


def _suppression(it):
    return '%s，范围 %s：%s（%s）' % (txt(it.get('time')), txt(it.get('scope')), brief(it.get('evidence')), txt(it.get('source')))


def _coexisting(it):
    return '%s，设备 %s，末次 %s，范围 %s：%s（%s）' % (txt(it.get('key')), txt(it.get('device')), txt(it.get('last_time')),
                                              txt(it.get('scope')), brief(it.get('evidence')), txt(it.get('source')))


def _protection(it):
    detail = txt(it.get('text')) if it.get('text') else brief(it.get('evidence'))
    return '%s %s：%s 次，%s ~ %s；%s（%s）' % (txt(it.get('key')), txt(it.get('title')), fmt(it.get('count')),
                                          txt(it.get('first_time')), txt(it.get('last_time')), detail,
                                          txt(it.get('first_source') or it.get('source')))


CATEGORY_SPECS = (
    ('warnings', '重要提示', 'W', lambda it: txt(it.get('text'))),
    ('launch_findings', '启动问题', 'L', _launch),
    ('partial_radar_suppressions', '部分雷达抑制', 'S', _suppression),
    ('coexisting_radar_conditions', '并存雷达状况', 'X', _coexisting),
    ('runtime_protections', '运行时保护', 'P', _protection),
)


def categories(report):
    out = []
    for field, title, prefix, describe in CATEGORY_SPECS:
        items = []
        for i, raw in enumerate(lst(report.get(field))):
            it = raw if isinstance(raw, dict) else {'text': raw}
            items.append({'ref': '%s%d' % (prefix, i + 1), 'text': describe(it), 'detail': clean(it)})
        out.append({'field': field, 'title': title, 'count': len(items), 'items': items})
    return out


def status_lines(report):
    rb = d(report.get('runtime_build'))
    comp = d(rb.get('compiled'))
    st = rb.get('status')
    label = {'matched': '运行构建与合同匹配', 'mismatch': '运行构建与合同不匹配：结论可能不适用于该构建',
             'unverified': '运行构建未验证：不能确认合同一致'}.get(st, '运行构建信息缺失或未知状态(' + txt(st) + ')')
    lines = ['构建：%s；源码版本 %s，dirty %s，构建于 %s；期望合同 %s' % (
        label, txt(comp.get('source_revision')), fmt(comp.get('source_dirty')), txt(comp.get('built_at')),
        short(rb.get('expected_contract')))]
    me = d(report.get('motion_evidence'))
    line = '运动证据：状态 ' + txt(me.get('status'))
    if str(me.get('status')).upper() == 'UNKNOWN' or truthy(me.get('unknown_advisory')):
        line += '（未知不等于通过）'
    if me.get('explanation'):
        line += '；' + txt(me.get('explanation'))
    lines.append(line)
    cp = d(report.get('compute_profile'))
    lines.append('计算剖析：状态 %s，启用 %s，路径 %s' % (txt(cp.get('status')), fmt(cp.get('enabled')), txt(cp.get('path'))))
    lines.append('采集器累计：%s；窗口内采集丢弃增量 %s' % (brief(d(report.get('collector_totals')), 8),
                                                short(report.get('window_collector_drop_delta'))))
    lines.append('启动证据 %s；匹配证据 %s' % (fmt(report.get('has_launch_evidence')), fmt(report.get('has_match_evidence'))))
    return lines


def build_summary(report):
    '''返回确定、可 JSON 序列化的摘要字典；不修改 report。'''
    report = d(report)
    thresholds = d(report.get('thresholds'))
    latest = d(report.get('latest'))
    current = [finding_view(f, thresholds, 'C%d' % (i + 1)) for i, f in enumerate(lst(report.get('current_findings')))]
    timing = global_timing(report)
    summary = {
        'version': SUMMARY_VERSION,
        'selection': {'selection': report.get('selection'), 'selected': report.get('selected'),
                      'resolved_run': report.get('resolved_run'), 'evidence_scope': report.get('evidence_scope')},
        'freshness': freshness(report),
        'latest_record': dict(runstate(latest, len(current)), time=latest.get('_time'), source=latest.get('_source'),
                              scope=latest.get('_scope'), policy_version=latest.get('policy_version')),
        'current_findings': current,
        'event': {'status': report.get('event_status'), 'trigger_time': report.get('trigger_time'),
                  'trigger_reasons': lst(report.get('trigger_reasons'))},
        'history': history(report, thresholds),
        'input_timing': timing,
        'collection_gaps': timing['gaps'],
        'match_failure': match_failure(report),
        'upstream_lidar': upstream(report),
        'categories': categories(report),
        'runtime_protection_scope': report.get('runtime_protection_scope'),
        'status_lines': status_lines(report),
        'map_overlap_policy': clean(d(d(report.get('overlap_observations')).get('policy'))),
        'overlap_observations': clean(d(report.get('overlap_observations'))),
        'details': {'thresholds': thresholds, 'collector_totals': report.get('collector_totals'),
                    'window_collector_drop_delta': report.get('window_collector_drop_delta'),
                    'image_collector_counters': d(d(d(report.get('match_failure_image')).get('collector')).get('counters'))},
    }
    return clean(summary)


def _event_key(e):
    when = e.get('observed_wall_ns')
    return (txt(e.get('category')), txt(e.get('stream')), txt(e.get('source')), txt(e.get('topic')),
            txt(e.get('time') if when is None else when))


def _refs_short(refs):
    refs = list(refs)
    if len(refs) <= 3:
        return '、'.join(refs)
    return '%s、%s…%s（共 %d）' % (refs[0], refs[1], refs[-1], len(refs))


def _aggregate(items):
    totals, order = {}, []
    for it in items:
        det = d(it.get('detail'))
        name = txt(det.get('key') or det.get('title') or det.get('scope') or det.get('text'))
        count = finite(det.get('count'))
        weight = int(count) if count is not None and count == int(count) and count > 0 else 1
        if name not in totals:
            totals[name] = [0, []]
            order.append(name)
        totals[name][0] += weight
        totals[name][1].append(it['ref'])
    return '；'.join('%s×%d [%s]' % (n, totals[n][0], _refs_short(totals[n][1])) for n in order)


def _finding_compact(v, out):
    tags = '，'.join(t for t in (v['key'], v['confidence'], PHASE_TEXT.get(v.get('phase'), '')) if t)
    measurements = v.get('measures') or []
    split_translation = v['key'] == 'translation_condition' and len(measurements) > 1
    body = measurements[0] if split_translation else '；'.join(measurements)
    if not body:
        body = '记录字段 ' + '，'.join(v['flags']) if v.get('flags') else '实测值未记录'
    out.append('  [%s] 码 %s %s（%s）时间 %s：%s' % (v['ref'], v['code'], txt(v.get('title')), tags,
                                             txt(v.get('time')), body))
    if split_translation:
        out.extend('      ' + line for line in measurements[1:])
    extra = []
    mr = v.get('matching_rejection')
    if mr:
        extra.append('匹配拒绝 ' + mr['text'] + '；' + mr['note'])
    extra.extend(v.get('notes') or [])
    if extra:
        out.append('      ' + '；'.join(extra))


def _episode_compact(v, refs, out):
    tags = [v['key']]
    for tag in (PHASE_TEXT.get(v.get('phase'), ''), ('设备 ' + txt(v.get('device'))) if v.get('device') else ''):
        if tag:
            tags.append(tag)
    measurements = v.get('onset_measures') or v.get('measures') or []
    split_translation = v['key'] == 'translation_condition' and len(measurements) > 1
    if measurements:
        values = ('触发值 ' if v.get('onset_measures') else '末样本 ') + (
            measurements[0] if split_translation else '；'.join(measurements))
    else:
        values = '数值未记录'
    line = '  [%s] 码 %s（%s）%s；开始 %s，末次 %s；%s；同类 %d 段' % (
        v['ref'], v['code'], '，'.join(tags), v['state_text'], txt(v.get('first_time')), txt(v.get('last_time')),
        values, len(refs))
    if v.get('coverage') != COVERAGE_TEXT['range_overlap']:
        line += '；' + txt(v.get('coverage'))
    out.append(line)
    if split_translation:
        out.extend('      ' + item for item in measurements[1:])
    if v.get('matching_rejection'):
        out.append('      历史匹配拒绝：' + v['matching_rejection']['text'])


def _timing_compact(eps, shown, timing, out):
    seen, events = set(), []
    for ref in shown:
        for e in d(eps[ref].get('timing')).get('events') or []:
            key = _event_key(e)
            if key not in seen:
                seen.add(key)
                events.append(e)
    if events:
        counts = {}
        for e in events:
            cat = txt(e.get('category'))
            counts[cat] = counts.get(cat, 0) + 1
        out.append('时序（所示段关联，仅时间相关，去重 %d 条）：%s' % (
            len(events), '，'.join('%s %d' % (k, counts[k]) for k in sorted(counts))))
        delays = [e for e in events if e.get('category') == 'delivery_delay' and e.get('stream') == 'lidar']
        if delays:
            top = max(delays, key=lambda e: _val(e.get('arrival_gap_sec')))
            out.append('  雷达延迟交付 %d 条；最突出同一对相邻帧 source_gap=%s / arrival_gap(回调到达间隔)=%s / 取证阈值=%s，时间 %s，话题 %s [%s]' % (
                len(delays), fmt(top.get('source_gap_sec'), 's'), fmt(top.get('arrival_gap_sec'), 's'),
                fmt(top.get('threshold_sec'), 's'), txt(top.get('time')), txt(top.get('topic')), txt(top.get('ref'))))
            out.append('  结论：已观察到延迟交付，环节未定；候选 %s；缺 %s，不能确认 DDS；相邻间隔比较不是绝对网络单向延迟，取证阈值不是 Safety 阈值' % (
                '、'.join(DELAY_CANDIDATES), '、'.join(MISSING_STAMPS)))
        gaps = [e for e in events if e.get('category') == 'source_gap']
        if gaps:
            top = max(gaps, key=lambda e: _val(e.get('source_gap_sec')))
            out.append('  源间隔超限 %d 条，最大 source_gap=%s（取证阈值 %s）[%s]：不能区分停发/中途丢弃/跳时，不等于硬件停发' % (
                len(gaps), fmt(top.get('source_gap_sec'), 's'), fmt(top.get('threshold_sec'), 's'), txt(top.get('ref'))))
        if 'callback_silence' in counts:
            out.append('  callback_silence 只表示回调暂未到达，不说明原因')
    elif shown:
        out.append('时序：所示段无关联输入时序证据，不建立时序链')
    g = d(timing)
    if g.get('status') == 'available' and g.get('lines'):
        out.append('  ' + g['lines'][0])
    else:
        limits = [line for line in g.get('lines') or [] if line.startswith('局限：')]
        out.append('时序：输入时序不可用（状态 %s）%s' % (txt(g.get('status')), ('；' + limits[0]) if limits else ''))


def _color_summary(lines, summary):
    '''只给摘要加终端样式；历史和未知状态不用当前正常/故障的颜色。'''
    current = d(summary.get('freshness')).get('is_current') is True
    severity = d(summary.get('latest_record')).get('severity')
    state_style = {'normal': '1;32', 'fault': '1;31'}.get(severity, '1;33') if current else '1;33'
    colored = []
    for line in lines:
        history_style = history_line_style(line)
        style = history_style[2:-1] if history_style else None
        if style:
            pass
        elif line.startswith(('【', '历史判据', '时序', '匹配图像/快照', '上游驱动', '完整：', '运动证据：')):
            style = '1;36'
        elif line.startswith('范围/新鲜度：'):
            style = '36' if current else '1;33'
        elif line.startswith(('最新记录：', '当前记录的判据', '所读最新记录（非当前）的判据')):
            style = state_style
        elif line.startswith('  [C'):
            style = '31' if current and severity == 'fault' else '33'
        elif line.startswith(('  码 ', '重要提示', '  结论：', '  雷达延迟交付')):
            style = '1;33'
        elif line.startswith(('  [H', '  [W', '  [U', '所选事件：', '  源间隔超限')):
            style = '33'
        elif line.startswith('      '):
            style = '1;33'
        elif line.startswith('构建：'):
            style = '36' if line.startswith('构建：运行构建与合同匹配；') else '1;33'
        if style:
            # 保留原始数值与阈值文字，仅强调每组比较；复位后继续使用该行颜色。
            parts = []
            for part in line.split('；'):
                before, marker, values = part.partition('实测 ')
                if not marker:
                    before, marker, values = part.partition('实测误差 ')
                if marker and ' / 阈值 ' in values:
                    part = before + '\x1b[1m' + marker + values + '\x1b[22m\x1b[' + style + 'm'
                parts.append(part)
            line = '\x1b[' + style + 'm' + '；'.join(parts) + '\x1b[0m'
        colored.append(line)
    return colored


def render_text(s, color=False):
    '''默认紧凑渲染：每类证据 1–2 行；完整条目见 render_index 与 JSON。'''
    s = d(s)
    out = []
    title = '【定位诊断摘要】'
    out.append(title)
    sel, fr, lr = d(s.get('selection')), d(s.get('freshness')), d(s.get('latest_record'))
    out.append('范围/新鲜度：' + txt(fr.get('text')))
    out.append('最新记录：' + '，'.join(lr.get('parts') or []) + ' → ' + txt(lr.get('verdict')))
    out.append('记录来源：时间 %s，%s，策略版本 %s；选择 %s，运行 %s' % (
        txt(lr.get('time')), txt(lr.get('source')), fmt(lr.get('policy_version')),
        txt(sel.get('selection')), txt(sel.get('resolved_run'))))
    label = '当前记录' if fr.get('is_current') else '所读最新记录（非当前）'
    current = s.get('current_findings') or []
    out.append('%s的判据（%d 项）%s' % (label, len(current), '' if current else '：未检出已知故障判据'))
    for v in current:
        _finding_compact(v, out)
    ev = d(s.get('event'))
    if ev.get('status') or ev.get('trigger_time') or ev.get('trigger_reasons'):
        reasons = '；'.join(short(r) for r in ev.get('trigger_reasons') or []) or '未记录'
        out.append('所选事件：状态 %s，触发时间 %s；触发原因：%s' % (txt(ev.get('status')), txt(ev.get('trigger_time')), reasons))
    hist = d(s.get('history'))
    eps = {v['ref']: v for v in hist.get('episodes') or []}
    groups = hist.get('groups') or []
    by_latest = {g['latest']: g for g in groups}
    shown = hist.get('shown') or []
    def format_history_episode(episode, refs):
        lines = []
        _episode_compact(episode, refs, lines)
        return lines
    out.extend(render_history_sections(hist, format_history_episode, txt, HISTORY_LIMIT))
    timing_refs = list(dict.fromkeys(shown + [f['ref'] for f in hist.get('first_onsets') or []]))
    _timing_compact(eps, timing_refs, s.get('input_timing'), out)
    for i, line in enumerate(d(s.get('match_failure')).get('compact') or []):
        out.append(('匹配图像/快照：' if i == 0 else '  ') + line)
    up = d(s.get('upstream_lidar'))
    found = up.get('findings') or []
    out.append('上游驱动（独立证据，不拼接为因果；active 或退出码 0 不等于健康）：%s；状态 %s%s' % (
        txt(up.get('kind')), txt(up.get('status')), '' if found else '；无异常条目，不能据此确认当前设备或数据流正常'))
    for f in found[:3]:
        out.append('  [%s] %s %s/%s：%s' % (f['ref'], txt(f.get('severity')), txt(f.get('role')),
                                         txt(f.get('stage')), txt(f.get('message'))))
    if len(found) > 3:
        out.append('  另 %d 条驱动异常见索引（%s）' % (len(found) - 3, _refs_short([f['ref'] for f in found[3:]])))
    out.extend('  [%s] %s' % (w['ref'], w['text']) for w in up.get('warnings') or [])
    cats = {c.get('field'): c for c in s.get('categories') or []}
    warns = d(cats.get('warnings')).get('items') or []
    if warns:
        out.append('重要提示（%d）：' % len(warns))
        out.extend('  [%s] %s' % (it['ref'], it['text']) for it in warns)
    empty = []
    for field in ('launch_findings', 'partial_radar_suppressions', 'coexisting_radar_conditions', 'runtime_protections'):
        c = d(cats.get(field))
        if not c.get('items'):
            empty.append(txt(c.get('title')))
            continue
        line = '%s %d：%s' % (c['title'], c['count'], _aggregate(c['items']))
        if field == 'runtime_protections' and s.get('runtime_protection_scope'):
            line += '（范围 %s）' % short(s.get('runtime_protection_scope'))
        out.append(line)
    if empty:
        out.append('无记录：' + '、'.join(empty))
    status = s.get('status_lines') or []
    gaps = s.get('collection_gaps') or []
    out.extend(status[:2])
    overlap = d(s.get('map_overlap_policy'))
    if isinstance(overlap.get('map_overlap_check_enabled'), bool):
        enabled = overlap['map_overlap_check_enabled']
        out.append('跟踪阶段重合率检查：%s（启动日志确认，仅作用于跟踪）；%s；实际阈值 %s' % (
            '开启' if enabled else '关闭',
            '参与匹配拒绝' if enabled else '只记录重合率，不作为拒绝项',
            fmt(overlap.get('map_min_overlap_ratio'))))
    else:
        out.append('跟踪阶段重合率检查：未知（该运行日志未记录，不按新默认值推断）')
    gravity_fields = ('tracking_map_max_gravity_angle_deg', 'tracking_map_min_gravity_cos',
                      'initialization_map_min_gravity_cos', 'relocalization_map_min_gravity_cos')
    if any(overlap.get(key) is not None for key in gravity_fields):
        out.append('重力一致性配置（启动日志）：跟踪夹角上限 %s° / 余弦下限 %s；初始化余弦下限 %s；重定位余弦下限 %s' %
                   tuple(fmt(overlap.get(key)) for key in gravity_fields))
    if overlap.get('map_consistency_enabled') in (0, False, '0', 'false'):
        out.append('地图一致性总开关已关闭；此时也没有重合率测量，按原总开关语义解释')
    observation = d(s.get('overlap_observations'))
    latest_observation = d(observation.get('latest'))
    if latest_observation:
        info = d(latest_observation.get('evidence'))
        out.append('最近低重合率观察 %s：实测 %s / 阈值 %s；%s；连续 %s 帧 / %s 秒' % (
            txt(latest_observation.get('time')), fmt(info.get('ratio')), fmt(info.get('threshold')),
            '已接受，非故障' if info.get('is_observation') else '证据需核对：' + txt(info.get('result')),
            fmt(info.get('streak_frames')), fmt(info.get('streak_seconds'))))
        out.append('  观察累计低重合率 %s 帧，因限频未请求采集 %s 帧；不是故障次数' % (
            fmt(info.get('total_low_frames')), fmt(info.get('capture_suppressed'))))
    collector = d(observation.get('collector'))
    if collector.get('status') == 'available':
        counts = d(collector.get('counters'))
        dropped = sum(counts.get(k, 0) for k in ('dropped_busy','dropped_queue','dropped_oversize','rejected_disk','rejected_map','writer_errors'))
        out.append('低重合率观察取证：已保存 %s，限频跳过 %s，丢弃/写入异常计数 %s；独立配额 %s MiB' % (
            fmt(counts.get('saved')), fmt(counts.get('suppressed')), dropped, fmt(collector.get('disk_mb'))))
    out.append('采集：缺口 %s；%s' % ('；'.join(gaps) if gaps else '输入时序未记录截断/丢失/序号缺口',
                                  status[3] if len(status) > 3 else '采集器统计未记录'))
    out.append('完整：--all（摘要 + 旧版全部判据段 + 证据索引）；结构化：--json 的 diagnostic_summary 字段')
    return NL.join(_color_summary(out, s) if color else out)


def render_summary(report, color=False):
    return render_text(build_summary(report), color)


def _plain(value):
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return fmt(value)


def _dump(out, title, obj, indent='    '):
    obj = d(obj)
    if obj:
        out.append(indent[:-2] + title + '：')
        out.extend(indent + str(k) + ' = ' + _plain(obj[k]) for k in sorted(obj, key=str))


def _index_view(out, v, seen=None):
    head = '[%s] 码 %s %s（%s）%s；时间 %s；来源 %s；范围 %s；阶段 %s' % (
        v['ref'], v['code'], txt(v.get('title')), v['key'], v['confidence'], txt(v.get('time')),
        txt(v.get('source')), txt(v.get('scope')), txt(v.get('phase')))
    if v.get('state_text'):
        head += '；%s；%s；%s ~ %s，结束 %s；触发来源 %s' % (
            v['state_text'], v['coverage'], txt(v.get('first_time')), txt(v.get('last_time')),
            txt(v.get('end_time')), txt(v.get('onset_source')))
    out.append(head)
    out.extend('  触发时 ' + m for m in v.get('onset_measures') or [])
    out.extend('  ' + m for m in v.get('measures') or [])
    out.extend('  注 ' + n for n in v.get('notes') or [])
    for label in ('advice', 'onset_advice'):
        if v.get(label):
            out.append('  ' + label + ' ' + txt(v[label]))
    _dump(out, '证据', v.get('evidence'))
    _dump(out, '触发证据', v.get('onset_evidence'))
    mr = v.get('matching_rejection')
    if mr:
        out.append('  匹配拒绝 ' + mr['text'] + '；' + mr['note'])
        _dump(out, '匹配拒绝证据', mr.get('evidence'))
    for key in ('event_paths', 'onset_context', 'context'):
        out.extend('  ' + key + ' ' + _plain(item) for item in v.get(key) or [])
    timing = d(v.get('timing'))
    out.extend('  时序 ' + line for line in timing.get('lines') or [])
    seen = {} if seen is None else seen
    for e in timing.get('events') or []:
        key = _event_key(e)
        if key in seen:
            out.append('  [%s] 同 [%s]' % (e.get('ref'), seen[key]))
            continue
        seen[key] = e.get('ref')
        out.append('  [%s] %s' % (e.get('ref'), _plain({k: x for k, x in e.items() if k != 'ref'})))


def render_index(summary):
    '''全部编号的逐项索引；每个字段一行，不把整个 JSON 塞进一行。'''
    s = d(summary)
    out = ['==== 证据索引（全部编号） ====', '-- 所读最新记录判据 --']
    for v in s.get('current_findings') or []:
        _index_view(out, v)
    hist = d(s.get('history'))
    out.append('-- 历史判据类型（%d 种，证据段 %s 条） --' % (len(hist.get('groups') or []), fmt(hist.get('record_count'))))
    for g in hist.get('groups') or []:
        out.append('类型 码 %s %s 设备 %s 阶段 %s：%s（最近 %s）' % (g['code'], g['key'], txt(g.get('device')),
                                                         txt(g.get('phase')), '、'.join(g['refs']), g['latest']))
    for f in hist.get('first_onsets') or []:
        out.append('首个触发 码 %s [%s] %s：%s' % (f['code'], f['ref'], f['key'], '；'.join(f.get('measures') or []) or '触发值未记录'))
    seen = {}
    for v in hist.get('episodes') or []:
        _index_view(out, v, seen)
    for c in s.get('categories') or []:
        out.append('-- %s（%d） --' % (c['title'], c['count']))
        for it in c['items']:
            out.append('[%s] %s' % (it['ref'], it['text']))
            _dump(out, '详情', it.get('detail'))
    up = d(s.get('upstream_lidar'))
    out.append('-- 上游雷达驱动（独立证据） --')
    out.extend(up.get('lines') or [])
    for f in up.get('findings') or []:
        out.append('[%s] 证据 %s' % (f.get('ref'), _plain(f.get('evidence'))))
    mf = d(s.get('match_failure'))
    out.append('-- 匹配失败图像与快照 --')
    out.extend(mf.get('lines') or [])
    for r in mf.get('snapshots') or []:
        out.append('[%s] %s' % (r.get('ref'), _plain({k: x for k, x in r.items() if k != 'ref'})))
    out.append('-- 输入时序、采集、构建与运动 --')
    out.extend(d(s.get('input_timing')).get('lines') or [])
    out.extend(s.get('collection_gaps') or [])
    out.extend(s.get('status_lines') or [])
    det = d(s.get('details'))
    _dump(out, '采集器累计', det.get('collector_totals'), '  ')
    _dump(out, '图像采集器计数', det.get('image_collector_counters'), '  ')
    _dump(out, '有效阈值配置', det.get('thresholds'), '  ')
    return NL.join(out)
