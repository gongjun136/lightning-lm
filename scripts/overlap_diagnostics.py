"""跟踪阶段地图重叠诊断（Python 3.10，纯标准库）。

只解释单行记录中的实际字段，不推断默认值，不做故障分类，不修改输入。
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

LOW_OVERLAP_EVENT = 'LOC_TRACKING_LOW_OVERLAP'
MAP_OVERLAP_REASON = 'map_overlap'

HELP_SUMMARY = (
    'LOC_TRACKING_LOW_OVERLAP 是观测记录，不是故障：该帧已被跟踪接受（accepted=1），'
    '地图检查已执行（map_checked=1），且运行中程序的跟踪重叠门控关闭'
    '（map_overlap_check_enabled=0），仅说明重叠率低于阈值。'
    '旧日志缺少门控字段时显示为未知，不按默认值推断。'
)

_TRUE_STRINGS = frozenset({'true', 'True', 'TRUE', '1'})
_FALSE_STRINGS = frozenset({'false', 'False', 'FALSE', '0'})
_EXACT_FLOAT_INT_LIMIT = 2 ** 53
_MISSING = object()

_GATE_LINES = {
    'enabled': '地图重叠门控：启用（重叠率低于阈值将拒绝跟踪帧）',
    'disabled': '地图重叠门控：关闭（仅观测，低重叠不拒绝跟踪帧）',
    'unknown': '地图重叠门控：未知（记录中无有效门控字段，不推断默认值）',
}

_RESULT_LINES = {
    'pass': '重叠结果：通过（重叠率不低于阈值，相等即通过）',
    'fail': '重叠结果：未通过（重叠门控拒绝）',
    'bypassed_low': '重叠结果：低于阈值，但门控关闭，仅作观测',
    'invalid': '重叠结果：证据无效（缺失、非有限或越界），不视为正常',
    'unavailable': '重叠结果：本帧未执行地图检查，无可用重叠证据',
    'unknown': '重叠结果：未知（证据不足以判定）',
    'inconsistent': '重叠结果：记录自相矛盾，请核对日志来源与版本',
}

_OBSERVATION_COUNTERS = (
    ('observation_streak_frames', 'streak_frames'),
    ('observation_total_low_frames', 'total_low_frames'),
    ('observation_capture_suppressed', 'capture_suppressed'),
)


def help_summary() -> str:
    """返回低重叠观测事件的帮助说明。"""
    return HELP_SUMMARY


def _parse_bool(value: Any) -> tuple[str, bool | None]:
    """严格解析布尔值；返回 (状态, 值)，状态为 ok/missing/malformed。"""
    if value is _MISSING or value is None:
        return 'missing', None
    if isinstance(value, bool):
        return 'ok', value
    if isinstance(value, int):
        if value in (0, 1):
            return 'ok', value == 1
        return 'malformed', None
    if isinstance(value, str):
        if value in _TRUE_STRINGS:
            return 'ok', True
        if value in _FALSE_STRINGS:
            return 'ok', False
    return 'malformed', None


def _parse_unit(value: Any) -> tuple[str, float | None]:
    """解析 [0,1] 比例；返回 (状态, 数值)，状态为 ok/missing/invalid。"""
    if value is _MISSING or value is None:
        return 'missing', None
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        return 'invalid', None
    try:
        number = float(value)
    except (ValueError, OverflowError):
        return 'invalid', None
    if math.isfinite(number) and 0.0 <= number <= 1.0:
        return 'ok', number
    # 非有限值不保留数值（置为 None），保证结果可用 json.dumps(allow_nan=False) 序列化；分类仍为 invalid。
    return 'invalid', (number if math.isfinite(number) else None)


def _parse_counter(value: Any) -> tuple[str, int | None]:
    """解析精确的非负整数计数器；不对超出精确范围的浮点数取整。"""
    if value is _MISSING or value is None:
        return 'missing', None
    if isinstance(value, bool):
        return 'malformed', None
    if isinstance(value, int):
        return ('ok', value) if value >= 0 else ('malformed', None)
    if isinstance(value, str):
        if value.isascii() and value.isdigit():
            try:
                return 'ok', int(value)
            except ValueError:
                return 'malformed', None
        return 'malformed', None
    if isinstance(value, float):
        if math.isfinite(value) and value.is_integer() and 0.0 <= value < _EXACT_FLOAT_INT_LIMIT:
            return 'ok', int(value)
        return 'malformed', None
    return 'malformed', None


def _parse_seconds(value: Any) -> tuple[str, float | None]:
    """解析有限非负秒数。"""
    if value is _MISSING or value is None:
        return 'missing', None
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        return 'malformed', None
    try:
        number = float(value)
    except (ValueError, OverflowError):
        return 'malformed', None
    if math.isfinite(number) and number >= 0.0:
        return 'ok', number
    return 'malformed', None


def _text(value: Any) -> str | None:
    return value if isinstance(value, str) else None


def _fmt(value: float | None) -> str:
    return '缺失' if value is None else f'{value:.4f}'


def _display(value: Any) -> str:
    return '未记录' if value is None else str(value)


def describe_overlap(row: Mapping[str, Any]) -> dict[str, Any]:
    """解释一行记录的地图重叠门控状态与结果；返回新字典，不修改输入，无副作用。"""
    warnings: list[str] = []
    lines: list[str] = []
    if not isinstance(row, Mapping):
        warnings.append('输入不是映射类型，所有字段按缺失处理')
        row = {}

    gate_state, gate_value = _parse_bool(row.get('map_overlap_check_enabled', _MISSING))
    if gate_state == 'ok':
        gate = 'enabled' if gate_value else 'disabled'
    else:
        gate = 'unknown'
        if gate_state == 'malformed':
            warnings.append('map_overlap_check_enabled 格式无法识别，门控状态按未知处理')

    checked_state, checked = _parse_bool(row.get('map_checked', _MISSING))
    if checked_state == 'malformed':
        warnings.append('map_checked 格式无法识别，按未知处理')
    accepted_state, accepted = _parse_bool(row.get('accepted', _MISSING))
    if accepted_state == 'malformed':
        warnings.append('accepted 格式无法识别，接受状态未知')

    ratio_state, ratio = _parse_unit(row.get('map_overlap', _MISSING))
    threshold_state, threshold = _parse_unit(row.get('map_min_overlap', _MISSING))
    reason = _text(row.get('rejection_reason', _MISSING))
    event = _text(row.get('event', _MISSING))

    contradictions: list[str] = []
    if reason == MAP_OVERLAP_REASON:
        if gate == 'disabled':
            contradictions.append('门控关闭却记录了 map_overlap 拒绝原因')
        if accepted is True:
            contradictions.append('帧已接受却记录了 map_overlap 拒绝原因')
        if checked is False:
            contradictions.append('未执行地图检查却记录了 map_overlap 拒绝原因')
    if event == LOW_OVERLAP_EVENT:
        if accepted is False:
            contradictions.append('低重叠观测事件要求 accepted=1')
        if checked is False:
            contradictions.append('低重叠观测事件要求 map_checked=1')
        if gate == 'enabled':
            contradictions.append('低重叠观测事件要求 map_overlap_check_enabled=0')

    if checked is None:
        base = 'unknown'
        if ratio_state != 'missing' or threshold_state != 'missing':
            warnings.append('map_checked 缺失，重叠数值未经确认')
    elif checked is False:
        base = 'unavailable'
        if ratio_state != 'missing' or threshold_state != 'missing':
            warnings.append('map_checked=0：忽略残留的重叠率/阈值数值')
        ratio = None
        threshold = None
    elif ratio_state != 'ok' or threshold_state != 'ok':
        base = 'invalid'
        warnings.append('重叠率或阈值缺失、非有限或超出 [0,1]')
        if event == LOW_OVERLAP_EVENT:
            contradictions.append('低重叠观测事件的重叠证据无效')
    elif ratio >= threshold:
        base = 'pass'
        if reason == MAP_OVERLAP_REASON:
            contradictions.append('重叠率达到阈值却记录了 map_overlap 拒绝原因')
        if event == LOW_OVERLAP_EVENT:
            contradictions.append('低重叠观测事件中重叠率并未低于阈值')
    elif gate == 'enabled':
        base = 'fail'
        if accepted is True:
            contradictions.append('门控启用且重叠率低于阈值，帧却被标记为已接受')
    elif gate == 'disabled':
        base = 'bypassed_low'
    else:
        base = 'unknown'

    result = 'inconsistent' if contradictions else base
    warnings.extend('矛盾：' + item for item in contradictions)

    lines.append(_GATE_LINES[gate])
    if ratio is not None or threshold is not None:
        lines.append(f'重叠率 {_fmt(ratio)} / 阈值 {_fmt(threshold)}')
    lines.append(_RESULT_LINES[result])
    if accepted is True:
        lines.append('跟踪接受：本帧已被接受')
    elif accepted is False:
        suffix = f'（拒绝原因：{reason}）' if reason else ''
        lines.append('跟踪接受：本帧未被接受' + suffix)
    else:
        lines.append('跟踪接受：记录未提供有效 accepted 字段，接受状态未知')

    is_observation = result == 'bypassed_low' and accepted is True
    if is_observation:
        lines.append('说明：这是已接受帧的低重叠观测记录，不是故障')
    elif result == 'bypassed_low' and accepted is False:
        lines.append('说明：重叠门控未拒绝本帧，但本帧因其他原因未被接受，不代表接受成功')
    elif result == 'bypassed_low':
        lines.append('说明：低重叠仅为观测，但无法确认本帧是否被接受')

    counters: dict[str, int | None] = {}
    for key, name in _OBSERVATION_COUNTERS:
        state, value = _parse_counter(row.get(key, _MISSING))
        if state == 'malformed':
            warnings.append(f'{key} 不是精确的非负整数，已忽略（不填充为 0）')
        counters[name] = value
    seconds_state, streak_seconds = _parse_seconds(row.get('observation_streak_sec', _MISSING))
    if seconds_state == 'malformed':
        warnings.append('observation_streak_sec 不是有限非负数，已忽略（不填充为 0）')
    if any(value is not None for value in counters.values()) or streak_seconds is not None:
        seconds_text = '未记录' if streak_seconds is None else f'{streak_seconds:.3f}'
        lines.append('观测统计：连续 {} 帧 / {} 秒，累计低重叠 {} 帧，未采集 {} 帧'.format(
            _display(counters['streak_frames']),
            seconds_text,
            _display(counters['total_low_frames']),
            _display(counters['capture_suppressed']),
        ))

    return {
        'gate': gate,
        'result': result,
        'lines': lines,
        'warnings': warnings,
        'ratio': ratio,
        'threshold': threshold,
        'accepted': accepted,
        'map_checked': checked,
        'is_observation': is_observation,
        'streak_frames': counters['streak_frames'],
        'streak_seconds': streak_seconds,
        'total_low_frames': counters['total_low_frames'],
        'capture_suppressed': counters['capture_suppressed'],
    }
