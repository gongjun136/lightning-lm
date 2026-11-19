'''平移一致性诊断证据:解析扁平字段契约并生成简洁中文说明。

仅用于诊断展示。本模块:
  * 不根据基础配置或旧的 error_x/error_y 计算阈值;
  * 不重新判定安全结论,只在浮点表示容差内核对记录是否自洽;
  * 不读取文件、不运行子进程、不修改输入字典。
缺失或无效证据不会被当作正常,也不会被编造为确切触发项。
'''

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any, Callable

SCHEMA_VERSION = 1
# 现有分类器中表示“平移失败”的条件编码;其他编码在此不作解读。
TRANSLATION_CONDITION_FAIL = 2

VERSION_KEY = 'translation_evidence_version'
SEQUENCE_KEY = 'translation_evidence_sequence'
STAMP_KEY = 'translation_evaluation_stamp'
SELECTED_KEY = 'translation_selected_window'
ORIGIN_KEY = 'translation_evidence_origin'
AVAILABLE_KEY = 'translation_evidence_available'
LATCH_KEY = 'translation_failure_latched'
CONDITION_KEY = 'translation_condition'

META_KEYS = (VERSION_KEY, SEQUENCE_KEY, STAMP_KEY, SELECTED_KEY, ORIGIN_KEY, AVAILABLE_KEY)

SLOT_SUFFIXES = (
    'available', 'raw_distinct', 'from_stamp', 'to_stamp', 'dt_sec', 'distance_m',
    'predicted_x_m', 'predicted_y_m', 'measured_x_m', 'measured_y_m',
    'residual_x_m', 'residual_y_m', 'error_x_m', 'error_y_m',
    'limit_x_m', 'limit_y_m', 'fail_x', 'fail_y',
)
FLAG_SUFFIXES = ('available', 'raw_distinct', 'fail_x', 'fail_y')
NUMERIC_SUFFIXES = tuple(s for s in SLOT_SUFFIXES if s not in FLAG_SUFFIXES)

# (槽位名, 区间, 分支),顺序与 C++ 槽位 0..3 一致。
SLOTS = (
    ('adjacent_final', 'adjacent', 'final'),
    ('adjacent_raw', 'adjacent', 'raw'),
    ('window_final', 'window', 'final'),
    ('window_raw', 'window', 'raw'),
)
SELECTED_VALUES = ('adjacent', 'window', 'none')
ORIGIN_VALUES = ('current', 'latched_failure')


def slot_prefix(slot: str) -> str:
    '''返回槽位字段前缀,例如 translation_adjacent_final_。'''
    return 'translation_' + slot + '_'


SCHEMA_KEYS = META_KEYS + tuple(
    slot_prefix(slot) + suffix for slot, _interval, _branch in SLOTS for suffix in SLOT_SUFFIXES
)

_LEGACY_AGGREGATE_KEYS = ('error_x', 'error_y', 'translation_error_x', 'translation_error_y')
_INTERVAL_LABEL = {'adjacent': '相邻区间', 'window': '窗口区间'}
_BRANCH_LABEL = {'final': '发布分支', 'raw': '未平滑分支'}
_TRUE_TEXT = frozenset(('1', 'true', 'yes', 'on'))
_FALSE_TEXT = frozenset(('0', 'false', 'no', 'off'))


# ---------- 基础解析 ----------

def _as_text(value: Any) -> str | None:
    if isinstance(value, bytes):
        try:
            return value.decode('ascii')
        except UnicodeDecodeError:
            return None
    if isinstance(value, str):
        return value
    return None


def _parse_bool(value: Any) -> bool | None:
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return (value == 1) if value in (0, 1) else None
    if isinstance(value, float):
        if value == 0.0:
            return False
        if value == 1.0:
            return True
        return None
    text = _as_text(value)
    if text is None:
        return None
    text = text.strip().lower()
    if text in _TRUE_TEXT:
        return True
    if text in _FALSE_TEXT:
        return False
    return None


def _parse_float(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        try:
            return float(value)
        except OverflowError:
            return None
    text = _as_text(value)
    if text is None:
        return None
    try:
        return float(text.strip())
    except ValueError:
        return None


def _parse_int(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if math.isfinite(value) and value.is_integer():
            return int(value)
        return None
    text = _as_text(value)
    if text is None:
        return None
    text = text.strip()
    try:
        return int(text, 10)
    except ValueError:
        pass
    number = _parse_float(text)
    if number is not None and math.isfinite(number) and number.is_integer():
        return int(number)
    return None


def _parse_enum(value: Any, allowed: tuple[str, ...]) -> str | None:
    text = _as_text(value)
    if text is None:
        return None
    text = text.strip()
    return text if text in allowed else None


# ---------- 容差与一致性 ----------

def _close(a: float, b: float, *refs: float) -> bool:
    # 浮点表示容差:按参与量级取 4 个 ULP。仅用于核对记录,不用于重新判定。
    scale = max([1.0, abs(a), abs(b)] + [abs(r) for r in refs])
    return abs(a - b) <= 4.0 * math.ulp(scale)


def _flag_consistent(error: float, limit: float, flag: bool) -> bool:
    if error == limit:
        return not flag  # 相等按通过处理
    if _close(error, limit):
        return True  # 舍入范围内接受记录的判定,不重新分类
    return flag == (error > limit)


def _recorded_fail(comp: dict[str, Any]) -> bool:
    return comp['available'] is True and (comp['fail_x'] is True or comp['fail_y'] is True)


def _parse_slot(row: Mapping[str, Any], slot: str, interval: str, branch: str) -> dict[str, Any]:
    prefix = slot_prefix(slot)
    comp: dict[str, Any] = {
        'slot': slot, 'interval': interval, 'branch': branch, 'status': 'ok', 'role': None,
        'raw_fallback': None, 'issues': [], 'missing_fields': [],
    }
    for suffix in FLAG_SUFFIXES:
        key = prefix + suffix
        comp[suffix] = None
        if key not in row:
            comp['missing_fields'].append(key)
            continue
        comp[suffix] = _parse_bool(row[key])
        if comp[suffix] is None:
            comp['issues'].append(suffix + ' 不是 0/1 标志')
    for suffix in NUMERIC_SUFFIXES:
        key = prefix + suffix
        comp[suffix] = None
        if key not in row:
            comp['missing_fields'].append(key)
            continue
        comp[suffix] = _parse_float(row[key])
        if comp[suffix] is None:
            comp['issues'].append(suffix + ' 不是数值')

    if comp['available'] is None:
        comp['status'] = 'inconsistent' if comp['issues'] else 'missing'
        return comp
    if comp['available'] is False:
        # 不可用槽位只能保持默认值:标志为 0,数值为 NaN。
        if any(comp[s] is True for s in ('raw_distinct', 'fail_x', 'fail_y')):
            comp['issues'].append('不可用的比较带有已置位标志')
        if any(isinstance(comp[s], float) and math.isfinite(comp[s]) for s in NUMERIC_SUFFIXES):
            comp['issues'].append('不可用的比较含有有限数值')
        comp['status'] = 'inconsistent' if comp['issues'] else 'unavailable'
        return comp

    if branch == 'raw' and comp['raw_distinct'] is not None:
        comp['raw_fallback'] = not comp['raw_distinct']
    if not comp['issues'] and not comp['missing_fields']:
        _check_values(comp)
    if comp['issues']:
        comp['status'] = 'inconsistent'
    elif comp['missing_fields']:
        comp['status'] = 'missing'
    return comp


def _check_values(comp: dict[str, Any]) -> None:
    issues = comp['issues']
    for suffix in NUMERIC_SUFFIXES:
        if not math.isfinite(comp[suffix]):
            issues.append(suffix + ' 非有限值')
    if issues:
        return
    start, end, dt = comp['from_stamp'], comp['to_stamp'], comp['dt_sec']
    if dt <= 0.0:
        issues.append('dt_sec 非正')
    if end <= start:
        issues.append('源时间没有向前推进')
    if not _close(dt, end - start, start, end):
        issues.append('dt_sec 与起止时间之差不一致')
    if comp['distance_m'] < 0.0:
        issues.append('D 为负')
    for axis in ('x', 'y'):
        predicted = comp[f'predicted_{axis}_m']
        measured = comp[f'measured_{axis}_m']
        residual = comp[f'residual_{axis}_m']
        error = comp[f'error_{axis}_m']
        limit = comp[f'limit_{axis}_m']
        flag = comp[f'fail_{axis}']
        if limit <= 0.0:
            issues.append(f'limit_{axis}_m 非正')
        if error < 0.0:
            issues.append(f'error_{axis}_m 为负')
        if not _close(residual, measured - predicted, measured, predicted):
            issues.append(f'residual_{axis}_m 与 实测-预测 不一致')
        if not _close(error, abs(residual), residual):
            issues.append(f'error_{axis}_m 与残差绝对值不一致')
        if limit > 0.0 and error >= 0.0 and not _flag_consistent(error, limit, flag):
            issues.append(f'fail_{axis} 与记录的误差/阈值矛盾')


# ---------- 公共接口 ----------

def extract_translation_evidence(row: Mapping[str, Any]) -> dict[str, Any]:
    '''返回新字典,只含本契约的 translation_* 字段及存在时的当前 translation_condition。'''
    if not isinstance(row, Mapping):
        return {}
    keys = SCHEMA_KEYS + (LATCH_KEY, CONDITION_KEY)
    return {key: row[key] for key in keys if key in row}


def translation_latch_attribution(row: Mapping[str, Any]) -> bool | None:
    '''返回明确记录的平移失败锁存;旧记录或字段无效时返回 None(无法归因)。'''
    if not isinstance(row, Mapping) or LATCH_KEY not in row:
        return None
    return _parse_bool(row[LATCH_KEY])


def needs_translation_finding(row: Mapping[str, Any]) -> bool:
    '''分类器使用:当前平移判定为失败,或明确记录了平移失败锁存时返回 True。

    偏航或协方差锁存不会使其返回 True;缺少 translation_failure_latched 的旧记录
    仍由现有通用锁存发现处理。
    '''
    if not isinstance(row, Mapping):
        return False
    condition = _parse_int(row[CONDITION_KEY]) if CONDITION_KEY in row else None
    return condition == TRANSLATION_CONDITION_FAIL or translation_latch_attribution(row) is True


def describe_translation(row: Mapping[str, Any], include_passing: bool = False) -> dict[str, Any]:
    '''解析平移证据,返回结构化细节与中文说明。

    status:complete(字段齐全且自洽)、partial(字段缺失或所选区间无可用比较)、
    missing(早于本格式的记录)、inconsistent(格式错误或字段相互矛盾)。
    triggers 只包含所选区间中明确记录为超限且自洽的比较。
    include_passing=True 时额外列出通过的比较,并明确标注为非故障。
    '''
    result: dict[str, Any] = {
        'status': 'missing', 'origin': None, 'selected_window': None, 'sequence': None,
        'evaluation_stamp': None, 'comparisons': [], 'triggers': [], 'lines': [], 'warnings': [],
    }
    if not isinstance(row, Mapping):
        result['lines'].append('没有可解析的记录,平移证据缺失,不视为正常。')
        return result
    if VERSION_KEY not in row:
        return _describe_legacy(row, result)
    if _parse_int(row[VERSION_KEY]) != SCHEMA_VERSION:
        result['status'] = 'inconsistent'
        result['warnings'].append('无法识别的平移证据版本:' + repr(row[VERSION_KEY]))
        result['lines'].append('平移证据版本无法识别,未解析分支、区间或阈值,不指认触发项。')
        return result

    issues: list[str] = []
    gaps: list[str] = []
    sequence = _meta(row, SEQUENCE_KEY, _parse_int, issues, gaps)
    if sequence is not None and sequence < 0:
        issues.append(SEQUENCE_KEY + ' 为负')
        sequence = None
    stamp = _meta(row, STAMP_KEY, _parse_float, issues, gaps)
    if stamp is not None and not math.isfinite(stamp):
        issues.append(STAMP_KEY + ' 非有限值')
        stamp = None
    selected = _meta(row, SELECTED_KEY, lambda v: _parse_enum(v, SELECTED_VALUES), issues, gaps)
    origin = _meta(row, ORIGIN_KEY, lambda v: _parse_enum(v, ORIGIN_VALUES), issues, gaps)
    meta_available = _meta(row, AVAILABLE_KEY, _parse_bool, issues, gaps)
    latch = _optional(row, LATCH_KEY, _parse_bool, issues)
    condition = _optional(row, CONDITION_KEY, _parse_int, issues)
    if condition is not None and condition not in (0, 1, 2):
        issues.append('translation_condition 不是已知判定编码')
    result.update(origin=origin, selected_window=selected, sequence=sequence,
                  evaluation_stamp=stamp)

    comparisons = [_parse_slot(row, slot, interval, branch) for slot, interval, branch in SLOTS]
    result['comparisons'] = comparisons
    for comp in comparisons:
        if comp['status'] == 'inconsistent':
            issues.append(comp['slot'] + ':' + ';'.join(comp['issues']))
        gaps.extend(comp['missing_fields'])

    interval = selected if selected in _INTERVAL_LABEL else None
    in_selected = [c for c in comparisons if interval is not None and c['interval'] == interval]
    selected_available = any(c['available'] is True for c in in_selected)
    selected_failed = any(_recorded_fail(c) for c in in_selected)
    adjacent_failed = any(_recorded_fail(c) for c in comparisons if c['interval'] == 'adjacent')

    # 与已验证调用方行为的一致性核对(只核对,不改判)。
    if selected is not None:
        if adjacent_failed and selected != 'adjacent':
            issues.append('相邻区间有超限记录却未被选为判定区间,与“相邻失败优先”矛盾')
        if meta_available is not None and meta_available != selected_available:
            issues.append(AVAILABLE_KEY + ' 与所选区间的比较可用性不一致')
    if origin == 'latched_failure':
        if not selected_failed:
            issues.append('保留的锁存失败证据中,所选区间没有超限记录')
        if latch is False:
            issues.append('证据来源为锁存失败,但 translation_failure_latched=0')
        if condition == TRANSLATION_CONDITION_FAIL:
            issues.append('当前判定为失败时不应展示保留的历史证据')
        if condition == 1:
            issues.append('当前判定已通过,不应保留平移失败锁存')
    elif origin == 'current' and condition is not None:
        failing_now = condition == TRANSLATION_CONDITION_FAIL
        if failing_now and not selected_failed:
            issues.append('当前判定为失败,但所选区间没有超限记录')
        if not failing_now and selected_failed:
            issues.append('所选区间有超限记录,但当前判定不是失败')
        if not failing_now and latch is True:
            issues.append('平移失败锁存仍有效,但未展示保留的失败证据')

    if issues:
        status = 'inconsistent'
    elif gaps or not selected_available:
        status = 'partial'
    else:
        status = 'complete'
    result['status'] = status
    _assign_roles(comparisons, interval, status)
    result['triggers'] = [c for c in comparisons if c['role'] == 'trigger']
    result['warnings'].extend(issues)
    if gaps:
        result['warnings'].append('缺少字段:' + '、'.join(gaps[:6]) + ('等' if len(gaps) > 6 else ''))
    _render_lines(result, comparisons, interval, origin, selected_available, selected_failed,
                  include_passing, sequence, stamp, status, condition)
    return result


# ---------- 内部辅助 ----------

def _meta(row: Mapping[str, Any], key: str, parser: Callable[[Any], Any],
          issues: list[str], gaps: list[str]) -> Any:
    if key not in row:
        gaps.append(key)
        return None
    value = parser(row[key])
    if value is None:
        issues.append(key + ' 格式非法:' + repr(row[key]))
    return value


def _optional(row: Mapping[str, Any], key: str, parser: Callable[[Any], Any],
              issues: list[str]) -> Any:
    if key not in row:
        return None
    value = parser(row[key])
    if value is None:
        issues.append(key + ' 格式非法:' + repr(row[key]))
    return value


def _assign_roles(comparisons: list[dict[str, Any]], interval: str | None, status: str) -> None:
    for comp in comparisons:
        if comp['available'] is not True:
            comp['role'] = comp['status'] if comp['status'] in ('missing', 'inconsistent') else 'unavailable'
            continue
        in_selected = comp['interval'] == interval
        failed = _recorded_fail(comp)
        if comp['status'] != 'ok' or (in_selected and failed and status == 'inconsistent'):
            comp['role'] = 'unconfirmed'
        elif failed:
            comp['role'] = 'trigger' if in_selected else 'other_interval_failure'
        else:
            comp['role'] = 'selected_pass' if in_selected else 'other_interval_pass'


def _describe_legacy(row: Mapping[str, Any], result: dict[str, Any]) -> dict[str, Any]:
    result['status'] = 'missing'
    result['lines'].append(
        '该记录早于平移证据格式:触发的分支(发布/未平滑)、区间(相邻/窗口)与实际动态阈值均未记录,'
        '无法从逐轴最大残差或基础容差重建;不能据此判定为正常,也不能指认触发项。'
    )
    if any(key in row for key in _LEGACY_AGGREGATE_KEYS):
        result['warnings'].append(
            '旧字段 error_x/error_y 是跨分支、跨区间的逐轴最大值:X 与 Y 可能来自不同比较,'
            '最大残差也不一定是相对阈值最严重的超限;本模块不使用它们推断触发项。'
        )
    if any(key in row for key in SCHEMA_KEYS):
        result['warnings'].append('存在部分平移证据字段但缺少版本号,未使用这些字段。')
    if CONDITION_KEY in row and _parse_int(row[CONDITION_KEY]) == TRANSLATION_CONDITION_FAIL:
        result['lines'].append('记录的当前平移判定为失败,但缺少可归因的比较证据。')
    return result


def _render_lines(result: dict[str, Any], comparisons: list[dict[str, Any]], interval: str | None,
                  origin: str | None, selected_available: bool, selected_failed: bool,
                  include_passing: bool, sequence: int | None, stamp: float | None,
                  status: str, condition: int | None) -> None:
    lines = result['lines']
    warnings = result['warnings']
    seq_text = str(sequence) if sequence is not None else '未知'
    stamp_text = _stamp(stamp)
    if status == 'inconsistent':
        lines.append('平移证据存在矛盾或格式错误,不据此指认触发项;详见警告。')
    if origin == 'latched_failure':
        label = '[锁存的历史失败证据]'
        current = {0: '未知(UNKNOWN)', 1: '通过(PASS)', 2: '失败(FAIL)'}.get(condition, '未记录或无效')
        lines.append(
            '当前平移判定:' + current + ';以下为记录中保留的历史失败证据'
            '(序号 ' + seq_text + ',评估时刻 ' + stamp_text + '),不是新的比较。'
        )
    else:
        label = '[当前比较]'
        if interval is None:
            lines.append('本次评估没有判定区间(无长区间锚点或证据缺失),平移结果为 UNKNOWN,不视为通过。')
        elif not selected_available:
            lines.append('所选' + _INTERVAL_LABEL[interval] + '没有可用比较,证据不完整,不视为通过。')
        elif selected_failed and status != 'inconsistent':
            lines.append('当前平移判定:失败(序号 ' + seq_text + ',评估时刻 ' + stamp_text
                         + ',判定区间:' + _INTERVAL_LABEL[interval] + ')。')
        elif selected_failed:
            lines.append('所选区间含超限标志,但与其他记录矛盾,不能确认本次触发原因。')
        else:
            lines.append('本次评估所选' + _INTERVAL_LABEL[interval] + '没有超限记录(序号 ' + seq_text + ')。')

    def by_role(role: str) -> list[dict[str, Any]]:
        return [c for c in comparisons if c['role'] == role]

    triggers = by_role('trigger')
    for comp in triggers:
        lines.append(_render(comp, label))
    for comp in by_role('unconfirmed'):
        if _recorded_fail(comp):
            lines.append(_render(comp, label + '[证据矛盾或不完整,未确认]'))
    others = by_role('other_interval_failure')
    for comp in others:
        lines.append(_render(comp, label + '[非判定区间,未作为判定依据]'))
    if include_passing:
        for comp in by_role('selected_pass') + by_role('other_interval_pass'):
            lines.append(_render(comp, label + '[通过,非故障]'))
    if {c['branch'] for c in triggers} == {'final', 'raw'}:
        warnings.append('发布与未平滑分支同属一个估计器且区间重叠;两者同时超限是相关事件,不是独立佐证。')
    if any(c['raw_fallback'] for c in triggers + others):
        warnings.append('未平滑分支至少一个端点回退为发布位姿,该比较不是独立的未平滑证据。')


def _render(comp: dict[str, Any], label: str) -> str:
    branch = _BRANCH_LABEL[comp['branch']]
    if comp['raw_fallback']:
        branch += '(端点回退为发布位姿,非独立证据)'
    parts = [
        label + ' ' + _INTERVAL_LABEL[comp['interval']] + '/' + branch,
        '源时间 ' + _stamp(comp['from_stamp']) + ' → ' + _stamp(comp['to_stamp'])
        + ',实际时长 ' + _duration(comp['dt_sec']) + ' s',
        'D=' + _metres(comp['distance_m']) + ' m',
    ]
    for axis in ('x', 'y'):
        parts.append(_axis_text(comp, axis))
    return ';'.join(parts)


def _axis_text(comp: dict[str, Any], axis: str) -> str:
    error_text, limit_text = _pair(comp[f'error_{axis}_m'], comp[f'limit_{axis}_m'])
    flag = comp[f'fail_{axis}']
    if flag is True:
        verdict = '超限'
    elif flag is False:
        verdict = '未超限'
    else:
        verdict = '标志未记录'
    return axis.upper() + ' 实测误差 ' + error_text + ' m / 阈值 ' + limit_text + ' m:' + verdict


def _stamp(value: float | None) -> str:
    # 源时间戳按最短可往返表示输出:精确且紧凑。
    return '未记录' if value is None else repr(float(value))


def _duration(value: float | None) -> str:
    if value is None:
        return '未记录'
    return format(value, '.6g') if math.isfinite(value) else repr(value)


def _metres(value: float | None) -> str:
    if value is None:
        return '未记录'
    return format(value, '.4f') if math.isfinite(value) else repr(value)


def _pair(a: float | None, b: float | None) -> tuple[str, str]:
    # 若四位小数无法区分误差与阈值,则改用精确表示,避免“相等却超限”的误读。
    if a is None or b is None or not (math.isfinite(a) and math.isfinite(b)):
        return _metres(a), _metres(b)
    text_a, text_b = format(a, '.4f'), format(b, '.4f')
    if text_a == text_b and a != b:
        return repr(a), repr(b)
    return text_a, text_b
