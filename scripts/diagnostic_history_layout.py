#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""诊断历史摘要分层渲染（紧凑版）。

只负责展示层：分区、子标题、缩进与可选着色。
- 不改变输入排序（first_onsets 按首次时间升序、shown 按每类末次时间降序，均原样使用）；
- 不改变稳定证据索引（H1 等）、JSON schema、归因与故障逻辑；
- 不修改输入对象；
- H 是稳定证据索引，不是严重度，也不是排序号。

最近记录的拆行依赖原回调第一行的分隔约定：
  “[H] 码X（key/阶段）状态；开始 ...，末次 ...；测量1；测量2；同类 N 段”
只把括号外的“；”替换为换行，其余原文逐字保留；
找不到“；开始 ”时，整行原样作为条目首行。
"""

import re

__all__ = [
    'render_history_sections',
    'history_line_style',
    'apply_line_style',
]

TITLE_LINE = '历史判据（仅供回溯，不代表当前状态）'
EMPTY_LINE = '历史判据：无'
COUNT_NOTE_TEMPLATE = ('  共 {groups} 类 / {episodes} 段；条数不代表物理故障次数'
                       '（同一实际故障可能有日志与时间线多份记录）')
REF_NOTE_LINE = '  [Hn] 为本报告内证据索引，非严重度、非排序；同一 H 出现在两个分区不代表两次故障'
FIRST_HEADING = '  首次触发（每个故障码最早一次，时间从早到晚）'
RECENT_HEADING = '  最近记录（每类原因最近一段，末次时间从近到远）'
NONE_ITEM_LINE = '    （无）'
ENTRY_PREFIX = '    - '
DETAIL_INDENT = '        '
START_SEPARATOR = '；开始 '
START_LABEL = '开始 '
SEGMENT_SEPARATOR = '；'
FIRST_TIME_LABEL = '首次 '
ORIGINAL_CONT_INDENT = '      '
EXTRA_CONT_INDENT = '  '
MISSING_EPISODE_TEXT = '（该证据记录缺失，无法展开）'
EMPTY_EPISODE_TEXT = '（无详细内容）'
DEFAULT_LIMIT = 6

STYLE_HEADING = '\033[96m'   # 亮青：总标题与两个分区子标题
STYLE_EVIDENCE = '\033[93m'  # 亮黄：历史证据条目首行
STYLE_RESET = '\033[0m'
STYLE_SUBHEADING = STYLE_HEADING  # 兼容旧名称

_HEADING_LINES = frozenset([TITLE_LINE, FIRST_HEADING, RECENT_HEADING])
# 只匹配本模块生成的条目首行：
#   首次："    - 码 X [Hn] ..."    最近："    - [Hn] ..."
_ENTRY_RE = re.compile(r'^    - (?:码 \S+ )?\[[^\[\]\s]+\]')
_OPEN_BRACKETS = '（([【'
_CLOSE_BRACKETS = '）)]】'


def _as_list(value):
    """list/tuple 复制为新列表；其它类型视为空，保证不完整输入安全退化。"""
    if isinstance(value, (list, tuple)):
        return list(value)
    return []


def _field(item, key):
    if isinstance(item, dict):
        return item.get(key)
    return None


def _hashable(value):
    try:
        hash(value)
    except TypeError:
        return False
    return True


def _txt(text, value):
    result = text(value)
    return '' if result is None else str(result)


def _list_len_or_none(value):
    if isinstance(value, (list, tuple)):
        return len(value)
    return None


def _count_text(count, text):
    return str(count) if count is not None else _txt(text, None)


def _effective_limit(limit, total):
    """limit=None 表示全部展开（对应 --all）。"""
    if limit is None:
        return total
    try:
        n = int(limit)
    except (TypeError, ValueError):
        n = DEFAULT_LIMIT
    return max(0, min(n, total))


def _lines_of(raw):
    if raw is None:
        return []
    if isinstance(raw, str):
        return [raw]
    try:
        return [str(x) for x in raw]
    except TypeError:
        return []


def _index_episodes(episodes):
    index = {}
    for episode in episodes:
        ref = _field(episode, 'ref')
        if ref is not None and _hashable(ref) and ref not in index:
            index[ref] = episode
    return index


def _index_group_refs(groups):
    index = {}
    for group in groups:
        latest = _field(group, 'latest')
        if latest is None or not _hashable(latest) or latest in index:
            continue
        index[latest] = _as_list(_field(group, 'refs'))
    return index


def _split_top_level(value, separator=SEGMENT_SEPARATOR):
    """只在括号外的分隔符处拆分，避免把“（证据；…）”之类拆碎。"""
    parts = []
    buf = []
    depth = 0
    for ch in value:
        if ch in _OPEN_BRACKETS:
            depth += 1
        elif ch in _CLOSE_BRACKETS and depth > 0:
            depth -= 1
        if ch == separator and depth == 0:
            parts.append(''.join(buf))
            buf = []
        else:
            buf.append(ch)
    parts.append(''.join(buf))
    return parts


def _render_first_onset(item, text):
    """首行“码 X [Hn] 原因（key）”，次行首次时间，其后每项判据一行。"""
    title = _field(item, 'title')
    key = _field(item, 'key')
    reason = _txt(text, title)
    if key is not None and key != title:
        reason = '{}（{}）'.format(reason, _txt(text, key))
    head = '码 {} [{}] {}'.format(_txt(text, _field(item, 'code')),
                                 _txt(text, _field(item, 'ref')), reason)
    lines = [
        ENTRY_PREFIX + head,
        DETAIL_INDENT + FIRST_TIME_LABEL + _txt(text, _field(item, 'time')),
    ]
    measures = _field(item, 'measures')
    if isinstance(measures, str):
        measures = [measures]
    for measure in _as_list(measures):
        value = measure if isinstance(measure, str) else _txt(text, measure)
        if value.strip():
            lines.append(DETAIL_INDENT + value)
    return lines


def _reindent_continuation(line):
    if line.startswith(ORIGINAL_CONT_INDENT):
        return EXTRA_CONT_INDENT + line
    return DETAIL_INDENT + line.lstrip(' ')


def _render_recent(ref, episodes_by_ref, group_refs, format_episode, text):
    """首行直接使用原回调的“[H] 码X（…）状态”，次行开始/末次，其后测量、计数各一行。"""
    ref_tag = '[{}]'.format(_txt(text, ref))
    episode = episodes_by_ref.get(ref) if _hashable(ref) else None
    if episode is None:
        return [ENTRY_PREFIX + ref_tag + ' ' + MISSING_EPISODE_TEXT]
    refs = list(group_refs.get(ref, [ref]))  # 传副本，回调修改不影响输入
    raw = _lines_of(format_episode(episode, refs))
    if not raw:
        return [ENTRY_PREFIX + ref_tag + ' ' + EMPTY_EPISODE_TEXT]
    first = raw[0].strip(' ')
    head, sep, tail = first.partition(START_SEPARATOR)
    head = head.rstrip(' ')
    if not head.startswith(ref_tag):
        # 回调首行没有 [H] 时补上，保证证据索引不丢。
        head = (ref_tag + ' ' + head) if head else ref_tag
    lines = [ENTRY_PREFIX + head]
    if sep:
        parts = _split_top_level(tail)
        lines.append(DETAIL_INDENT + START_LABEL + parts[0].strip(' '))
        for part in parts[1:]:
            part = part.strip(' ')
            if part:
                lines.append(DETAIL_INDENT + part)
    for extra in raw[1:]:
        if extra.strip():
            lines.append(_reindent_continuation(extra))
    return lines


def render_history_sections(history, format_episode, text, limit=DEFAULT_LIMIT):
    """渲染历史判据分区，返回字符串行列表（不含颜色）。

    history: 已构建好的字典（episodes/groups/shown/first_onsets）。
    format_episode(episode, refs): 原有单条格式化，第一行两空格缩进，续行六空格。
    text(value): 原有缺失值中文转换器。
    limit: 每个分区最多展开条数；None 表示全部（--all）。
    """
    if not isinstance(history, dict):
        return [EMPTY_LINE]
    episodes = _as_list(history.get('episodes'))
    groups = _as_list(history.get('groups'))
    shown = _as_list(history.get('shown'))
    raw_first = _as_list(history.get('first_onsets'))
    first_onsets = [x for x in raw_first if isinstance(x, dict)]
    skipped_first = len(raw_first) - len(first_onsets)
    if not first_onsets and not shown and not groups:
        return [EMPTY_LINE]

    episodes_by_ref = _index_episodes(episodes)
    group_refs = _index_group_refs(groups)

    count_line = COUNT_NOTE_TEMPLATE.format(
        groups=_count_text(_list_len_or_none(history.get('groups')), text),
        episodes=_count_text(_list_len_or_none(history.get('episodes')), text))
    lines = [TITLE_LINE, count_line, REF_NOTE_LINE, '', FIRST_HEADING]

    first_count = _effective_limit(limit, len(first_onsets))
    if not first_onsets:
        lines.append(NONE_ITEM_LINE)
    for item in first_onsets[:first_count]:
        lines.extend(_render_first_onset(item, text))
    hidden_first = len(first_onsets) - first_count
    if hidden_first > 0:
        lines.append('    另有 {} 个故障码的首次记录未展开（使用 --all 查看全部）'.format(hidden_first))
    if skipped_first > 0:
        lines.append('    另有 {} 条首次记录格式不完整，已略过'.format(skipped_first))

    lines.append('')
    lines.append(RECENT_HEADING)
    recent_count = _effective_limit(limit, len(shown))
    displayed = shown[:recent_count]
    if not displayed:
        lines.append(NONE_ITEM_LINE)
    for ref in displayed:
        lines.extend(_render_recent(ref, episodes_by_ref, group_refs, format_episode, text))
    total = len(groups) if groups else len(shown)
    hidden_recent = max(0, total - len(displayed))
    if hidden_recent > 0:
        lines.append('    另有 {} 类原因未展开（使用 --all 查看全部）'.format(hidden_recent))
    return lines


def history_line_style(line):
    """返回历史分区专用的完整 ANSI 序列；非本模块特定标签的行返回 None，交给原有逻辑。

    只匹配精确标题行和条目首行（“    - 码 X [Hn]” / “    - [Hn]”），不按缩进笼统上色；
    历史行从不使用红色。
    """
    if not isinstance(line, str):
        return None
    if line in _HEADING_LINES:
        return STYLE_HEADING
    if _ENTRY_RE.match(line):
        return STYLE_EVIDENCE
    return None


def apply_line_style(line, style):
    """style 为空时原样返回，保证关闭颜色时文本完全一致。"""
    if not style:
        return line
    return style + line + STYLE_RESET
