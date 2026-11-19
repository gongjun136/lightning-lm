#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""diagnostic_history_layout 的合成场景单元测试（不使用真实日志、路径、日期或地址）。"""

import copy
import os
import re
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import diagnostic_history_layout as dhl  # noqa: E402

ANSI_RE = re.compile(r'\x1b\[[0-9;]*m')


def fake_text(value):
    """模拟现有缺失值中文转换器。"""
    if value is None or value == '':
        return '未知'
    return str(value)


def fake_format_episode(episode, refs):
    first = '  [{}] 码 {}（{}）{}；开始 {}，末次 {}；触发值 {}；同类 {} 段'.format(
        fake_text(episode.get('ref')), fake_text(episode.get('code')),
        fake_text(episode.get('key')), episode.get('state', '判据发生变化（未必恢复）'),
        fake_text(episode.get('start')), fake_text(episode.get('end')),
        episode.get('measure', '帧年龄 实测 0.31 s / 阈值 0.3 s（证据；示例）'),
        len(refs))
    return [first, '      附注：续行示例']


def make_history():
    episodes = [
        {'ref': 'H4', 'code': 4, 'key': 'overlap_drop', 'title': '重合率下降', 'start': 'T+700s', 'end': 'T+720s'},
        {'ref': 'H5', 'code': 2, 'key': 'bound_far', 'title': '越界（远端）', 'start': 'T+900s', 'end': 'T+905s',
         'state': '已观测到原始码回 0'},
        {'ref': 'H6', 'code': 4, 'key': 'overlap_low', 'title': '重合率偏低', 'start': 'T+800s', 'end': None},
        {'ref': 'H7', 'code': 2, 'key': 'bound_near', 'title': '越界（近端）', 'start': 'T+12s', 'end': 'T+60s'},
        {'ref': 'H10', 'code': 4, 'key': 'overlap_drop', 'title': '重合率下降', 'start': 'T+340s', 'end': 'T+350s'},
    ]
    groups = [
        {'latest': 'H5', 'refs': ['H5'], 'code': 2, 'key': 'bound_far'},
        {'latest': 'H6', 'refs': ['H6'], 'code': 4, 'key': 'overlap_low'},
        {'latest': 'H4', 'refs': ['H10', 'H4'], 'code': 4, 'key': 'overlap_drop'},
        {'latest': 'H7', 'refs': ['H7'], 'code': 2, 'key': 'bound_near'},
    ]
    first_onsets = [
        {'code': 2, 'ref': 'H7', 'title': '越界（近端）', 'key': 'bound_near',
         'time': 'T+12s', 'measures': ['阈值 0.95', '实测 0.90']},
        {'code': 4, 'ref': 'H10', 'title': '重合率下降', 'key': 'overlap_drop',
         'time': 'T+340s', 'measures': ['重合率 0.41 < 0.55']},
    ]
    return {
        'episodes': episodes,
        'groups': groups,
        'shown': ['H5', 'H6', 'H4', 'H7'],
        'first_onsets': first_onsets,
    }


def episode_by_ref(history, ref):
    for episode in history['episodes']:
        if episode.get('ref') == ref:
            return episode
    raise KeyError(ref)


def render(history, fmt=fake_format_episode, limit=dhl.DEFAULT_LIMIT):
    return dhl.render_history_sections(history, fmt, fake_text, limit=limit)


def section(lines, heading):
    start = lines.index(heading) + 1
    end = start
    while end < len(lines) and lines[end] != '':
        end += 1
    return lines[start:end]


def recent_refs(lines):
    heads = [l for l in section(lines, dhl.RECENT_HEADING) if l.startswith('    - ')]
    return [re.match(r'    - \[([^\]]+)\]', h).group(1) for h in heads]


class RenderHistorySectionsTest(unittest.TestCase):

    def test_sections_and_non_contiguous_order(self):
        lines = render(make_history())
        i_first = lines.index(dhl.FIRST_HEADING)
        i_recent = lines.index(dhl.RECENT_HEADING)
        self.assertLess(i_first, i_recent)
        self.assertEqual(lines[i_recent - 1], '')
        first_part = section(lines, dhl.FIRST_HEADING)
        self.assertEqual([l for l in first_part if l.startswith('    - ')],
                         ['    - 码 2 [H7] 越界（近端）（bound_near）',
                          '    - 码 4 [H10] 重合率下降（overlap_drop）'])
        self.assertEqual(recent_refs(lines), ['H5', 'H6', 'H4', 'H7'])
        self.assertFalse(any('证据：' in l for l in lines))  # 不再有重复的“证据 H”行

    def test_first_onset_compact(self):
        first_part = section(render(make_history()), dhl.FIRST_HEADING)
        self.assertEqual(first_part[:4], [
            '    - 码 2 [H7] 越界（近端）（bound_near）',
            '        首次 T+12s',
            '        阈值 0.95',
            '        实测 0.90',
        ])

    def test_hierarchy_code2_not_parent_of_later_code4(self):
        lines = render(make_history())
        i_recent = lines.index(dhl.RECENT_HEADING)
        self.assertTrue(dhl.RECENT_HEADING.startswith('  ') and not dhl.RECENT_HEADING.startswith('   '))
        self.assertEqual(lines[i_recent - 1], '')
        for part in (section(lines, dhl.FIRST_HEADING), section(lines, dhl.RECENT_HEADING)):
            for line in part:
                indent = len(line) - len(line.lstrip(' '))
                self.assertIn(indent, (4, 8), line)

    def test_overview_counts_and_ref_note(self):
        lines = render(make_history())
        self.assertEqual(lines[0], dhl.TITLE_LINE)
        self.assertIn('不代表当前', lines[0])
        self.assertEqual(lines[1], '  共 4 类 / 5 段；条数不代表物理故障次数（同一实际故障可能有日志与时间线多份记录）')
        self.assertEqual(lines[2], dhl.REF_NOTE_LINE)
        self.assertIn('同一 H 出现在两个分区不代表两次故障', dhl.REF_NOTE_LINE)
        self.assertIn('非严重度、非排序', dhl.REF_NOTE_LINE)

    def test_recent_split_preserves_text(self):
        history = make_history()
        lines = render(history)
        original = fake_format_episode(episode_by_ref(history, 'H6'), ['H6'])[0].strip()
        idx = next(i for i, l in enumerate(lines) if l.startswith('    - [H6]'))
        head = lines[idx][len(dhl.ENTRY_PREFIX):]
        time_line, measure_line, count_line = lines[idx + 1:idx + 4]
        self.assertEqual(time_line, '        开始 T+800s，末次 未知')
        self.assertTrue(measure_line.startswith('        触发值 '))
        self.assertIn('（证据；示例）', measure_line)  # 括号内分号不拆
        self.assertEqual(count_line, '        同类 1 段')
        self.assertEqual(lines[idx + 4], '        附注：续行示例')
        rebuilt = '；'.join([head, time_line.strip(), measure_line.strip(), count_line.strip()])
        self.assertEqual(rebuilt, original)
        self.assertIn('判据发生变化（未必恢复）', head)
        self.assertFalse(any('开始 ' in l and '触发值' in l for l in lines))
        self.assertIn('    - [H5] 码 2（bound_far）已观测到原始码回 0', lines)

    def test_refs_passed_to_format_episode(self):
        calls = []

        def recording(episode, refs):
            calls.append((episode['ref'], list(refs)))
            return fake_format_episode(episode, refs)

        render(make_history(), recording)
        self.assertEqual(calls, [('H5', ['H5']), ('H6', ['H6']),
                                 ('H4', ['H10', 'H4']), ('H7', ['H7'])])

    def test_unknown_time(self):
        history = make_history()
        history['first_onsets'][0]['time'] = None
        episode_by_ref(history, 'H5')['start'] = None
        lines = render(history)
        self.assertIn('        首次 未知', lines)
        self.assertIn('        开始 未知，末次 T+905s', lines)

    def test_empty_history(self):
        for history in (None, {}, [], 'bad', {'episodes': []},
                        {'episodes': [], 'groups': [], 'shown': [], 'first_onsets': []}):
            self.assertEqual(render(history), [dhl.EMPTY_LINE])

    def test_incomplete_input_degrades_safely(self):
        history = {
            'episodes': [{'ref': 'H5', 'code': 2}, 'junk', {'code': 9}],
            'shown': ['H5', 'H99', None],
            'first_onsets': [None, {'ref': 'H7', 'measures': None}, 'junk'],
        }
        lines = render(history)
        self.assertTrue(lines[1].startswith('  共 未知 类 / 3 段'))
        self.assertIn('    - 码 未知 [H7] 未知', lines)
        self.assertIn('        首次 未知', lines)
        self.assertIn('    - [H99] ' + dhl.MISSING_EPISODE_TEXT, lines)
        self.assertIn('    - [未知] ' + dhl.MISSING_EPISODE_TEXT, lines)
        self.assertIn('    另有 2 条首次记录格式不完整，已略过', lines)

    def test_empty_episode_format(self):
        lines = render(make_history(), lambda e, r: [])
        self.assertEqual(sum(l.endswith(dhl.EMPTY_EPISODE_TEXT) for l in lines), 4)
        self.assertIn('    - [H5] ' + dhl.EMPTY_EPISODE_TEXT, lines)

    def test_first_only_shows_none_in_recent(self):
        history = make_history()
        history['shown'] = []
        history['groups'] = []
        self.assertEqual(section(render(history), dhl.RECENT_HEADING), [dhl.NONE_ITEM_LINE])

    def test_input_not_mutated(self):
        history = make_history()
        snapshot = copy.deepcopy(history)

        def mutating(episode, refs):
            refs.append('HX')
            return fake_format_episode(episode, refs)

        render(history, mutating)
        render(history, limit=None)
        self.assertEqual(history, snapshot)

    def test_limit_and_all_hint(self):
        lines = render(make_history(), limit=2)
        self.assertEqual(recent_refs(lines), ['H5', 'H6'])
        self.assertIn('    另有 2 类原因未展开（使用 --all 查看全部）', lines)
        lines = render(make_history(), limit=1)
        self.assertIn('    另有 1 个故障码的首次记录未展开（使用 --all 查看全部）', lines)
        lines = render(make_history(), limit=None)
        self.assertFalse(any('--all' in l for l in lines))
        history = make_history()
        history['shown'] = ['H5', 'H6', 'H4']
        self.assertIn('    另有 1 类原因未展开（使用 --all 查看全部）', render(history, limit=None))


    def test_ref_tag_prepended_when_callback_lacks_it(self):
        lines = render(make_history(), lambda e, r: ['  码 9 无标签；开始 T+1s，末次 T+2s'])
        self.assertIn('    - [H5] 码 9 无标签', lines)
        self.assertIn('        开始 T+1s，末次 T+2s', lines)

    def test_no_start_separator_keeps_whole_line(self):
        lines = render(make_history(), lambda e, r: ['  [{}] 码 9 无时间；字段'.format(e['ref'])])
        self.assertIn('    - [H5] 码 9 无时间；字段', lines)
        self.assertFalse(any(l.startswith('        开始') for l in lines))

    def test_blank_continuation_dropped(self):
        lines = render(make_history(), lambda e, r: ['  [{}] 码 1；开始 T+1s，末次 T+2s'.format(e['ref']), '', '      续'])
        self.assertEqual(recent_refs(lines), ['H5', 'H6', 'H4', 'H7'])
        self.assertIn('        续', lines)


class HistoryLineStyleTest(unittest.TestCase):

    def test_styles(self):
        for heading in (dhl.TITLE_LINE, dhl.FIRST_HEADING, dhl.RECENT_HEADING):
            self.assertEqual(dhl.history_line_style(heading), '\033[96m')
        self.assertEqual(dhl.history_line_style('    - 码 2 [H7] 越界（bound_near）'), '\033[93m')
        self.assertEqual(dhl.history_line_style('    - [H14] 码 4（x，运行阶段）判据发生变化（未必恢复）'), '\033[93m')
        self.assertEqual(dhl.history_line_style('    - [H99] ' + dhl.MISSING_EPISODE_TEXT), '\033[93m')

    def test_unrelated_lines_not_styled(self):
        lines = render(make_history())
        for line in (lines[1], dhl.REF_NOTE_LINE, dhl.EMPTY_LINE, '        开始 T+1s，末次 T+2s',
                     '        触发值 实测 1 / 阈值 2', '        同类 4 段', '    - 当前状态：异常',
                     '    - 码 2 当前', '     - [H1] x', '      - [H1] x', '    -[H1] x',
                     '', None, 42):
            self.assertIsNone(dhl.history_line_style(line), repr(line))

    def test_never_red_and_color_off_identical(self):
        lines = render(make_history())
        colored = []
        for line in lines:
            style = dhl.history_line_style(line)
            if style:
                self.assertNotIn('[31', style)
                self.assertNotIn('[91', style)
            colored.append(dhl.apply_line_style(line, style))
        self.assertEqual([ANSI_RE.sub('', l) for l in colored], lines)
        self.assertEqual([dhl.apply_line_style(l, None) for l in lines], lines)


if __name__ == '__main__':
    unittest.main()
