#!/usr/bin/env python3
"""Guard Chinese production debug descriptions and their machine-readable fields."""
import ast
from collections import defaultdict
from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
LITERAL = re.compile(r'"(?:\\.|[^"\\])*"')
TOKEN = re.compile(r'"(?:\\.|[^"\\])*"|//[^\n]*|/\*[\s\S]*?\*/')
CALL = re.compile(r'\bdebug_event::(EmitThrottled|ReportState|Emit)\s*\(')
CHINESE = re.compile(r'[\u4e00-\u9fff]')


def arguments(text, start):
    result, begin, depth, index = [], start, 0, start
    while index < len(text):
        token = TOKEN.match(text, index)
        if token:
            index = token.end()
            continue
        char = text[index]
        if char in '([{':
            depth += 1
        elif char in ')]}':
            if char == ')' and depth == 0:
                result.append(text[begin:index].strip())
                return result
            depth -= 1
        elif char == ',' and depth == 0:
            result.append(text[begin:index].strip())
            begin = index + 1
        index += 1
    raise AssertionError('Unclosed debug-event call')


def descriptions(text):
    # Ignore comments while preserving quoted text and line positions.
    text = TOKEN.sub(lambda match: match[0] if match[0].startswith('"')
                     else re.sub(r'[^\n]', ' ', match[0]), text)
    for match in CALL.finditer(text):
        args = arguments(text, match.end())
        kind = match[1]
        key = args[0] if kind != 'Emit' else None
        bodies = ([args[0]] if kind == 'Emit' else
                  [args[1]] if kind == 'EmitThrottled' else args[2:4])
        for body in bodies:
            if '.str()' in body:
                variable = re.fullmatch(r'(\w+)\.str\(\)', body)
                assert variable, 'Unsupported debug description: ' + body
                begin = text.rfind('std::ostringstream ' + variable[1] + ';', 0, match.start())
                assert begin >= 0, 'Cannot find debug description builder: ' + body
                body = text[begin:match.start()]
            literals = [ast.literal_eval(value) for value in LITERAL.findall(body)]
            assert literals, 'Missing debug description: ' + body
            yield text.count('\n', 0, match.start()) + 1, key, literals


class DebugEventLanguageTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.entries = []
        for source in (ROOT / 'src').rglob('*'):
            if source.suffix not in ('.cc', '.cpp', '.cxx', '.h', '.hpp'):
                continue
            if source.relative_to(ROOT / 'src').parts[0] == 'test':
                continue
            text = source.read_text(encoding='utf-8')
            if not CALL.search(text):
                continue
            for line, key, literals in descriptions(text):
                cls.entries.append((str(source.relative_to(ROOT)) + ':' + str(line), key, literals))

    def test_every_production_description_is_chinese(self):
        self.assertTrue(self.entries, 'No production debug events were inspected')
        for location, _, literals in self.entries:
            with self.subTest(location=location):
                self.assertRegex(''.join(literals), CHINESE)
                for literal in literals:
                    if CHINESE.search(literal):
                        continue
                    # Punctuation, units, and existing key=value fragments stay ASCII.
                    machine_text = re.sub(r'\b[A-Za-z_]\w*\s*=', '', literal)
                    machine_text = re.sub(r'\b(?:sec|ms|m|rad|deg)\b', '', machine_text)
                    self.assertNotRegex(machine_text, r'[A-Za-z]', literal)

    def test_internal_event_keys_remain_machine_identifiers(self):
        for location, key, _ in self.entries:
            if key is None:
                continue
            with self.subTest(location=location):
                literals = LITERAL.findall(key)
                self.assertTrue(literals)
                for literal in literals:
                    self.assertRegex(ast.literal_eval(literal), r'^[a-z][a-z0-9_]*$')

    def test_numeric_and_enum_field_names_are_preserved(self):
        fields = defaultdict(set)
        for _, key, literals in self.entries:
            name = ast.literal_eval(LITERAL.findall(key)[0]) if key else 'Emit'
            fields[name].update(re.findall(r'\b([A-Za-z_]\w*)=', ''.join(literals)))
        expected = {
            'long_imu_propagation_interval': {'interval_sec', 'warning_threshold_sec'},
            'stale_lidar_correction': {'age_sec', 'hard_deadline_sec'},
            'adaptive_lidar_load_step': {'processing_ms', 'latency_ms', 'tracking_healthy'},
            'fused_lidar_timestamp_rollback': {'timestamp', 'last_timestamp'},
            'smoother_motion_jump': {'distance_m', 'allowed_distance_m'},
            'pgo_dr_timestamp_rollback': {'rollback_sec'},
            'pgo_lidar_odom_timestamp_rollback': {'rollback_sec'},
            'pgo_lidar_loc_parked_epoch_drop': {'age_sec'},
            'pgo_lidar_loc_timestamp_rollback': {'rollback_sec'},
            'processed_sensor_timestamp_rollback': {'rollback_sec', 'timestamp'},
            'imu_timestamp_rollback': {'rollback_sec'},
            'live_output_timestamp_rollback': {'event', 'timestamp', 'last_accepted', 'rollback_sec'},
            'pose_vel_timestamp_rollback': {'event', 'timestamp', 'last_accepted', 'rollback_sec'},
            'stale_lidar_input_': {'event', 'lidar_id', 'timestamp', 'newest_timestamp',
                                  'lag_sec', 'configured_max_lag_sec'},
            'Emit': {'reason', 'candidate', 'score', 'search_ms', 'confirmations',
                     'retrieval_score', 'refinement_backend', 'refinement_score', 'overlap', 'count'},
        }
        for name, required in expected.items():
            with self.subTest(event=name):
                self.assertTrue(required <= fields[name], (name, fields[name], required))


if __name__ == '__main__':
    unittest.main()
