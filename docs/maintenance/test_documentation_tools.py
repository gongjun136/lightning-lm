"""Regression checks for source excerpts, inline navigation and the version notice."""
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from check_links import check_snippets
from prepare_html import manual_html_names, prepare_html, source_version, stamp_source, strip_fragment_links


class SnippetChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for directory in ('docs', 'src/core', 'scripts'):
            (self.root / directory).mkdir(parents=True)
        self.source = self.root / 'src/core/queue.h'
        self.source.write_text('// [capacity]\npop_front();\n// [capacity]\n', encoding='utf-8')
        self.page = self.root / 'docs/queue.md'
        self.page.write_text('@snippet{lineno} core/queue.h capacity\n', encoding='utf-8')

    def test_real_source_and_api_comment_are_checked(self):
        (self.root / 'src/api.h').write_text(
            '/**\n * @snippet core/queue.h capacity\n */', encoding='utf-8')
        self.assertEqual(check_snippets(self.root), [])
        self.source.write_text('// [capacity]\npop_front();\n', encoding='utf-8')
        self.assertEqual(len(check_snippets(self.root)), 2)

    def test_renamed_file_ambiguous_basename_and_empty_excerpt(self):
        self.page.write_text('@snippet missing.h capacity\n', encoding='utf-8')
        self.assertTrue(check_snippets(self.root))
        self.page.write_text('@snippet queue.h capacity\n', encoding='utf-8')
        (self.root / 'src/queue.h').write_text(self.source.read_text(), encoding='utf-8')
        self.assertTrue(check_snippets(self.root))
        self.page.write_text('@snippet core/queue.h capacity\n', encoding='utf-8')
        self.source.write_text('// [capacity]\n// [capacity]\n', encoding='utf-8')
        self.assertTrue(check_snippets(self.root))

    def test_fenced_instruction_is_not_a_live_reference(self):
        self.page.write_text('```text\n@snippet nonexistent.h example\n```\n', encoding='utf-8')
        self.assertEqual(check_snippets(self.root), [])


class SourceVersionChecks(unittest.TestCase):
    def test_dirty_checkout_is_visible_and_reprocessing_replaces_notice(self):
        page = '<div id="top">title</div><!-- top --><div id="doc-content">code</div>'
        version = {'commit': 'abcdef1234567890', 'modified': True,
                   'generated_at_utc': '2026-10-24T00:00:00+00:00'}
        stamped = stamp_source(page, version)
        self.assertIn('含未提交修改', stamped)
        self.assertEqual(stamped.count('id="lightning-source-version"'), 1)
        self.assertEqual(stamp_source(stamped, version), stamped)
        version['modified'] = False
        self.assertNotIn(' · 含未提交修改', stamp_source(stamped, version))
        self.assertIn('<div id="doc-content">code</div>', stamped)

    def test_no_git_information_never_claims_a_clean_commit(self):
        with patch('prepare_html.subprocess.check_output', side_effect=FileNotFoundError):
            version = source_version(Path('.'))
        self.assertIsNone(version['commit'])
        self.assertIsNone(version['modified'])
        self.assertIn('源码版本未知', stamp_source('</div><!-- top -->', version))

    def test_status_failure_preserves_commit_but_marks_state_unknown(self):
        with patch('prepare_html.subprocess.check_output', side_effect=[
                'abcdef1234567890\n', subprocess.CalledProcessError(1, 'git')]):
            version = source_version(Path('.'))
        self.assertEqual(version['commit'], 'abcdef1234567890')
        self.assertIsNone(version['modified'])
        self.assertIn('工作区状态未知', stamp_source('</div><!-- top -->', version))


class FragmentLinkChecks(unittest.TestCase):
    def test_only_fragment_links_and_nested_tooltips_are_removed(self):
        prefix = '<p><a href="api.html#method">API 方法</a></p>'
        line = ('<div class="line"><a name="l00524"></a><span class="lineno">524</span>&#160;'
                '<span class="keyword">return</span> '
                '<a class="code" href="wrong.html#exp"><span class="name">std::exp</span></a>'
                '(x &lt; 2); &quot;中文 &amp; &lt;a href=example&gt;&quot;</div>')
        tooltip = ('<div class="ttc" id="wrong"><div class="ttname"><a href="wrong.html">exp</a>'
                   '</div><div class="ttdef"><b>Definition:</b> wrong.cc</div></div>')
        outside = '<div class="ttc"><a href="kept.html">outside tooltip</a></div>'
        page = prefix + '<div class="fragment">' + line + tooltip + '</div>' + outside
        expected_line = line.replace('<a class="code" href="wrong.html#exp">', '').replace('</span></a>', '</span>')
        expected = prefix + '<div class="fragment">' + expected_line + '</div>' + outside
        self.assertEqual(strip_fragment_links(page), expected)
        self.assertEqual(strip_fragment_links(expected), expected)

    def test_named_link_targets_and_multiple_fragments_survive(self):
        page = ('<div class="fragment"><div class="line">'
                '<a href="wrong.html" name="line_target" id="line_id">frame</a>'
                '</div></div><p><a href="source_source.html#l00001">完整源码</a></p>'
                '<div class="fragment"><div class="line"><a href="wrong.html">count</a></div></div>')
        expected = page.replace('<a href="wrong.html" name="line_target" id="line_id">',
                                '<a id="line_id" name="line_target">').replace(
                                    '<a href="wrong.html">count</a>', 'count')
        self.assertEqual(strip_fragment_links(page), expected)

    def test_preparation_discovers_manual_pages_and_preserves_full_source_navigation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            docs = root / 'docs'
            docs.mkdir()
            (docs / 'Doxyfile.in').write_text('INPUT = docs\n', encoding='utf-8')
            (docs / 'index.md').write_text('# Manual\n', encoding='utf-8')
            (docs / 'new.md').write_text('@page new_topic 新主题\n', encoding='utf-8')
            (docs / 'imported.dox').write_text('@page imported_topic 导入主题\n', encoding='utf-8')
            (docs / 'shell_sources.dox').write_text('@example{lineno} run.sh\n', encoding='utf-8')
            self.assertEqual(manual_html_names(root), {'index.html', 'new_topic.html', 'imported_topic.html'})
            folder = root / 'html'
            folder.mkdir()
            (folder / 'navtree.js').write_text("script.src = scriptName+'.js';\n", encoding='utf-8')
            base_style = 'body { margin: 0; }\n'
            (folder / 'doxygen.css').write_text(base_style, encoding='utf-8')
            fragment = ('<div class="fragment"><div class="line"><a name="l00001"></a>'
                        '<span class="lineno">1</span><a class="code" href="wrong.html">x</a></div>'
                        '<div class="ttc"><div class="ttname">tooltip</div></div></div>')
            page = ('<html><head><link href="doxygen.css?v=old" rel="stylesheet"/></head>'
                    '<body><div id="top"></div><!-- top -->' + fragment + '</body></html>')
            names = ['index.html', 'new_topic.html', 'imported_topic.html',
                     'queue_8h_source.html', 'run_8sh-example.html', 'classQueue.html']
            for name in names:
                (folder / name).write_text(page, encoding='utf-8')
            inventory = manual_html_names(root)
            version = {'commit': None, 'modified': None, 'generated_at_utc': '2026-10-25T00:00:00+00:00'}
            with patch('prepare_html.manual_html_names', return_value=inventory), \
                    patch('prepare_html.source_version', return_value=version):
                prepare_html(folder)
                first = {name: (folder / name).read_text(encoding='utf-8') for name in names}
                prepare_html(folder)
            self.assertTrue((folder / 'lightning-dark.css').read_text(encoding='utf-8').startswith(base_style))
            for name in names:
                result = (folder / name).read_text(encoding='utf-8')
                self.assertEqual(result, first[name], name)
                self.assertEqual(result.count('id="lightning-theme"'), 1, name)
                self.assertNotIn('href="doxygen.css', result, name)
                if name in inventory:
                    self.assertIn(strip_fragment_links(fragment), result)
                    self.assertNotIn('href="wrong.html"', result)
                    self.assertNotIn('class="ttc"', result)
                else:
                    self.assertIn(fragment, result, name)


if __name__ == '__main__':
    unittest.main()
