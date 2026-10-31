#!/usr/bin/env python3
"""Prepare readable manual excerpts, source provenance and versioned Doxygen assets."""
import argparse
from datetime import datetime, timezone
import hashlib
from html import escape
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import subprocess
from urllib.parse import urlsplit

from check_links import PAGE, manual_pages


LOADER = re.compile(r"""script\.src\s*=\s*scriptName\s*\+\s*(['"])\.js(?:\?v=[0-9a-f]+)?\1""")
ASSET = re.compile(r"""(<[A-Za-z][^>]*?\b(?:src|href)=["'])([^"']+\.(?:js|css)(?:\?[^"']*)?)(["'])""")
VERSION_META = re.compile(r'<meta name="lightning-navigation-version" content="[^"]*" */?>')
THEME_LINK = re.compile(r'<link id="lightning-theme"[^>]*>\n?')
BASE_STYLESHEET = re.compile(r'''<link\b[^>]*\bhref=["'](?:\.\./)*doxygen\.css(?:\?[^"']*)?["'][^>]*>\n?''')
SOURCE_BANNER = re.compile(r'<!-- lightning-source-version:start -->.*?<!-- lightning-source-version:end -->\n?', re.S)


class FragmentLinks(HTMLParser):
    """Locate generated links without reserializing code, entities or whitespace."""
    def __init__(self, text):
        super().__init__(convert_charrefs=False)
        self.text = text
        self.line_offsets = [0, *(match.end() for match in re.finditer('\n', text))]
        self.edits = []
        self.div_depth = 0
        self.fragment_depth = None
        self.tooltip = None
        self.anchors = []

    def source_offset(self):
        line, column = self.getpos()
        return self.line_offsets[line - 1] + column

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == 'div':
            self.div_depth += 1
            classes = attrs.get('class', '').split()
            if 'fragment' in classes and self.fragment_depth is None:
                self.fragment_depth = self.div_depth
            if self.fragment_depth is not None and 'ttc' in classes and self.tooltip is None:
                self.tooltip = (self.div_depth, self.source_offset())
        if tag == 'a':
            unlink = self.fragment_depth is not None and self.tooltip is None and 'href' in attrs
            # Preserve source-line/named targets even if Doxygen puts a link on the same tag.
            named = ''.join(f' {key}="{escape(attrs[key], quote=True)}"'
                            for key in ('id', 'name') if attrs.get(key))
            self.anchors.append(unlink and not named)
            if unlink:
                start = self.source_offset()
                self.edits.append((start, start + len(self.get_starttag_text()),
                                   '<a' + named + '>' if named else ''))

    def handle_endtag(self, tag):
        if tag == 'a' and self.anchors:
            if self.anchors.pop():
                start = self.source_offset()
                self.edits.append((start, self.text.index('>', start) + 1, ''))
        if tag == 'div':
            if self.tooltip is not None and self.tooltip[0] == self.div_depth:
                self.edits.append((self.tooltip[1], self.text.index('>', self.source_offset()) + 1, ''))
                self.tooltip = None
            if self.fragment_depth == self.div_depth:
                self.fragment_depth = None
            self.div_depth -= 1


def strip_fragment_links(text):
    """Unwrap inline code links and remove their generated tooltip subtrees."""
    parser = FragmentLinks(text)
    parser.feed(text)
    parser.close()
    for start, end, replacement in sorted(parser.edits, reverse=True):
        text = text[:start] + replacement + text[end:]
    return text


def manual_html_names(root):
    """Use the same hand-written page inventory as the navigation checker."""
    names = {'index.html'}
    for source in manual_pages(root):
        names.update(name + '.html' for name, _ in PAGE.findall(source.read_text(encoding='utf-8')))
    # Full source/example pages retain Doxygen's symbol navigation.
    return {name for name in names if not name.endswith(('_source.html', '-example.html'))}


def source_version(root):
    """Read the checkout being documented, not the installed binary's version."""
    version = {'commit': None, 'modified': None,
               'generated_at_utc': datetime.now(timezone.utc).isoformat(timespec='seconds')}
    try:
        version['commit'] = subprocess.check_output(
            ['git', '-C', str(root), 'rev-parse', 'HEAD'], stderr=subprocess.DEVNULL,
            text=True).strip()
        status = subprocess.check_output(
            ['git', '-C', str(root), 'status', '--porcelain', '--untracked-files=normal',
             '--', 'docs', 'src', 'scripts', 'config', 'cmake', 'CMakeLists.txt', 'package.xml'],
            stderr=subprocess.DEVNULL, text=True)
        version['modified'] = bool(status.strip())
    except (OSError, subprocess.CalledProcessError):
        # An exported source tree may have no Git metadata. Do not claim it is clean.
        pass
    return version


def stamp_source(text, version):
    text = SOURCE_BANNER.sub('', text)
    commit = version['commit']
    label = f'源码版本 {commit[:12]}' if commit else '源码版本未知（无 Git 信息）'
    if version['modified'] is True:
        label += ' · 含未提交修改'
    elif version['modified'] is None and commit:
        label += ' · 工作区状态未知'
    label += ' · 文档生成于 ' + version['generated_at_utc'] + '（UTC）'
    banner = ('<!-- lightning-source-version:start -->\n'
              '<div id="lightning-source-version" role="note" '
              'title="片段来自生成时的工作区；含未提交修改时，不等同于该提交的纯净源码。">'
              + escape(label) + '</div>\n<!-- lightning-source-version:end -->\n')
    return text.replace('</div><!-- top -->', banner + '</div><!-- top -->', 1)


def prepare_html(folder):
    root = Path(__file__).resolve().parents[2]
    # Deliver the layout and theme together, including in search-result iframes.
    theme_name = 'lightning-dark.css'
    base_style = (folder / 'doxygen.css').read_bytes()
    (folder / theme_name).write_bytes(base_style + b'\n' + Path(__file__).with_name(theme_name).read_bytes())
    provenance = source_version(root)
    manual_targets = {folder / name for name in manual_html_names(root)}
    loader_file = folder / 'navtree.js'
    loader = loader_file.read_text(encoding='utf-8')
    if len(LOADER.findall(loader)) != 1:
        raise ValueError('Unrecognized Doxygen navigation loader in navtree.js')
    digest = hashlib.sha256(Path(__file__).read_bytes())
    assets = sorted([*folder.rglob('*.js'), *folder.rglob('*.css')])
    for asset in assets:
        if asset == loader_file:
            continue
        digest.update(asset.relative_to(folder).as_posix().encode('utf-8'))
        digest.update(asset.read_bytes())
    version = digest.hexdigest()[:16]
    loader = LOADER.sub(lambda match: "script.src = scriptName+'.js?v=" + version + "'", loader)
    loader_file.write_bytes(loader.encode('utf-8'))

    def version_url(match):
        url = urlsplit(match[2])
        if url.scheme or url.netloc:
            return match[0]
        return match[1] + url._replace(query='v=' + version).geturl() + match[3]

    meta = f'<meta name="lightning-navigation-version" content="{version}"/>'
    for page in folder.rglob('*.html'):
        text = page.read_text(encoding='utf-8')
        header, separator, body = text.partition('</head>')
        if not separator:
            continue
        if VERSION_META.search(header):
            header = VERSION_META.sub(meta, header)
        else:
            header += meta + '\n'
        theme_url = '../' * (len(page.relative_to(folder).parts) - 1) + theme_name
        header = THEME_LINK.sub('', BASE_STYLESHEET.sub('', header))
        header += f'<link id="lightning-theme" href="{theme_url}" rel="stylesheet" type="text/css"/>\n'
        updated = header + separator + body
        if page in manual_targets:
            updated = strip_fragment_links(updated)
        updated = stamp_source(ASSET.sub(version_url, updated), provenance)
        if updated != text:
            page.write_bytes(updated.encode('utf-8'))
    (folder / 'navigation_version.txt').write_bytes((version + '\n').encode('ascii'))
    (folder / 'source_version.json').write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    return version


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--html', type=Path, required=True)
    args = parser.parse_args()
    print('Navigation asset version: ' + prepare_html(args.html.resolve()))
