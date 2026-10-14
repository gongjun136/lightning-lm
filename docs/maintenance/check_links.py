#!/usr/bin/env python3
"""Validate manual navigation/relative links and, optionally, generated HTML targets."""
import argparse
from collections import Counter
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import shlex
from urllib.parse import parse_qs, unquote, urlsplit


PAGE = re.compile(r'^\s*(?:\*\s*)?@page\s+(\w+)[ \t]*(.*)$', re.M)
TOP_LEVEL_MARKDOWN = {'index.md', 'architecture.md'}
TOPIC_DIRECTORIES = {'guide', 'getting_started', 'flows', 'modules',
                     'algorithms', 'reference', 'operations', 'maintenance'}


def manual_pages(root):
    """Include imported Markdown and .dox pages, not just docs/*.md."""
    docs = root / 'docs'
    candidates = set(docs.rglob('*.md')) | set(docs.rglob('*.dox'))
    config = re.sub(r'\\\r?\n', ' ', (docs / 'Doxyfile.in').read_text(encoding='utf-8'))
    inputs = re.search(r'^INPUT\s*=\s*(.*)$', config, re.M)
    for entry in shlex.split(inputs.group(1), comments=True) if inputs else []:
        source = Path(entry.replace('@LIGHTNING_SOURCE_DIR@', root.as_posix()))
        if not source.is_absolute():
            source = root / source
        if source.is_dir():
            candidates.update(source.rglob('*.md'))
            candidates.update(source.rglob('*.dox'))
        elif source.suffix in {'.md', '.dox'} and source.is_file():
            candidates.add(source)
    return sorted(page for page in candidates
                  if page.suffix == '.md' or PAGE.search(page.read_text(encoding='utf-8')))


def plain_title(title):
    title = re.sub(r'<[^>]+>', '', title)
    title = re.sub(r'[*_\x60]', '', title)
    return re.sub(r'\s+', ' ', title).strip()


def subpages(text):
    prose = re.sub(r'```.*?```|`[^`\n]*`', '', text, flags=re.S)
    return re.findall(r'@subpage\s+(\w+)', prose)


def navigation_parents(root):
    parents = {}
    for page in manual_pages(root):
        text = page.read_text(encoding='utf-8')
        declaration = PAGE.search(text)
        parent = 'index' if page == root / 'docs/index.md' else (
            declaration[1] if declaration else None)
        for child in subpages(text):
            parents[child] = parent
    return parents


def check_sources(root):
    errors = []
    pages = manual_pages(root)
    ids = []
    titles = []
    references = []
    navigation = {}
    guide_ids = set()
    index = root / 'docs/index.md'
    for page in pages:
        text = page.read_text(encoding='utf-8')
        declarations = PAGE.findall(text)
        page_ids = [name for name, title in declarations]
        ids.extend(page_ids)
        if not page.is_relative_to(root / 'docs'):
            errors.append(f'{page}: manual page outside docs; keep README as an entry link')
        elif page.parent == root / 'docs':
            if page.name not in TOP_LEVEL_MARKDOWN:
                errors.append(f'{page}: top-level topic must move into a documentation subdirectory')
        elif page.relative_to(root / 'docs').parts[0] not in TOPIC_DIRECTORIES:
            errors.append(f'{page}: unknown documentation category')
        for name, title in declarations:
            normalized = plain_title(title)
            titles.append(normalized)
            if not normalized or title != normalized:
                errors.append(f'{page}: page title must be nonempty plain text: {name}')
        if page.parent == root / 'docs/guide':
            guide_ids.update(page_ids)
        if page != index and len(page_ids) != 1:
            errors.append(f'{page}: expected one stable page id')
        # Fenced and inline code are examples, not navigation instructions.
        children = subpages(text)
        if children and page != index and page.parent != root / 'docs/guide':
            errors.append(f'{page}: only index.md and guide pages may declare subpages')
        if page == index or page_ids:
            navigation['index' if page == index else page_ids[0]] = children
        references.extend(children)
    duplicate = [name for name, count in Counter(ids).items() if count > 1]
    errors.extend(f'duplicate page id: {name}' for name in duplicate)
    errors.extend(f'duplicate page title: {title}'
                  for title, count in Counter(titles).items() if count > 1)
    if 'index' in ids:
        errors.append('reserved main page id: index')
    errors.extend(f'index subpage must be a guide page: {name}'
                  for name in navigation.get('index', []) if name not in guide_ids)
    errors.extend(f'unknown subpage: {name}' for name in references if name not in ids)
    errors.extend(f'page missing from navigation: {name}' for name in ids if name not in references)
    errors.extend(f'repeated subpage: {name}' for name,count in Counter(references).items() if count > 1)

    reached = set()
    def visit(name, ancestors):
        if name in ancestors:
            errors.append('navigation cycle: ' + ' -> '.join((*ancestors, name)))
            return
        if name in reached:
            return
        reached.add(name)
        for child in navigation.get(name, []):
            visit(child, (*ancestors, name))
    visit('index', ())
    unreachable = sorted(set(ids) - reached)
    errors.extend(f'page unreachable from index: {name}' for name in unreachable)
    # Inspect disconnected components too, so an orphaned cycle is diagnosed.
    for name in unreachable:
        visit(name, ())

    # Shell examples include the real files at build time, without copying source.
    registry = root / 'docs/shell_sources.dox'
    if not registry.is_file():
        errors.append('missing shell source registry: docs/shell_sources.dox')
    else:
        examples = re.findall(r'@example(?:\{[^}]*\})?\s+(\S+)',
                              registry.read_text(encoding='utf-8'))
        scripts = {path.relative_to(root / 'scripts').as_posix()
                   for path in (root / 'scripts').rglob('*.sh')}
        errors.extend(f'duplicate shell source: {name}'
                      for name, count in Counter(examples).items() if count > 1)
        errors.extend(f'shell source missing from registry: {name}'
                      for name in sorted(scripts - set(examples)))
        errors.extend(f'unknown shell source: {name}'
                      for name in sorted(set(examples) - scripts))

    for page in pages + [root/'README.md', root/'README_CN.md',
                         root/'scripts/README.md', root/'config/sany/README.md']:
        text = re.sub(r'```.*?```', '', page.read_text(encoding='utf-8'), flags=re.S)
        for target in re.findall(r'!?\[[^\]]*\]\(([^\s)]+)\)', text):
            url = urlsplit(target)
            if url.scheme or url.netloc or not url.path:
                continue
            if not (page.parent / unquote(url.path)).exists():
                errors.append(f'{page.relative_to(root)}: missing {target}')
    return errors, len(pages)


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()
        self.targets = []
        self.navigation_version = None

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if attrs.get('id'): self.ids.add(attrs['id'])
        if tag == 'a' and attrs.get('name'): self.ids.add(attrs['name'])
        if tag == 'meta' and attrs.get('name') == 'lightning-navigation-version':
            self.navigation_version = attrs.get('content')
        for key in ('href', 'src', 'data'):
            if key in attrs: self.targets.append(attrs[key])


def check_navigation(folder, page_ids, parents):
    errors, seen = [], Counter()

    def read_array(name):
        try:
            if not re.fullmatch(r'\w+', name):
                raise ValueError('invalid script name')
            text = (folder / (name + '.js')).read_text(encoding='utf-8')
            match = re.search(r'\bvar\s+\w+\s*=\s*(\[.*?\]);', text, re.S)
            data = json.loads(match[1]) if match else None
            if not isinstance(data, list):
                raise ValueError('missing array')
            return data
        except (OSError, ValueError):
            errors.append(f'invalid or missing rendered navigation: {name}.js')
            return []

    def walk(nodes, parent, ancestors):
        for node in nodes:
            if not isinstance(node, list) or len(node) != 3 or not isinstance(node[1], str):
                errors.append('invalid rendered navigation node')
                continue
            title, target, children = node
            url = urlsplit(target)
            name = Path(url.path).stem
            child_parent = parent
            if not url.fragment:
                if name == 'index' or name in page_ids:
                    seen[name] += 1
                    expected = parents.get(name)
                    if parent != expected:
                        errors.append(f'rendered parent mismatch: {name}: expected {expected}, got {parent}')
                    child_parent = name
                elif name.startswith('md_'):
                    errors.append(f'unknown rendered manual topic: {title} ({target})')
                    continue
                else:
                    # API indexes are separate from the hand-written manual.
                    continue
            if isinstance(children, str):
                if children in ancestors:
                    errors.append(f'rendered navigation cycle: {children}')
                else:
                    walk(read_array(children), child_parent, (*ancestors, children))
            elif isinstance(children, list):
                walk(children, child_parent, ancestors)

    tree = read_array('navtreedata')
    if tree and len(tree[0]) == 3 and isinstance(tree[0][2], list):
        walk(tree[0][2], None, ('navtreedata',))
    else:
        errors.append('missing rendered manual navigation')
    for name in ['index', *page_ids]:
        if seen[name] != 1:
            errors.append(f'rendered manual page must occur once: {name} ({seen[name]} occurrences)')
    return errors


def check_html(folder, page_ids, parents):
    errors, cache = check_navigation(folder, page_ids, parents), {}
    version_file = folder / 'navigation_version.txt'
    version = version_file.read_text(encoding='ascii').strip() if version_file.is_file() else ''
    if not re.fullmatch(r'[0-9a-f]{16}', version):
        errors.append('missing navigation asset version; rebuild the docs target')
    else:
        loader_file = folder / 'navtree.js'
        loader = loader_file.read_text(encoding='utf-8') if loader_file.is_file() else ''
        if ".js?v=" + version not in loader:
            errors.append('navigation loader uses an outdated asset version')
    # Doxygen renders unresolved Markdown references as plain text, so they
    # disappear from the HTML link list. Inspect the build log as well.
    warning_log = folder.parent.parent / 'doxygen-warnings.log'
    if warning_log.is_file():
        errors.extend(line for line in warning_log.read_text(encoding='utf-8').splitlines()
                      if re.search(r'[/\\]docs[/\\].*: warning:', line))
    def parse(path):
        if path not in cache:
            parser = Links()
            parser.feed(path.read_text(encoding='utf-8'))
            cache[path] = parser
        return cache[path]
    # Audit manual and shell pages, plus their API/source/resource targets.
    source_pages = sorted(path.stem for path in folder.glob('*-example.html'))
    for name in ['index', 'examples', *page_ids, *source_pages]:
        page = folder / (name + '.html')
        if not page.is_file():
            errors.append(f'missing rendered page: {page.name}')
            continue
        if parse(page).navigation_version != version:
            errors.append(f'{page.name}: navigation asset version mismatch')
        for target in parse(page).targets:
            url = urlsplit(target)
            if url.scheme or url.netloc: continue
            if Path(url.path).suffix in {'.js', '.css'} and parse_qs(url.query).get('v') != [version]:
                errors.append(f'{page.name}: unversioned or outdated asset {target}')
            path = (page.parent/unquote(url.path)).resolve() if url.path else page
            if not path.is_file():
                errors.append(f'{page.name}: missing {target}')
            elif url.fragment and path.suffix == '.html' and unquote(url.fragment) not in parse(path).ids:
                errors.append(f'{page.name}: missing anchor {target}')
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--html', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    # maintenance -> docs -> repository
    errors, count = check_sources(root)
    if args.html:
        ids = []
        for page in manual_pages(root):
            ids.extend(name for name, title in PAGE.findall(page.read_text(encoding='utf-8')))
        errors.extend(check_html(args.html.resolve(), ids, navigation_parents(root)))
    for error in errors: print(error)
    print(f'{count} manual pages; {len(errors)} link/navigation errors')
    return bool(errors)


if __name__ == '__main__':
    raise SystemExit(main())
