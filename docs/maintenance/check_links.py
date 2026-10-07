#!/usr/bin/env python3
"""Validate manual navigation/relative links and, optionally, generated HTML targets."""
import argparse
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit


def check_sources(root):
    errors = []
    pages = sorted((root / 'docs').rglob('*.md'))
    ids = []
    references = []
    navigation = {}
    guide_ids = set()
    index = root / 'docs/index.md'
    for page in pages:
        text = page.read_text(encoding='utf-8')
        page_ids = re.findall(r'^@page\s+(\w+)', text, re.M)
        ids.extend(page_ids)
        if page.parent == root / 'docs/guide':
            guide_ids.update(page_ids)
        if page != index and len(page_ids) != 1:
            errors.append(f'{page}: expected one stable page id')
        # Fenced and inline code are examples, not navigation instructions.
        prose = re.sub(r'```.*?```|`[^`\n]*`', '', text, flags=re.S)
        children = re.findall(r'@subpage\s+(\w+)', prose)
        if children and page != index and page.parent != root / 'docs/guide':
            errors.append(f'{page}: only index.md and guide pages may declare subpages')
        if page == index or page_ids:
            navigation['index' if page == index else page_ids[0]] = children
        references.extend(children)
    duplicate = [name for name, count in Counter(ids).items() if count > 1]
    errors.extend(f'duplicate page id: {name}' for name in duplicate)
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

    for page in pages + [root/'README.md', root/'README_CN.md']:
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

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if attrs.get('id'): self.ids.add(attrs['id'])
        if tag == 'a' and attrs.get('name'): self.ids.add(attrs['name'])
        for key in ('href', 'src', 'data'):
            if key in attrs: self.targets.append(attrs[key])


def check_html(folder, page_ids):
    errors, cache = [], {}
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
        for target in parse(page).targets:
            url = urlsplit(target)
            if url.scheme or url.netloc: continue
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
        for page in (root/'docs').rglob('*.md'):
            ids.extend(re.findall(r'^@page\s+(\w+)',page.read_text(encoding='utf-8'),re.M))
        errors.extend(check_html(args.html.resolve(), ids))
    for error in errors: print(error)
    print(f'{count} manual pages; {len(errors)} link/navigation errors')
    return bool(errors)


if __name__ == '__main__':
    raise SystemExit(main())
