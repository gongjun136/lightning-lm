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
    pages = list((root / 'docs').rglob('*.md'))
    ids = []
    references = []
    for page in pages:
        text = page.read_text(encoding='utf-8')
        ids.extend(re.findall(r'^@page\s+(\w+)', text, re.M))
        # Fenced and inline code are examples, not navigation instructions.
        prose = re.sub(r'```.*?```|`[^`\n]*`', '', text, flags=re.S)
        children = re.findall(r'@subpage\s+(\w+)', prose)
        if children and page != root / 'docs/index.md':
            errors.append(f'{page}: only index.md may declare subpages')
        references.extend(children)
    duplicate = [name for name, count in Counter(ids).items() if count > 1]
    errors.extend(f'duplicate page id: {name}' for name in duplicate)
    errors.extend(f'unknown subpage: {name}' for name in references if name not in ids)
    errors.extend(f'page missing from index: {name}' for name in ids if name not in references)
    errors.extend(f'repeated subpage: {name}' for name,count in Counter(references).items() if count > 1)
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
    # Audit all manual pages and every directly linked API/source/resource target.
    for name in ['index', *page_ids]:
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
