#!/usr/bin/env python3
"""Version Doxygen assets, including lazily loaded navigation scripts."""
import argparse
import hashlib
from pathlib import Path
import re
from urllib.parse import urlsplit


LOADER = re.compile(r"""script\.src\s*=\s*scriptName\s*\+\s*(['"])\.js(?:\?v=[0-9a-f]+)?\1""")
ASSET = re.compile(r"""(<[A-Za-z][^>]*?\b(?:src|href)=["'])([^"']+\.(?:js|css)(?:\?[^"']*)?)(["'])""")
VERSION_META = re.compile(r'<meta name="lightning-navigation-version" content="[^"]*" */?>')


def prepare_html(folder):
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
        updated = ASSET.sub(version_url, header + separator + body)
        if updated != text:
            page.write_bytes(updated.encode('utf-8'))
    (folder / 'navigation_version.txt').write_bytes((version + '\n').encode('ascii'))
    return version


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--html', type=Path, required=True)
    args = parser.parse_args()
    print('Navigation asset version: ' + prepare_html(args.html.resolve()))
