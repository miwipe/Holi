#!/usr/bin/env python3
"""Export proposed wiki pages without publishing or changing the live wiki."""
import argparse
from pathlib import Path
import re
from urllib.parse import quote

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('destination', type=Path)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
source = root / 'docs/wiki'
destination = args.destination.resolve()
if destination.exists() and any(destination.iterdir()):
    parser.error('Use an empty export directory to avoid overwriting an existing wiki checkout.')
destination.mkdir(parents=True, exist_ok=True)
for page in source.glob('*.md'):
    text = page.read_text()
    # Main-repository archive filenames contain balanced parentheses; handle those
    # first, then the simple local wiki page links.
    text = re.sub(r'\]\((\.\./wiki-archive/.*?\.md)\)', lambda m: '](' + 'https://github.com/miwipe/Holi/blob/holi_dev/docs/wiki-archive/' + quote(m[1].split('/')[-1]) + ')', text)
    text = text.replace('](../DEVELOPMENT.md)', '](https://github.com/miwipe/Holi/blob/holi_dev/docs/DEVELOPMENT.md)')
    text = re.sub(r'\]\(([A-Za-z0-9_-]+)\.md(#[^)]*)?\)', lambda m: '](' + m[1] + (m[2] or '') + ')', text)
    (destination / page.name).write_text(text)
print(f'Prepared {len(list(source.glob("*.md")))} wiki pages in {destination}; nothing published.')
