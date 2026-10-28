#!/usr/bin/env python3
"""Check local links in maintained Markdown and compatibility paths."""
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit

root = Path(__file__).resolve().parents[1]
pages = [root/'README.md', root/'docs/README.md', root/'docs/DEVELOPMENT.md']
pages += list((root/'docs/wiki').glob('*.md'))
pages += [root/'config/README.md', root/'literature/README.md', root/'data/literature/README.md', root/'dashboard/vendor/README.md']
errors = []
links = 0
for page in pages:
    text = page.read_text()
    text = re.sub(r'```.*?```', '', text, flags=re.S)
    for match in re.finditer(r'\]\(', text):
        i = match.end(); start = i; depth = 1
        while i < len(text) and depth:
            if text[i] == '(':
                depth += 1
            elif text[i] == ')':
                depth -= 1
            i += 1
        target = text[start:i-1].strip().strip('<>')
        if not target or target.startswith('#') or urlsplit(target).scheme:
            continue
        target = unquote(target.split('#')[0])
        links += 1
        if not (page.parent/target).exists():
            errors.append(f'{page.relative_to(root)}: {target}')
for target in ['bin/holi.yml','bin/butteracid.yml','sedaDNA_Shotgun_CaptureEnrichment_litterature/sedaDNA_lit_map/generate_figures.R']:
    if not (root/target).exists():
        errors.append(f'Broken compatibility path: {target}')
if errors:
    raise SystemExit('\n'.join(errors))
print(f'Passed: {len(pages)} maintained Markdown files, {links} local links, and compatibility paths.')
