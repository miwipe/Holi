#!/usr/bin/env python3
"""Validate the explicit Nextflow manifests before biological tasks start."""
import argparse
import csv
import hashlib
import json
import re
import shutil
from pathlib import Path

INDEX_PARTS = ('1', '2', '3', '4', 'rev.1', 'rev.2')


def reference(root, value):
    if not isinstance(value, str) or not value or Path(value).is_absolute() or '..' in Path(value).parts:
        raise ValueError(f'Reference paths must be relative to --references: {value!r}')
    return root / value


def validate_settings(settings, root):
    allowed = {'microbial_assignment', 'microbial_indexes', 'eukaryotic_indexes', 'taxonomy', 'tools', 'parameters'}
    if set(settings) - allowed:
        raise ValueError(f'Unknown settings: {sorted(set(settings) - allowed)}')
    mode = settings.get('microbial_assignment', 'getrtax')
    if mode not in ('getrtax', 'lca'):
        raise ValueError('microbial_assignment must be getrtax or lca')
    for key in ('microbial_indexes', 'eukaryotic_indexes'):
        indexes = settings.get(key)
        if not isinstance(indexes, list) or not indexes or len(set(indexes)) != len(indexes):
            raise ValueError(f'{key} must be a nonempty list of distinct Bowtie2 prefixes')
        for value in indexes:
            prefix = reference(root, value)
            if not any(all(Path(f'{prefix}.{part}.{ext}').is_file() for part in INDEX_PARTS)
                       for ext in ('bt2', 'bt2l')):
                raise ValueError(f'Incomplete Bowtie2 index (need all six .bt2 or .bt2l files): {value}')
    tax = settings.get('taxonomy', {})
    keys = ['names', 'nodes', 'acc2tax']
    keys += ['microbial_taxid_map'] if mode == 'getrtax' else ['microbial_names', 'microbial_nodes', 'microbial_acc2tax']
    for key in keys:
        values = tax.get(key)
        if key == 'acc2tax':
            if not isinstance(values, list) or not values:
                raise ValueError('taxonomy.acc2tax must be a nonempty list')
        else:
            values = [values]
        for value in values:
            if not reference(root, value).is_file():
                raise ValueError(f'Missing taxonomy.{key}: {value}')
    defaults = json.loads(Path(__file__).with_name('defaults.json').read_text())
    overrides = settings.get('parameters', {})
    if set(overrides) - set(defaults):
        raise ValueError(f'Unknown scientific parameters: {sorted(set(overrides) - set(defaults))}')
    for key, value in overrides.items():
        if isinstance(defaults[key], (int, float)):
            if isinstance(value, bool) or not isinstance(value, (int, float)) or value < 0:
                raise ValueError(f'{key} must be a nonnegative number')
            if isinstance(defaults[key], int) and not isinstance(value, int):
                raise ValueError(f'{key} must be an integer')
        elif not isinstance(value, str) or not value:
            raise ValueError(f'{key} must be a nonempty string')
    params = {**defaults, **overrides}
    for prefix in ('ALNFILT',):
        if not 0 <= params[f'{prefix}_MINANI'] <= params[f'{prefix}_MAXANI'] <= 100:
            raise ValueError('ALNFILT ANI bounds must satisfy 0 <= min <= max <= 100')
    for prefix in ('METADMG_BAC', 'METADMG_EUK'):
        if not 0 <= params[f'{prefix}_SIM_SCORE_LOW'] <= params[f'{prefix}_SIM_SCORE_HIGH'] <= 1:
            raise ValueError(f'{prefix} similarity bounds must satisfy 0 <= low <= high <= 1')
    tools = settings.get('tools', {})
    expected = {'fastp', 'vsearch', 'sga', 'bowtie2', 'samtools', 'seqtk', 'getRTax', 'unicorn', 'metadmg', 'compressbam'}
    if set(tools) - expected or any(not isinstance(v, str) or not v for v in tools.values()):
        raise ValueError('tools must map known tool names to executable paths/names, without arguments')


def validate_samples(path):
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ['sample', 'fastq_1', 'fastq_2']:
            raise ValueError('CSV header must be sample,fastq_1,fastq_2')
        rows = list(reader)
    if not rows:
        raise ValueError('Sample sheet is empty')
    seen = set()
    for row in rows:
        sample = row['sample']
        if None in row or any(v is None for v in row.values()):
            raise ValueError('Each sample row must have exactly three CSV fields')
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sample) or sample in seen:
            raise ValueError(f'Invalid or duplicate sample ID: {sample!r}')
        seen.add(sample)
        reads = [row['fastq_1']] + ([row['fastq_2']] if row['fastq_2'] else [])
        if not row['fastq_1'] or len({Path(p).name for p in reads}) != len(reads):
            raise ValueError(f'{sample}: reads need distinct basenames and a nonempty fastq_1')
        for value in reads:
            file = Path(value)
            if not file.is_absolute() or not file.is_file() or file.stat().st_size == 0:
                raise ValueError(f'{sample}: FASTQ must be an existing, nonempty absolute path: {value}')
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--settings', type=Path, required=True)
    parser.add_argument('--references', type=Path, required=True)
    parser.add_argument('--check-tools', action='store_true')
    args = parser.parse_args()
    settings = json.loads(args.settings.read_text())
    validate_settings(settings, args.references)
    rows = validate_samples(args.input)
    tools = {}
    if args.check_tools:
        expected = ['fastp', 'vsearch', 'sga', 'bowtie2', 'samtools', 'seqtk', 'unicorn', 'metadmg', 'compressbam', 'sort', 'comm', 'gzip']
        if settings.get('microbial_assignment', 'getrtax') == 'getrtax':
            expected.append('getRTax')
        for key in expected:
            executable = settings.get('tools', {}).get(key, 'metaDMG-cpp' if key == 'metadmg' else key)
            found = shutil.which(executable)
            if not found:
                raise ValueError(f'Missing executable {key}: {executable}. Activate your environment or configure tools.{key}.')
            checksum = hashlib.sha256()
            with open(found, 'rb') as source:
                for chunk in iter(lambda: source.read(1024 * 1024), b''):
                    checksum.update(chunk)
            tools[key] = {'path': found, 'sha256': checksum.hexdigest()}
    with open('validated.csv', 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['sample', 'fastq_1', 'fastq_2'])
        writer.writeheader()
        writer.writerows(rows)
    Path('validation.json').write_text(json.dumps({'samples': len(rows), 'settings': settings, 'tools': tools,
        'validation': 'Manifest paths and index components only; not biological reference or tool compatibility validation.'}, indent=2) + '\n')


if __name__ == '__main__':
    main()
