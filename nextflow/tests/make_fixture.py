#!/usr/bin/env python3
"""Tiny synthetic fixtures; fake indexes are suitable ONLY for stub/mock tests."""
import csv
import gzip
import json
import sys
from pathlib import Path


def create(destination):
    root = destination.resolve()
    root.mkdir(parents=True, exist_ok=True)
    refs = root / 'references'
    refs.mkdir(exist_ok=True)
    for prefix in ['microbe', 'nuclear', 'organelle']:
        for part in ['1', '2', '3', '4', 'rev.1', 'rev.2']:
            (refs / f'{prefix}.{part}.bt2').write_text('FAKE INDEX FOR WIRING TEST ONLY\n')
    for name in ['names.dmp', 'nodes.dmp', 'taxid.map']:
        (refs / name).write_text('FAKE TAXONOMY FOR WIRING TEST ONLY\n')
    with gzip.open(refs / 'acc2tax.gz', 'wt') as target:
        target.write('accession\taccession.version\ttaxid\tgi\nREF\tREF.1\t2\t0\n')
    rows = []
    for sample, paired in [('paired', True), ('single', False)]:
        reads = []
        for mate in range(1, 3 if paired else 2):
            path = root / f"{sample} sample's R{mate}.fq.gz"
            with gzip.open(path, 'wt') as target:
                target.write('@bacterial description\n' + 'ACGT' * 10 + '\n+\n' + 'I' * 40 + '\n')
                target.write('@eukaryotic description\n' + 'TGCA' * 10 + '\n+\n' + 'I' * 40 + '\n')
            reads.append(str(path))
        rows.append([sample, reads[0], reads[1] if paired else ''])
    with (root / 'samples.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['sample', 'fastq_1', 'fastq_2'])
        writer.writerows(rows)
    settings = {'microbial_assignment': 'getrtax', 'microbial_indexes': ['microbe'],
        'eukaryotic_indexes': ['nuclear', 'organelle'], 'taxonomy': {
            'microbial_taxid_map': 'taxid.map', 'microbial_names': 'names.dmp', 'microbial_nodes': 'nodes.dmp',
            'microbial_acc2tax': 'acc2tax.gz', 'names': 'names.dmp', 'nodes': 'nodes.dmp', 'acc2tax': ['acc2tax.gz']},
        'parameters': {}}
    (root / 'settings.json').write_text(json.dumps(settings, indent=2))
    return root


if __name__ == '__main__':
    print(create(Path(sys.argv[1])))
