#!/usr/bin/env python3
"""Protocol mocks. These DO NOT implement biological tools or create valid BAMs."""
import gzip
import os
import shutil
import sys
from pathlib import Path

tool = Path(sys.argv[0]).name
args = sys.argv[1:]


def value(flag):
    return args[args.index(flag) + 1]


def put(path, content='MOCK OUTPUT - NOT SCIENTIFIC DATA\n'):
    if path.endswith('.gz'):
        with gzip.open(path, 'wt') as handle:
            handle.write(content)
    else:
        Path(path).write_text(content)


if os.environ.get('HOLI_MOCK_FAIL') == tool:
    sys.exit(42)
if tool == 'fastp':
    with gzip.open(value('-i'), 'rt') as source:
        put(value('--merged_out') if '--merged_out' in args else value('-o'), source.read())
    put(value('-h'))
    put(value('-j'), '{}')
elif tool == 'vsearch':
    shutil.copyfile(value('--fastx_uniques'), value('--fastqout'))
elif tool == 'sga':
    shutil.copyfile(args[args.index('-o') - 1], value('-o'))
elif tool == 'bowtie2':
    print('MOCK SAM')
elif tool == 'samtools':
    if args[0] == 'view':
        sys.stdin.buffer.read()
        print('MOCK BAM')
    elif args[0] == 'sort':
        shutil.copyfile(args[-1], value('-o'))
    elif args[0] == 'merge':
        put(value('-f'))
    elif args[0] == 'quickcheck':
        assert Path(args[-1]).stat().st_size > 0
    else:
        raise ValueError(args)
elif tool == 'getRTax':
    put(value('-p') + '.bacteria.gz', 'bacterial\n')
elif tool == 'seqtk':
    ids = set(Path(args[-1]).read_text().splitlines())
    with gzip.open(args[-2], 'rt') as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            record = header + ''.join(handle.readline() for _ in range(3))
            if header[1:].split()[0] in ids:
                sys.stdout.write(record)
elif tool == 'compressbam':
    shutil.copyfile(value('--input'), value('--output'))
elif tool == 'unicorn':
    if args[0] in ('alnfilt', 'refstats'):
        shutil.copyfile(value('-b'), value('--outbam'))
        put(value('--outstat') + '.gz')
    else:
        shutil.copyfile(value('-b'), value('-o'))
        put(value('--outstat') + '.gz')
elif tool == 'metaDMG-cpp':
    prefix = value('--out_prefix')
    if args[0] == 'lca':
        for suffix in ['.stat.gz', '.bdamage.gz']:
            put(prefix + suffix)
        put(prefix + '.lca.gz', 'bacterial\tBacteria\n')
        if '--filtered_acc2tax' in args:
            put(value('--filtered_acc2tax'))
    elif args[0] == 'dfit':
        put(prefix + '.dfit.gz')
    elif args[0] == 'aggregate':
        put(prefix + '.stat.gz')
    else:
        raise ValueError(args)
else:
    raise ValueError(f'Unsupported mock tool: {tool}')
