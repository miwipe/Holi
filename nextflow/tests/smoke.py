#!/usr/bin/env python3
"""Real Nextflow orchestration tests with synthetic data and protocol mocks."""
import argparse
import csv
import gzip
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from make_fixture import create


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nextflow', default='nextflow', help='Installed executable, version 25.10.2 or newer')
    parser.add_argument('--keep', type=Path, help='Keep artifacts in this new/empty directory')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    executable = shutil.which(args.nextflow)
    if not executable:
        parser.error('Nextflow executable was not found')
    root = args.keep.resolve() if args.keep else Path(tempfile.mkdtemp(prefix='holi-nf-test-'))
    root.mkdir(exist_ok=True, parents=True)
    fixture = create(root / 'fixture')
    mock_bin = root / 'mock-bin'
    mock_bin.mkdir(exist_ok=True)
    mock = mock_bin / 'mock-tool'
    shutil.copyfile(repo / 'nextflow/tests/mock_tool.py', mock)
    mock.chmod(0o755)
    for tool in ['fastp', 'vsearch', 'sga', 'samtools', 'bowtie2', 'seqtk', 'getRTax', 'unicorn', 'metaDMG-cpp', 'compressbam']:
        (mock_bin / tool).symlink_to(mock.name)
    env = {**os.environ, 'NXF_SYNTAX_PARSER': 'v2', 'NXF_OFFLINE': 'true', 'PATH': str(mock_bin) + os.pathsep + os.environ['PATH']}
    results = []
    try:
        for mode in ['stub', 'getrtax', 'lca']:
            run = root / mode
            run.mkdir()
            settings = json.loads((fixture / 'settings.json').read_text())
            settings['microbial_assignment'] = 'lca' if mode == 'lca' else 'getrtax'
            settings_path = run / 'settings.json'
            settings_path.write_text(json.dumps(settings))
            command = [executable, 'run', str(repo / 'nextflow/main.nf'), '--input', str(fixture / 'samples.csv'),
                '--settings', str(settings_path), '--references', str(fixture / 'references'), '--outdir', str(run / 'results'),
                '--cpus', '1', '--memory', '1 GB', '-work-dir', str(run / 'work'), '-ansi-log', 'false']
            if mode == 'stub':
                command.append('-stub-run')
            with (run / 'console.log').open('w') as log:
                subprocess.run(command, cwd=run, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            for sample in ['paired', 'single']:
                for stage in range(1, 9):
                    output = run / f'results/stage{stage}/{sample}.stage{stage}'
                    provenance = json.loads((output / 'provenance.json').read_text())
                    assert provenance['stub'] == (mode == 'stub')
                if mode != 'stub':
                    with gzip.open(run / f'results/stage3/{sample}.stage3/{sample}.euk.fastq.gz', 'rt') as handle:
                        reads = handle.read()
                    assert '@eukaryotic' in reads and '@bacterial' not in reads, reads
            # Resume must reuse all 17 tasks, including validation.
            with (run / 'resume.log').open('w') as log:
                subprocess.run(command + ['-resume'], cwd=run, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            traces = sorted((run / 'results/execution').glob('*/trace.tsv'))
            with traces[-1].open() as handle:
                trace = list(csv.DictReader(handle, delimiter='\t'))
            assert len(trace) == 17 and all(row['status'] == 'CACHED' for row in trace), trace
            results.append({'mode': mode, 'tasks': 17, 'resume_cached': 17, 'samples': 2})
            print(f'PASS {mode}: paired + single, 17 tasks, 17 cached on resume', flush=True)
        (root / 'summary.json').write_text(json.dumps(results, indent=2) + '\n')
        print('Protocol mocks verify orchestration and read partitioning; they do not validate biological results.')
    except Exception:
        print(f'Failure artifacts: {root}', flush=True)
        raise
    else:
        if not args.keep:
            shutil.rmtree(root)


if __name__ == '__main__':
    main()
