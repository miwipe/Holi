#!/usr/bin/env python3
"""Holi stage commands. Arguments are passed as argv, never evaluated by a shell."""
import argparse
import gzip
import json
import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path


def open_text(path):
    with open(path, 'rb') as handle:
        compressed = handle.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if compressed else open(path)


class Stage:
    def __init__(self, args):
        self.args = args
        self.sample = args.sample
        self.config = json.loads(args.settings.read_text())
        self.params = json.loads(Path(__file__).with_name('defaults.json').read_text())
        self.params.update(self.config.get('parameters', {}))
        self.tools = {'metadmg': 'metaDMG-cpp', **self.config.get('tools', {})}
        self.root = args.references.resolve()
        self.up = args.upstream.resolve() if args.upstream else None
        self.reads = [p.resolve() for p in args.reads]
        self.cpus = args.cpus
        # samtools -m is per thread; leave room for its main thread and overhead.
        self.sort_mem = f'{max(1, args.memory_mb // (self.cpus + 2))}M'
        self.out = Path(f'{self.sample}.stage{args.stage}')
        self.out.mkdir()
        self.commands = []

    def tool(self, name):
        value = self.tools.get(name, name)
        found = shutil.which(value)
        if not found:
            raise ValueError(f'Missing executable {name}: {value}. Activate the analysis environment or set tools.{name}.')
        return found

    def argv(self, command):
        return [self.tool(command[0]), *map(str, command[1:])]

    def run(self, command, stdout=None):
        argv = self.argv(command)
        self.commands.append(shlex.join(argv))
        with open('commands.log', 'a') as log:
            print('$ ' + shlex.join(argv), file=log, flush=True)
            subprocess.run(argv, check=True, stdout=stdout or log, stderr=log)

    def pipe(self, first, second, output):
        a, b = self.argv(first), self.argv(second)
        self.commands.append(shlex.join(a) + ' | ' + shlex.join(b))
        with open('commands.log', 'a') as log, open(output, 'wb') as target:
            print('$ ' + self.commands[-1], file=log, flush=True)
            with subprocess.Popen(a, stdout=subprocess.PIPE, stderr=log) as producer:
                try:
                    consumer = subprocess.run(b, stdin=producer.stdout, stdout=target, stderr=log)
                finally:
                    producer.stdout.close()
                producer_rc = producer.wait()
                if producer_rc:
                    raise subprocess.CalledProcessError(producer_rc, a)
                consumer.check_returncode()

    def p(self, key):
        return str(self.params[key])

    def taxonomy(self, key):
        return self.root / self.config['taxonomy'][key]

    def tax_args(self, microbial=False):
        prefix = 'microbial_' if microbial else ''
        return ['--names', self.taxonomy(prefix + 'names'), '--nodes', self.taxonomy(prefix + 'nodes')]

    def sort(self, source, target, name=False):
        self.run(['samtools', 'sort', *(['-n'] if name else []), '-@', max(0, self.cpus - 1),
                  '-m', self.sort_mem, '-o', target, source])

    def check_bam(self, path):
        self.run(['samtools', 'quickcheck', '-v', path])

    def merge(self, inputs, output, name=False):
        self.run(['samtools', 'merge', '-@', max(0, self.cpus - 1), *(['-n'] if name else []), '-f', output, *inputs])
        self.check_bam(output)

    def compress(self, source, target):
        self.run(['compressbam', '--threads', self.cpus, '--input', source, '--output', target])
        self.check_bam(target)

    def preprocess(self):
        s = self.sample
        command = ['fastp', '-i', self.reads[0]]
        command += ['-I', self.reads[1], '-m', '--merged_out', f'{s}.ppm.fq', '--detect_adapter_for_pe'] if len(self.reads) == 2 else ['-o', f'{s}.ppm.fq']
        command += ['-V', '-D', '--dup_calc_accuracy', self.p('FASTP_DUP_ACCURACY'), '-g', '-x',
                    '-q', self.p('FASTP_QUAL'), '-e', self.p('FASTP_MIN_AVG_QUAL'), '-l', self.p('FASTP_MIN_LEN'),
                    '-y', '-c', '-p', '-h', f'{s}.fastp.report.html', '-j', f'{s}.fastp.report.json', '-w', self.cpus]
        self.run(command)
        self.run(['vsearch', '--fastx_uniques', f'{s}.ppm.fq', '--fastqout', f'{s}.ppm.vs.fq',
                  '--minseqlength', self.p('VSEARCH_MIN_SEQLEN'), '--strand', 'both', '--threads', self.cpus])
        self.run(['sga', 'preprocess', '--dust-threshold=' + self.p('SGA_DUST_THRESHOLD'), '-m', self.p('SGA_MIN_LEN'),
                  f'{s}.ppm.vs.fq', '-o', f'{s}.ppm.vs.d4.fq'])
        with open(f'{s}.ppm.vs.d4.fq', 'rb') as source, gzip.open(f'{s}.ppm.vs.d4.fq.gz', 'wb') as target:
            shutil.copyfileobj(source, target)
        for path in Path('.').glob('*.fq'):
            path.unlink()

    def mapping(self, microbial):
        s = self.sample
        input_fastq = self.up / (f'{s}.ppm.vs.d4.fq.gz' if microbial else f'{s}.euk.fastq.gz')
        bams = []
        for i, prefix in enumerate(self.config['microbial_indexes' if microbial else 'eukaryotic_indexes'], 1):
            bam = f'{s}.{"gtdb" if microbial else "euk"}.{i}.bam'
            command = ['bowtie2', '--threads', self.cpus, '-x', self.root / prefix, '-U', input_fastq,
                       '-k', self.p('BT2_GTDB_k' if microbial else 'BT2_EUK_k'), '--mm', '--no-unal']
            if microbial:
                for flag in ('D', 'R', 'N', 'L', 'i', 'np', 'mp', 'rdg', 'rfg', 'score_min'):
                    command += [('-' if len(flag) == 1 else '--') + flag.replace('_', '-'), self.p('BT2_GTDB_' + flag)]
            else:
                command += ['-t']
            self.pipe(command, ['samtools', 'view', '-bS', '-'], bam)
            self.check_bam(bam)
            if microbial:
                sorted_bam = bam.replace('.bam', '.sorted.bam')
                self.sort(bam, sorted_bam)
                Path(bam).unlink()
                bams.append(sorted_bam)
        if microbial:
            self.merge(bams, f'{s}.gtdb.merged.bam')
            for path in bams:
                Path(path).unlink()

    def split(self):
        s = self.sample
        bam = self.up / f'{s}.gtdb.merged.bam'
        ids = f'{s}.bact_reads.raw.txt'
        if self.config.get('microbial_assignment', 'getrtax') == 'lca':
            self.sort(bam, f'{s}.gtdb.namesorted.bam', name=True)
            self.lca(f'{s}.gtdb.namesorted.bam', f'{s}.microbial', microbial=True)
            with open_text(f'{s}.microbial.lca.gz') as source, open(ids, 'w') as target:
                for line in source:
                    if re.search(r'Archaea|virus|bacteria', line, re.I):
                        target.write(line.split('\t')[0].split()[0] + '\n')
        else:
            prefix = f'{s}.getrtax'
            self.run(['getRTax', '--bam', bam, '-T', self.taxonomy('microbial_taxid_map'),
                      '-r', json.dumps({'domain': ['d__Bacteria', 'd__Archaea', 'd__Viruses']}),
                      '--threads', self.cpus, '--unique', '--only-read-ids', '-p', prefix])
            files = sorted(Path('.').glob(prefix + '*'))
            if not files:
                raise ValueError('getRTax did not write read-ID files; cannot treat a failed classification as zero hits')
            with open(ids, 'w') as target:
                for path in files:
                    with open_text(path) as source:
                        for line in source:
                            if line.strip():
                                target.write(line.split()[0] + '\n')
        # External sort/comm keep memory bounded for large read sets.
        with open(f'{s}.bact_reads.txt', 'w') as target:
            self.run(['sort', '-u', ids], stdout=target)
        with open_text(self.reads[0]) as source, open('all_reads.raw.txt', 'w') as target:
            while True:
                header = source.readline()
                if not header:
                    break
                sequence, plus, quality = (source.readline() for _ in range(3))
                if not header.startswith('@') or not plus.startswith('+') or len(sequence.rstrip()) != len(quality.rstrip()) or not quality:
                    raise ValueError('Malformed four-line FASTQ after preprocessing')
                target.write(header[1:].split()[0] + '\n')
        with open('all_reads.txt', 'w') as target:
            self.run(['sort', '-u', 'all_reads.raw.txt'], stdout=target)
        with open(f'{s}.euk_reads.txt', 'w') as target:
            self.run(['comm', '-23', 'all_reads.txt', f'{s}.bact_reads.txt'], stdout=target)
        for ids_file, output in [(f'{s}.bact_reads.txt', f'{s}.bact_reads.fq.gz'), (f'{s}.euk_reads.txt', f'{s}.euk.fastq.gz')]:
            self.pipe(['seqtk', 'subseq', self.reads[0], ids_file], ['gzip', '-c'], output)

    def merge_euk(self):
        bams = []
        for i in range(1, len(self.config['eukaryotic_indexes']) + 1):
            prefix = f'{self.sample}.euk.{i}'
            self.compress(self.up / f'{prefix}.bam', f'{prefix}.comp.bam')
            self.sort(f'{prefix}.comp.bam', f'{prefix}.comp.sorted.bam', name=True)
            bams.append(f'{prefix}.comp.sorted.bam')
        merged = f'{self.sample}.merged.bam'
        self.merge(bams, merged, name=True)
        self.compress(merged, f'{self.sample}.comp.bam')
        for path in [merged, *bams, *[f'{self.sample}.euk.{i}.comp.bam' for i in range(1, len(bams) + 1)]]:
            Path(path).unlink()

    def filter(self):
        s = self.sample
        self.run(['unicorn', 'alnfilt', '-b', self.up / f'{s}.comp.bam', '-t', self.cpus,
                  '--mode', self.p('ALNFILT_MODE'), '--outbam', f'{s}.alnfilt.bam', '--outstat', f'{s}.alnfilt.refstats',
                  '--minani', self.p('ALNFILT_MINANI'), '--maxani', self.p('ALNFILT_MAXANI')])
        self.check_bam(f'{s}.alnfilt.bam')
        self.run(['unicorn', 'refstats', '-b', f'{s}.alnfilt.bam', '-t', self.cpus,
                  '--minreads', self.p('REFSTATS_MINREADS'), '--outbam', f'{s}.alnfilt.unicorn.bam',
                  '--outstat', f'{s}.comp.unicorn.refstats', *self.tax_args()])
        self.check_bam(f'{s}.alnfilt.unicorn.bam')

    def combined_taxonomy(self):
        target = Path('combined.acc2taxid')
        with target.open('w') as output:
            for value in self.config['taxonomy']['acc2tax']:
                with open_text(self.root / value) as source:
                    shutil.copyfileobj(source, output)
        return target

    def lca(self, bam, prefix, microbial=False):
        scope = 'METADMG_BAC' if microbial else 'METADMG_EUK'
        command = ['metadmg', 'lca', *self.tax_args(microbial), '--acc2tax',
                   self.taxonomy('microbial_acc2tax') if microbial else self.combined_taxonomy(),
                   '--fix_ncbi', 0, '--threads', self.cpus, '--bam', bam, '--out_prefix', prefix]
        for flag in ('sim_score_low', 'sim_score_high', 'how_many', 'weight_type'):
            command += ['--' + flag, self.p(scope + '_' + flag.upper())]
        if not microbial:
            command += ['--filtered_acc2tax', f'{self.sample}.acc2tax']
        self.run(command)

    def damage(self):
        prefix = f'{self.sample}.sort.comp.filtered'
        # Match Bash: both final analyses consume alnfilt, not the refstats BAM.
        self.sort(self.up / f'{self.sample}.alnfilt.bam', prefix + '.bam', name=True)
        self.lca(prefix + '.bam', prefix)
        self.run(['metadmg', 'dfit', prefix + '.bdamage.gz', '--threads', self.cpus, *self.tax_args(),
                  '--showfits', 2, '--nopt', self.p('METADMG_DFIT_NOPT'), '--nbootstrap', self.p('METADMG_DFIT_NBOOTSTRAP'),
                  '--seed', self.p('METADMG_DFIT_SEED'), '--lib', self.p('METADMG_DFIT_LIB'), '--out_prefix', prefix])
        self.run(['metadmg', 'aggregate', prefix + '.bdamage.gz', *self.tax_args(), '--lcastat', prefix + '.stat.gz',
                  '--dfit', prefix + '.dfit.gz', '--out_prefix', prefix + '.agg'])
        Path('combined.acc2taxid').unlink()

    def taxstats(self):
        table = self.combined_taxonomy()
        for rank in ('family', 'genus', 'species'):
            prefix = f'{self.sample}.comp.filtered.{rank}.taxstats'
            self.run(['unicorn', 'taxstats', '-b', self.up / f'{self.sample}.alnfilt.bam', '-t', self.cpus,
                      '-o', prefix + '.bam', *self.tax_args(), '--acc2tax', table, '-k', self.p('TAXSTATS_k'),
                      '--outstat', prefix, '--minreads', self.p('TAXSTATS_MINREADS'), '--rank', rank])
            self.check_bam(prefix + '.bam')
        table.unlink()

    def execute(self):
        os.chdir(self.out)
        os.environ['LC_ALL'] = 'C'
        try:
            if self.args.stub:
                Path('STUB_ONLY.txt').write_text('Synthetic workflow wiring test; NO biological results.\n')
                # Downstream stubs consume directories, not their contents.
            else:
                {1: self.preprocess, 2: lambda: self.mapping(True), 3: self.split, 4: lambda: self.mapping(False),
                 5: self.merge_euk, 6: self.filter, 7: self.damage, 8: self.taxstats}[self.args.stage]()
            for name in required_outputs(self.args.stage, self.sample, len(self.config['eukaryotic_indexes'])) if not self.args.stub else []:
                choices = (name,) if isinstance(name, str) else name
                if not any(Path(choice).is_file() and Path(choice).stat().st_size > 0 for choice in choices):
                    raise ValueError(f'Stage {self.args.stage} did not produce required output: {name}')
        finally:
            Path('provenance.json').write_text(json.dumps({'sample': self.sample, 'stage': self.args.stage,
                'stub': self.args.stub, 'cpus': self.cpus, 'settings': self.config, 'resolved_parameters': self.params,
                'commands': self.commands}, indent=2) + '\n')


def required_outputs(stage, s, count):
    return {1: [f'{s}.ppm.vs.d4.fq.gz', f'{s}.fastp.report.html', f'{s}.fastp.report.json'],
            2: [f'{s}.gtdb.merged.bam'], 3: [f'{s}.euk.fastq.gz', f'{s}.bact_reads.fq.gz'],
            4: [f'{s}.euk.{i}.bam' for i in range(1, count + 1)], 5: [f'{s}.comp.bam'],
            6: [f'{s}.alnfilt.bam', f'{s}.alnfilt.unicorn.bam',
                (f'{s}.alnfilt.refstats', f'{s}.alnfilt.refstats.gz'),
                (f'{s}.comp.unicorn.refstats', f'{s}.comp.unicorn.refstats.gz')],
            7: [f'{s}.sort.comp.filtered.bdamage.gz', f'{s}.sort.comp.filtered.stat.gz', f'{s}.sort.comp.filtered.dfit.gz',
                f'{s}.sort.comp.filtered.agg.stat.gz'],
            8: [item for rank in ('family', 'genus', 'species') for item in
                [f'{s}.comp.filtered.{rank}.taxstats.bam',
                 (f'{s}.comp.filtered.{rank}.taxstats', f'{s}.comp.filtered.{rank}.taxstats.gz')]]}[stage]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', type=int, choices=range(1, 9), required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--settings', type=Path, required=True)
    parser.add_argument('--references', type=Path, required=True)
    parser.add_argument('--upstream', type=Path)
    parser.add_argument('--reads', type=Path, nargs='*', default=[])
    parser.add_argument('--cpus', type=int, required=True)
    parser.add_argument('--memory-mb', type=int, required=True)
    parser.add_argument('--stub', action='store_true')
    args = parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', args.sample) or args.cpus < 1 or args.memory_mb < 1:
        parser.error('Invalid sample ID or resources')
    Stage(args).execute()


if __name__ == '__main__':
    main()
