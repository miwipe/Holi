import csv
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'runner'))
from prepare import validate_samples, validate_settings
from stage import Stage
from make_fixture import create


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = create(Path(self.temp.name))
        self.settings = json.loads((self.root / 'settings.json').read_text())

    def tearDown(self):
        self.temp.cleanup()

    def test_paired_and_single_inputs(self):
        rows = validate_samples(self.root / 'samples.csv')
        self.assertEqual(len(rows), 2)
        self.assertTrue(rows[0]['fastq_2'])
        self.assertEqual(rows[1]['fastq_2'], '')
        validate_settings(self.settings, self.root / 'references')

    def test_scientific_defaults_match_bash_example(self):
        repo = Path(__file__).resolve().parents[2]
        defaults = json.loads((repo / 'nextflow/runner/defaults.json').read_text())
        expected = {}
        for line in (repo / 'config/holi.example.yml').read_text().splitlines():
            if not line or line.startswith('#') or ':' not in line:
                continue
            key, value = line.split(':', 1)
            if key.startswith(('FASTP_', 'VSEARCH_', 'SGA_', 'BT2_', 'METADMG_', 'ALNFILT_', 'REFSTATS_', 'TAXSTATS_')):
                try:
                    value = json.loads(value.strip())
                except ValueError:
                    value = value.strip()
                expected[key] = value
        self.assertEqual(defaults, expected)

    def test_duplicate_or_unsafe_sample_ids(self):
        for identifier in ['paired', 'bad;touch pwned', '../escape']:
            path = self.root / 'samples.csv'
            create(self.root)
            with path.open() as source:
                rows = list(csv.reader(source))
            rows[2][0] = identifier
            with path.open('w', newline='') as target:
                csv.writer(target).writerows(rows)
            with self.assertRaises(ValueError):
                validate_samples(path)

    def test_missing_and_partial_indexes(self):
        (self.root / 'references/microbe.rev.2.bt2').unlink()
        with self.assertRaisesRegex(ValueError, 'Incomplete Bowtie2'):
            validate_settings(self.settings, self.root / 'references')

    def test_unknown_parameters_and_bad_ranges(self):
        for params in [{'BT2_TYPO': 1}, {'ALNFILT_MINANI': 101}, {'METADMG_EUK_SIM_SCORE_LOW': 2}, {'BT2_EUK_k': 2.5}]:
            self.settings['parameters'] = params
            with self.assertRaises(ValueError):
                validate_settings(self.settings, self.root / 'references')

    def test_reference_escape_rejected(self):
        self.settings['taxonomy']['names'] = '../outside'
        with self.assertRaisesRegex(ValueError, 'relative'):
            validate_settings(self.settings, self.root / 'references')

    def test_lca_requires_its_taxonomy(self):
        self.settings['microbial_assignment'] = 'lca'
        del self.settings['taxonomy']['microbial_nodes']
        with self.assertRaises(ValueError):
            validate_settings(self.settings, self.root / 'references')

    def test_mapper_failure_is_not_hidden_by_consumer_success(self):
        stage = Stage.__new__(Stage)
        stage.commands = []
        stage.argv = lambda command: command
        previous = Path.cwd()
        try:
            os.chdir(self.root)
            with self.assertRaises(subprocess.CalledProcessError) as raised:
                stage.pipe([sys.executable, '-c', 'import sys; sys.exit(42)'],
                           [sys.executable, '-c', 'import sys; sys.stdin.read()'], 'output.bam')
            self.assertEqual(raised.exception.returncode, 42)
        finally:
            os.chdir(previous)

    @unittest.skipUnless(all(shutil.which(tool) for tool in ('bowtie2', 'bowtie2-build', 'samtools')),
                         'Optional real Bowtie2/samtools test needs installed executables')
    def test_real_mapping_creates_readable_bam(self):
        refs = self.root / 'real-refs'
        refs.mkdir()
        fasta = refs / 'tiny.fa'
        fasta.write_text('>REF\n' + 'A' * 20 + 'ACGT' * 10 + 'C' * 20 + '\n')
        subprocess.run(['bowtie2-build', str(fasta), str(refs / 'tiny')],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        self.settings['microbial_indexes'] = ['tiny']
        settings = self.root / 'real-settings.json'
        settings.write_text(json.dumps(self.settings))
        upstream = self.root / 'preprocessed'
        upstream.mkdir()
        shutil.copyfile(next(self.root.glob('single*.fq.gz')), upstream / 'single.ppm.vs.d4.fq.gz')
        args = argparse.Namespace(sample='single', stage=2, settings=settings, references=refs,
            upstream=upstream, reads=[], cpus=1, memory_mb=256, stub=False)
        previous = Path.cwd()
        try:
            os.chdir(self.root)
            Stage(args).execute()
            bam = self.root / 'single.stage2/single.gtdb.merged.bam'
            count = subprocess.check_output(['samtools', 'view', '-c', str(bam)], text=True)
            self.assertGreater(int(count.strip()), 0)
        finally:
            os.chdir(previous)


if __name__ == '__main__':
    unittest.main()
