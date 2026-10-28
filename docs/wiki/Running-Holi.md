# Running Holi

Holi supports the existing Bash runner described below and an [experimental Nextflow alternative](../../nextflow/README.md). Both cover the eight main stages. Nextflow uses explicit CSV/JSON manifests and `-resume`; its execution tests do not yet establish scientific equivalence on real data.

Run commands below from the repository root. `bin/holi.sh` is the current entry point; names such as `holi_prefilter.sh` and `test_tools_safe.sh` in the historical wiki do not identify files in this checkout.

## Prepare the input and configuration

```bash
cp config/holi.example.yml config/holi.local.yml
cp config/samples.example.list config/samples.local.list
```

Replace all `/path/to/` placeholders and the example sample identifiers. The list has one sample identifier per line. For paired-end inputs, the script searches under `INPATH` using `<sample>*R1*.fastq.gz` and `<sample>*R2*.fastq.gz`. Single-end mode uses a different pattern, `*_<sample>_*R1*_001.fastq.gz`. Inspect the preprocessing block before preparing inputs with different names or lane structures.

The config parser reads flat `KEY: value` lines. Keys must begin at column one, with no indentation or inline comments on value lines. It is not a general YAML parser. Preserve the existing `DB_PATH_Norwary` spelling: the implementation expects that key. Keep output and input paths free of whitespace for the analysis scripts.

Important settings:

| Settings | Meaning |
| --- | --- |
| `INPATH` | Directory containing input FASTQs |
| `RESULT_PATH`, `MICROB_OUT`, `EUK_OUT` | Output directories |
| `LOG_FILE`, `LOGS`, `TMP` | Master log, task logs, temporary data |
| `THREADSP`, `THREADS` | Parallel jobs and threads per job; budget their product against allocated CPUs |
| `DB_PATH*`, `TAX_PATH*` | Indexed reference and taxonomy locations |
| `FASTP_*`, `BT2_*`, `ALNFILT_*`, `METADMG_*`, `TAXSTATS_*` | Stage parameters; see `config/holi.example.yml` |

The mapping loops assume the existing catalogue layout, including fixed reference shards and organelle targets. A path change does not adapt the workflow to arbitrary new index naming or shard counts.

## Execute and resume

```bash
bash bin/holi.sh config/holi.local.yml config/samples.local.list
bash bin/holi.sh config/holi.local.yml config/samples.local.list --from-stage 4 --to-stage 6
bash bin/holi.sh config/holi.local.yml config/samples.local.list --single-end
```

The script uses existing output files to infer progress. Inspect logs and files before trusting an automatically resumed run; a file's presence alone does not establish scientific validity. `--force-stage N` overrides the inferred starting point. Options come after the configuration and sample-list arguments.

Cleanup flags control whether intermediates remain available:

```bash
bash bin/holi.sh config/holi.local.yml config/samples.local.list \
  --skip-preprocessing-cleanup --skip-gtdb-cleanup --skip-euk-cleanup
```

`--storage-friendly` requests deletion of stage intermediates. Choose cleanup intentionally when planning restarts or troubleshooting. The usage text in `bin/holi.sh` lists the complete skip options and stage controls; calling it with no arguments prints that usage and exits with status 1.

For a cluster job, use your institution's scheduler configuration, activate the environment in the job shell, and allocate CPUs and memory for the chosen reference sets. No generic Slurm submission script is included or implied to be validated here.
