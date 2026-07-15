# Holi — Snakemake workflow

This is a Snakemake port of `bin/holi.sh`, the 8-stage shotgun metagenomics
pipeline. `bin/holi.sh` (bash + GNU `parallel`, hardcoded cluster paths,
manual stage/resume flags) is still present and untouched; this workflow is
an alternative way to run the same stages that's easier to troubleshoot:

- **Every step is its own rule** with explicit inputs/outputs, so a failure
  points at exactly which file couldn't be produced, instead of a wall of
  interleaved `parallel` output.
- **Resuming is automatic.** Snakemake only reruns rules whose outputs are
  missing or older than their inputs — there's no equivalent of
  `--from-stage`/`--force-stage`/auto-detection to reason about.
- **Cleanup is automatic.** Intermediate files are marked `temp()` and get
  deleted once nothing downstream needs them anymore — this replaces all of
  `holi.sh`'s `--skip-*-cleanup` / `--storage-friendly` flags. Pass
  `--notemp` to Snakemake if you want to keep everything for inspection
  while debugging.
- **All tool/database paths live in one place:** [`config/config.yaml`](../config/config.yaml).
- Stage 1 preprocessing outputs now land under `results/preprocessing/`
  instead of the working directory `holi.sh` wrote them to.

Microbial classification (`bin/Holi_microbe_classify.sh`) and QC reporting
(`bin/butteracid.sh`) were **not** converted — this workflow only covers the
main `holi.sh` pipeline (stages 1–8).

## Layout

```
config/
  config.yaml   # all paths and parameters — edit this before running
  samples.txt   # one sample name per line
workflow/
  Snakefile
  rules/*.smk   # one file per pipeline stage
  envs/holi.yaml  # conda env for the common bioinformatics tools
```

## Setup

1. Install Snakemake (e.g. `conda create -n snakemake -c conda-forge -c bioconda snakemake` or `pip install snakemake`).
2. Install the tools in `workflow/envs/holi.yaml` (fastp, vsearch, sga,
   bowtie2, samtools, seqtk, ...), e.g. `conda env create -f workflow/envs/holi.yaml`.
3. `metaDMG-cpp`, `unicorn`, `getRTax`, and `compressbam` are not on
   conda/bioconda — install them separately and point `config.yaml`'s
   `tools:` section at their binaries (absolute paths are fine if they're
   not on `PATH`).
4. Edit `config/config.yaml`: set `inpath`, the output paths, all `db_path_*`
   / `tax_path_*` entries, and `single_end` / `lca_assignment` as needed.
5. List your samples in `config/samples.txt`.

## Running

```bash
# See what would run without doing anything:
snakemake -n --cores 16

# Visualize the DAG for a sanity check:
snakemake --dag | dot -Tpng > dag.png

# Actually run:
snakemake --cores 16

# Run only up to a given stage, e.g. just through GTDB mapping:
snakemake --cores 16 results/microbial/mysample.gtdb.merged.bam

# Keep intermediates instead of auto-deleting them (useful while debugging):
snakemake --cores 16 --notemp
```

`--cores` controls how many rule instances (samples/chunks) run
concurrently; each rule's own thread count comes from
`threads_per_job` in `config.yaml`.
