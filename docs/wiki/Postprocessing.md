# Butteracid postprocessing

Butteracid is under active development. Its Bash wrapper collects QC metrics and invokes the Unicorn R postprocessing script, which joins rank-level Unicorn and metaDMG evidence with taxonomy.

## Files and environment

Create `envs/butteracid.yml` with Conda, then inspect `config/butteracid.example.config`. This is a Bash-style `KEY="value"` configuration, not the flat YAML-like format used by `holi.sh`.

The wrapper and R script resolve neighbouring R files relative to the working directory. Run from `bin/`:

```bash
conda activate Renv
cd bin
bash butteracid.sh ../config/butteracid.local.config
```

Prepare `butteracid.local.config` from the example first, with absolute paths to your data and output locations. The wrapper currently calls `module load samtools`; on systems without environment modules, adapt that site-specific line in your own branch before running.

`MODE="unicorn"` is the supported checked-in R path. The wrapper names a `butteracid_bamfilter.R` alternative, but that file is absent from this checkout. Do not select `MODE="bamfilter"` without providing and validating that implementation.

## Outputs and reruns

The configured `QC_OUTPUT` receives the QC table. The R script writes merged species/genus TSVs and caches under `processed_data_ss/` relative to the working directory, and creates a local taxonomy database. Existing final merged tables cause an early exit. `--force-qc` rebuilds the wrapper's QC table; it does not invalidate all R caches or final merged outputs.

Review read-stage definitions before interpreting percentages. A prefilter-retention measure is not a mapping-success measure. The checked-in scripts are the authoritative implementation for current filenames, required config variables, and joins.

The detailed [original QC page](../wiki-archive/6.1-Postprocessing-pipeline:-Butteracid-‐-QC-(under-active-development).md) and [original merging page](../wiki-archive/6.2-Postprocessing-pipeline:-Butteracid-‐-Merging-(Under-active-development).md) preserve the methods and earlier examples for comparison.
