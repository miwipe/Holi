# Configuration examples

- `holi.example.yml`: flat `KEY: value` input for `bin/holi.sh`. Keep keys at column one and replace every `/path/to/` placeholder.
- `samples.example.list`: one sample identifier per line; replace the two examples with identifiers matching your FASTQ filenames.
- `nextflow.example.json` and `nextflow.samples.example.csv`: explicit reference/settings and sample manifests for the [Nextflow runner](../nextflow/README.md). Copy to `*.local.json` and `*.local.csv`.
- `butteracid.example.config`: Bash-style literal values for the postprocessing wrapper and its R parser. Use absolute paths and `MODE="unicorn"`.

Copy examples to `*.local.yml`, `*.local.list`, or `*.local.config`; those files are ignored. The existing project-specific configs in `bin/` are retained for compatibility and are not portable examples. See [Running Holi](../docs/wiki/Running-Holi.md) and [Postprocessing](../docs/wiki/Postprocessing.md).
