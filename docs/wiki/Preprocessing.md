# Preprocessing

Stage 1 of `bin/holi.sh` combines fastp, vsearch, and sga. The checked-in `FASTP_*`, `VSEARCH_*`, and `SGA_*` settings in the example configuration document the current parameters.

1. **fastp:** trim and quality-filter reads; merge paired ends when running paired-end mode. Single-end mode follows its own input-discovery pattern and skips paired-read merging.
2. **vsearch:** remove duplicate sequences from the prepared reads.
3. **sga:** filter low-complexity sequences to produce the `.ppm.vs.d4.fq.gz` preprocessing output.

Review library preparation, read layout, adapter behaviour, and the input filenames before selecting parameters. Report the read set used at each step; raw-read counts and merged-fragment counts are not interchangeable.

Use `--skip-preprocessing-cleanup` if intermediate outputs are needed for diagnostics. The historical [preprocessing tutorial](../wiki-archive/3.-Data-Preprocessing-QC.md) preserves manual command examples; current automated settings are defined by the script and configuration in this branch.
