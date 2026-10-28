# Pipeline workflow

![Holi workflow](https://raw.githubusercontent.com/miwipe/Holi/main/figures/Holi_Workflow2025.png)

The current `bin/holi.sh` implementation runs the following stages in order. Paths below are representative outputs under the configured directories; consult the script for complete intermediate filenames.

| Stage | Operation | Representative output |
| --- | --- | --- |
| 1 | fastp trimming/merging, vsearch duplicate removal, sga complexity filtering | `<sample>.ppm.vs.d4.fq.gz` |
| 2 | Mapping to the configured microbial reference shards, sorting and merging | `<sample>.gtdb.merged.bam` |
| 3 | Read assignment and microbial/eukaryotic splitting | `<sample>.euk.fastq.gz` |
| 4 | Mapping to nuclear, mitochondrial, plastid, and additional configured catalogues | Reference-set alignment BAMs |
| 5 | Competitive sorting and merging | `<sample>.comp.bam` |
| 6 | Unicorn alignment filtering, reference statistics, BAM statistics | `<sample>.comp.filtered.bam` |
| 7 | metaDMG LCA, damage fitting, and aggregation | `<sample>.sort.comp.filtered.agg.stat.gz` |
| 8 | Unicorn rank-level taxonomic statistics | Family, genus, and species taxstats |

The default microbial split uses getRTax. `--lca-assignment` selects the metaDMG LCA alternative. The reference sets and shard counts are tied to the current indexed-catalogue layout.

Stage 3 retention, stage 4 mapping, and stage 6 filtering describe different read sets. Keep their denominators separate when reporting retention, mapping rates, and taxonomic support. See [Authentication](Authentication.md) for interpretation and [Running Holi](Running-Holi.md) for execution.
