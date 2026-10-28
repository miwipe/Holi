# Microbial prefiltering

Stages 2 and 3 first map the prepared reads to configured microbial references, then use assignments to split reads for subsequent analysis.

The default path uses getRTax. `--lca-assignment` selects metaDMG LCA for microbial assignment instead. The implementation expects the microbial reference shards and their corresponding taxonomy/accession mappings to be in the layout used by the existing catalogue.

Reads retained for eukaryotic analysis are a candidate set, not confirmed eukaryotic identifications. They can include reads absent from the microbial catalogue or reads without a confident microbial assignment. Eukaryotic mapping, competitive filtering, and authenticity assessment are subsequent steps.

Record the microbial catalogue build, taxonomy version, similarity settings, and read counts before and after splitting. Keep prefilter retention separate from eukaryotic mapping success when writing reports.

The [original tutorial](../wiki-archive/4.-Microbial-prefiltering.md) retains detailed manual examples and parameter explanations. For the executed commands and current defaults, use `bin/holi.sh` and `config/holi.example.yml`.
