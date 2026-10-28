# Taxonomic classification

Stages 4–8 combine the reference sets before filtering and reporting assignments. This competitive setup lets a read's support be compared across the configured nuclear and organelle references.

## Mapping and filtering

Stage 4 maps against the reference sets described in the script. Stage 5 combines the alignments into a competitive BAM. Stage 6 uses Unicorn `alnfilt`, `refstats`, and `bamstats`; parameters such as `ALNFILT_MODE`, `ALNFILT_MINANI`, `ALNFILT_MAXANI`, and `REFSTATS_MINREADS` are read from the config.

Read support depends on catalogue composition, alignment settings, reference redundancy, and the filtering strategy. Do not interpret the presence of a species label as independent proof of species-level specificity.

## Damage and rank-level summaries

Stage 7 runs metaDMG assignment and damage fitting. Stage 8 computes Unicorn family-, genus-, and species-level statistics. [Butteracid](Postprocessing.md) subsequently joins evidence for downstream review.

Inspect fragment lengths, identity distributions, coverage, damage uncertainty, and plausible competing assignments. Keep raw, prefiltered, mapped, and filtered denominators separate. The [historical classification page](../wiki-archive/5.-Taxonomic-Classification.md) includes older manual workflows; its alternative bam-filter instructions do not describe the current stage-6 implementation.
