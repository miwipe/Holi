# Ancient DNA authentication

Holi provides evidence for authenticity assessment; it does not replace interpretation of the sample, controls, reference catalogue, and preservation environment.

Review these lines of evidence together:

- **Fragment lengths:** examine distributions for the assigned reads and consider library preparation and mapping selection.
- **Damage patterns:** inspect terminal misincorporation estimates, uncertainty, library treatment, and read support.
- **Reference coverage:** check breadth, evenness, local pile-ups, and whether the supported regions distinguish plausible alternatives.
- **Taxonomic specificity:** inspect competitive alignments and catalogue redundancy before drawing species-level conclusions.
- **Controls and context:** consider laboratory controls, depositional history, contamination routes, and ecological plausibility.

The automated damage analysis is in stage 7 (`metaDMG-cpp`). Fragment-length and coverage information comes from the upstream alignments and statistics. Reads retained after microbial filtering are not a measure of eukaryotic mapping success.

The archived authentication page contains exploratory PMD tools and references to `scripts/readlengthplot_from_bam.sh` and `scripts/split_bam_by_PMD.sh`. Those helper scripts are not distributed in this checkout. They are not required steps in the documented automatic workflow.
