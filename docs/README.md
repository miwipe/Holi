# Holi documentation

Start with the workflow that matches your task. The command examples refer to files in this checkout; detailed historical methods remain available in the original wiki snapshot.

| Topic | Guide |
| --- | --- |
| Installation and dependencies | [Getting started](wiki/Getting-started.md) |
| Analysis stages and outputs | [Pipeline workflow](wiki/Pipeline-workflow.md) |
| Inputs, configuration, and execution | [Running Holi](wiki/Running-Holi.md) |
| Alternative Nextflow runner and validation limits | [Nextflow](../nextflow/README.md) and [branch review](NEXTFLOW_REVIEW.md) |
| Constructing reference sets | [Reference catalogue](wiki/Reference-catalogue.md) |
| Trimming, merging, and filtering | [Preprocessing](wiki/Preprocessing.md) |
| Microbial assignment and read splitting | [Microbial prefiltering](wiki/Microbial-prefiltering.md) |
| Competitive mapping and classification | [Taxonomic classification](wiki/Taxonomic-classification.md) |
| Interpretation of ancient DNA evidence | [Authentication](wiki/Authentication.md) |
| QC and merged taxonomic tables | [Postprocessing](wiki/Postprocessing.md) |
| Literature data and original plots | [Literature and figures](wiki/Literature-and-figures.md) |
| Dashboard development and wiki publication | [Development](DEVELOPMENT.md) |

## Version and provenance

These pages were reorganised against the scripts on `holi_dev` on 27 October 2026. They describe the checked-in implementation; they do not certify a successful full analysis run. The original public wiki was imported at commit `1a476eaa2bbb75f5dc391cc74d8a8f159bd8e69e` and is preserved in [wiki-archive](wiki-archive/README.md).

GitHub's public wiki is a separate Git repository. The files in [wiki](wiki/Home.md) are the proposed replacement pages and can be reviewed with the code before publishing. Changes here do not automatically update the live wiki.
