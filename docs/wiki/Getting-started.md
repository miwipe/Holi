# Getting started

The supplied `holi` environment targets Linux: it includes Linux compiler packages and is intended for a suitably resourced workstation or cluster. Resource needs depend on the input size, catalogue size, and number of simultaneous mapping jobs. The earlier wiki's example of a large-memory node is not a universal minimum.

## Install from the repository

```bash
git clone https://github.com/miwipe/Holi.git
cd Holi
conda env create -f envs/holi.yml
conda activate holi
```

The environment is named `holi`, not `holi.yml`. The environment definition already contains a pinned bam-filter source dependency; do not replace it with an unpinned installation as part of routine setup.

`bin/holi.yml` remains a compatibility link to `envs/holi.yml`. Linux compiler dependencies prevent treating this environment file as a validated native macOS installation recipe.

## Check external executables

`bin/holi.sh` checks these commands before starting:

```bash
for tool in parallel fastp vsearch sga bowtie2 samtools conda seqtk unicorn metaDMG-cpp; do
  command -v "$tool" || echo "Missing: $tool"
done
```

The default microbial split also invokes `getRTax`; check that command if using the default mode. Verify the actual executables and their versions after environment creation, including Unicorn and metaDMG-cpp. An environment solve alone does not verify these tools or a successful pipeline run.

Reference genomes, Bowtie2 indexes, and taxonomy files are not bundled. Continue with [reference catalogue construction](Reference-catalogue.md) and [input configuration](Running-Holi.md).

## Butteracid environment

```bash
conda env create -f envs/butteracid.yml
conda activate Renv
```

Read the environment's `name:` field if you customise it. [Postprocessing](Postprocessing.md) describes required configuration and working-directory assumptions.
