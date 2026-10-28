# Run Holi with Nextflow

This experimental DSL2 runner provides an alternative to `bin/holi.sh`. Both are maintained in `holi_dev`; the Bash entry point and existing plots are unchanged. The older `nextflow_porting` branch contains no Nextflow implementation or unique commits relative to the inspected development history. It is not a second maintained pipeline.

## Choose a runner

| | Bash | Nextflow |
| --- | --- | --- |
| Entry point | `bash bin/holi.sh CONFIG SAMPLE_LIST` | `nextflow run nextflow/main.nf ...` |
| Inputs | Sample identifiers and filename patterns | Explicit CSV; paired and single samples can coexist |
| References | Existing fixed catalogue/shard layout | Explicit index-prefix lists under one reference directory |
| Scheduling | GNU Parallel inside the script | One task per sample per stage; reference shards run sequentially inside each mapping task |
| Resume | Existing files and stage flags | `-resume` with retained Nextflow cache and work directory |
| Readiness | Existing implementation, with site-specific paths | Execution tested; real-data scientific equivalence is not yet established |

Use Bash for continuity with existing runs. Use Nextflow for an explicitly configured pilot, task isolation, and cached restarts. Nextflow is not yet a validated replacement for production results.

## Requirements

The runner was tested with Nextflow **25.10.2**, Java 17, Python 3.9, and the v2 syntax parser. Nextflow 25.10.2 or newer is required. The biological environment remains the Linux environment in [envs/holi.yml](../envs/holi.yml); a successful orchestration test on macOS does not validate that environment on macOS.

Activate the analysis environment before launching. `fastp`, `vsearch`, `sga`, `bowtie2`, `samtools`, `seqtk`, `getRTax` (default assignment), `unicorn`, `metaDMG-cpp`, and metaDMG's separate `compressbam` executable must be available. Python's `metadmg` package does not substitute for the C++ executable. `tools` in the JSON can override executable names with absolute paths, without extra arguments. Validation checks executables and records their paths and SHA-256 hashes. It does not prove API compatibility; record and retain your tested tool builds, particularly Unicorn and metaDMG-cpp.

No tools, Conda environments, reference databases, or containers are downloaded automatically. There are no automatic environment activations inside tasks.

## Prepare manifests

```bash
cp config/nextflow.example.json config/nextflow.local.json
cp config/nextflow.samples.example.csv config/nextflow.samples.local.csv
```

The CSV has exactly `sample,fastq_1,fastq_2` columns. Use absolute FASTQ paths and quote CSV fields containing commas. Leave `fastq_2` empty for single-end input. Each row represents one sample; concatenate lanes deliberately before preparing the sheet. Sample identifiers must be unique and contain only letters, digits, underscores, dots, or hyphens, starting with a letter or digit. Paired FASTQs must have distinct basenames. Spaces and apostrophes in FASTQ paths are supported.

Set `--references` to a directory containing the indexed references and taxonomy. JSON reference paths are relative to this directory, with no `..` components. Every Bowtie2 prefix needs all six `.bt2` or all six `.bt2l` components. List **every** microbial and eukaryotic shard to be included: the short example is a layout illustration, not the complete Holi catalogue. Include nuclear, mitochondrial, plastid, core NT and phylo-Norway references as appropriate to the existing analysis. Missing entries alter the competing reference set and therefore the scientific result.

`taxonomy.acc2tax` lists the accession tables for the complete eukaryotic reference set, including phylo-Norway when used. Their uncompressed concatenation is created task-locally for metaDMG and Unicorn. Allocate scratch space for this additional table. Taxonomy identifiers, accession versions, and reference sequence names must agree with the indexes; existence checks cannot establish that agreement. References should be immutable, versioned builds. Staging uses symlinks, not database copies. A bundle can contain symlinks to existing files rather than moving large catalogues.

Choose `microbial_assignment` as `getrtax` (default, using `microbial_taxid_map`) or `lca` (using `microbial_names`, `microbial_nodes`, and `microbial_acc2tax`). Values in `parameters` override the defaults in [runner/defaults.json](runner/defaults.json), which were copied from [config/holi.example.yml](../config/holi.example.yml). Unknown parameter names and invalid ANI/similarity ranges are rejected.

## Execute and resume

From the repository root:

```bash
export NXF_SYNTAX_PARSER=v2
nextflow run nextflow/main.nf \
  --input config/nextflow.samples.local.csv \
  --settings config/nextflow.local.json \
  --references /absolute/path/to/reference_bundle \
  --outdir /absolute/path/to/run01/results \
  -work-dir /absolute/path/to/run01/work \
  --cpus 8 --memory '64 GB'
```

Repeat the same command with **`-resume`** to reuse successful tasks. Retain both the work directory and the launching directory's `.nextflow/` state. A new work directory alone does not replace the launching directory's cache. Do not delete work files while planning a restart. Nextflow does not use Bash's `--from-stage`, `--to-stage`, or skip flags. There is no automatic intermediate cleanup. Published stage outputs are copied in addition to retained work files, so budget disk space for both.

`--cpus` and `--memory` describe each task. `--max_forks` defaults to 1 per stage, not 1 for the entire workflow: different stages and the two final branches can run concurrently. Set local executor limits in a separate config to fit your machine, for example:

```groovy
executor {
    cpus = 16
    memory = '96 GB'
}
```

Pass that config with `-c /path/to/site.config`. `-profile slurm` selects Slurm, but you must supply your institution's queue/account/resource settings and a shared filesystem. This profile has not been tested on a cluster. There are no cloud profiles.

Output directories are `stage1/<sample>.stage1/` through `stage8/<sample>.stage8/`, plus `validation/` and timestamped `execution/` reports, timelines, and TSV traces. Each stage records resolved scientific settings and commands in `provenance.json`, with tool output in `commands.log`. Task `.command.*` files remain in `work/`. Final stages run independently from the stage-6 output. Files in a chosen output directory can be replaced when tasks rerun; use a separate output directory for each scientific experiment, and inspect exit status and trace rather than trusting old published files.

Nextflow hashes task inputs, settings, and runner source files for caching. Reference directories and executable names do not provide a complete content identity for a database build or software environment. **Do not edit reference bundles or replace tool builds in place and rely on `-resume` to detect the change.** Use a new versioned reference directory and a new run/cache location when either changes. The validation record's executable hashes describe that validation task; a resumed validation task will retain its earlier record.

## Scientific compatibility and deliberate differences

The eight stages use the same main commands and scientific parameter defaults as `bin/holi.sh`. Stages 7 and 8 intentionally consume `<sample>.alnfilt.bam`, matching Bash; stage 6 also retains the refstats-filtered BAM. Despite the Bash stage heading, the current script does not invoke `unicorn bamstats`, so this runner does not invent a bamstats step. `--filtered_acc2tax` is an output from metaDMG LCA, not an input to stage 6.

Changes in execution include explicit FASTQ paths, configurable shard counts, task-scoped files, recorded fastp JSON reports, bounded sort-memory allocations, and threads taken from the task allocation rather than hard-coded 6/8/10/12-thread values. Failed producers in a mapping pipe abort the task. Unicorn statistics are accepted as plain files or with the `.gz` suffix used by newer builds. Produced BAMs are checked with `samtools quickcheck`; required final outputs must exist and be nonempty. Those checks detect some failures but do not establish valid read counts or taxonomy.

Read IDs are consistently reduced to their first whitespace-delimited token, sorted uniquely, and partitioned with external `sort`/`comm`, avoiding whole-sample Python sets. This fixes a difference between the original LCA and default splitting branches. getRTax read-ID output may be plain or gzip; missing output files cause failure rather than silently treating classification failure as zero assignments. Compressed outputs and a first-token read partition preserve valid empty microbial/eukaryotic read sets; downstream tools can still reject empty biological inputs, in which case the runner stops rather than fabricating results.

The temporary merge in stage 5 is explicitly BAM rather than the misleading `.sam.gz` suffix used by Bash. No literature PDFs are added to this runner, the repository, or the wiki.

## Validation

```bash
python3 -m unittest discover -s nextflow/tests -p 'test_runner.py' -v
python3 nextflow/tests/smoke.py --nextflow nextflow
```

The unit suite also runs a tiny real Bowtie2/samtools mapping check when those tools are installed (otherwise it reports a skip). The smoke test uses a real Nextflow engine, tiny paired/single FASTQs, fake reference files, and clearly labelled **protocol mocks** for biological executables. It exercises all eight stages, both assignment modes, read partitioning, output contracts, and cache reuse. Stub results contain `STUB_ONLY.txt`; neither stubs nor mocked BAMs are scientific outputs. Test files are created in a temporary directory and removed on success; `--keep /absolute/new/test-directory` retains logs.

Before adopting this runner for production, run the same small real dataset, complete reference build, taxonomy snapshot and tool builds through both runners. Compare preprocessing read counts, microbial/retained read IDs, per-reference mapping counts, alignment filtering, taxonomic assignments and damage estimates, and investigate differences in thread-dependent or numerical output. That real-data comparison has not yet been performed.

See [Nextflow's 25.10 migration notes](https://docs.seqera.io/nextflow/migrations/25-10) and [metaDMG-cpp's command documentation](https://github.com/metaDMG-dev/metaDMG-cpp) for upstream interfaces.

Nextflow outputs are separated by stage and sample. Existing Butteracid wrappers expect the Bash output layout; direct Butteracid invocation against the new directory tree is not supported by this initial port. Gather/adapt the required final tables deliberately before postprocessing.
