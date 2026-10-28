# Nextflow branch review and integration

Reviewed 28 October 2026 against `holi_dev` at `bb88afb2663c744a3ed2b981e12480d7abcc32a6` and `origin/nextflow_porting` at `ce4713cb47f4f04f471a4228f3351819c80394de`.

## What was present

The inspected `nextflow_porting` tip is an ancestor of `holi_dev`, with **zero unique commits** relative to it. Its tracked files and available Git history contain no `.nf` workflow or `nextflow.config`. The biological scripts match the development branch. There was therefore no existing Nextflow implementation to grade or merge; the branch name describes an intention rather than an implemented port.

Recommendation: maintain **both execution options together in `holi_dev`**, with shared scientific defaults and documented differences. A permanently separate Nextflow branch would duplicate maintenance and make drift harder to detect. The old branch is left untouched. The new port was prepared on `codex/nextflow-dual-runner` before integration.

## Findings affecting a port

| Finding in current Bash implementation | Consequence | Treatment in the alternative runner |
| --- | --- | --- |
| Absolute institution-specific Unicorn/metaDMG paths | A portable config alone cannot run the full pipeline elsewhere | Executables resolved from the active environment or explicit JSON paths |
| Fixed microbial/eukaryotic shard counts and catalogue naming | Substituting a directory does not define a new competing reference set | Explicit complete index-prefix lists, checked for six index components |
| Mapper output piped into samtools without a pipeline-wide producer check | A mapper failure can be obscured by downstream command status | Both subprocess exit codes checked; produced BAMs checked with samtools |
| File-presence restart heuristics | Partial/stale outputs can be mistaken for completed stages | Nextflow task isolation and cached task completion; published files alone are not success evidence |
| Literal thread counts for several tools, independent of general config | Resource requests may not match actual concurrency | Task CPU allocation propagated; samtools sort memory divided per thread |
| LCA and default split branches normalise read headers differently | Whitespace-bearing read names can partition differently | Consistent first-token read IDs, external sorting and set difference |
| Stage 6 heading mentions bamstats, but no bamstats command is executed | Documentation can overstate produced analyses | No invented bamstats stage; documentation describes actual commands |
| Final metaDMG and taxstats inputs are alnfilt BAMs | Accidentally selecting the refstats BAM would change scientific behaviour | Preserve the actual Bash inputs explicitly |

The original Bash implementation is retained unchanged to preserve existing operation. This review is not a claim that existing results are invalid; the listed issues are portability, failure-detection, reproducibility, and interpretation concerns requiring attention when comparing runs.

## Scope and quality of the new runner

The alternative contains eight separately cached DSL2 processes plus manifest validation. It supports paired and single samples in one CSV, both getRTax and microbial LCA assignment, local execution, an untested opt-in Slurm profile, explicit scientific defaults, command/resource provenance, executable hashes, trace/report/timeline outputs, and resumable work directories. Mapping shards currently run sequentially within each sample's mapping stage; there is no per-shard scheduling or retry granularity.

Validation covers malformed manifests, missing index components, reference path escapes, unknown parameters and invalid ranges, safe sample IDs, mapper failure propagation, and shared scientific-default consistency. An optional tiny real Bowtie2/samtools test checks that the mapping stage produces a readable BAM. Real Nextflow stub and mocked-tool runs cover both sample types, both assignment modes, read partitioning, all eight stages, output contracts including compressed Unicorn statistics, and caching on resume. Mocks and stubs are clearly identified and are not scientific results.

Local validation on 28 October 2026 passed all nine unit/integration checks, including the real mapping test. With Nextflow 25.10.2 and the v2 parser, the stub, getRTax mock, and LCA mock scenarios each completed 17 tasks (one validation plus eight stages for two samples); each subsequent resume reused all 17 tasks. Documentation links, wiki export, dashboard validation and original-figure checks also passed. BAM read-count checking in the real mapping test establishes a nonempty readable mapping result, not end-to-end biological equivalence.

The runner is **experimental**. No complete real-data Bash/Nextflow comparison has been performed, no cluster profile has been validated, and custom Unicorn/metaDMG versions are not installed or pinned automatically. References must be immutable/versioned; changing an installed executable or editing files inside a reference bundle is not a reliable cache invalidation strategy. A resumed validation task retains its earlier executable hashes. Butteracid's existing wrappers expect the Bash output layout and are not directly adapted to the stage/sample directory tree in this initial port.

Before production adoption, compare both runners using the same real reads, complete indexes, taxonomy snapshot, parameters and executable builds. Check intermediate read IDs and counts, per-reference alignment counts, filtered BAM content, taxonomic tables, and damage estimates. Record any permitted numerical/thread-dependent differences rather than declaring equivalence from exit status.

See the [Nextflow guide](../nextflow/README.md) for commands, manifests, output layout and test instructions. Existing literature plots, dashboard data, Bash scripts and PDF storage remain unchanged by this integration.
