import glob
import os
import re

# ------------------------------------------------------------------
# Samples
# ------------------------------------------------------------------
with open(config["samples"]) as fh:
    SAMPLES = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

if not SAMPLES:
    raise WorkflowError(f"No samples found in {config['samples']}")

wildcard_constraints:
    sample = "|".join(re.escape(s) for s in SAMPLES),

# ------------------------------------------------------------------
# Shorthand paths (mirrors the variable names used in bin/holi.sh)
# ------------------------------------------------------------------
INPATH = config["inpath"]
RESULT_PATH = config["result_path"]
PREPROC_OUT = f"{RESULT_PATH}/preprocessing"
MICROB_OUT = config["microb_out"]
EUK_OUT = config["euk_out"]
LOGS = config["logs"]
TMP = config["tmp"]

SINGLE_END = bool(config.get("single_end", False))
LCA_ASSIGNMENT = bool(config.get("lca_assignment", False))

os.makedirs(LOGS, exist_ok=True)
os.makedirs(TMP, exist_ok=True)
os.environ["TMPDIR"] = TMP
os.environ["TEMP"] = TMP
os.environ["TMP"] = TMP


# ------------------------------------------------------------------
# FASTQ discovery helpers
#
# holi.sh locates raw FASTQs by globbing INPATH for the sample name.
# Resolving the glob here (at DAG-build time) instead of inside the
# shell command means Snakemake fails fast with a clear error if a
# sample's FASTQ is missing or ambiguous, rather than fastp silently
# receiving a literal, unexpanded glob pattern.
# ------------------------------------------------------------------

def _one_match(pattern, sample, label):
    matches = sorted(glob.glob(pattern))
    if len(matches) == 0:
        raise WorkflowError(f"No {label} FASTQ found for sample '{sample}' matching: {pattern}")
    if len(matches) > 1:
        raise WorkflowError(
            f"Multiple {label} FASTQs matched for sample '{sample}': {matches}. "
            "Refine the glob pattern or merge lanes before running."
        )
    return matches[0]


def get_se_fastq(sample):
    return _one_match(f"{INPATH}/*_{sample}_*R1*_001.fastq.gz", sample, "R1 (single-end)")


def get_r1(sample):
    return _one_match(f"{INPATH}/{sample}*R1*.fastq.gz", sample, "R1")


def get_r2(sample):
    return _one_match(f"{INPATH}/{sample}*R2*.fastq.gz", sample, "R2")


# ------------------------------------------------------------------
# Stage 4 mapping "part" names shared between euk_mapping.smk and
# comp_merge.smk (every bowtie2 target that stage 5 has to
# compress/sort/merge into the final per-sample comp.bam).
# ------------------------------------------------------------------

def euk_mapping_parts():
    parts = [f"euk.{i}" for i in range(1, config["n_euk_chunks"] + 1)]
    parts.append("mito")
    parts += [f"phyNor.{i}" for i in range(1, config["n_norway_chunks"] + 1)]
    parts.append("core_nt")
    parts.append("pla")
    return parts


def acc2tax_glob():
    """Space-separated glob of *.acc2taxid.gz files merged for eukaryotic
    metaDMG lca / unicorn taxstats (NCBI taxdump + phylo-Norway)."""
    return f"{config['tax_path_ncbi']}/*.acc2taxid.gz {config['tax_path_norway_acc']}/*.acc2taxid.gz"
