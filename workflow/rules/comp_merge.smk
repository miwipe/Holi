# ============================================================
# Stage 5 — Compress, sort, merge BAMs
#
# Each Stage 4 mapping BAM is compressed with metaDMG's compressbam and
# name-sorted, then all of them are merged into one name-sorted SAM.gz
# per sample, which is itself compressed once more into the final
# per-sample comp.bam used by every later stage.
# ============================================================

wildcard_constraints:
    part = r"(euk\.\d+|mito|phyNor\.\d+|core_nt|pla)",


rule compress_and_sort_bam:
    input:
        f"{EUK_OUT}/{{sample}}.{{part}}.bam",
    output:
        comp=temp(f"{EUK_OUT}/{{sample}}.{{part}}.comp.bam"),
        sorted_=temp(f"{EUK_OUT}/{{sample}}.{{part}}.comp.sorted.bam"),
    log:
        compress=f"{LOGS}/{{sample}}.{{part}}__compressbam.log",
        sort=f"{LOGS}/{{sample}}.{{part}}__sort_comp.log",
    threads: config["threads_per_job"]
    params:
        compressbam=config["tools"]["compressbam"],
        samtools=config["tools"]["samtools"],
    shell:
        """
        {params.compressbam} --threads {threads} --input {input} --output {output.comp} > {log.compress} 2>&1
        {params.samtools} sort -n -@ {threads} -m 4G -o {output.sorted_} {output.comp} > {log.sort} 2>&1
        """


rule merge_comp_sam:
    input:
        lambda wc: expand(
            f"{EUK_OUT}/{{sample}}.{{part}}.comp.sorted.bam",
            sample=wc.sample,
            part=euk_mapping_parts(),
        ),
    output:
        temp(f"{EUK_OUT}/{{sample}}.comp.sam.gz"),
    threads: config["threads_per_job"]
    params:
        samtools=config["tools"]["samtools"],
    shell:
        "{params.samtools} merge -@ {threads} -n -f {output} {input}"


rule compress_merged_comp:
    input:
        f"{EUK_OUT}/{{sample}}.comp.sam.gz",
    output:
        f"{EUK_OUT}/{{sample}}.comp.bam",
    log:
        f"{LOGS}/{{sample}}__compress_merged.log",
    threads: config["threads_per_job"]
    params:
        compressbam=config["tools"]["compressbam"],
    shell:
        "{params.compressbam} --threads {threads} --input {input} --output {output} > {log} 2>&1"
