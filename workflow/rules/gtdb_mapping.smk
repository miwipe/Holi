# ============================================================
# Stage 2 — GTDB mapping
# bowtie2 against each of the N GTDB chunks -> sort each -> merge
# ============================================================

wildcard_constraints:
    chunk = r"\d+",


rule bowtie2_gtdb_chunk:
    input:
        fq=f"{PREPROC_OUT}/{{sample}}.ppm.vs.d4.fq.gz",
    output:
        bam=temp(f"{MICROB_OUT}/{{sample}}.gtdb.{{chunk}}.bam"),
    log:
        f"{LOGS}/{{sample}}.gtdb.{{chunk}}.bowtie2.log",
    params:
        db=lambda wc: f"{config['db_path_bac']}.{wc.chunk}.fas.gz",
        k=config["bt2_gtdb"]["k"],
        D=config["bt2_gtdb"]["D"],
        R=config["bt2_gtdb"]["R"],
        N=config["bt2_gtdb"]["N"],
        L=config["bt2_gtdb"]["L"],
        i=config["bt2_gtdb"]["i"],
        np=config["bt2_gtdb"]["np"],
        mp=config["bt2_gtdb"]["mp"],
        rdg=config["bt2_gtdb"]["rdg"],
        rfg=config["bt2_gtdb"]["rfg"],
        score_min=config["bt2_gtdb"]["score_min"],
        bt2=config["tools"]["bowtie2"],
        samtools=config["tools"]["samtools"],
    threads: config["threads_per_job"]
    shell:
        """
        {params.bt2} --threads {threads} \
          -x {params.db} \
          -U {input.fq} \
          -k {params.k} -D {params.D} -R {params.R} -N {params.N} -L {params.L} \
          -i {params.i} --np {params.np} --mp '{params.mp}' \
          --rdg '{params.rdg}' --rfg '{params.rfg}' \
          --score-min '{params.score_min}' --mm --no-unal \
          2> {log} \
          | {params.samtools} view -bS - > {output.bam}
        """


rule sort_gtdb_chunk:
    input:
        f"{MICROB_OUT}/{{sample}}.gtdb.{{chunk}}.bam",
    output:
        temp(f"{MICROB_OUT}/{{sample}}.gtdb.{{chunk}}.sorted.bam"),
    threads: config["threads_per_job"]
    params:
        samtools=config["tools"]["samtools"],
    shell:
        "{params.samtools} sort -@ {threads} -m 4G -o {output} {input}"


rule merge_gtdb:
    input:
        lambda wc: expand(
            f"{MICROB_OUT}/{{sample}}.gtdb.{{chunk}}.sorted.bam",
            sample=wc.sample,
            chunk=range(1, config["n_gtdb_chunks"] + 1),
        ),
    output:
        f"{MICROB_OUT}/{{sample}}.gtdb.merged.bam",
    threads: config["threads_per_job"]
    params:
        samtools=config["tools"]["samtools"],
    shell:
        "{params.samtools} merge -@ {threads} -f {output} {input}"
