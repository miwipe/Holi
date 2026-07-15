# ============================================================
# Stage 4 — Eukaryotic mapping
# bowtie2 against: N eukaryote WGS chunks, mitochondrion, N phylo-Norway
# shards, core_nt, and plastid. Each produces one BAM per sample that
# Stage 5 will compress/sort/merge into a single comp.bam.
# ============================================================

wildcard_constraints:
    euk_chunk = r"\d+",
    norway_chunk = r"\d+",


rule bowtie2_euk_chunk:
    input:
        fq=f"{EUK_OUT}/{{sample}}.euk.fastq.gz",
    output:
        temp(f"{EUK_OUT}/{{sample}}.euk.{{euk_chunk}}.bam"),
    log:
        f"{LOGS}/{{sample}}__eukmap_part_{{euk_chunk}}.log",
    params:
        db=lambda wc: f"{config['db_path']}.{wc.euk_chunk}.fas.gz",
        k=config["bt2_euk"]["k"],
        bt2=config["tools"]["bowtie2"],
        samtools=config["tools"]["samtools"],
    threads: config["threads_per_job"]
    shell:
        """
        {params.bt2} --threads {threads} -k {params.k} -t \
          -x {params.db} \
          -U {input.fq} \
          --no-unal --mm -t 2> {log} \
          | {params.samtools} view -bS - > {output}
        """


rule bowtie2_mito:
    input:
        fq=f"{EUK_OUT}/{{sample}}.euk.fastq.gz",
    output:
        temp(f"{EUK_OUT}/{{sample}}.mito.bam"),
    log:
        f"{LOGS}/{{sample}}__mitochondrion.log",
    params:
        db=f"{config['db_path_clean']}/refseq_mitochondrion.genomic.fas.gz",
        k=config["bt2_euk"]["k"],
        bt2=config["tools"]["bowtie2"],
        samtools=config["tools"]["samtools"],
    threads: config["threads_per_job"]
    shell:
        """
        {params.bt2} --threads {threads} -k {params.k} -t \
          -x {params.db} \
          -U {input.fq} \
          --no-unal --mm -t 2> {log} \
          | {params.samtools} view -bS - > {output}
        """


rule bowtie2_phynor_chunk:
    input:
        fq=f"{EUK_OUT}/{{sample}}.euk.fastq.gz",
    output:
        temp(f"{EUK_OUT}/{{sample}}.phyNor.{{norway_chunk}}.bam"),
    log:
        f"{LOGS}/{{sample}}__phynor_part_{{norway_chunk}}.log",
    params:
        db=lambda wc: f"{config['db_path_norway']}.{wc.norway_chunk}-of-{config['n_norway_chunks']}",
        k=config["bt2_euk"]["k"],
        bt2=config["tools"]["bowtie2"],
        samtools=config["tools"]["samtools"],
    threads: config["threads_per_job"]
    shell:
        """
        {params.bt2} --threads {threads} -k {params.k} -t \
          -x {params.db} \
          -U {input.fq} \
          --no-unal --mm -t 2> {log} \
          | {params.samtools} view -bS - > {output}
        """


rule bowtie2_core_nt:
    input:
        fq=f"{EUK_OUT}/{{sample}}.euk.fastq.gz",
    output:
        temp(f"{EUK_OUT}/{{sample}}.core_nt.bam"),
    log:
        f"{LOGS}/{{sample}}__core_nt.log",
    params:
        db=f"{config['db_path_clean']}/core_nt.fas.gz",
        k=config["bt2_euk"]["k"],
        bt2=config["tools"]["bowtie2"],
        samtools=config["tools"]["samtools"],
    threads: config["threads_per_job"]
    shell:
        """
        {params.bt2} --threads {threads} -k {params.k} -t \
          -x {params.db} \
          -U {input.fq} \
          --no-unal --mm -t 2> {log} \
          | {params.samtools} view -bS - > {output}
        """


rule bowtie2_plastid:
    input:
        fq=f"{EUK_OUT}/{{sample}}.euk.fastq.gz",
    output:
        temp(f"{EUK_OUT}/{{sample}}.pla.bam"),
    log:
        f"{LOGS}/{{sample}}__plastid.log",
    params:
        db=f"{config['db_path_clean']}/refseq_plastid.genomic.fas.gz",
        k=config["bt2_euk"]["k"],
        bt2=config["tools"]["bowtie2"],
        samtools=config["tools"]["samtools"],
    threads: config["threads_per_job"]
    shell:
        """
        {params.bt2} --threads {threads} -k {params.k} -t \
          -x {params.db} \
          -U {input.fq} \
          --no-unal --mm -t 2> {log} \
          | {params.samtools} view -bS - > {output}
        """
