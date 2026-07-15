# ============================================================
# Stage 6 — BAM filtering (Unicorn: alnfilt -> refstats)
# ============================================================

rule unicorn_alnfilt:
    input:
        f"{EUK_OUT}/{{sample}}.comp.bam",
    output:
        bam=f"{EUK_OUT}/{{sample}}.alnfilt.bam",
        stat=f"{EUK_OUT}/{{sample}}.alnfilt.refstats",
    log:
        f"{LOGS}/{{sample}}__unicorn_alnfilt.log",
    threads: config["threads_per_job"]
    params:
        unicorn=config["tools"]["unicorn"],
        mode=config["alnfilt"]["mode"],
        minani=config["alnfilt"]["minani"],
        maxani=config["alnfilt"]["maxani"],
    shell:
        """
        {params.unicorn} alnfilt \
          -b {input} \
          -t {threads} --mode {params.mode} \
          --outbam {output.bam} --outstat {output.stat} \
          --minani {params.minani} --maxani {params.maxani} \
          > {log} 2>&1
        """


rule unicorn_refstats:
    input:
        f"{EUK_OUT}/{{sample}}.alnfilt.bam",
    output:
        bam=f"{EUK_OUT}/{{sample}}.alnfilt.unicorn.bam",
        stat=f"{EUK_OUT}/{{sample}}.comp.unicorn.refstats",
    log:
        f"{LOGS}/{{sample}}__unicorn_refstats.log",
    threads: config["threads_per_job"]
    params:
        unicorn=config["tools"]["unicorn"],
        minreads=config["refstats"]["minreads"],
        names=f"{config['tax_path_ncbi']}/taxdump/names.dmp",
        nodes=f"{config['tax_path_ncbi']}/taxdump/nodes.dmp",
    shell:
        """
        {params.unicorn} refstats \
          -b {input} \
          -t {threads} --minreads {params.minreads} \
          --outbam {output.bam} --outstat {output.stat} \
          --names {params.names} --nodes {params.nodes} \
          > {log} 2>&1
        """
