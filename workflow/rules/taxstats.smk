# ============================================================
# Stage 8 — Unicorn taxstats (family, genus, species)
# ============================================================

wildcard_constraints:
    rank = "family|genus|species",


rule unicorn_taxstats:
    input:
        f"{EUK_OUT}/{{sample}}.alnfilt.bam",
    output:
        bam=f"{EUK_OUT}/{{sample}}.comp.filtered.{{rank}}.taxstats.bam",
        stat=f"{EUK_OUT}/{{sample}}.comp.filtered.{{rank}}.taxstats",
    log:
        f"{LOGS}/{{sample}}__unicorn_taxstats_{{rank}}.log",
    threads: config["threads_per_job"]
    params:
        unicorn=config["tools"]["unicorn"],
        names=f"{config['tax_path_ncbi']}/taxdump/names.dmp",
        nodes=f"{config['tax_path_ncbi']}/taxdump/nodes.dmp",
        acc2tax_glob=lambda wc: acc2tax_glob(),
        k=config["taxstats"]["k"],
        minreads=config["taxstats"]["minreads"],
    shell:
        """
        {params.unicorn} taxstats \
          -b {input} \
          -t {threads} \
          -o {output.bam} \
          --names {params.names} --nodes {params.nodes} \
          --acc2tax <(zcat {params.acc2tax_glob}) \
          -k {params.k} --outstat {output.stat} \
          --minreads {params.minreads} --rank {wildcards.rank} \
          > {log} 2>&1
        """
