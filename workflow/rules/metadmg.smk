# ============================================================
# Stage 7 — metaDMG (LCA, dfit, aggregate)
# ============================================================

rule sort_for_metadmg:
    input:
        f"{EUK_OUT}/{{sample}}.alnfilt.bam",
    output:
        f"{EUK_OUT}/{{sample}}.sort.comp.filtered.bam",
    log:
        f"{LOGS}/{{sample}}__sort_comp_filtered.log",
    threads: config["threads_per_job"]
    params:
        samtools=config["tools"]["samtools"],
    shell:
        "{params.samtools} sort -n -@ {threads} -m 10G -o {output} {input} > {log} 2>&1"


rule metadmg_lca_euk:
    input:
        f"{EUK_OUT}/{{sample}}.sort.comp.filtered.bam",
    output:
        acc2tax=f"{EUK_OUT}/{{sample}}.acc2tax",
        bdamage=f"{EUK_OUT}/{{sample}}.sort.comp.filtered.bdamage.gz",
        stat=f"{EUK_OUT}/{{sample}}.sort.comp.filtered.stat.gz",
    log:
        f"{LOGS}/{{sample}}__metadmg_lca.log",
    threads: 10
    params:
        metadmg=config["tools"]["metadmg"],
        names=f"{config['tax_path_ncbi']}/taxdump/names.dmp",
        nodes=f"{config['tax_path_ncbi']}/taxdump/nodes.dmp",
        acc2tax_glob=lambda wc: acc2tax_glob(),
        sim_low=config["metadmg_euk"]["sim_score_low"],
        sim_high=config["metadmg_euk"]["sim_score_high"],
        how_many=config["metadmg_euk"]["how_many"],
        weight_type=config["metadmg_euk"]["weight_type"],
        prefix=f"{EUK_OUT}/{{sample}}.sort.comp.filtered",
    shell:
        """
        {params.metadmg} lca \
          --names {params.names} \
          --nodes {params.nodes} \
          --acc2tax <(zcat {params.acc2tax_glob}) \
          --sim_score_low {params.sim_low} \
          --sim_score_high {params.sim_high} \
          --how_many {params.how_many} \
          --weight_type {params.weight_type} \
          --fix_ncbi 0 --threads {threads} \
          --filtered_acc2tax {output.acc2tax} \
          --bam {input} \
          --out_prefix {params.prefix} \
          > {log} 2>&1
        """


rule metadmg_dfit:
    input:
        f"{EUK_OUT}/{{sample}}.sort.comp.filtered.bdamage.gz",
    output:
        f"{EUK_OUT}/{{sample}}.sort.comp.filtered.dfit.gz",
    log:
        f"{LOGS}/{{sample}}__metadmg_dfit.log",
    threads: 6
    params:
        metadmg=config["tools"]["metadmg"],
        names=f"{config['tax_path_ncbi']}/taxdump/names.dmp",
        nodes=f"{config['tax_path_ncbi']}/taxdump/nodes.dmp",
        nopt=config["metadmg_dfit"]["nopt"],
        nbootstrap=config["metadmg_dfit"]["nbootstrap"],
        seed=config["metadmg_dfit"]["seed"],
        lib=config["metadmg_dfit"]["lib"],
        prefix=f"{EUK_OUT}/{{sample}}.sort.comp.filtered",
    shell:
        """
        {params.metadmg} dfit {input} --threads {threads} \
          --names {params.names} --nodes {params.nodes} \
          --showfits 2 --nopt {params.nopt} \
          --nbootstrap {params.nbootstrap} \
          --seed {params.seed} --lib {params.lib} \
          --out_prefix {params.prefix} \
          > {log} 2>&1
        """


rule metadmg_aggregate:
    input:
        bdamage=f"{EUK_OUT}/{{sample}}.sort.comp.filtered.bdamage.gz",
        stat=f"{EUK_OUT}/{{sample}}.sort.comp.filtered.stat.gz",
        dfit=f"{EUK_OUT}/{{sample}}.sort.comp.filtered.dfit.gz",
    output:
        f"{EUK_OUT}/{{sample}}.sort.comp.filtered.agg.stat.gz",
    log:
        f"{LOGS}/{{sample}}__metadmg_agg.log",
    params:
        metadmg=config["tools"]["metadmg"],
        names=f"{config['tax_path_ncbi']}/taxdump/names.dmp",
        nodes=f"{config['tax_path_ncbi']}/taxdump/nodes.dmp",
        prefix=f"{EUK_OUT}/{{sample}}.sort.comp.filtered.agg",
    shell:
        """
        {params.metadmg} aggregate {input.bdamage} \
          --names {params.names} --nodes {params.nodes} \
          --lcastat {input.stat} --dfit {input.dfit} \
          --out_prefix {params.prefix} \
          > {log} 2>&1
        """
