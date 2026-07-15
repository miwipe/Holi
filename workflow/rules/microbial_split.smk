# ============================================================
# Stage 3 — Microbial / eukaryotic read splitting
#
# Two mutually-exclusive branches, chosen by config["lca_assignment"]:
#   - getRTax (default)
#   - metaDMG LCA (--lca-assignment in the old bash script)
# Both converge on the same output: {EUK_OUT}/{sample}.euk.fastq.gz,
# so Stage 4 doesn't need to know which branch produced it.
# ============================================================

if not LCA_ASSIGNMENT:

    rule getrtax_classify:
        # getRTax writes one or more "<prefix>.txt*" files whose exact
        # naming is tool-version-dependent, so the classification call and
        # the zcat merge are kept in a single rule rather than declared as
        # separate Snakemake outputs.
        input:
            bam=f"{MICROB_OUT}/{{sample}}.gtdb.merged.bam",
        output:
            f"{MICROB_OUT}/{{sample}}.bact_reads_all.txt",
        log:
            f"{LOGS}/{{sample}}__getRTax.log",
        params:
            taxmap=f"{config['tax_path_bac']}/taxid.map",
            bin=config["tools"]["getRTax"],
            prefix=f"{MICROB_OUT}/{{sample}}.bact_reads.txt",
        threads: 8
        shell:
            """
            {params.bin} \
              --bam {input.bam} \
              -T {params.taxmap} \
              -r '{{"domain":["d__Bacteria", "d__Archaea", "d__Viruses"]}}' \
              --threads {threads} \
              --unique \
              --only-read-ids \
              -p {params.prefix} > {log} 2>&1
            zcat {params.prefix}* > {output}
            """

    rule subset_bact_fastq:
        input:
            fq=f"{PREPROC_OUT}/{{sample}}.ppm.vs.d4.fq.gz",
            ids=f"{MICROB_OUT}/{{sample}}.bact_reads_all.txt",
        output:
            f"{MICROB_OUT}/{{sample}}.bact_reads.fq.gz",
        log:
            f"{LOGS}/{{sample}}__seqtk_bact.log",
        params:
            bin=config["tools"]["seqtk"],
        shell:
            "{params.bin} subseq {input.fq} {input.ids} 2> {log} | gzip > {output}"

    rule all_read_ids:
        input:
            f"{PREPROC_OUT}/{{sample}}.ppm.vs.d4.fq.gz",
        output:
            temp(f"{MICROB_OUT}/{{sample}}.all_reads.txt"),
        log:
            f"{LOGS}/{{sample}}__allreads.log",
        shell:
            """
            zcat {input} | awk 'NR%4==1 {{split(substr($0, 2), a, " "); print a[1]}}' > {output} 2> {log}
            """

    rule sort_read_lists:
        input:
            all=f"{MICROB_OUT}/{{sample}}.all_reads.txt",
            bact=f"{MICROB_OUT}/{{sample}}.bact_reads_all.txt",
        output:
            all_sorted=temp(f"{MICROB_OUT}/{{sample}}.all_reads.sorted.txt"),
            bact_sorted=temp(f"{MICROB_OUT}/{{sample}}.bact_reads_all.sorted.txt"),
        shell:
            """
            sort {input.all} > {output.all_sorted}
            sort {input.bact} > {output.bact_sorted}
            """

    rule euk_read_list:
        input:
            all_sorted=f"{MICROB_OUT}/{{sample}}.all_reads.sorted.txt",
            bact_sorted=f"{MICROB_OUT}/{{sample}}.bact_reads_all.sorted.txt",
        output:
            temp(f"{EUK_OUT}/{{sample}}.euk_reads.txt"),
        log:
            f"{LOGS}/{{sample}}__comm_euk.log",
        shell:
            "comm -23 {input.all_sorted} {input.bact_sorted} > {output} 2> {log}"

    rule extract_euk_fastq:
        input:
            fq=f"{PREPROC_OUT}/{{sample}}.ppm.vs.d4.fq.gz",
            ids=f"{EUK_OUT}/{{sample}}.euk_reads.txt",
        output:
            f"{EUK_OUT}/{{sample}}.euk.fastq.gz",
        log:
            f"{LOGS}/{{sample}}__seqtk_euk.log",
        params:
            bin=config["tools"]["seqtk"],
        shell:
            "{params.bin} subseq {input.fq} {input.ids} 2> {log} | gzip > {output}"

else:

    rule sort_gtdb_for_lca:
        input:
            f"{MICROB_OUT}/{{sample}}.gtdb.merged.bam",
        output:
            temp(f"{MICROB_OUT}/{{sample}}.gtdb.merged.sorted.bam"),
        log:
            f"{LOGS}/{{sample}}__sortbam.log",
        threads: config["threads_per_job"]
        params:
            samtools=config["tools"]["samtools"],
        shell:
            "{params.samtools} sort -n -@ {threads} -m 10G -o {output} {input} > {log} 2>&1"

    rule metadmg_lca_microbial:
        input:
            f"{MICROB_OUT}/{{sample}}.gtdb.merged.sorted.bam",
        output:
            f"{MICROB_OUT}/{{sample}}.lca.gz",
        log:
            f"{LOGS}/{{sample}}__metadmg.log",
        params:
            metadmg=config["tools"]["metadmg"],
            names=f"{config['tax_path_bac']}/names.dmp",
            nodes=f"{config['tax_path_bac']}/nodes.dmp",
            acc2tax=f"{config['tax_path_bac_acc']}/hires-organelles-viruses-smags.acc2taxid.gz",
            sim_low=config["metadmg_bac"]["sim_score_low"],
            sim_high=config["metadmg_bac"]["sim_score_high"],
            how_many=config["metadmg_bac"]["how_many"],
            weight_type=config["metadmg_bac"]["weight_type"],
            prefix=f"{MICROB_OUT}/{{sample}}",
        threads: 10
        shell:
            """
            {params.metadmg} lca \
              --names {params.names} \
              --nodes {params.nodes} \
              --acc2tax {params.acc2tax} \
              --sim_score_low {params.sim_low} \
              --sim_score_high {params.sim_high} \
              --how_many {params.how_many} \
              --weight_type {params.weight_type} \
              --fix_ncbi 0 \
              --threads {threads} \
              --bam {input} \
              --out_prefix {params.prefix} > {log} 2>&1
            """

    rule extract_bact_reads_lca:
        input:
            f"{MICROB_OUT}/{{sample}}.lca.gz",
        output:
            temp(f"{MICROB_OUT}/{{sample}}.bact_reads.txt"),
        log:
            f"{LOGS}/{{sample}}__extract_bact.log",
        shell:
            "zgrep -i -E 'Archaea|virus|bacteria' {input} | cut -f1 > {output} 2> {log}"

    rule subset_bact_fastq_lca:
        input:
            fq=f"{PREPROC_OUT}/{{sample}}.ppm.vs.d4.fq.gz",
            ids=f"{MICROB_OUT}/{{sample}}.bact_reads.txt",
        output:
            f"{MICROB_OUT}/{{sample}}.bact_reads.fq.gz",
        log:
            f"{LOGS}/{{sample}}__seqtk_bact.log",
        params:
            bin=config["tools"]["seqtk"],
        shell:
            "{params.bin} subseq {input.fq} {input.ids} 2> {log} | gzip > {output}"

    rule all_read_ids_lca:
        input:
            f"{PREPROC_OUT}/{{sample}}.ppm.vs.d4.fq.gz",
        output:
            temp(f"{MICROB_OUT}/{{sample}}.all_reads.txt"),
        log:
            f"{LOGS}/{{sample}}__allreads.log",
        shell:
            "zcat {input} | awk 'NR%4==1 {{print substr($0, 2)}}' > {output} 2> {log}"

    rule euk_read_list_lca:
        input:
            all=f"{MICROB_OUT}/{{sample}}.all_reads.txt",
            bact=f"{MICROB_OUT}/{{sample}}.bact_reads.txt",
        output:
            temp(f"{EUK_OUT}/{{sample}}.euk_reads.txt"),
        log:
            f"{LOGS}/{{sample}}__comm_euk.log",
        shell:
            "comm -23 <(sort {input.all}) <(sort {input.bact}) > {output} 2> {log}"

    rule extract_euk_fastq_lca:
        input:
            fq=f"{PREPROC_OUT}/{{sample}}.ppm.vs.d4.fq.gz",
            ids=f"{EUK_OUT}/{{sample}}.euk_reads.txt",
        output:
            f"{EUK_OUT}/{{sample}}.euk.fastq.gz",
        log:
            f"{LOGS}/{{sample}}__seqtk_euk.log",
        params:
            bin=config["tools"]["seqtk"],
        shell:
            "{params.bin} subseq {input.fq} {input.ids} 2> {log} | gzip > {output}"
