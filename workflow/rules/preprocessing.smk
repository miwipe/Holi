# ============================================================
# Stage 1 — Preprocessing
# fastp (trim/merge) -> vsearch (dedup) -> sga (low-complexity filter) -> gzip
#
# Outputs live under PREPROC_OUT instead of the working directory
# holi.sh used to write them to — keeps per-sample intermediates
# out of the repo root and makes them easy to find.
# ============================================================

if SINGLE_END:

    rule fastp_se:
        input:
            r1=lambda wc: get_se_fastq(wc.sample),
        output:
            merged=temp(f"{PREPROC_OUT}/{{sample}}.ppm.fq"),
            html=f"{PREPROC_OUT}/{{sample}}.fastp.report.html",
        log:
            f"{LOGS}/{{sample}}__fastp.log",
        params:
            qual=config["fastp"]["qual"],
            min_avg_qual=config["fastp"]["min_avg_qual"],
            min_len=config["fastp"]["min_len"],
            dup_accuracy=config["fastp"]["dup_accuracy"],
            bin=config["tools"]["fastp"],
        threads: 1
        shell:
            """
            {params.bin} \
              -i {input.r1} \
              -o {output.merged} \
              -V \
              -D --dup_calc_accuracy {params.dup_accuracy} \
              -g -x -q {params.qual} -e {params.min_avg_qual} -l {params.min_len} -y -c -p \
              -h {output.html} \
              -w 1 > {log} 2>&1
            """

else:

    rule fastp_pe:
        input:
            r1=lambda wc: get_r1(wc.sample),
            r2=lambda wc: get_r2(wc.sample),
        output:
            merged=temp(f"{PREPROC_OUT}/{{sample}}.ppm.fq"),
            html=f"{PREPROC_OUT}/{{sample}}.fastp.report.html",
        log:
            f"{LOGS}/{{sample}}__fastp.log",
        params:
            qual=config["fastp"]["qual"],
            min_avg_qual=config["fastp"]["min_avg_qual"],
            min_len=config["fastp"]["min_len"],
            dup_accuracy=config["fastp"]["dup_accuracy"],
            bin=config["tools"]["fastp"],
        threads: 1
        shell:
            """
            {params.bin} \
              -i {input.r1} -I {input.r2} \
              -m --merged_out {output.merged} \
              -V --detect_adapter_for_pe \
              -D --dup_calc_accuracy {params.dup_accuracy} \
              -g -x -q {params.qual} -e {params.min_avg_qual} -l {params.min_len} -y -c -p \
              -h {output.html} \
              -w 1 > {log} 2>&1
            """


rule vsearch_dedup:
    input:
        f"{PREPROC_OUT}/{{sample}}.ppm.fq",
    output:
        temp(f"{PREPROC_OUT}/{{sample}}.ppm.vs.fq"),
    log:
        f"{LOGS}/{{sample}}__vsearch.log",
    params:
        min_seqlen=config["vsearch"]["min_seqlen"],
        bin=config["tools"]["vsearch"],
    threads: 1
    shell:
        """
        {params.bin} \
          --fastx_uniques {input} \
          --fastqout {output} \
          --minseqlength {params.min_seqlen} \
          --strand both > {log} 2>&1
        """


rule sga_lowcomplexity:
    input:
        f"{PREPROC_OUT}/{{sample}}.ppm.vs.fq",
    output:
        temp(f"{PREPROC_OUT}/{{sample}}.ppm.vs.d4.fq"),
    log:
        f"{LOGS}/{{sample}}__sga.log",
    params:
        dust=config["sga"]["dust_threshold"],
        minlen=config["sga"]["min_len"],
        bin=config["tools"]["sga"],
    threads: 1
    shell:
        """
        {params.bin} preprocess --dust-threshold={params.dust} -m {params.minlen} \
          {input} -o {output} > {log} 2>&1
        """


rule gzip_preprocessed:
    input:
        f"{PREPROC_OUT}/{{sample}}.ppm.vs.d4.fq",
    output:
        f"{PREPROC_OUT}/{{sample}}.ppm.vs.d4.fq.gz",
    shell:
        "gzip -c {input} > {output}"
