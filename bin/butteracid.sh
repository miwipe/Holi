#!/bin/bash
set -euo pipefail

# ---------------------------
# Optional flag: --force-qc
# ---------------------------
FORCE_QC=false
if [[ "${1:-}" == "--force-qc" || "${1:-}" == "--force-QC" ]]; then
    FORCE_QC=true
    shift
fi

# ---------------------------
# Load configuration
# ---------------------------
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 [--force-qc] CONFIG_FILE"
    exit 1
fi

CONFIG_FILE="$1"
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "[ERROR] Config file '$CONFIG_FILE' not found!"
    exit 1
fi
source "$CONFIG_FILE"

# Ensure required config variables are defined
: "${MODE:?MODE not set in config file}"
: "${LIBRARY_LIST:?LIBRARY_LIST not set in config file}"
: "${QC_OUTPUT:?QC_OUTPUT not set in config file}"
: "${FASTQ_LANE_BASE:?FASTQ_LANE_BASE not set in config file}"
: "${FASTQ_LOW_COMPLEXITY_BASE:?FASTQ_LOW_COMPLEXITY_BASE not set in config file}"
: "${FASTQ_PREFILTER_BASE:?FASTQ_PREFILTER_BASE not set in config file}"
: "${BAM_BASE:?BAM_BASE not set in config file}"
: "${PREFILTER_BAM_BASE:?PREFILTER_BAM_BASE not set in config file}"
SCRATCH_DIR="${SCRATCH_DIR:-/tmp}"   # optional with default

module load samtools


# ---------------------------
# Helper: get 16 seqkit stats for one FASTQ
# ---------------------------
get_seqkit_stats_16() {
    local fq="$1"
    local lib_id="$2"

    if [[ -f "$fq" ]]; then
        # Take columns 4..NF (num_seqs .. sum_n) as a single space-separated line
        local line
        line=$(seqkit stats -a -T "$fq" \
               | awk 'NR==2 {for(i=4;i<=NF;i++) printf "%s%s", $i,(i==NF?ORS:OFS)}')

        # Split into bash array
        local stats=()
        read -r -a stats <<< "$line"

        local nf=${#stats[@]}
        if (( nf != 16 )); then
            echo "[WARN] $lib_id: $fq produced $nf fields instead of 16 from seqkit. Padding/truncating." >&2
        fi

        # Pad or truncate to exactly 16 fields
        while (( nf < 16 )); do
            stats+=("NA")
            ((nf++))
        done
        if (( nf > 16 )); then
            stats=("${stats[@]:0:16}")
        fi

        echo "${stats[*]}"
    else
        # File missing → log + 16 NAs
        echo "[WARN] $lib_id: FASTQ not found: $fq. Filling with NA." >&2
        local na_fields=()
        for _ in {1..16}; do
            na_fields+=("NA")
        done
        echo "${na_fields[*]}"
    fi
}

export -f get_seqkit_stats_16

# ---------------------------
# Function to process a single library
# ---------------------------
process_library() {
    local LIB_ID="$1"

    # mapped reads
    local MAPPED_READS_UNIQ="$BAM_BASE/${LIB_ID}.gtdb.merged.bam"
    local STAT_G="NA"
    if [[ -f "$MAPPED_READS_UNIQ" ]]; then
        STAT_G=$(samtools view "$MAPPED_READS_UNIQ" \
                 | cut -f1 \
                 | sort -u --temporary-directory="$SCRATCH_DIR" \
                 | wc -l)
    fi

    local MAPPED_READS_UNIQ_P="$PREFILTER_BAM_BASE/${LIB_ID}.comp.bam"
    local STAT_H="NA"
    if [[ -f "$MAPPED_READS_UNIQ_P" ]]; then
        STAT_H=$(samtools view "$MAPPED_READS_UNIQ_P" \
                 | cut -f1 \
                 | sort -u --temporary-directory="$SCRATCH_DIR" \
                 | wc -l)
    fi

    # FASTQ files (adjust patterns if you need globs)
    local FASTQS=(
        $FASTQ_LANE_BASE/${LIB_ID}_R1_001.fastq.gz
        $FASTQ_LANE_BASE/${LIB_ID}_R2_001.fastq.gz
        $FASTQ_LANE_BASE/${LIB_ID}_L003_R1_001.fastq.gz
        $FASTQ_LANE_BASE/${LIB_ID}_L004_R1_001.fastq.gz
        "$FASTQ_LOW_COMPLEXITY_BASE/${LIB_ID}.ppm.vs.d4.fq.gz"
        "$FASTQ_PREFILTER_BASE/${LIB_ID}.euk.fastq.gz"
    )

    # Collect all 6×16 stats into one array
    local SEQKIT_LINE=()
    local fq
    for fq in "${FASTQS[@]}"; do
        read -r -a stats <<< "$(get_seqkit_stats_16 "$fq" "$LIB_ID")"
        SEQKIT_LINE+=("${stats[@]}")
    done

    local total_stats=${#SEQKIT_LINE[@]}
    if (( total_stats != 96 )); then
        echo "[WARN] $LIB_ID: collected $total_stats seqkit fields (expected 96). Check FASTQs above." >&2
    fi

    # Final output line: plain space-separated, no quotes
    echo "${LIB_ID} ${STAT_G} ${STAT_H} ${SEQKIT_LINE[*]}"
}

export -f process_library
export FASTQ_LANE_BASE FASTQ_LOW_COMPLEXITY_BASE FASTQ_PREFILTER_BASE \
       BAM_BASE PREFILTER_BAM_BASE SCRATCH_DIR

# ---------------------------
# QC generation block
# ---------------------------
if [[ -s "$QC_OUTPUT" && "${FORCE_QC:-false}" == false ]]; then
    echo "[INFO] Found existing non-empty QC table at: $QC_OUTPUT"
    echo "[INFO] Skipping QC generation and using existing file."
else
    echo "[INFO] Generating QC table at: $QC_OUTPUT"

    # Header as a single logical line (space-separated), no "\" written to the file
    echo "library_id total_no_reads_mapped_uniq total_no_reads_mapped_uniq_prefilter \
R1_num_seqs R1_sum_len R1_min_len R1_avg_len R1_max_len R1_Q1 R1_Q2 R1_Q3 R1_sum_gap R1_N50 R1_N50_num R1_Q20(%) R1_Q30(%) R1_AvgQual R1_GC(%) R1_sum_n \
R2_num_seqs R2_sum_len R2_min_len R2_avg_len R2_max_len R2_Q1 R2_Q2 R2_Q3 R2_sum_gap R2_N50 R2_N50_num R2_Q20(%) R2_Q30(%) R2_AvgQual R2_GC(%) R2_sum_n \
R3_num_seqs R3_sum_len R3_min_len R3_avg_len R3_max_len R3_Q1 R3_Q2 R3_Q3 R3_sum_gap R3_N50 R3_N50_num R3_Q20(%) R3_Q30(%) R3_AvgQual R3_GC(%) R3_sum_n \
R4_num_seqs R4_sum_len R4_min_len R4_avg_len R4_max_len R4_Q1 R4_Q2 R4_Q3 R4_sum_gap R4_N50 R4_N50_num R4_Q20(%) R4_Q30(%) R4_AvgQual R4_GC(%) R4_sum_n \
microb_in_num_seqs microb_in_sum_len microb_in_min_len microb_in_avg_len microb_in_max_len microb_in_Q1 microb_in_Q2 microb_in_Q3 microb_in_sum_gap microb_in_N50 microb_in_N50_num microb_in_Q20(%) microb_in_Q30(%) microb_in_AvgQual microb_in_GC(%) microb_in_sum_n \
prefiltered_num_seqs prefiltered_sum_len prefiltered_min_len prefiltered_avg_len prefiltered_max_len prefiltered_Q1 prefiltered_Q2 prefiltered_Q3 prefiltered_sum_gap prefiltered_N50 prefiltered_N50_num prefiltered_Q20(%) prefiltered_Q30(%) prefiltered_AvgQual prefiltered_GC(%) prefiltered_sum_n" > "$QC_OUTPUT"

    # Run all libraries in parallel and append to QC_OUTPUT
    parallel -j 10 process_library :::: "$LIBRARY_LIST" >> "$QC_OUTPUT"

    echo "[INFO] QC table with mapped reads + full seqkit stats written to: $QC_OUTPUT"
fi
# ---------------------------
# Getting tsv files: R-SCRIPT
# ---------------------------
case "$MODE" in
    bamfilter)
        Rscript butteracid_bamfilter.R "$CONFIG_FILE"
        ;;
    unicorn)
        Rscript butteracid_unicorn.R "$CONFIG_FILE"
        ;;
    *)
        echo "[ERROR] Unknown MODE: '$MODE'. Expected 'bamfilter' or 'unicorn'." >&2
        exit 1
        ;;
esac
