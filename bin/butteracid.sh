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

# fastp integration config
: "${FASTP_REPORT_BASE:?FASTP_REPORT_BASE not set in config file}"
: "${FASTP_REPORT_PATTERN:?FASTP_REPORT_PATTERN not set in config file}"

SCRATCH_DIR="${SCRATCH_DIR:-/tmp}"

module load samtools
TAB=$'\t'

# ---------------------------
# Helper: get 16 seqkit stats for one FASTQ
# ---------------------------
get_seqkit_stats_16() {
  local fq="$1"
  local lib_id="$2"

  if [[ -f "$fq" ]]; then
    local line
    # seqkit stats -T outputs TSV; we take columns 4..end from row 2
	line=$(seqkit stats -a -T "$fq" \
	  | awk -F $'\t' 'NR==2 {for(i=4;i<=NF;i++) printf "%s%s", $i,(i==NF?ORS:OFS)}' OFS=$'\t')

    local stats=()
    # Read as TSV (not whitespace)
    IFS=$'\t' read -r -a stats <<< "$line"

    local nf=${#stats[@]}
    if (( nf != 16 )); then
      echo "[WARN] $lib_id: $fq produced $nf fields instead of 16 from seqkit. Padding/truncating." >&2
    fi

    while (( nf < 16 )); do stats+=("NA"); ((nf++)); done
    if (( nf > 16 )); then stats=("${stats[@]:0:16}"); fi

    printf "%s" "${stats[0]}"
    for ((i=1;i<16;i++)); do printf "\t%s" "${stats[i]}"; done
    printf "\n"
  else
    echo "[WARN] $lib_id: FASTQ not found: $fq. Filling with NA." >&2
    printf "NA"
    for ((i=2;i<=16;i++)); do printf "\tNA"; done
    printf "\n"
  fi
}
export -f get_seqkit_stats_16


# ---------------------------
# fastp HTML extraction helpers (per file)
# ---------------------------
extract_cell_text() {
  awk -F'>' '{print $5}' | awk -F'<' '{print $1}'
}

nth_match() {
  local pattern="$1"; local n="$2"; local file="$3"
  grep -Fi "$pattern" "$file" 2>/dev/null | extract_cell_text | awk -v n="$n" 'NR==n {print; exit}'
}

nth_paren_percent() {
  local pattern="$1"; local n="$2"; local file="$3"
  grep -Fi "$pattern" "$file" 2>/dev/null \
    | extract_cell_text \
    | awk -v n="$n" 'NR==n {print; exit}' \
    | awk -F'(' '{print $2}' \
    | awk -F'%' '{print $1}'
}

normalize_count() {
  local s="${1:-}"
  s="$(echo "$s" | tr -d ' ,\t\r\n')"
  [[ -z "$s" ]] && { echo ""; return; }

  echo "$s" | awk '
    BEGIN { IGNORECASE=1 }
    {
      v=$0
      suffix=""
      if (v ~ /[kKmM]$/) {
        suffix=substr(v, length(v), 1)
        v=substr(v, 1, length(v)-1)
      }
      if (v == "" || v ~ /[^0-9.]/) { print $0; exit }
      mult=1
      if (suffix=="K" || suffix=="k") mult=1000
      if (suffix=="M" || suffix=="m") mult=1000000
      n = v + 0
      out = n * mult
      printf "%.0f\n", out
    }'
}

# Find the fastp report for a given library id, using a config pattern.
# Supports globs in FASTP_REPORT_PATTERN.
find_fastp_html() {
  local lib_id="$1"
  lib_id="${lib_id//$'\r'/}"
  local f="$FASTP_REPORT_BASE/${lib_id}.fastp.report.html"
  [[ -f "$f" ]] && echo "$f" || echo ""
}
export -f find_fastp_html


# Return 9 fastp fields for one library (tab-separated, but we’ll store as array)
# Fields:
# sequencing_mode duplication_rate_pct total_reads_before total_reads_after
# reads_too_short_pct low_complexity_pct low_quality_pct gc_content insert_size_peak
get_fastp_fields_9() {
  local lib_id="$1"
  local html
  html="$(find_fastp_html "$lib_id")"

  if [[ -z "$html" || ! -f "$html" ]]; then
    echo "[WARN] $lib_id: fastp HTML not found (in $FASTP_REPORT_BASE). Filling with NA." >&2
    printf "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"
    return 0
  fi

  local sequencing_mode duplication_rate total_before_raw total_after_raw total_before total_after
  local reads_too_short low_complexity low_quality gc_content insert_size_peak

  sequencing_mode="$(nth_match "sequencing" 1 "$html")"

  duplication_rate="$(
    grep -Fi "duplication rate" "$html" 2>/dev/null \
      | extract_cell_text \
      | awk -F'%' 'NR==1 {print $1; exit}'
  )"

  total_before_raw="$(nth_match "total reads" 1 "$html")"
  total_after_raw="$(nth_match "total reads" 2 "$html")"
  total_before="$(normalize_count "$total_before_raw")"
  total_after="$(normalize_count "$total_after_raw")"

  reads_too_short="$(nth_paren_percent "reads too short" 1 "$html")"
  low_complexity="$(nth_paren_percent "low complexity" 1 "$html")"

  low_quality="$(
    grep -Fi "low quality" "$html" 2>/dev/null \
      | extract_cell_text \
      | awk 'NR==1 {print; exit}' \
      | awk -F'(' '{print $2}' \
      | awk -F'%' '{print $1}'
  )"

  gc_content="$(nth_match "GC content" 2 "$html")"
  insert_size_peak="$(nth_match "Insert size peak" 1 "$html")"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${sequencing_mode:-NA}" \
    "${duplication_rate:-NA}" \
    "${total_before:-NA}" \
    "${total_after:-NA}" \
    "${reads_too_short:-NA}" \
    "${low_complexity:-NA}" \
    "${low_quality:-NA}" \
    "${gc_content:-NA}" \
    "${insert_size_peak:-NA}"
}

export -f extract_cell_text nth_match nth_paren_percent normalize_count find_fastp_html get_fastp_fields_9

# ---------------------------
# Function to process a single library
# ---------------------------
process_library() {
    local LIB_ID="$1"
	LIB_ID="${LIB_ID//$'\r'/}"

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

    # fastp fields (9)
    local FASTP_FIELDS=()
	IFS=$'\t' read -r -a FASTP_FIELDS <<< "$(get_fastp_fields_9 "$LIB_ID")"
    if (( ${#FASTP_FIELDS[@]} != 9 )); then
      echo "[WARN] $LIB_ID: fastp fields count ${#FASTP_FIELDS[@]} (expected 9). Padding/truncating." >&2
      while (( ${#FASTP_FIELDS[@]} < 9 )); do FASTP_FIELDS+=("NA"); done
      FASTP_FIELDS=("${FASTP_FIELDS[@]:0:9}")
    fi

    # FASTQ files
    local FASTQS=(
        $FASTQ_LANE_BASE/${LIB_ID}_L001_R1_001.fastq.gz
        $FASTQ_LANE_BASE/${LIB_ID}_L002_R2_001.fastq.gz
        $FASTQ_LANE_BASE/${LIB_ID}_L003_R1_001.fastq.gz
        $FASTQ_LANE_BASE/${LIB_ID}_L004_R1_001.fastq.gz
        "$FASTQ_LOW_COMPLEXITY_BASE/${LIB_ID}.ppm.vs.d4.fq.gz"
        "$FASTQ_PREFILTER_BASE/${LIB_ID}.euk.fastq.gz"
    )

    local SEQKIT_LINE=()
    local fq
    for fq in "${FASTQS[@]}"; do
		IFS=$'\t' read -r -a stats <<< "$(get_seqkit_stats_16 "$fq" "$LIB_ID")"
        SEQKIT_LINE+=("${stats[@]}")
    done

    local total_stats=${#SEQKIT_LINE[@]}
    if (( total_stats != 96 )); then
        echo "[WARN] $LIB_ID: collected $total_stats seqkit fields (expected 96). Check FASTQs above." >&2
    fi

	printf "%s\t%s\t%s" "$LIB_ID" "$STAT_G" "$STAT_H"
	for v in "${FASTP_FIELDS[@]}"; do printf "\t%s" "$v"; done
	for v in "${SEQKIT_LINE[@]}"; do printf "\t%s" "$v"; done
	printf "\n"
    
}

export -f process_library
export FASTQ_LANE_BASE FASTQ_LOW_COMPLEXITY_BASE FASTQ_PREFILTER_BASE \
       BAM_BASE PREFILTER_BAM_BASE SCRATCH_DIR \
       FASTP_REPORT_BASE FASTP_REPORT_PATTERN

# ---------------------------
# QC generation block
# ---------------------------
if [[ -s "$QC_OUTPUT" && "${FORCE_QC:-false}" == false ]]; then
    echo "[INFO] Found existing non-empty QC table at: $QC_OUTPUT"
    echo "[INFO] Skipping QC generation and using existing file."
else
    echo "[INFO] Generating QC table at: $QC_OUTPUT"

	HEADER="library_id total_no_reads_mapped_uniq total_no_reads_mapped_uniq_prefilter \
	fastp_sequencing_mode fastp_duplication_rate_pct fastp_total_reads_before_filtering fastp_total_reads_after_filtering fastp_reads_too_short_pct fastp_low_complexity_pct fastp_low_quality_pct fastp_gc_content fastp_insert_size_peak \
	R1_num_seqs R1_sum_len R1_min_len R1_avg_len R1_max_len R1_Q1 R1_Q2 R1_Q3 R1_sum_gap R1_N50 R1_N50_num R1_Q20(%) R1_Q30(%) R1_AvgQual R1_GC(%) R1_sum_n \
	R2_num_seqs R2_sum_len R2_min_len R2_avg_len R2_max_len R2_Q1 R2_Q2 R2_Q3 R2_sum_gap R2_N50 R2_N50_num R2_Q20(%) R2_Q30(%) R2_AvgQual R2_GC(%) R2_sum_n \
	R3_num_seqs R3_sum_len R3_min_len R3_avg_len R3_max_len R3_Q1 R3_Q2 R3_Q3 R3_sum_gap R3_N50 R3_N50_num R3_Q20(%) R3_Q30(%) R3_AvgQual R3_GC(%) R3_sum_n \
	R4_num_seqs R4_sum_len R4_min_len R4_avg_len R4_max_len R4_Q1 R4_Q2 R4_Q3 R4_sum_gap R4_N50 R4_N50_num R4_Q20(%) R4_Q30(%) R4_AvgQual R4_GC(%) R4_sum_n \
	microb_in_num_seqs microb_in_sum_len microb_in_min_len microb_in_avg_len microb_in_max_len microb_in_Q1 microb_in_Q2 microb_in_Q3 microb_in_sum_gap microb_in_N50 microb_in_N50_num microb_in_Q20(%) microb_in_Q30(%) microb_in_AvgQual microb_in_GC(%) microb_in_sum_n \
	prefiltered_num_seqs prefiltered_sum_len prefiltered_min_len prefiltered_avg_len prefiltered_max_len prefiltered_Q1 prefiltered_Q2 prefiltered_Q3 prefiltered_sum_gap prefiltered_N50 prefiltered_N50_num prefiltered_Q20(%) prefiltered_Q30(%) prefiltered_AvgQual prefiltered_GC(%) prefiltered_sum_n"

	echo "$HEADER" | tr ' ' '\t' > "$QC_OUTPUT"
    

    parallel -j 10 process_library :::: "$LIBRARY_LIST" >> "$QC_OUTPUT"
    echo "[INFO] QC table with mapped reads + fastp + full seqkit stats written to: $QC_OUTPUT"
fi

# ---------------------------
# Getting tsv files: R-SCRIPT
# ---------------------------
case "$MODE" in
    bamfilter)
        Rscript butteracid_bamfilter.R "$CONFIG_FILE"
        ;;
    unicorn)
        #Rscript butteracid_unicorn.R "$CONFIG_FILE"
        ;;
    *)
        echo "[ERROR] Unknown MODE: '$MODE'. Expected 'bamfilter' or 'unicorn'." >&2
        exit 1
        ;;
esac
