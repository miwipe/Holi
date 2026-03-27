#!/bin/bash
# Script: holi_prefilter.sh
# Description: Comprehensive automated pipeline for preprocessing, mapping, filtering, and taxonomic classification of FASTQ files.

# -------------------------------
# Usage information
# -------------------------------
usage() {
    echo "Usage: $0 <config.yml> <sample.list> [OPTIONS]"
    echo
    echo "Arguments:"
    echo "  <config.yml>            Path to YAML configuration file."
    echo "  <sample.list>           File containing list of sample names."
    echo
    echo "Stage control (stages run in order 1-8):"
    echo "  Stage 1  preprocessing"
    echo "  Stage 2  gtdb-mapping"
    echo "  Stage 3  microbial-split"
    echo "  Stage 4  euk-mapping"
    echo "  Stage 5  comp-merge"
    echo "  Stage 6  bam-filtering"
    echo "  Stage 7  metadmg"
    echo "  Stage 8  unicorn-tidstats"
    echo
    echo "  --from-stage <N>        Skip all stages before N (inclusive shortcut)."
    echo "  --to-stage <N>          Skip all stages after N (inclusive shortcut)."
    echo
    echo "Fine-grained skip options:"
    echo "  --skip-preprocessing              Skip trimming, merging, duplicate removal, and complexity filtering."
    echo "  --skip-preprocessing-cleanup      Do NOT delete preprocessing intermediates after final .ppm.vs.d4.fq.gz exists."
    echo "  --skip-gtdb-mapping               Skip GTDB mapping + merge/sort steps."
    echo "  --skip-gtdb-cleanup               Do NOT delete GTDB chunk BAMs after .gtdb.merged.bam exists."
    echo "  --skip-microbial-split            Skip bacterial/eukaryotic read ID extraction + FASTQ subsetting."
    echo "  --skip-euk-mapping                Skip all eukaryotic mapping steps (129 parts + mito + phyNor + core_nt + plastid)."
    echo "  --skip-euk-cleanup                Do NOT delete euk/core_nt/phyNor/mito/pla intermediates after .comp.bam exists."
    echo "  --skip-comp-merge                 Skip compress+sort+merge that produces *.comp.bam (assumes it already exists)."
    echo "  --skip-bam-filtering              Skip Unicorn filtering step (alnfilt → refstats → bamstats)."
    echo "  --skip-metadmg                    Skip metaDMG lca/dfit/aggregate steps."
    echo "  --skip-unicorn-tidstats           Skip final unicorn tidstats step."
    echo
    echo "Mode options:"
    echo "  --single-end            Treat input as single-end FASTQ (no R2, no merging)."
    echo "                          Looks for files matching: *_<sample>_*R1*_001.fastq.gz"
    echo "  --storage-friendly      Delete intermediate files, keeping only final outputs per stage."
    echo "  --lca-assignment        Use metaDMG LCA for microbial read splitting (default: getRTax)."
    echo "  --force-stage <N>       Force pipeline to start from stage N, bypassing auto-detection."
    echo
    echo "Examples:"
    echo "  $0 config.yml sample.list"
    echo "  $0 config.yml sample.list --from-stage 4"
    echo "  $0 config.yml sample.list --force-stage 4"
    echo "  $0 config.yml sample.list --from-stage 4 --to-stage 6"
    echo "  $0 config.yml sample.list --from-stage 4 --skip-euk-cleanup"
    echo "  $0 config.yml sample.list --skip-preprocessing --skip-gtdb-mapping"
    exit 1
}

# Show usage if no arguments are provided
if [ $# -lt 2 ]; then
    usage
fi

# -------------------------------
# Check for required tools
# -------------------------------

REQUIRED_TOOLS=(
  parallel
  fastp
  vsearch
  sga
  bowtie2
  samtools
  conda
  seqtk
  unicorn
  metaDMG-cpp
)

echo "Checking required tools..."

missing_tools=()

for tool in "${REQUIRED_TOOLS[@]}"; do
  if ! command -v "$tool" &> /dev/null; then
    missing_tools+=("$tool")
  fi
done

if [ ${#missing_tools[@]} -ne 0 ]; then
  echo "The following required tools are missing:"
  for t in "${missing_tools[@]}"; do
    echo "  - $t"
  done
  echo "Please install them and re-run the script."
  exit 1
else
  echo "All required tools found."
fi

# -------------------------------
# Load config file
# -------------------------------

CONFIG="${1:-config.yml}"
SAMPLE_LIST="$2"

load_config() {
    CONFIG_FILE="$1"
    get_value() {
        grep "^$1:" "$CONFIG_FILE" | sed 's/^.*:[[:space:]]*//' | sed 's/"//g'
    }

    # Infrastructure
    LOG_FILE=$(get_value "LOG_FILE")
    THREADSP=$(get_value "THREADSP")
    THREADS=$(get_value "THREADS")
    INPATH=$(get_value "INPATH")
    RESULT_PATH=$(get_value "RESULT_PATH")
    MICROB_OUT=$(get_value "MICROB_OUT")
    EUK_OUT=$(get_value "EUK_OUT")
    LOGS=$(get_value "LOGS")
    TMP=$(get_value "TMP")

    # Databases
    DB_PATH=$(get_value "DB_PATH")
    DB_PATH_clean=$(get_value "DB_PATH_clean")
    DB_PATH_Norwary=$(get_value "DB_PATH_Norwary")
    DB_PATH_bac=$(get_value "DB_PATH_bac")

    # Taxonomy
    TAX_PATH_NCBI=$(get_value "TAX_PATH_NCBI")
    TAX_PATH_BAC=$(get_value "TAX_PATH_BAC")
    TAX_PATH_BAC_ACC=$(get_value "TAX_PATH_BAC_ACC")

    # Stage 1 — fastp
    FASTP_QUAL=$(get_value "FASTP_QUAL")
    FASTP_MIN_AVG_QUAL=$(get_value "FASTP_MIN_AVG_QUAL")
    FASTP_MIN_LEN=$(get_value "FASTP_MIN_LEN")
    FASTP_DUP_ACCURACY=$(get_value "FASTP_DUP_ACCURACY")

    # Stage 1 — vsearch
    VSEARCH_MIN_SEQLEN=$(get_value "VSEARCH_MIN_SEQLEN")

    # Stage 1 — sga
    SGA_DUST_THRESHOLD=$(get_value "SGA_DUST_THRESHOLD")
    SGA_MIN_LEN=$(get_value "SGA_MIN_LEN")

    # Stage 2 — bowtie2 GTDB
    BT2_GTDB_k=$(get_value "BT2_GTDB_k")
    BT2_GTDB_D=$(get_value "BT2_GTDB_D")
    BT2_GTDB_R=$(get_value "BT2_GTDB_R")
    BT2_GTDB_N=$(get_value "BT2_GTDB_N")
    BT2_GTDB_L=$(get_value "BT2_GTDB_L")
    BT2_GTDB_i=$(get_value "BT2_GTDB_i")
    BT2_GTDB_np=$(get_value "BT2_GTDB_np")
    BT2_GTDB_mp=$(get_value "BT2_GTDB_mp")
    BT2_GTDB_rdg=$(get_value "BT2_GTDB_rdg")
    BT2_GTDB_rfg=$(get_value "BT2_GTDB_rfg")
    BT2_GTDB_score_min=$(get_value "BT2_GTDB_score_min")

    # Stage 3 — metaDMG LCA microbial
    METADMG_BAC_SIM_SCORE_LOW=$(get_value "METADMG_BAC_SIM_SCORE_LOW")
    METADMG_BAC_SIM_SCORE_HIGH=$(get_value "METADMG_BAC_SIM_SCORE_HIGH")
    METADMG_BAC_HOW_MANY=$(get_value "METADMG_BAC_HOW_MANY")
    METADMG_BAC_WEIGHT_TYPE=$(get_value "METADMG_BAC_WEIGHT_TYPE")

    # Stage 4 — bowtie2 euk
    BT2_EUK_k=$(get_value "BT2_EUK_k")

    # Stage 6 — unicorn alnfilt
    ALNFILT_MODE=$(get_value "ALNFILT_MODE")
    ALNFILT_MINANI=$(get_value "ALNFILT_MINANI")
    ALNFILT_MAXANI=$(get_value "ALNFILT_MAXANI")

    # Stage 6 — unicorn refstats
    REFSTATS_MINREADS=$(get_value "REFSTATS_MINREADS")

    # Stage 7 — metaDMG LCA euk
    METADMG_EUK_SIM_SCORE_LOW=$(get_value "METADMG_EUK_SIM_SCORE_LOW")
    METADMG_EUK_SIM_SCORE_HIGH=$(get_value "METADMG_EUK_SIM_SCORE_HIGH")
    METADMG_EUK_HOW_MANY=$(get_value "METADMG_EUK_HOW_MANY")
    METADMG_EUK_WEIGHT_TYPE=$(get_value "METADMG_EUK_WEIGHT_TYPE")

    # Stage 7 — metaDMG dfit
    METADMG_DFIT_NOPT=$(get_value "METADMG_DFIT_NOPT")
    METADMG_DFIT_NBOOTSTRAP=$(get_value "METADMG_DFIT_NBOOTSTRAP")
    METADMG_DFIT_SEED=$(get_value "METADMG_DFIT_SEED")
    METADMG_DFIT_LIB=$(get_value "METADMG_DFIT_LIB")

    # Stage 8 — unicorn taxstats
    TAXSTATS_k=$(get_value "TAXSTATS_k")
    TAXSTATS_MINREADS=$(get_value "TAXSTATS_MINREADS")
}

load_config "$CONFIG"
echo "[INFO] Loaded config from $CONFIG"

# Ensure TMP exists and is used by Python tools (filterBAM, etc.)
: "${TMP:=/tmp}"
mkdir -p "$TMP" || { echo "[ERROR] Failed to create TMP dir: $TMP"; exit 1; }
export TMPDIR="$TMP"
export TEMP="$TMP"
export TMP="$TMP"

# -----------------------------
# Logging helpers
# -----------------------------
GREEN="\033[0;32m"
YELLOW="\033[1;33m"
RED="\033[0;31m"
RESET="\033[0m"

timestamp() { date +"%Y-%m-%d %H:%M:%S"; }

log_step()  { printf "${YELLOW}[%s][STEP]${RESET} %s\n"  "$(timestamp)" "$1"; }
log_info()  { printf "${GREEN}[%s][INFO]${RESET} %s\n"   "$(timestamp)" "$1"; }
log_error() { printf "${RED}[%s][ERROR]${RESET} %s\n"    "$(timestamp)" "$1"; }

check_success() {
    if [ $? -eq 0 ]; then
        log_info "$1 completed"
    else
        log_error "$1 failed. Stopping."
        exit 1
    fi
}

# -----------------------------
# Ensure essential output directories exist
# -----------------------------
REQUIRED_DIRS=("$RESULT_PATH" "$EUK_OUT" "$MICROB_OUT" "$LOGS" "$TMP")

for dir in "${REQUIRED_DIRS[@]}"; do
    if [ ! -d "$dir" ]; then
        log_info "Directory $dir does not exist. Creating it..."
        mkdir -p "$dir" || { log_error "Failed to create directory $dir. Exiting."; exit 1; }
    else
        log_info "Directory $dir exists."
    fi
done

# ---------------------------
# Parse flags
# ---------------------------

SKIP_PREPROCESSING=false
SKIP_PREPROCESSING_CLEANUP=false
SKIP_GTDB_MAPPING=false
SKIP_GTDB_CLEANUP=false
SKIP_MICROBIAL_SPLIT=false
SKIP_EUK_MAPPING=false
SKIP_EUK_CLEANUP=false
SKIP_COMP_MERGE=false
SKIP_BAM_FILTERING=false
SKIP_METADMG=false
SKIP_UNICORN_TIDSTATS=false

LCA_ASSIGN=false
SINGLE_END=false
STORAGE_FRIENDLY=false

FROM_STAGE=1
TO_STAGE=8
FORCE_STAGE=""
USER_SET_FROM_STAGE=false

shift 2   # consume positional args (config + sample list)

while [[ $# -gt 0 ]]; do
  case "$1" in
    --from-stage)                 FROM_STAGE="$2"; USER_SET_FROM_STAGE=true; shift ;;
    --to-stage)                   TO_STAGE="$2";             shift ;;
    --force-stage)                FORCE_STAGE="$2";          shift ;;
    --single-end)                 SINGLE_END=true ;;
    --storage-friendly)           STORAGE_FRIENDLY=true ;;
    --skip-preprocessing)         SKIP_PREPROCESSING=true ;;
    --skip-preprocessing-cleanup) SKIP_PREPROCESSING_CLEANUP=true ;;
    --skip-gtdb-mapping)          SKIP_GTDB_MAPPING=true ;;
    --skip-gtdb-cleanup)          SKIP_GTDB_CLEANUP=true ;;
    --skip-microbial-split)       SKIP_MICROBIAL_SPLIT=true ;;
    --skip-euk-mapping)           SKIP_EUK_MAPPING=true ;;
    --skip-euk-cleanup)           SKIP_EUK_CLEANUP=true ;;
    --skip-comp-merge)            SKIP_COMP_MERGE=true ;;
    --skip-bam-filtering)         SKIP_BAM_FILTERING=true ;;
    --skip-metadmg)               SKIP_METADMG=true ;;
    --skip-unicorn-tidstats)      SKIP_UNICORN_TIDSTATS=true ;;
    --lca-assignment)             LCA_ASSIGN=true ;;
    *) log_error "Unknown option: $1"; usage ;;
  esac
  shift
done

# ---------------------------
# Auto-detection: check output integrity for each stage across all samples
# Returns the first stage where any sample is incomplete/corrupt.
# ---------------------------

check_bam()  { samtools quickcheck "$1" 2>/dev/null; }
check_gz() {
    # Check file is non-empty and starts with gzip magic bytes (1f 8b)
    # then verify it is not truncated by confirming the EOF marker exists
    [[ -s "$1" ]] || return 1
    local magic
    magic=$(xxd -l 2 -p "$1" 2>/dev/null)
    [[ "$magic" == "1f8b" ]] || return 1
    # Fast truncation check: last 2 bytes of a valid gzip must be non-zero
    # (they encode the original file size mod 2^32)
    local tail_bytes
    tail_bytes=$(tail -c 8 "$1" | wc -c)
    [[ "$tail_bytes" -ge 8 ]] || return 1
    return 0
}
check_file() { [[ -s "$1" ]]; }

# Returns 0 (complete) or 1 (incomplete/missing/corrupt) for a given stage
stage_complete() {
    local stage=$1

    while IFS= read -r sample; do
        local ok=true
        case "$stage" in
            1) check_gz  "${sample}.ppm.vs.d4.fq.gz"                               || ok=false ;;
            2) check_bam "$MICROB_OUT/${sample}.gtdb.merged.bam"                    || ok=false ;;
            3) check_gz  "$EUK_OUT/${sample}.euk.fastq.gz"                          || ok=false ;;
            4) # pla.bam may have been cleaned up — accept comp.bam as evidence stage 4+5 completed
               { check_bam "$EUK_OUT/${sample}.pla.bam" || check_bam "$EUK_OUT/${sample}.comp.bam"; } || ok=false ;;
            5) check_bam "$EUK_OUT/${sample}.comp.bam"                              || ok=false ;;
            6) check_bam "$EUK_OUT/${sample}.comp.filtered.bam"                     || ok=false ;;
            7) check_gz  "$EUK_OUT/${sample}.sort.comp.filtered.agg.stat.gz"        || ok=false ;;
            8) check_file "$EUK_OUT/${sample}.comp.filtered.species.taxstats"       || ok=false ;;
        esac

        if [ "$ok" = false ]; then
            log_info "Stage $stage incomplete for sample: $sample" >&2
            return 1
        fi
    done < "$SAMPLE_LIST"

    return 0
}

detect_resume_stage() {
    log_info "Auto-detecting resume stage..." >&2
    for stage in 1 2 3 4 5 6 7 8; do
        if ! stage_complete "$stage"; then
            log_info "Auto-detected start stage: $stage" >&2
            echo "$stage"
            return
        fi
    done
    echo "done"
}

# ---------------------------
# Resolve FROM_STAGE: force > explicit > auto-detect
# ---------------------------

if [ -n "$FORCE_STAGE" ]; then
    if ! [[ "$FORCE_STAGE" =~ ^[1-8]$ ]]; then
        log_error "--force-stage must be an integer between 1 and 8."
        exit 1
    fi
    FROM_STAGE="$FORCE_STAGE"
    log_info "Forced start stage: $FROM_STAGE (auto-detection skipped)"

elif [ "$USER_SET_FROM_STAGE" = false ]; then
    DETECTED=$(detect_resume_stage)
    if [ "$DETECTED" = "done" ]; then
        log_info "All stages appear complete for all samples. Nothing to do."
        log_info "Use --force-stage <N> to rerun from a specific stage."
        exit 0
    fi
    FROM_STAGE="$DETECTED"
    log_info "Resuming from stage: $FROM_STAGE"

else
    log_info "Manual --from-stage $FROM_STAGE provided. Skipping auto-detection."
fi

# ---------------------------
# Apply --from-stage / --to-stage
# These set skip flags but never override an already-set skip flag,
# so fine-grained --skip-X options always win.
# ---------------------------

# Validate stage range
if ! [[ "$FROM_STAGE" =~ ^[1-8]$ ]] || ! [[ "$TO_STAGE" =~ ^[1-8]$ ]]; then
    log_error "--from-stage and --to-stage must be integers between 1 and 8."
    exit 1
fi
if [ "$FROM_STAGE" -gt "$TO_STAGE" ]; then
    log_error "--from-stage ($FROM_STAGE) cannot be greater than --to-stage ($TO_STAGE)."
    exit 1
fi

# Skip stages before FROM_STAGE
[ "$FROM_STAGE" -gt 1 ] && SKIP_PREPROCESSING=true
[ "$FROM_STAGE" -gt 2 ] && SKIP_GTDB_MAPPING=true
[ "$FROM_STAGE" -gt 3 ] && SKIP_MICROBIAL_SPLIT=true
[ "$FROM_STAGE" -gt 4 ] && SKIP_EUK_MAPPING=true
[ "$FROM_STAGE" -gt 5 ] && SKIP_COMP_MERGE=true
[ "$FROM_STAGE" -gt 6 ] && SKIP_BAM_FILTERING=true
[ "$FROM_STAGE" -gt 7 ] && SKIP_METADMG=true

# Skip stages after TO_STAGE
[ "$TO_STAGE" -lt 8 ] && SKIP_UNICORN_TIDSTATS=true
[ "$TO_STAGE" -lt 7 ] && SKIP_METADMG=true
[ "$TO_STAGE" -lt 6 ] && SKIP_BAM_FILTERING=true
[ "$TO_STAGE" -lt 5 ] && SKIP_COMP_MERGE=true
[ "$TO_STAGE" -lt 4 ] && SKIP_EUK_MAPPING=true
[ "$TO_STAGE" -lt 3 ] && SKIP_MICROBIAL_SPLIT=true
[ "$TO_STAGE" -lt 2 ] && SKIP_GTDB_MAPPING=true

# ---------------------------
# Print resolved stage plan
# ---------------------------
log_info "Stage plan:"
log_info "  Stage 1  preprocessing      : $([ "$SKIP_PREPROCESSING"    = true ] && echo SKIP || echo RUN)"
log_info "  Stage 2  gtdb-mapping       : $([ "$SKIP_GTDB_MAPPING"     = true ] && echo SKIP || echo RUN)"
log_info "  Stage 3  microbial-split    : $([ "$SKIP_MICROBIAL_SPLIT"  = true ] && echo SKIP || echo RUN)"
log_info "  Stage 4  euk-mapping        : $([ "$SKIP_EUK_MAPPING"      = true ] && echo SKIP || echo RUN)"
log_info "  Stage 5  comp-merge         : $([ "$SKIP_COMP_MERGE"       = true ] && echo SKIP || echo RUN)"
log_info "  Stage 6  bam-filtering      : $([ "$SKIP_BAM_FILTERING"    = true ] && echo SKIP || echo RUN)"
log_info "  Stage 7  metadmg            : $([ "$SKIP_METADMG"          = true ] && echo SKIP || echo RUN)"
log_info "  Stage 8  unicorn-tidstats   : $([ "$SKIP_UNICORN_TIDSTATS" = true ] && echo SKIP || echo RUN)"

# ---------------------------
# Validate config and sample list
# ---------------------------

if [ ! -s "$CONFIG" ]; then
    echo "[ERROR] Config file ($CONFIG) is empty or does not exist." | tee -a "$LOG_FILE"
    exit 1
fi

if [ -z "$SAMPLE_LIST" ]; then
    echo "[ERROR] No sample.list provided." | tee -a "$LOG_FILE"
    exit 1
fi

if [ ! -s "$SAMPLE_LIST" ]; then
    echo "[ERROR] Sample list file ($SAMPLE_LIST) is empty or does not exist." | tee -a "$LOG_FILE"
    exit 1
fi

# ===========================================================
# STAGE 1 — Preprocessing
# ===========================================================

if [ "$SKIP_PREPROCESSING" = false ]; then
    ####################### QC ##############################

    if [ "$SINGLE_END" = true ]; then
        log_step "Trimming single-end reads with fastp..."

        export INPATH LOGS
        cat "$SAMPLE_LIST" | parallel -j "$THREADSP" --env INPATH --env LOGS '
          sample={}
          logfile="$LOGS/${sample}__fastp.log"

          matches=( ${INPATH}/*_${sample}_*R1*_001.fastq.gz )

          if [[ "${matches[0]}" == "${INPATH}/*_${sample}_*R1*_001.fastq.gz" ]]; then
            echo "ERROR: No FASTQ matched pattern: ${INPATH}/*_${sample}_*R1*_001.fastq.gz" > "$logfile"
            exit 1
          fi

          if (( ${#matches[@]} > 1 )); then
            echo "ERROR: Multiple FASTQs matched for ${sample}:" > "$logfile"
            printf "%s\n" "${matches[@]}" >> "$logfile"
            echo "Refine pattern or merge lanes before running." >> "$logfile"
            exit 1
          fi

          in_fq="${matches[0]}"

          fastp \
            -i "$in_fq" \
            -o "${sample}.ppm.fq" \
            -V \
            -D --dup_calc_accuracy $FASTP_DUP_ACCURACY \
            -g -x -q $FASTP_QUAL -e $FASTP_MIN_AVG_QUAL -l $FASTP_MIN_LEN -y -c -p \
            -h "${sample}.fastp.report.html" \
            -w 1 > "$logfile" 2>&1
        '
        check_success "Trimming single-end reads"

    else
        log_step "Trimming and merging reads with fastp..."
        cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "fastp \
          -i \"${INPATH}/{}\"*R1*.fastq.gz \
          -I \"${INPATH}/{}\"*R2*.fastq.gz \
          -m --merged_out '{}.ppm.fq' \
          -V --detect_adapter_for_pe \
          -D --dup_calc_accuracy $FASTP_DUP_ACCURACY \
          -g -x -q $FASTP_QUAL -e $FASTP_MIN_AVG_QUAL -l $FASTP_MIN_LEN -y -c -p \
          -h '{}.fastp.report.html' -w 1 > $LOGS/{}__fastp.log 2>&1"
        check_success "Trimming and merging reads"
    fi

    log_step "Removing duplicates with vsearch..."
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "vsearch \
      --fastx_uniques '{}.ppm.fq' \
      --fastqout '{}.ppm.vs.fq' \
      --minseqlength $VSEARCH_MIN_SEQLEN \
      --strand both > $LOGS/{}__vsearch.log 2>&1"
    check_success "Duplicate removal"

    log_step "Filtering low-complexity reads with SGA..."
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "sga preprocess --dust-threshold=$SGA_DUST_THRESHOLD -m $SGA_MIN_LEN '{}.ppm.vs.fq' -o '{}.ppm.vs.d4.fq' > $LOGS/{}__sga.log 2>&1"
    check_success "Low-complexity filtering"

    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "gzip '{}.ppm.vs.d4.fq'"
    check_success "Compressing filtered files"

else
    log_step "Skipping Stage 1 (preprocessing) as requested."
fi

if [ "$SKIP_PREPROCESSING_CLEANUP" = false ]; then
  log_step "Cleaning preprocessing intermediates where final *.ppm.vs.d4.fq.gz exists..."

  export LOGS

  parallel -j "$THREADSP" --env LOGS '
    sample={}
    logfile="$LOGS/${sample}__cleanup_preprocessing.log"
    final="${sample}.ppm.vs.d4.fq.gz"

    if [ ! -s "$final" ]; then
      echo "[SKIP] Final preprocessed file missing/empty: $final" >> "$logfile"
      exit 0
    fi

    candidates=""
    for f in \
      "${sample}.ppm.fq" \
      "${sample}.ppm.vs.fq" \
      "${sample}.ppm.vs.d4.fq" \
      "${sample}.ppm.fq.gz" \
      "${sample}.ppm.vs.fq.gz"
    do
      [ -e "$f" ] && candidates="$candidates $f"
    done

    if [ -z "$candidates" ]; then
      echo "[INFO] Nothing to delete for $sample" >> "$logfile"
      exit 0
    fi

    echo "[INFO] Deleting:$candidates" >> "$logfile"
    rm -f $candidates 2>> "$logfile"
    echo "[INFO] Cleanup done for $sample" >> "$logfile"
  ' :::: "$SAMPLE_LIST"

  check_success "Preprocessing cleanup"
else
  log_step "Skipping preprocessing cleanup (--skip-preprocessing-cleanup)."
fi

# ===========================================================
# STAGE 2 — GTDB mapping
# ===========================================================

if [ "$SKIP_GTDB_MAPPING" = false ]; then
  log_step "Starting mapping and taxonomic classification against GTDB..."

  for db in {1..7}; do
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "bowtie2 --threads $THREADS \
        -x $DB_PATH_bac.$db.fas.gz \
        -U {}.ppm.vs.d4.fq.gz \
        -k $BT2_GTDB_k -D $BT2_GTDB_D -R $BT2_GTDB_R -N $BT2_GTDB_N -L $BT2_GTDB_L \
        -i $BT2_GTDB_i --np $BT2_GTDB_np --mp '$BT2_GTDB_mp' \
        --rdg '$BT2_GTDB_rdg' --rfg '$BT2_GTDB_rfg' \
        --score-min '$BT2_GTDB_score_min' --mm --no-unal \
        2> $LOGS/{}.gtdb.$db.bowtie2.log \
        | samtools view -bS - > $MICROB_OUT/{}.gtdb.$db.bam"
    check_success "Mapping to GTDB chunk $db"
  done

  log_step "Sorting GTDB BAM files..."
  cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "for bam in $MICROB_OUT/{}.gtdb.*.bam; do
      sorted_bam=\$MICROB_OUT/\$(basename \"\$bam\" .bam).sorted.bam
      samtools sort -@ $THREADS -m 4G -o \"\$sorted_bam\" \"\$bam\"
    done"
  check_success "Sorting GTDB BAM files"

  log_step "Merging GTDB BAM files..."
  cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
    "samtools merge -@ $THREADS -f $MICROB_OUT/{}.gtdb.merged.bam $MICROB_OUT/{}.gtdb.*.sorted.bam"
  check_success "Merging GTDB BAM files"
else
  log_step "Skipping Stage 2 (GTDB mapping) as requested."
fi

if [ "$SKIP_GTDB_CLEANUP" = false ]; then
    log_step "Cleaning GTDB chunk BAMs now that merged BAM exists..."

    export MICROB_OUT LOGS

    parallel -j "$THREADSP" --env MICROB_OUT --env LOGS '
      sample={};
      logfile="$LOGS/${sample}__cleanup_gtdb_bams.log";
      merged="$MICROB_OUT/${sample}.gtdb.merged.bam";

      if [ ! -s "$merged" ]; then
        echo "[SKIP] Missing/empty merged BAM: $merged" > "$logfile";
        exit 0
      fi

      candidates=$(find "$MICROB_OUT" -maxdepth 1 -type f \
        \( -name "${sample}.gtdb.[1-7].bam" \
           -o -name "${sample}.gtdb.[1-7].sorted.bam" \
        \) \
        ! -name "${sample}.gtdb.merged.bam" \
      )

      if [ -z "$candidates" ]; then
        echo "[INFO] No GTDB chunk BAMs found to delete for $sample" > "$logfile"
        exit 0
      fi

      echo "$candidates" | wc -l | awk "{print \"[INFO] Deleting \" \$1 \" files for $sample\"}" > "$logfile"
      echo "$candidates" >> "$logfile"
      echo "$candidates" | xargs -r rm -f
      echo "[INFO] Cleanup done for $sample" >> "$logfile"
    ' :::: "$SAMPLE_LIST"

    check_success "GTDB chunk BAM cleanup"
else
    log_step "Skipping GTDB BAM cleanup (--skip-gtdb-cleanup)."
fi

# ===========================================================
# STAGE 3 — Microbial / eukaryotic read splitting
# ===========================================================

if [ "$SKIP_MICROBIAL_SPLIT" = false ]; then

    if [ "$LCA_ASSIGN" = true ]; then

        log_step "Sorting merged BAM file for metaDMG..."
        parallel -j "$THREADSP" "\
            samtools sort -n -@ $THREADS -m 10G \
                -o $MICROB_OUT/{}.gtdb.merged.sorted.bam \
                   $MICROB_OUT/{}.gtdb.merged.bam \
            > $LOGS/{}__sortbam.log 2>&1" \
            :::: "$SAMPLE_LIST"
        check_success "Sorting BAM file"

        log_step "Running metaDMG LCA assignment..."
        parallel -j "$THREADSP" "\
            /projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp lca \
                --names $TAX_PATH_BAC/names.dmp \
                --nodes $TAX_PATH_BAC/nodes.dmp \
                --acc2tax $TAX_PATH_BAC_ACC/hires-organelles-viruses-smags.acc2taxid.gz \
                --sim_score_low $METADMG_BAC_SIM_SCORE_LOW \
                --sim_score_high $METADMG_BAC_SIM_SCORE_HIGH \
                --how_many $METADMG_BAC_HOW_MANY \
                --weight_type $METADMG_BAC_WEIGHT_TYPE \
                --fix_ncbi 0 \
                --threads 10 \
                --bam $MICROB_OUT/{}.gtdb.merged.sorted.bam \
                --out_prefix $MICROB_OUT/{} \
            > $LOGS/{}__metadmg.log 2>&1" \
            :::: "$SAMPLE_LIST"
        check_success "metaDMG LCA assignment"

        log_step "Extracting bacterial/archaeal/viral reads..."
        parallel -j "$THREADSP" "\
            zgrep -i -E 'Archaea|virus|bacteria' \
                $MICROB_OUT/{}.lca.gz \
                | cut -f1 \
                > $MICROB_OUT/{}.bact_reads.txt \
            2> $LOGS/{}__extract_bact.log" \
            :::: "$SAMPLE_LIST"
        check_success "Extracting reads"

        log_step "Subsetting FASTQ for bacterial reads..."
        parallel -j "$THREADSP" "\
            seqtk subseq {}.ppm.vs.d4.fq.gz \
                $MICROB_OUT/{}.bact_reads.txt \
                | gzip > $MICROB_OUT/{}.bact_reads.fq.gz \
            2> $LOGS/{}__seqtk_bact.log" \
            :::: "$SAMPLE_LIST"
        check_success "Subsetting bacterial reads"

        log_step "Extracting all read IDs..."
        parallel -j "$THREADSP" "\
            zcat {}.ppm.vs.d4.fq.gz \
                | awk 'NR%4==1 {print substr(\$0, 2)}' \
                > $MICROB_OUT/{}.all_reads.txt \
            2> $LOGS/{}__allreads.log" \
            :::: "$SAMPLE_LIST"
        check_success "Extracting all read IDs"

        log_step "Computing eukaryotic read set..."
        parallel --shell /bin/bash -j "$THREADSP" "\
            comm -23 <(sort $MICROB_OUT/{}.all_reads.txt) \
                     <(sort $MICROB_OUT/{}.bact_reads.txt) \
                > $EUK_OUT/{}.euk_reads.txt \
            2> $LOGS/{}__comm_euk.log" \
            :::: "$SAMPLE_LIST"
        check_success "Generating eukaryotic read list"

        log_step "Extracting eukaryotic reads..."
        parallel -j "$THREADSP" "\
            seqtk subseq {}.ppm.vs.d4.fq.gz \
                $EUK_OUT/{}.euk_reads.txt \
                | gzip > $EUK_OUT/{}.euk.fastq.gz \
            2> $LOGS/{}__seqtk_euk.log" \
            :::: "$SAMPLE_LIST"
        check_success "Extracting eukaryotic FASTQ"

        if [ "$STORAGE_FRIENDLY" = true ]; then
            log_step "Cleaning up intermediate files..."
            parallel -j "$THREADSP" "\
                rm -f \
                    $MICROB_OUT/{}.lca.gz \
                    $MICROB_OUT/{}.bact_reads.txt \
                    $MICROB_OUT/{}.all_reads.txt \
                    $EUK_OUT/{}.euk_reads.txt \
                2> $LOGS/{}__cleanup.log" \
                :::: "$SAMPLE_LIST"
            check_success "Cleanup intermediate files"
        fi

    else

        log_step "Running getRTax classification..."
        parallel -j "$THREADSP" "\
            getRTax \
                --bam $MICROB_OUT/{}.gtdb.merged.bam \
                -T $TAX_PATH_BAC/taxid.map \
                -r '{\"domain\":[\"d__Bacteria\", \"d__Archaea\", \"d__Viruses\"]}' \
                --threads 8 \
                --unique \
                --only-read-ids \
                -p $MICROB_OUT/{}.bact_reads.txt \
            > $LOGS/{}__getRTax.log 2>&1" \
            :::: "$SAMPLE_LIST"
        check_success "getRTax classification"

        log_step "Merging bacterial read ID files..."
        parallel -j "$THREADSP" "\
            zcat $MICROB_OUT/{}.bact_reads.txt* \
                > $MICROB_OUT/{}.bact_reads_all.txt \
            2> $LOGS/{}__merge_bact.log" \
            :::: "$SAMPLE_LIST"
        check_success "Merging bacterial ID files"

        log_step "Subsetting bacterial reads..."
        parallel -j "$THREADSP" "\
            seqtk subseq {}.ppm.vs.d4.fq.gz \
                $MICROB_OUT/{}.bact_reads_all.txt \
                | gzip > $MICROB_OUT/{}.bact_reads.fq.gz \
            2> $LOGS/{}__seqtk_bact.log" \
            :::: "$SAMPLE_LIST"
        check_success "Extracting bacterial FASTQ"

        log_step "Extracting all read IDs..."
        parallel -j "$THREADSP" "\
            zcat {}.ppm.vs.d4.fq.gz \
                | awk 'NR%4==1 {split(substr(\$0, 2), a, \" \"); print a[1]}' \
                > $MICROB_OUT/{}.all_reads.txt \
            2> $LOGS/{}__allreads.log" \
            :::: "$SAMPLE_LIST"
        check_success "Extracting all read IDs"

        log_step "Sorting read lists..."
        parallel -j "$THREADSP" "sort $MICROB_OUT/{}.all_reads.txt > $MICROB_OUT/{}.all_reads.sorted.txt" :::: "$SAMPLE_LIST"
        parallel -j "$THREADSP" "sort $MICROB_OUT/{}.bact_reads_all.txt > $MICROB_OUT/{}.bact_reads_all.sorted.txt" :::: "$SAMPLE_LIST"
        check_success "Sorting"

        log_step "Computing euk read list..."
        parallel -j "$THREADSP" "\
            comm -23 \
                $MICROB_OUT/{}.all_reads.sorted.txt \
                $MICROB_OUT/{}.bact_reads_all.sorted.txt \
                > $EUK_OUT/{}.euk_reads.txt \
            2> $LOGS/{}__comm_euk.log" \
            :::: "$SAMPLE_LIST"
        check_success "Eukaryotic read list"

        log_step "Extracting euk FASTQ..."
        parallel -j "$THREADSP" "\
            seqtk subseq {}.ppm.vs.d4.fq.gz \
                $EUK_OUT/{}.euk_reads.txt \
                | gzip > $EUK_OUT/{}.euk.fastq.gz \
            2> $LOGS/{}__seqtk_euk.log" \
            :::: "$SAMPLE_LIST"
        check_success "Extracting eukaryotic FASTQ"

        if [ "$STORAGE_FRIENDLY" = true ]; then
            log_step "Cleaning up intermediate files..."
            parallel -j "$THREADSP" "\
                rm -f \
                    $MICROB_OUT/{}.bact_reads.txt* \
                    $MICROB_OUT/{}.bact_reads_all.txt \
                    $MICROB_OUT/{}.all_reads.txt \
                    $MICROB_OUT/{}.all_reads.sorted.txt \
                    $MICROB_OUT/{}.bact_reads_all.sorted.txt \
                    $EUK_OUT/{}.euk_reads.txt \
                2> $LOGS/{}__cleanup.log" \
                :::: "$SAMPLE_LIST"
            check_success "Cleanup intermediate files"
        fi

    fi
else
  log_step "Skipping Stage 3 (microbial read splitting) as requested."
fi

# ===========================================================
# STAGE 4 — Eukaryotic mapping
# ===========================================================

if [ "$SKIP_EUK_MAPPING" = false ]; then

    log_step "Mapping reads to eukaryote database (129 parts) with bowtie2..."
    for db in {1..129}; do
        parallel -j "$THREADSP" "\
            bowtie2 --threads $THREADS -k $BT2_EUK_k -t \
                -x $DB_PATH.$db.fas.gz \
                -U $EUK_OUT/{}.euk.fastq.gz \
                --no-unal --mm -t 2> $LOGS/{}__eukmap_part_${db}.log \
            | samtools view -bS - \
                > $EUK_OUT/{}.euk.$db.bam" \
            :::: "$SAMPLE_LIST"
        check_success "Mapping to eukaryote database part $db"
    done

    log_step "Mapping reads to mitochondrion database with bowtie2..."
    parallel -j "$THREADSP" "\
        bowtie2 --threads $THREADS -k $BT2_EUK_k -t \
            -x $DB_PATH_clean/refseq_mitochondrion.genomic.fas.gz \
            -U $EUK_OUT/{}.euk.fastq.gz \
            --no-unal --mm -t 2> $LOGS/{}__mitochondrion.log \
        | samtools view -bS - \
            > $EUK_OUT/{}.mito.bam" \
        :::: "$SAMPLE_LIST"
    check_success "Mapping to mitochondrion database"

    log_step "Mapping reads to phylonorwary database (10 parts) with bowtie2..."
    for db in {1..10}; do
        parallel -j "$THREADSP" "\
            bowtie2 --threads $THREADS -k $BT2_EUK_k -t \
                -x $DB_PATH_Norwary.$db-of-10 \
                -U $EUK_OUT/{}.euk.fastq.gz \
                --no-unal --mm -t 2> $LOGS/{}__phynor_part_${db}.log \
            | samtools view -bS - \
                > $EUK_OUT/{}.phyNor.$db.bam" \
            :::: "$SAMPLE_LIST"
        check_success "Mapping to phylonorwary database part $db"
    done

    log_step "Mapping reads to core NT database with bowtie2..."
    parallel -j "$THREADSP" "\
        bowtie2 --threads $THREADS -k $BT2_EUK_k -t \
            -x $DB_PATH_clean/core_nt.fas.gz \
            -U $EUK_OUT/{}.euk.fastq.gz \
            --no-unal --mm -t 2> $LOGS/{}__core_nt.log \
        | samtools view -bS - \
            > $EUK_OUT/{}.core_nt.bam" \
        :::: "$SAMPLE_LIST"
    check_success "Mapping to core NT database"

    log_step "Mapping reads to plastid database with bowtie2..."
    parallel -j "$THREADSP" "\
        bowtie2 --threads $THREADS -k $BT2_EUK_k -t \
            -x $DB_PATH_clean/refseq_plastid.genomic.fas.gz \
            -U $EUK_OUT/{}.euk.fastq.gz \
            --no-unal --mm -t 2> $LOGS/{}__plastid.log \
        | samtools view -bS - \
            > $EUK_OUT/{}.pla.bam" \
        :::: "$SAMPLE_LIST"
    check_success "Mapping to plastid database"

    log_step "Mapping finished. Continuing with merging..."
else
  log_step "Skipping Stage 4 (eukaryotic mapping) as requested."
fi

# ===========================================================
# STAGE 5 — Compress, sort, merge BAMs
# ===========================================================

if [ "$SKIP_COMP_MERGE" = false ]; then
    log_step "Compressing BAM files using metaDMG..."

    export EUK_OUT THREADS THREADSP LOGS

    parallel -j "$THREADSP" --env EUK_OUT --env LOGS '
      sample={};

      find "$EUK_OUT" -type f -name "${sample}*.bam" \
        | grep -E "${sample}\.(euk|mito|pla|phyNor|core_nt)(\.[0-9]+)?\.bam$" \
        | grep -vE "sorted|comp|merged" \
        | while read -r bam; do

            base=$(basename "$bam" .bam)
            dir=$(dirname "$bam")
            outpath="$dir/${base}.comp.bam"

            /projects/wintherpedersen/apps/metaDMG_14jun24/metaDMG-cpp/misc/compressbam \
              --threads 12 \
              --input "$bam" \
              --output "$outpath" \
              > "$LOGS/${sample}__compressbam.log" 2>&1

          done
    ' :::: "$SAMPLE_LIST"
    check_success "Compressing BAM files"

    log_step "Sorting each BAM file before merging..."

    parallel -j "$THREADSP" --env EUK_OUT --env THREADS --env LOGS '
      sample={};

      find "$EUK_OUT" -type f -name "${sample}*.comp.bam" \
        | while read -r bam; do

            base=$(basename "$bam" .bam)
            dir=$(dirname "$bam")
            sorted_bam="$dir/${base}.sorted.bam"

            samtools sort -n -@ "$THREADS" -m 4G -o "$sorted_bam" "$bam" \
              > "$LOGS/${sample}__sort_comp.log" 2>&1

          done
    ' :::: "$SAMPLE_LIST"
    check_success "BAM files sorted"

    log_step "Merging all sorted BAM files..."
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
      "samtools merge -@ $THREADS -n -f \
         $EUK_OUT/{}.comp.sam.gz $EUK_OUT/{}*.comp.sorted.bam"
    check_success "Merging BAM files to sam.gz"

    log_step "Compressing merged SAM..."
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
      "/projects/wintherpedersen/apps/metaDMG_14jun24/metaDMG-cpp/misc/compressbam \
         --threads 12 \
         --input $EUK_OUT/{}.comp.sam.gz \
         --output $EUK_OUT/{}.comp.bam \
         > $LOGS/{}__compress_merged.log 2>&1"
    check_success "Merged sam.gz compressed to comp.bam"
else
  log_step "Skipping Stage 5 (compress/sort/merge). Expecting $EUK_OUT/<sample>.comp.bam to already exist."
fi

log_step "Checking that merged comp.bam exists for all samples..."
missing_comp=0
while read -r s; do
  if [ ! -s "$EUK_OUT/${s}.comp.bam" ]; then
    log_error "Missing or empty: $EUK_OUT/${s}.comp.bam (needed for downstream steps)."
    missing_comp=1
  fi
done < "$SAMPLE_LIST"

if [ "$missing_comp" -ne 0 ]; then
  log_error "One or more .comp.bam files missing. Either run without --skip-comp-merge or provide the files."
  exit 1
fi

if [ "$SKIP_EUK_CLEANUP" = false ]; then
  log_step "Cleaning up per-database mapping BAMs now that merged *.comp.bam exists..."

  export EUK_OUT LOGS

  parallel -j "$THREADSP" --env EUK_OUT --env LOGS '
    sample={};
    logfile="$LOGS/${sample}__cleanup_mapping_bams.log";
    merged="$EUK_OUT/${sample}.comp.bam";

    if [ ! -s "$merged" ]; then
      echo "[SKIP] merged comp.bam missing: $merged" > "$logfile";
      exit 0
    fi

    candidates=$(find "$EUK_OUT" -maxdepth 1 -type f \
      \( -name "${sample}.euk.*.bam" \
         -o -name "${sample}.core_nt.bam" \
         -o -name "${sample}.phyNor.*.bam" \
         -o -name "${sample}.mito.bam" \
         -o -name "${sample}.pla.bam" \
         -o -name "${sample}.euk.*.comp.bam" \
         -o -name "${sample}.core_nt.comp.bam" \
         -o -name "${sample}.phyNor.*.comp.bam" \
         -o -name "${sample}.mito.comp.bam" \
         -o -name "${sample}.pla.comp.bam" \
         -o -name "${sample}.euk.*.comp.sorted.bam" \
         -o -name "${sample}.core_nt.comp.sorted.bam" \
         -o -name "${sample}.phyNor.*.comp.sorted.bam" \
         -o -name "${sample}.mito.comp.sorted.bam" \
         -o -name "${sample}.pla.comp.sorted.bam" \
         -o -name "${sample}.comp.sam.gz" \
      \) \
      ! -name "${sample}.comp.bam" \
      ! -name "${sample}.comp.filtered*.bam" \
      ! -name "${sample}.sort.comp*.bam" \
    )

    if [ -z "$candidates" ]; then
      echo "[INFO] No mapping intermediates found to delete for $sample" > "$logfile"
      exit 0
    fi

    echo "$candidates" | wc -l | awk "{print \"[INFO] Deleting \" \$1 \" files for $sample\"}" > "$logfile"
    echo "$candidates" >> "$logfile"
    echo "$candidates" | xargs -r rm -f
    echo "[INFO] Cleanup done for $sample" >> "$logfile"
  ' :::: "$SAMPLE_LIST"

  check_success "Cleanup of mapping intermediates"
else
  log_step "Skipping mapping-intermediate cleanup (--skip-euk-cleanup)."
fi

# ===========================================================
# STAGE 6 — BAM filtering (Unicorn: alnfilt → refstats → bamstats)
# ===========================================================

if [ "$SKIP_BAM_FILTERING" = false ]; then

    log_step "Running unicorn alnfilt..."
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
      "/projects/wintherpedersen/apps/unicorn/unicorn alnfilt \
         -b $EUK_OUT/{}.comp.bam \
         -t $THREADS --mode $ALNFILT_MODE \
         --outbam $EUK_OUT/{}.alnfilt.bam \
         --outstat $EUK_OUT/{}.alnfilt.refstats \
         --minani $ALNFILT_MINANI --maxani $ALNFILT_MAXANI \
         > $LOGS/{}__unicorn_alnfilt.log 2>&1"
    check_success "Unicorn alnfilt"

    log_step "Running unicorn refstats..."

    export EUK_OUT LOGS TAX_PATH_NCBI THREADS

    cat "$SAMPLE_LIST" | parallel -j 1 '
      sample={};
      inbam="$EUK_OUT/${sample}.alnfilt.bam";
      outbam="$EUK_OUT/${sample}.comp.unicorn.bam";
      outstat="$EUK_OUT/${sample}.comp.unicorn.refstats";
      logfile="$LOGS/${sample}__unicorn_refstats.log";

      if [[ -s "$outbam" && -s "$outstat" ]]; then
        echo "[SKIP] $sample: outputs exist and are non-empty" > "$logfile"
      else
        /projects/wintherpedersen/apps/unicorn/unicorn refstats \
          -b "$inbam" \
          -t "$THREADS" --minreads $REFSTATS_MINREADS \
          --outbam "$outbam" \
          --outstat "$outstat" \
          --names "$TAX_PATH_NCBI/taxdump/names.dmp" \
          --nodes "$TAX_PATH_NCBI/taxdump/nodes.dmp" \
          > "$logfile" 2>&1
      fi
    '
    check_success "Unicorn refstats"

    log_step "Running unicorn bamstats..."
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
      "/projects/wintherpedersen/apps/unicorn/unicorn bamstats \
         -b $EUK_OUT/{}.comp.unicorn.bam \
         -t $THREADS \
         --outbam $EUK_OUT/{}.comp.filtered.bam \
         --outstat $EUK_OUT/{}.comp.filtered.unicorn.bamstats \
         --printdists $EUK_OUT/{}.comp.filtered.unicorn \
         > $LOGS/{}__unicorn_bamstats.log 2>&1"
    check_success "Unicorn bamstats"

    if [ "$STORAGE_FRIENDLY" = true ]; then
        log_step "Cleaning up intermediate Unicorn BAM files..."
        cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
          "rm -f \
             $EUK_OUT/{}.comp.bam \
             $EUK_OUT/{}.alnfilt.bam \
             $EUK_OUT/{}.comp.unicorn.bam \
             > /dev/null 2> $LOGS/{}__cleanup_unicorn.log"
        check_success "Cleanup intermediate Unicorn BAM files"
    fi

else
    log_step "Skipping Stage 6 (BAM filtering) as requested."
fi

# ===========================================================
# STAGE 7 — metaDMG (LCA, dfit, aggregate)
# ===========================================================

if [ "$SKIP_METADMG" = false ]; then

    log_step "Sorting merged BAM file for metaDMG..."
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
      "samtools sort -n -@ $THREADS -m 10G \
         -o $EUK_OUT/{}.sort.comp.filtered.bam \
            $EUK_OUT/{}.comp.filtered.bam"
    check_success "Sorting BAM file"

    log_step "Running taxonomic classification with metaDMG (LCA)..."
    cat "$SAMPLE_LIST" | parallel --shell /bin/bash -j "$THREADSP" \
      "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp lca \
         --names $TAX_PATH_NCBI/taxdump/names.dmp \
         --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp \
         --acc2tax <(zcat $TAX_PATH_NCBI/*.acc2taxid.gz /datasets/caeg_dataset/references/phylo_norway/20250127/*.acc2taxid.gz) \
         --sim_score_low $METADMG_EUK_SIM_SCORE_LOW \
         --sim_score_high $METADMG_EUK_SIM_SCORE_HIGH \
         --how_many $METADMG_EUK_HOW_MANY \
         --weight_type $METADMG_EUK_WEIGHT_TYPE \
         --fix_ncbi 0 --threads 10 --filtered_acc2tax $EUK_OUT/{}.acc2tax \
         --bam $EUK_OUT/{}.sort.comp.filtered.bam \
         --out_prefix $EUK_OUT/{}.sort.comp.filtered \
         > $LOGS/{}__metadmg_lca.log 2>&1"
    check_success "Taxonomic classification (LCA)"

    log_step "Running damage estimation with metaDMG (dfit)..."
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
      "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp dfit \
         $EUK_OUT/{}.sort.comp.filtered.bdamage.gz --threads 6 \
         --names $TAX_PATH_NCBI/taxdump/names.dmp \
         --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp \
         --showfits 2 --nopt $METADMG_DFIT_NOPT \
         --nbootstrap $METADMG_DFIT_NBOOTSTRAP \
         --doboot 1 --seed $METADMG_DFIT_SEED --lib $METADMG_DFIT_LIB \
         --out_prefix $EUK_OUT/{}.sort.comp.filtered \
         > $LOGS/{}__metadmg_dfit.log 2>&1"
    check_success "Damage estimation (dfit)"

    log_step "Aggregating metaDMG LCA and dfit results..."
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
      "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp aggregate \
         $EUK_OUT/{}.sort.comp.filtered.bdamage.gz \
         --names $TAX_PATH_NCBI/taxdump/names.dmp \
         --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp \
         --lcastat $EUK_OUT/{}.sort.comp.filtered.stat.gz \
         --dfit $EUK_OUT/{}.sort.comp.filtered.dfit.gz \
         --out_prefix $EUK_OUT/{}.sort.comp.filtered.agg \
         > $LOGS/{}__metadmg_agg.log 2>&1"
    check_success "metaDMG aggregation"
else
  log_step "Skipping Stage 7 (metaDMG) as requested."
fi

# ===========================================================
# STAGE 8 — Unicorn taxstats (family, genus, species)
# ===========================================================

if [ "$SKIP_UNICORN_TIDSTATS" = false ]; then

    for rank in family genus species; do
        log_step "Running unicorn taxstats at rank: ${rank}..."
        cat "$SAMPLE_LIST" | parallel --shell /bin/bash -j "$THREADSP" \
          "/projects/wintherpedersen/apps/unicorn/unicorn taxstats \
             -b $EUK_OUT/{}.comp.filtered.bam \
             -t $THREADS \
             -o $EUK_OUT/{}.comp.filtered.${rank}.taxstats.bam \
             --names $TAX_PATH_NCBI/taxdump/names.dmp \
             --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp \
             --acc2tax <(zcat $TAX_PATH_NCBI/*.acc2taxid.gz /datasets/caeg_dataset/references/phylo_norway/20250127/*.acc2taxid.gz) \
             -k $TAXSTATS_k --outstat $EUK_OUT/{}.comp.filtered.${rank}.taxstats \
             --minreads $TAXSTATS_MINREADS --rank ${rank} \
             > $LOGS/{}__unicorn_taxstats_${rank}.log 2>&1"
        check_success "Unicorn taxstats (${rank})"
    done

else
  log_step "Skipping Stage 8 (unicorn taxstats) as requested."
fi

echo "Pipeline completed successfully." | tee -a "$LOG_FILE"