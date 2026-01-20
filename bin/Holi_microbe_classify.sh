#!/bin/bash
# Script: Holi_microbe_classify.sh
# Description: Microbial workflow with two modes:
#   - reference: filterBAM filter (name-sorted output) + metaDMG getdamage + metaDMG print_ugly (--bam)
#   - both: reference steps + metaDMG lca + dfit + aggregate + unicorn tidstats
#
# Naming convention (per sample):
#   - reference products:  <sample>.ref.*
#   - LCA products:        <sample>.lca.*
#
# Usage:
#   ./Holi_microbe_classify.sh <config.yml> <sample.list> --mode reference
#   ./Holi_microbe_classify.sh <config.yml> <sample.list> --mode both

set -euo pipefail

usage() {
    echo "Usage: $0 <config.yml> <sample.list> --mode <reference|both>"
    echo
    echo "Arguments:"
    echo "  <config.yml>            Path to YAML configuration file."
    echo "  <sample.list>           File containing list of sample names."
    echo
    echo "Options:"
    echo "  --mode reference        Run filterBAM filter + metaDMG getdamage + print_ugly."
    echo "  --mode both             Run reference steps then LCA workflow."
    echo
    exit 1
}

if [ $# -lt 4 ]; then
    usage
fi

CONFIG="$1"
SAMPLE_LIST="$2"
shift 2

MODE=""

while [ $# -gt 0 ]; do
  case "$1" in
    --mode)
      MODE="${2:-}"; shift 2;;
    -h|--help)
      usage;;
    *)
      echo "Unknown option: $1"
      usage;;
  esac
done

if [ -z "$MODE" ]; then
  echo "Missing --mode"
  usage
fi

if [ "$MODE" != "reference" ] && [ "$MODE" != "both" ]; then
  echo "Invalid --mode: $MODE"
  usage
fi

if [ ! -f "$SAMPLE_LIST" ]; then
  echo "Sample list not found: $SAMPLE_LIST"
  exit 1
fi

# -------------------------------
# Tool paths (pin to your install)
# -------------------------------
META_DMG_BIN="metaDMG-cpp"
UNICORN_BIN="unicorn"
BAMFILTER_BIN="filterBAM"

# -------------------------------
# Load config file (simple YAML KEY: VALUE)
# -------------------------------
load_config() {
    CONFIG_FILE="$1"
    get_value() {
        grep "^$1:" "$CONFIG_FILE" \
          | sed 's/^.*:[[:space:]]*//' \
          | sed 's/"//g' \
          | sed 's/^[[:space:]]*//; s/[[:space:]]*$//'
    }

    THREADSP=$(get_value "THREADSP")
    THREADS=$(get_value "THREADS")

    TAX_PATH_BAC=$(get_value "TAX_PATH_BAC")
    TAX_PATH_BAC_ACC=$(get_value "TAX_PATH_BAC_ACC")

    MICROB_OUT=$(get_value "MICROB_OUT")
    RESULT_PATH=$(get_value "RESULT_PATH")
    LOGS=$(get_value "LOGS")
}

load_config "$CONFIG"
echo "[INFO] Loaded config from $CONFIG"

# -----------------------------
# Logging functions
# -----------------------------
GREEN="\033[0;32m"
YELLOW="\033[1;33m"
RED="\033[0;31m"
RESET="\033[0m"

timestamp() { date +"%Y-%m-%d %H:%M:%S"; }
log_step() { printf "${YELLOW}[%s][STEP]${RESET} %s\n" "$(timestamp)" "$1"; }
log_info() { printf "${GREEN}[%s][INFO]${RESET} %s\n" "$(timestamp)" "$1"; }
log_error(){ printf "${RED}[%s][ERROR]${RESET} %s\n" "$(timestamp)" "$1"; }

check_success() {
    if [ $? -eq 0 ]; then
        log_info "$1 completed"
    else
        log_error "$1 failed. Stopping."
        exit 1
    fi
}

tool_exists() {
  local t="$1"
  if [[ "$t" == */* ]]; then
    [ -x "$t" ]
  else
    command -v "$t" >/dev/null 2>&1
  fi
}

# -----------------------------
# Ensure essential output directories exist
# -----------------------------
REQUIRED_DIRS=("$RESULT_PATH" "$MICROB_OUT" "$LOGS")
for dir in "${REQUIRED_DIRS[@]}"; do
    if [ ! -d "$dir" ]; then
        log_info "Directory $dir does not exist. Creating it..."
        mkdir -p "$dir"
    else
        log_info "Directory $dir exists."
    fi
done

# -------------------------------
# File naming
# -------------------------------
BAM_SUFFIX=".gtdb.merged.bam"
FILTERED_BAM_SUFFIX=".filtered.qname.bam"

STATS_ALL_SUFFIX=".stats-all.tsv.gz"
STATS_FILT_SUFFIX=".stats-filtered.tsv.gz"

# Prefix tags
REF_TAG=".ref"
LCA_TAG=".lca"

# Common suffixes
BDAMAGE_SUFFIX=".bdamage.gz"
STAT_SUFFIX=".stat.gz"
DFIT_SUFFIX=".dfit.gz"

# -------------------------------
# Required tools check
# -------------------------------
REQUIRED_TOOLS=( "parallel" "$META_DMG_BIN" "$BAMFILTER_BIN" )
if [ "$MODE" = "both" ]; then
  REQUIRED_TOOLS+=( "$UNICORN_BIN" )
fi

echo "Checking required tools..."
missing_tools=()
for tool in "${REQUIRED_TOOLS[@]}"; do
  if ! tool_exists "$tool"; then
    missing_tools+=("$tool")
  fi
done

if [ ${#missing_tools[@]} -ne 0 ]; then
  echo "The following required tools are missing or not executable:"
  for t in "${missing_tools[@]}"; do
    echo "  - $t"
  done
  exit 1
else
  echo "All required tools found."
fi

# -------------------------------
# Steps
# -------------------------------
run_filter_and_getdamage() {
  GENOMES_LEN_MAP="/datasets/caeg_dataset/references/hires-organelles-viruses-smags/20240313/misc/genomes.len.map"

  log_step "Running filterBAM filter (output queryname sorted via -N)..."
  parallel -j "$THREADSP" "\
      INBAM=\"$MICROB_OUT/{}${BAM_SUFFIX}\"; \
      OUTBAM=\"$MICROB_OUT/{}${FILTERED_BAM_SUFFIX}\"; \
      STATS_ALL=\"$MICROB_OUT/{}${STATS_ALL_SUFFIX}\"; \
      STATS_FILT=\"$MICROB_OUT/{}${STATS_FILT_SUFFIX}\"; \
      $BAMFILTER_BIN filter \
        --bam \$INBAM \
        --stats \$STATS_ALL \
        --stats-filtered \$STATS_FILT \
        --bam-filtered \$OUTBAM \
        -r $GENOMES_LEN_MAP \
        -m 5G \
        -t $THREADS \
        -N \
        -n 3 \
        -A 92 \
        -a 94 \
        -l 30 \
      > \"$LOGS/{}__filterBAM.log\" 2>&1" \
      :::: "$SAMPLE_LIST"
  check_success "filterBAM filter"

  log_step "Running metaDMG getdamage (per reference) + print_ugly (per reference name via --bam)..."
  parallel -j "$THREADSP" "\
      SAMPLE=\"{}\"; \
      QBAM=\"$MICROB_OUT/\${SAMPLE}${FILTERED_BAM_SUFFIX}\"; \
      REF_PREFIX=\"$MICROB_OUT/\${SAMPLE}${REF_TAG}\"; \
      REF_BDAMAGE=\"\${REF_PREFIX}${BDAMAGE_SUFFIX}\"; \
      GETDAMAGE_STDOUT=\"\${REF_PREFIX}.getdamage.stdout.tsv\"; \
      \
      # 1) getdamage: write stdout to a parseable file, stderr to log \
      $META_DMG_BIN getdamage \
        -r 1 \
        -p 15 \
        -l 30 \
        -n $THREADS \
        -o \"\$REF_PREFIX\" \
        \"\$QBAM\" \
      > \"\$GETDAMAGE_STDOUT\" 2> \"$LOGS/\${SAMPLE}__metadmg_getdamage.ref.log\"; \
      \
      gzip -f \"\$GETDAMAGE_STDOUT\"; \
      \
      if [ ! -s \"\$REF_BDAMAGE\" ]; then \
        echo \"[ERROR] Missing bdamage output: \$REF_BDAMAGE\" >&2; \
        exit 1; \
      fi; \
      \
      # 2) print_ugly: produces .ref.uglyprint.* (mismatch + stat) \
      $META_DMG_BIN print_ugly \"\$REF_BDAMAGE\" \
        --bam \"\$QBAM\" \
        --out_prefix \"\$REF_PREFIX\" \
      > \"$LOGS/\${SAMPLE}__metadmg_print_ugly.ref.log\" 2>&1" \
      :::: "$SAMPLE_LIST"
  check_success "metaDMG getdamage + print_ugly"
}

run_lca_workflow() {
  log_step "Running metaDMG LCA assignment (on filtered, queryname sorted BAM)..."
  parallel -j "$THREADSP" "\
      SAMPLE=\"{}\"; \
      QBAM=\"$MICROB_OUT/\${SAMPLE}${FILTERED_BAM_SUFFIX}\"; \
      LCA_PREFIX=\"$MICROB_OUT/\${SAMPLE}${LCA_TAG}\"; \
      \
      $META_DMG_BIN lca \
          --names $TAX_PATH_BAC/names.dmp \
          --nodes $TAX_PATH_BAC/nodes.dmp \
          --acc2tax $TAX_PATH_BAC_ACC/hires-organelles-viruses-smags.acc2taxid.gz \
          --sim_score_low 0.95 \
          --sim_score_high 1.0 \
          --how_many 15 \
          --weight_type 1 \
          --fix_ncbi 0 \
          --threads $THREADS \
          --lca_rank species \
          --filtered_acc2tax $MICROB_OUT/\${SAMPLE}.lca.acc2tax.gz \
          --bam \"\$QBAM\" \
          --out_prefix \"\$LCA_PREFIX\" \
      > \"$LOGS/\${SAMPLE}__metadmg_lca.lca.log\" 2>&1" \
      :::: "$SAMPLE_LIST"
  check_success "metaDMG lca"

  log_step "Running damage estimation with metaDMG (dfit) on LCA bdamage..."
  cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
    "SAMPLE={}; \
     LCA_PREFIX=\"$MICROB_OUT/\${SAMPLE}${LCA_TAG}\"; \
     $META_DMG_BIN dfit \
       \"\${LCA_PREFIX}${BDAMAGE_SUFFIX}\" --threads $THREADS \
       --names $TAX_PATH_BAC/names.dmp \
       --nodes $TAX_PATH_BAC/nodes.dmp \
       --showfits 2 --nopt 10 \
       --nbootstrap 20 --doboot 1 --seed 1234 --lib ds \
       --out_prefix \"\$LCA_PREFIX\" \
       > \"$LOGS/\${SAMPLE}__metadmg_dfit.lca.log\" 2>&1"
  check_success "metaDMG dfit"

  log_step "Aggregating LCA + dfit outputs with metaDMG (aggregate)..."
  cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
    "SAMPLE={}; \
     LCA_PREFIX=\"$MICROB_OUT/\${SAMPLE}${LCA_TAG}\"; \
     $META_DMG_BIN aggregate \
       \"\${LCA_PREFIX}${BDAMAGE_SUFFIX}\" \
       --names $TAX_PATH_BAC/names.dmp \
       --nodes $TAX_PATH_BAC/nodes.dmp \
       --lcastat \"\${LCA_PREFIX}${STAT_SUFFIX}\" \
       --dfit \"\${LCA_PREFIX}${DFIT_SUFFIX}\" \
       --out_prefix \"$MICROB_OUT/\${SAMPLE}.lca.agg\" \
       > \"$LOGS/\${SAMPLE}__metadmg_agg.lca.log\" 2>&1"
  check_success "metaDMG aggregate"

  log_step "Unicorn per taxID statistics (on filtered, queryname sorted BAM)..."
  cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
    "$UNICORN_BIN tidstats \
       -b $MICROB_OUT/{}${FILTERED_BAM_SUFFIX} \
       -t $THREADS \
       -o $MICROB_OUT/{}.lca.unicorn.tidstats \
       --names $TAX_PATH_BAC/names.dmp \
       --nodes $TAX_PATH_BAC/nodes.dmp \
       --acc2tax $TAX_PATH_BAC_ACC/hires-organelles-viruses-smags.acc2taxid.gz \
       > $LOGS/{}__unicorn_tidstats.lca.log 2>&1"
  check_success "unicorn tidstats"
}

# -------------------------------
# Dispatch
# -------------------------------
if [ "$MODE" = "reference" ]; then
  run_filter_and_getdamage
  log_info "Reference mode finished."
  exit 0
fi

if [ "$MODE" = "both" ]; then
  run_filter_and_getdamage
  run_lca_workflow
  log_info "Reference + LCA mode finished."
  exit 0
fi
