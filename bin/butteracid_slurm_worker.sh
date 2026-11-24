#!/usr/bin/env bash
set -Eeuo pipefail
if (( $# != 2 )); then
  echo "Usage: $0 SETTINGS STAGE" >&2
  exit 2
fi
source "$1"
stage=$2
lib=GLOBAL
trap 'status=$?; printf "[%s] STAGE=%s LIBRARY_ID=%s JOB_ID=%s EXIT=%s\n" "$(date -Is)" "$stage" "$lib" "${SLURM_JOB_ID:-NA}" "$status" >> "$RUN/failures.log"; exit "$status"' ERR
if [[ "$stage" == prepare ]]; then
  printf 'prepare\t%s\n' "${SLURM_JOB_ID:?Missing Slurm job ID}" >> "$RUN/jobs.tsv"
fi
cd "$LAUNCH_DIR"
source "$CONFIG"
# Optional file that loads modules/activates an environment on every compute job.
[[ -z ${SLURM_ENV_SETUP:-} ]] || source "$SLURM_ENV_SETUP"
export STEP5_MODE=rerun
: "${OUT_PREFIX:=butteracid}"
: "${MERGED_BAM_DIR:=$WORKDIR/merged_bams}"
: "${FLOWCELL_SUMMARY:=$WORKDIR/${OUT_PREFIX}_flowcell_summary.tsv}"
: "${LIBRARY_BAM_LIST:=$WORKDIR/${OUT_PREFIX}_library_bams.tsv}"
: "${MERGED_QC_STATS_OUTPUT:=$WORKDIR/${OUT_PREFIX}_qc_merged_stats.tsv}"
: "${PARQUET_DIR:=$WORKDIR/${OUT_PREFIX}_parquet}"
: "${RSCRIPT_PARQUET:=$(dirname "$RSCRIPT_RERUN")/butteracid_parquet_export.R}"
: "${MISSING_DATA_LOG:=$WORKDIR/missing_data.log}"
: "${FAILED_JOBS_LOG:=$WORKDIR/failed_jobs.log}"
: "${LOG_DIR:=$WORKDIR/logs}"
: "${MERGE_SORT_MEM:=1G}"
UNICORN_BIN=${UNICORN_BIN:-unicorn}
METADMG_BIN=${METADMG_BIN:-metaDMG-cpp}
ONEBAM_BIN=${ONEBAM_BIN:-onebam}
PROCESSED_DIR="$LAUNCH_DIR/${OUT_PREFIX}_processed_data"
mkdir -p "$EUK_DIR" "$LOG_DIR" "$SORT_TMPDIR"

submit_job() {
  local name=$1 dep=$2 array=$3 cpus=$4 mem=$5 time=$6
  local args=(--parsable --ntasks=1 --cpus-per-task="$cpus" --mem="$mem" --time="$time"
              --chdir="$LAUNCH_DIR" --export=ALL --job-name="ba_${name}"
              --kill-on-invalid-dep=yes)
  [[ -z "$dep" ]] || args+=(--dependency="$dep")
  if [[ -n "$array" ]]; then
    args+=(--array="$array" --output="$RUN/logs/${name}_%A_%a.out" --error="$RUN/logs/${name}_%A_%a.err")
  else
    args+=(--output="$RUN/logs/${name}_%j.out" --error="$RUN/logs/${name}_%j.err")
  fi
  [[ -z ${SLURM_ACCOUNT:-} ]] || args+=(--account="$SLURM_ACCOUNT")
  [[ -z ${SLURM_PARTITION:-} ]] || args+=(--partition="$SLURM_PARTITION")
  [[ -z ${SLURM_QOS:-} ]] || args+=(--qos="$SLURM_QOS")
  [[ -z ${SLURM_RESERVATION:-} ]] || args+=(--reservation="$SLURM_RESERVATION")
  local result id
  result=$(sbatch "${args[@]}" "$RUN/worker.sh" "$RUN/settings.sh" "$name") || return 1
  id=${result%%;*}
  [[ "$id" =~ ^[0-9]+$ ]] || { echo "Invalid sbatch response: $result" >&2; return 1; }
  printf '%s\t%s\n' "$name" "$id" >> "$RUN/jobs.tsv"
  printf '%s\n' "$id"
}

if [[ "$stage" == prepare ]]; then
  # The original preparation runs on this allocation, including its GNU parallel jobs.
  export LOG_DIR FAILED_JOBS_LOG
  export SHELL=$BASH SHELLOPTS
  export JOBLOG="$RUN/prepare.parallel.joblog"
  before=0
  [[ ! -f "$FAILED_JOBS_LOG" ]] || before=$(wc -l < "$FAILED_JOBS_LOG")
  flags=(--rerun)
  [[ "$ONLY" != 1 ]] || flags+=(--only-step5)
  [[ "$REBUILD" != 1 ]] || flags+=(--rebuild-paths)
  bash "$RUN/prepare.sh" "$CONFIG" "${flags[@]}"
  after=0
  [[ ! -f "$FAILED_JOBS_LOG" ]] || after=$(wc -l < "$FAILED_JOBS_LOG")
  if (( after > before )); then
    echo "Preparation recorded failed parallel jobs; fix them before submitting downstream jobs." >&2
    exit 1
  fi
  [[ -s "$SAMPLE_LIST" && -s "$LIBRARY_BAM_LIST" ]]
  cap=${SLURM_MAX_CONCURRENT:-10}
  [[ "$cap" =~ ^[1-9][0-9]*$ ]] || { echo 'SLURM_MAX_CONCURRENT must be positive' >&2; exit 1; }
  cp "$LIBRARY_BAM_LIST" "$RUN/manifest.tsv"
  # Avoid concurrent tasks modifying the shared BAM list. It is updated at finalization.
  awk -F'\t' '
    NF != 4 || $1 !~ /^[[:alnum:]_][[:alnum:]_.-]*$/ || $2 == "" ||
    $3 !~ /^(single|merged|single_usable)$/ || $4 == "" || seen[$1]++ {
      print "Invalid manifest row " NR > "/dev/stderr"
      bad = 1
    }
    END { exit bad }
  ' "$RUN/manifest.tsv"
  # Materialize acc2tax once. Unlike process substitution, a failed zcat is fatal.
  shopt -s nullglob
  maps=("$TAX_PATH_NCBI"/*.acc2taxid.gz)
  shopt -u nullglob
  if (( ${#maps[@]} == 0 )); then
    echo "No .acc2taxid.gz files found in $TAX_PATH_NCBI" >&2
    exit 1
  fi
  [[ -z ${EXTRA_ACC2TAX_GZ:-} ]] || maps+=("$EXTRA_ACC2TAX_GZ")
  [[ -s "$TAX_PATH_NCBI/taxdump/names.dmp" && -s "$TAX_PATH_NCBI/taxdump/nodes.dmp" ]]
  gzip -dc -- "${maps[@]}" > "$RUN/acc2tax.tsv"
  [[ -s "$RUN/acc2tax.tsv" ]]
  n=$(awk 'END { print NR }' "$RUN/manifest.tsv")
  (( n > 0 ))
  array="1-${n}%${cap}"
  prev=''
  for name in precheck alnfilt refstats nsort lca dfit aggregate xrsort family genus species; do
    dep="afterok:${SLURM_JOB_ID}"
    [[ -z "$prev" ]] || dep="aftercorr:$prev"
    cpus=$THREADS
    case "$name" in
      alnfilt) cpus=1 ;;
      lca|dfit) cpus=$METADMG_THREADS ;;
      aggregate) cpus=1 ;;
    esac
    # Optional overrides: SLURM_NSORT_MEM, SLURM_LCA_TIME, etc.
    stage_key=$(printf '%s' "$name" | tr '[:lower:]' '[:upper:]')
    key="SLURM_${stage_key}_MEM"
    mem=${!key:-${SLURM_SAMPLE_MEM:-32G}}
    key="SLURM_${stage_key}_TIME"
    time=${!key:-${SLURM_SAMPLE_TIME:-1-00:00:00}}
    prev=$(submit_job "$name" "$dep" "$array" "$cpus" "$mem" "$time")
  done
  final=$(submit_job finalize "afterok:$prev" '' "${SLURM_FINAL_CPUS:-$THREADS}" \
    "${SLURM_FINAL_MEM:-64G}" "${SLURM_FINAL_TIME:-1-00:00:00}")
  echo "All sample stages submitted. Final R/export job: $final"
  exit 0
fi

if [[ "$stage" == finalize ]]; then
  : > "$RUN/library_bams.final.tsv"
  n=$(awk 'END { print NR }' "$RUN/manifest.tsv")
  for (( i=1; i<=n; i++ )); do
    [[ -s "$RUN/input_maps/$i.tsv" ]]
    cat "$RUN/input_maps/$i.tsv" >> "$RUN/library_bams.final.tsv"
  done
  manifest_tmp=$(mktemp "${LIBRARY_BAM_LIST}.slurm.XXXXXX")
  trap 'rm -f -- "$manifest_tmp"' EXIT
  cp "$RUN/library_bams.final.tsv" "$manifest_tmp"
  mv "$manifest_tmp" "$LIBRARY_BAM_LIST"
  Rscript "$RSCRIPT_RERUN" "$CONFIG"
  # The R scripts now write Parquet themselves (Step 6 is commented out in v3).
  # Set RUN_PARQUET_EXPORT=1 in the config to also run the old separate export.
  if [[ "${RUN_PARQUET_EXPORT:-0}" != 1 ]]; then
    echo "Butteracid R post-processing finished."
    exit 0
  fi
  [[ -f "$RSCRIPT_PARQUET" ]]
  mkdir -p "$PARQUET_DIR"
  Rscript "$RSCRIPT_PARQUET" "$CONFIG" \
    --input-dir "$EUK_DIR" --output-dir "$PARQUET_DIR" --prefix "$OUT_PREFIX" \
    --step5-mode "$STEP5_MODE" --metadata "$WORKDIR/metadata" \
    --read-qc "$QC_DISREGARDED_SUMMARY" --lib-qc "$LIB_QC_STATS_OUTPUT" \
    --seq-qc "$QC_STATS_OUTPUT" --merged-qc "$MERGED_QC_STATS_OUTPUT" \
    --flowcell-summary "$FLOWCELL_SUMMARY" --library-bams "$LIBRARY_BAM_LIST" \
    --processed-dir "$PROCESSED_DIR" 2>&1 | tee -a "$LOG_FILE"
  echo "Butteracid finished: $PARQUET_DIR"
  exit 0
fi

idx=${SLURM_ARRAY_TASK_ID:?Missing Slurm array index}
[[ "$idx" =~ ^[1-9][0-9]*$ ]]
line=$(sed -n "${idx}p" "$RUN/manifest.tsv")
[[ -n "$line" ]]
IFS=$'\t' read -r lib bam source_type flowcells <<< "$line"
p="$EUK_DIR/$lib"
names="$TAX_PATH_NCBI/taxdump/names.dmp"
nodes="$TAX_PATH_NCBI/taxdump/nodes.dmp"
acc="$RUN/acc2tax.tsv"
if [[ "$stage" != precheck ]]; then
  IFS=$'\t' read -r lib bam source_type flowcells < "$RUN/input_maps/$idx.tsv"
fi
# Per-job scratch directory; trap cleanup never touches a production input.
samtools_threads=$(( THREADS - 1 ))
tmp=$(mktemp -d "$SORT_TMPDIR/ba_${lib}_${stage}.XXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
valid_bam() {
  samtools quickcheck -q "$1" || return $?
  samtools view -c "$1" > /dev/null  # full decode, not merely a header/EOF check
}
case "$stage" in
  precheck)
    valid_bam "$bam"
    hd=$(samtools view -H "$bam" | awk '$1=="@HD" {print}')
    if [[ "$hd" != *'SO:queryname'* ]]; then
      samtools sort -n -@ "$samtools_threads" -m "$MERGE_SORT_MEM" -T "$tmp/sort" -o "$tmp/sorted.bam" "$bam"
      valid_bam "$tmp/sorted.bam"
      # Never modify an input BAM, including an existing merged BAM.
      bam="$p.input.nsort.bam"
      mv "$tmp/sorted.bam" "$bam"
    fi
    printf '%s\t%s\t%s\t%s\n' "$lib" "$bam" "$source_type" "$flowcells" > "$RUN/input_maps/$idx.tsv"
    ;;
  alnfilt)
    "$UNICORN_BIN" alnfilt -b "$bam" --mode "$ALNFILT_MODE" --outbam "$tmp/out.bam" \
      --minani "$ALNFILT_MINANI" --maxani "$ALNFILT_MAXANI"
    valid_bam "$tmp/out.bam"
    mv "$tmp/out.bam" "$p.alnfilt.bam"
    ;;
  refstats)
    "$UNICORN_BIN" refstats -b "$p.alnfilt.bam" -o "$tmp/out.bam" \
      --outstat "$tmp/refstats.tsv" -t "$THREADS" --withtid \
      --names "$names" --nodes "$nodes" --acc2tax "$acc" --verbose
    valid_bam "$tmp/out.bam"
    [[ -s "$tmp/refstats.tsv" ]]
    mv "$tmp/refstats.tsv" "$p.alnfilt.refstats"
    mv "$tmp/out.bam" "$p.refstats.bam"
    ;;
  nsort)
    status=0
    if timeout --signal=TERM --kill-after=5m 1h "$ONEBAM_BIN" bamsort \
      -T "$THREADS" -m "${ONEBAM_SORT_MEM:-10G}" -o "$tmp/sorted.bam" "$p.refstats.bam"; then
      if ! valid_bam "$tmp/sorted.bam"; then status=1; fi
    else
      status=$?
    fi
    if (( status != 0 )); then
      # One log per library prevents concurrent append/header races.
      printf 'timestamp\tlibrary_id\tonebam_exit_code\n%s\t%s\t%s\n' "$(date -Is)" "$lib" "$status" \
        > "$RUN/logs/onebam_failure_${lib}.tsv"
      rm -f "$tmp/sorted.bam"
      samtools sort -n -@ "$samtools_threads" -m "$MERGE_SORT_MEM" -T "$tmp/sort" \
        -o "$tmp/sorted.bam" "$p.refstats.bam"
      valid_bam "$tmp/sorted.bam"
    fi
    hd=$(samtools view -H "$tmp/sorted.bam" | awk '$1=="@HD" {print}')
    [[ "$hd" == *'SO:queryname'* ]] || { echo 'Sorted BAM lacks SO:queryname' >&2; exit 1; }
    mv "$tmp/sorted.bam" "$p.sort.filtered.bam"
    ;;
  lca)
    "$METADMG_BIN" lca --names "$names" --nodes "$nodes" --acc2tax "$acc" \
      --sim_score_low "$METADMG_SIM_SCORE_LOW" --sim_score_high "$METADMG_SIM_SCORE_HIGH" \
      --how_many "$METADMG_HOW_MANY" --weight_type "$METADMG_WEIGHT_TYPE" \
      --lca_rank "$METADMG_LCA_RANK" --fix_ncbi "$METADMG_FIX_NCBI" --threads "$METADMG_THREADS" \
      --filtered_acc2tax "$tmp/$lib.acc2tax" --bam "$p.sort.filtered.bam" --out_prefix "$tmp/$lib.sort.filtered"
    gzip -t "$tmp/$lib.sort.filtered.bdamage.gz" "$tmp/$lib.sort.filtered.stat.gz"
    mv "$tmp/$lib."* "$EUK_DIR/"
    ;;
  dfit)
    "$METADMG_BIN" dfit "$p.sort.filtered.bdamage.gz" --threads "$METADMG_THREADS" \
      --names "$names" --nodes "$nodes" --showfits "$METADMG_SHOWFITS" --nopt "$METADMG_NOPT" \
      --nbootstrap "$METADMG_NBOOTSTRAP" --seed "$METADMG_SEED" --lib "$METADMG_LIB" \
      --out_prefix "$tmp/$lib.sort.filtered"
    gzip -t "$tmp/$lib.sort.filtered.dfit.gz"
    mv "$tmp/$lib."* "$EUK_DIR/"
    ;;
  aggregate)
    "$METADMG_BIN" aggregate "$p.sort.filtered.bdamage.gz" --names "$names" --nodes "$nodes" \
      --lcastat "$p.sort.filtered.stat.gz" --dfit "$p.sort.filtered.dfit.gz" --out_prefix "$tmp/$lib.sort.filtered.agg"
    gzip -t "$tmp/$lib.sort.filtered.agg.stat.gz"
    mv "$tmp/$lib."* "$EUK_DIR/"
    ;;
  xrsort)
    samtools sort -t XR -@ "$samtools_threads" -m "$MERGE_SORT_MEM" -T "$tmp/sort" \
      -o "$tmp/sorted.bam" "$p.refstats.bam"
    valid_bam "$tmp/sorted.bam"
    mv "$tmp/sorted.bam" "$p.XRsorted.bam"
    ;;
  family|genus|species)
    tag=$stage; [[ "$stage" != species ]] || tag=spec
    "$UNICORN_BIN" taxstats -b "$p.XRsorted.bam" -o "$tmp/taxstats.bam" -t "$THREADS" \
      --names "$names" --nodes "$nodes" --acc2tax "$acc" -k 17 \
      --outstat "$tmp/taxstats.tsv" --minreads 3 --rank "$stage"
    [[ -s "$tmp/taxstats.tsv" ]]
    valid_bam "$tmp/taxstats.bam"
    mv "$tmp/taxstats.tsv" "$p.filtered.$tag.taxstats"
    mv "$tmp/taxstats.bam" "$p.filtered.$tag.taxstats.bam"
    ;;
  *) echo "Unknown stage: $stage" >&2; exit 2 ;;
esac
printf '[%s] Completed %s for library %s\n' "$(date -Is)" "$stage" "$lib"
