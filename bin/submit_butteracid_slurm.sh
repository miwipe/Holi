#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 CONFIG PRODUCTION_SCRIPT [--rerun] [--rebuild-paths] [--only-step5]"
}

if [[ ${1:-} == --help || ${1:-} == -h ]]; then
  usage
  exit 0
fi
if (( $# < 2 )); then
  usage >&2
  exit 2
fi

config=$(realpath "$1")
original=$(realpath "$2")
shift 2
launch_dir=$PWD
worker=$(dirname "$(realpath "${BASH_SOURCE[0]}")")/butteracid_slurm_worker.sh
for file in "$config" "$original" "$worker"; do
  [[ -f "$file" && -r "$file" ]] || { echo "Cannot read $file" >&2; exit 2; }
done

source "$config"
only=${ONLY_STEP5:-0}
rebuild=${REBUILD_SAMPLE_LIST:-0}
for arg in "$@"; do
  case "$arg" in
    --rerun) ;;
    --prod) echo "--prod has been removed; Step 5 always recomputes results." >&2; exit 2 ;;
    --only-step5) only=1 ;;
    --rebuild-paths) rebuild=1 ;;
    *) echo "Unknown flag: $arg" >&2; usage >&2; exit 2 ;;
  esac
done

command -v sbatch >/dev/null
command -v python3 >/dev/null
: "${WORKDIR:?}" "${THREADS:?}" "${THREADSP:?}" "${METADMG_THREADS:?}"
for value in "$THREADS" "$THREADSP" "$METADMG_THREADS"; do
  [[ "$value" =~ ^[1-9][0-9]*$ ]] || { echo "Thread counts must be positive integers." >&2; exit 2; }
done

umask 077
mkdir -p "$WORKDIR/slurm_runs"
run=$(mktemp -d "$WORKDIR/slurm_runs/run_$(date +%Y%m%d_%H%M%S).XXXXXX")
run=$(realpath "$run")
mkdir -p "$run/logs" "$run/input_maps"
cp "$worker" "$run/worker.sh"

# Keep the production input discovery, merging and QC code.
python3 - "$original" "$run/prepare.sh" <<'PY'
import sys
from pathlib import Path

source = Path(sys.argv[1]).read_text().splitlines()
marker = "fi   # end of Steps 1-3 block"
if source.count(marker) != 1:
    raise SystemExit("Expected the Butteracid production wrapper with its Steps 1–3 end marker.")
prefix = source[:source.index(marker) + 1]
fixes = (
    "No need to comment this out anymore: if SAMPLE_LIST already exists it is",
    "reused. Use --rebuild-paths (or REBUILD_SAMPLE_LIST=1 in config) to force.",
)
lines = [f"# {line}" if line in fixes else line for line in prefix]
Path(sys.argv[2]).write_text("\n".join(lines) + "\n")
PY

bash -n "$run/prepare.sh"
bash -n "$run/worker.sh"
{
  printf 'CONFIG=%q\nRUN=%q\nLAUNCH_DIR=%q\nONLY=%q\nREBUILD=%q\n' \
    "$config" "$run" "$launch_dir" "$only" "$rebuild"
} > "$run/settings.sh"
printf 'stage\tjob_id\n' > "$run/jobs.tsv"

args=(--parsable --ntasks=1 --cpus-per-task="${SLURM_PREP_CPUS:-$(( THREADS * THREADSP ))}"
      --mem="${SLURM_PREP_MEM:-64G}" --time="${SLURM_PREP_TIME:-2-00:00:00}"
      --job-name=butteracid_prepare --chdir="$launch_dir" --export=ALL
      --output="$run/logs/prepare_%j.out" --error="$run/logs/prepare_%j.err")
[[ -z ${SLURM_ACCOUNT:-} ]] || args+=(--account="$SLURM_ACCOUNT")
[[ -z ${SLURM_PARTITION:-} ]] || args+=(--partition="$SLURM_PARTITION")
[[ -z ${SLURM_QOS:-} ]] || args+=(--qos="$SLURM_QOS")
[[ -z ${SLURM_RESERVATION:-} ]] || args+=(--reservation="$SLURM_RESERVATION")
if ! result=$(sbatch "${args[@]}" "$run/worker.sh" "$run/settings.sh" prepare); then
  echo "Preparation submission failed. Run directory: $run" >&2
  exit 1
fi
id=${result%%;*}
[[ "$id" =~ ^[0-9]+$ ]] || { echo "Unexpected sbatch response: $result (run directory: $run)" >&2; exit 1; }
printf 'Preparation job: %s\nRun directory: %s\n' "$id" "$run"
echo "Preparation submits the library arrays and final R job."
echo 'Monitor with squeue -u "$USER".'
