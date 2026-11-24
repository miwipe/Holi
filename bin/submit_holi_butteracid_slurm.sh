#!/usr/bin/env bash
set -euo pipefail

die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

# The Holi YAML reader supports flat keys without indentation or inline comments.
config_value() {
    awk -v key="$1" 'index($0, key ":") == 1 {
        sub(/^[^:]+:[[:space:]]*/, "")
        gsub(/"/, "")
        sub(/[[:space:]]+$/, "")
        print
        exit
    }' "$holi_config"
}

absolute_path() {
    case "$1" in
        /*) printf '%s\n' "$1" ;;
        *)  printf '%s/%s\n' "$launch_dir" "$1" ;;
    esac
}

config_path() {
    local value
    value=$(config_value "$1")
    [[ -n "$value" ]] || die "$1 is missing from $holi_config."
    absolute_path "$value"
}

# These values are read by both Bash and Butteracid's plain-text R config parser.
write_setting() {
    local key=$1 value=$2
    case "$value" in
        *[\"\$\`\\]*|*$'\n'*|*$'\r'*) die "$key contains characters unsupported by the Butteracid config reader." ;;
    esac
    printf '%s="%s"\n' "$key" "$value"
}

link_table() {
    local source=$1 target=$2
    if [[ -s "$source.gz" ]]; then
        ln -sfn -- "$source.gz" "$target.gz"
    elif [[ -s "$source" ]]; then
        ln -sfn -- "$source" "$target"
    else
        printf 'Missing or empty table: %s[.gz]\n' "$source" >&2
        return 1
    fi
}

run_butteracid() {
    local sample name source refstats
    cd "$run_dir/butteracid"
    mkdir -p inputs scratch

    # Limit the recursive R scan to this run's sample list and adapt the refstats name.
    while IFS= read -r sample; do
        for source in "$microb_out/$sample.gtdb.merged.bam" \
            "$euk_out/$sample.comp.bam" "$euk_out/$sample.euk.fastq.gz" \
            "$launch_dir/$sample.ppm.vs.d4.fq.gz" "$launch_dir/$sample.fastp.report.html"; do
            [[ -s "$source" ]] || die "Missing Holi QC input: $source"
        done
        for name in "$sample.sort.comp.filtered.agg.stat" "$sample.acc2tax"; do
            link_table "$euk_out/$name" "$PWD/inputs/$name" || exit 1
        done
        refstats="$euk_out/$sample.comp.unicorn.refstats"
        if [[ ! -s "$refstats" && ! -s "$refstats.gz" ]]; then
            refstats="$euk_out/$sample.comp.filtered.unicorn.refstats"
        fi
        link_table "$refstats" "$PWD/inputs/$sample.comp.filtered.unicorn.refstats" || exit 1
    done < "$run_dir/sample.list"

    ln -sfn -- "$bin_dir/butteracid_unicorn.R" "$bin_dir/postprocessing_functions.R" .
    {
        write_setting MODE unicorn
        write_setting LIBRARY_LIST "$run_dir/sample.list"
        write_setting QC_OUTPUT "$PWD/butteracid_qc.tsv"
        write_setting FASTQ_LANE_BASE "$lane_base"
        write_setting FASTQ_LOW_COMPLEXITY_BASE "$launch_dir"
        write_setting FASTQ_PREFILTER_BASE "$euk_out"
        write_setting BAM_BASE "$microb_out"
        write_setting PREFILTER_BAM_BASE "$euk_out"
        write_setting FASTP_REPORT_BASE "$launch_dir"
        write_setting FASTP_REPORT_PATTERN '*.fastp.report.html'
        write_setting SCRATCH_DIR "$PWD/scratch"
        write_setting LOG_FILE "$PWD/butteracid.log"
        write_setting TAX_PATH_NCBI "$tax_path"
        write_setting EUK_DIR "$PWD/inputs"
    } > butteracid.config
    local cpus=${SLURM_CPUS_PER_TASK:-10}
    (( cpus <= 10 )) || cpus=10
    BUTTERACID_QC_JOBS="$cpus" bash "$bin_dir/butteracid.sh" --force-qc "$PWD/butteracid.config"
}

load_environment() {
    [[ -n "$1" ]] || return 0
    # Conda and module initialization may inspect unset shell variables.
    set +u
    source "$1"
    set -u
}

if [[ ${1:-} == --worker ]]; then
    shift
    [[ $# -eq 2 ]] || die 'Expected a run settings file and a stage.'
    source "$1"
    cd "$launch_dir"
    case "$2" in
        holi)
            load_environment "$holi_env_setup"
            exec bash "$bin_dir/holi.sh" "$holi_config" "$run_dir/sample.list" ${holi_options[@]+"${holi_options[@]}"}
            ;;
        butteracid)
            load_environment "$butteracid_env_setup"
            run_butteracid
            ;;
        *) die "Unknown stage: $2" ;;
    esac
    exit
fi

usage() {
    echo "Usage: $0 HOLI_CONFIG SAMPLE_LIST [HOLI_OPTIONS...]"
}

if [[ ${1:-} == --help || ${1:-} == -h ]]; then
    usage
    exit 0
fi
if (( $# < 2 )); then
    usage >&2
    exit 2
fi

launch_dir=$PWD
script=$(realpath "${BASH_SOURCE[0]}")
bin_dir=$(dirname "$script")
holi_config=$(absolute_path "$1")
sample_list=$(absolute_path "$2")
shift 2
holi_options=("$@")
for option in "$@"; do
    [[ "$option" != --storage-friendly ]] || die 'Butteracid QC requires the merged BAMs removed by --storage-friendly.'
done
[[ -r "$holi_config" && -f "$holi_config" ]] || die "Cannot read $holi_config."
[[ -r "$sample_list" && -f "$sample_list" ]] || die "Cannot read $sample_list."
command -v sbatch >/dev/null || die 'sbatch is required.'

result_path=$(config_path RESULT_PATH)
euk_out=$(config_path EUK_OUT)
microb_out=$(config_path MICROB_OUT)
tax_path=$(config_path TAX_PATH_NCBI)
lane_base=$(absolute_path "${BUTTERACID_FASTQ_LANE_BASE:-$(config_path INPATH)}")
holi_env_setup=${HOLI_ENV_SETUP:+$(absolute_path "$HOLI_ENV_SETUP")}
butteracid_env_setup=${BUTTERACID_ENV_SETUP:+$(absolute_path "$BUTTERACID_ENV_SETUP")}
for setup in "$holi_env_setup" "$butteracid_env_setup"; do
    [[ -z "$setup" || -r "$setup" ]] || die "Cannot read environment setup: $setup"
done

threads=$(config_value THREADS)
jobs=$(config_value THREADSP)
for n in "$threads" "$jobs"; do
    [[ "$n" =~ ^[1-9][0-9]*$ ]] || die 'THREADS and THREADSP must be positive integers.'
done
threads_per_job=$(( threads + 1 ))
(( threads_per_job >= 12 )) || threads_per_job=12

umask 077
mkdir -p "$result_path/slurm_runs"
run_dir=$(mktemp -d "$result_path/slurm_runs/holi_butteracid.XXXXXX")
mkdir -p "$run_dir/logs" "$run_dir/butteracid"
cp -- "$holi_config" "$run_dir/holi.yml"
holi_config="$run_dir/holi.yml"
awk '
    { sub(/\r$/, "") }
    !NF { next }
    !/^[[:alnum:]_][[:alnum:]_.-]*$/ || seen[$0]++ {
        print "Invalid or repeated sample ID on line " NR > "/dev/stderr"
        bad = 1
        next
    }
    { print; count++ }
    END { exit (bad || !count) }
' "$sample_list" > "$run_dir/sample.list" || die "Expected a nonempty list of unique sample IDs: $sample_list"
{
    for key in launch_dir bin_dir holi_config euk_out microb_out tax_path lane_base holi_env_setup butteracid_env_setup run_dir; do
        printf '%s=%q\n' "$key" "${!key}"
    done
    # Use an explicit empty array for invocations without additional Holi options.
    printf 'holi_options=('
    if (( $# > 0 )); then
        printf '%q ' "${holi_options[@]}"
    fi
    printf ')\n'
} > "$run_dir/settings.sh"
printf 'stage\tjob_id\n' > "$run_dir/jobs.tsv"

submit() {
    local stage=$1 cpus=$2 mem=$3 time=$4 dependency=$5 result id
    local args=(--parsable --ntasks=1 --cpus-per-task="$cpus" --mem="$mem" --time="$time"
                --chdir="$launch_dir" --export=ALL --job-name="$stage"
                --output="$run_dir/logs/${stage}_%j.out" --error="$run_dir/logs/${stage}_%j.err")
    [[ -z "${SLURM_ACCOUNT:-}" ]] || args+=(--account="$SLURM_ACCOUNT")
    [[ -z "${SLURM_PARTITION:-}" ]] || args+=(--partition="$SLURM_PARTITION")
    [[ -z "${SLURM_QOS:-}" ]] || args+=(--qos="$SLURM_QOS")
    [[ -z "${SLURM_RESERVATION:-}" ]] || args+=(--reservation="$SLURM_RESERVATION")
    [[ -z "$dependency" ]] || args+=(--dependency="afterok:$dependency" --kill-on-invalid-dep=yes)
    result=$(sbatch "${args[@]}" "$script" --worker "$run_dir/settings.sh" "$stage") || return 1
    id=${result%%;*}
    [[ "$id" =~ ^[0-9]+$ ]] || { printf 'Unexpected sbatch response: %s\n' "$result" >&2; return 1; }
    printf '%s\t%s\n' "$stage" "$id" >> "$run_dir/jobs.tsv"
    printf '%s\n' "$id"
}

printf 'Run directory: %s\n' "$run_dir"
holi_job=$(submit holi "${SLURM_HOLI_CPUS:-$(( threads_per_job * jobs ))}" \
    "${SLURM_HOLI_MEM:-512G}" "${SLURM_HOLI_TIME:-3-00:00:00}" "") || die 'Holi submission failed.'
printf 'Holi job: %s\n' "$holi_job"
butteracid_job=$(submit butteracid "${SLURM_BUTTERACID_CPUS:-14}" \
    "${SLURM_BUTTERACID_MEM:-64G}" "${SLURM_BUTTERACID_TIME:-1-00:00:00}" "$holi_job") ||
    die "Butteracid submission failed. Holi job $holi_job is still submitted."
printf 'Butteracid job: %s (afterok:%s)\n' "$butteracid_job" "$holi_job"
