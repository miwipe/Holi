#!/bin/bash
# Script: Holi.sh
# Description: Comprehensive automated pipeline for preprocessing, mapping, filtering, and taxonomic classification of FASTQ files.
# Requirements: GNU Parallel, fastp, vsearch, sga, bowtie2, samtools, filterBAM, conda, seqtk, getRtax

# -------------------------------
# Usage information
# -------------------------------
usage() {
    echo "Usage: $0 <config.yml> <sample.list> [--skip-preprocessing] [--lca-assignment] [--unicorn]"
    echo
    echo "Arguments:"
    echo "  <config.yml>            Path to YAML configuration file."
    echo "  <sample.list>           File containing list of sample names."
    echo
    echo "Options:"
    echo "  --skip-preprocessing    Skip trimming, merging, duplicate removal, and complexity filtering."
    echo
    echo "Example:"
    echo "  $0 config.yml sample.list"
    echo "  $0 config.yml sample.list --skip-preprocessing"
	echo "  $0 config.yml sample.list --lca-assignment"
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
  filterBAM
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

CONFIG=${1:-config.yml}
SAMPLE_LIST="$2"

load_config() {
    CONFIG_FILE="$1"
    get_value() {
        grep "^$1:" "$CONFIG_FILE" | sed 's/^.*:[[:space:]]*//' | sed 's/"//g'
    }

    LOG_FILE=$(get_value "LOG_FILE")
    THREADSP=$(get_value "THREADSP")
    DB_PATH=$(get_value "DB_PATH")
    DB_PATH_clean=$(get_value "DB_PATH_clean")
	TAX_PATH_NCBI=$(get_value "TAX_PATH_NCBI")
    DB_PATH_Norwary=$(get_value "DB_PATH_Norwary")
    DB_PATH_bac=$(get_value "DB_PATH_bac")
    TAX_PATH_BAC=$(get_value "TAX_PATH_BAC")
    TAX_PATH_BAC_ACC=$(get_value "TAX_PATH_BAC_ACC")
    THREADS=$(get_value "THREADS")
    MICROB_OUT=$(get_value "MICROB_OUT")
    RESULT_PATH=$(get_value "RESULT_PATH")
    INPATH=$(get_value "INPATH")
	EUK_OUT=$(get_value "EUK_OUT")
	LOGS=$(get_value "LOGS")
}

load_config "$CONFIG"
echo "[INFO] Loaded config from $CONFIG"


# -----------------------------
# Logging function and success function - forces the pipeline to exit if something goes wrong
# -----------------------------
GREEN="\033[0;32m"
YELLOW="\033[1;33m"
RED="\033[0;31m"
RESET="\033[0m"

timestamp() {
    date +"%Y-%m-%d %H:%M:%S"
}

log_step() {
    printf "${YELLOW}[%s][STEP]${RESET} %s\n" "$(timestamp)" "$1"
}

log_info() {
    printf "${GREEN}[%s][INFO]${RESET} %s\n" "$(timestamp)" "$1"
}

log_error() {
    printf "${RED}[%s][ERROR]${RESET} %s\n" "$(timestamp)" "$1"
}

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
REQUIRED_DIRS=("$RESULT_PATH" "$EUK_OUT" "$MICROB_OUT" "$LOGS")

for dir in "${REQUIRED_DIRS[@]}"; do
    if [ ! -d "$dir" ]; then
        log_info "Directory $dir does not exist. Creating it..."
        mkdir -p "$dir"
        if [ $? -ne 0 ]; then
            log_error "Failed to create directory $dir. Exiting."
            exit 1
        fi
    else
        log_info "Directory $dir exists."
    fi
done
# ---------------------------
## Holi pipeline:
# ---------------------------

# Check for command-line arguments
SKIP_PREPROCESSING=false
LCA_ASSIGN=false
UNICORN=false


while [[ "$1" != "" ]]; do
    case $1 in
        --skip-preprocessing ) 
            SKIP_PREPROCESSING=true 
            ;;
        --lca-assignment )
            LCA_ASSIGN=true
            ;;
        --unicorn )
            UNICORN=true
            ;;
        * ) 
            SAMPLE_LIST="$1"
            ;;
    esac
    shift
done

# Check if config file exists and is not empty
if [ ! -s "$CONFIG" ]; then
    echo "[ERROR] Config file ($CONFIG) is empty or does not exist." | tee -a "$LOG_FILE"
    exit 1
fi

# Check if sample list file exists and is not empty
if [ -z "$SAMPLE_LIST" ]; then
    echo "[ERROR] No sample.list provided. Please provide a sample list file as an argument." | tee -a "$LOG_FILE"
    exit 1
fi

if [ ! -s "$SAMPLE_LIST" ]; then
    echo "[ERROR] Sample list file ($SAMPLE_LIST) is empty or does not exist." | tee -a "$LOG_FILE"
    exit 1
fi


if [ "$SKIP_PREPROCESSING" = false ]; then
    ####################### QC ##############################
    log_step "Trimming and merging reads with fastp..."
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "fastp \
      -i \"${INPATH}/{}\"*R1*.fastq.gz \
      -I \"${INPATH}/{}\"*R2*.fastq.gz \
      -m --merged_out '{}.ppm.fq' \
      -V --detect_adapter_for_pe \
      -D --dup_calc_accuracy 5 \
      -g -x -q 30 -e 25 -l 30 -y -c -p \
      -h '{}.fastp.report.html' -w 1 > $LOGS/{}__fastp.log 2>&1"
    check_success "Trimming and merging reads"

    log_step "Removing duplicates with vsearch..."
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "vsearch \
      --fastx_uniques '{}.ppm.fq' \
      --fastqout '{}.ppm.vs.fq' \
      --minseqlength 30 \
      --strand both > $LOGS/{}__vsearch.log 2>&1"
    check_success "Duplicate removal"

    log_step "Filtering low-complexity reads with SGA..."
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "sga preprocess --dust-threshold=1 -m 30 '{}.ppm.vs.fq' -o '{}.ppm.vs.d4.fq' > $LOGS/{}__sga.log 2>&1"
    check_success "Low-complexity filtering"

    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "gzip '{}.ppm.vs.d4.fq'"
    check_success "Compressing filtered files"

    log_step "Cleaning up intermediate files..."
    rm -f *.ppm.fq *.ppm.vs.fq
    log_step "Intermediate files removed."
else
    log_step "Skipping preprocessing as requested."
fi


# -------------------------------
# Step 1: Mapping against GTDB (7 chunks)
# -------------------------------
log_step "Starting mapping and taxonomic classification against GTDB..."

# Map against each of the 7 GTDB chunks
for db in {1..7}; do
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "bowtie2 --threads $THREADS -x $DB_PATH_bac.$db.fas.gz -U {}.ppm.vs.d4.fq.gz -k 1000 -D 15 -R 2 -N 1 -L 22 -i S,1,1.15 --np 1 --mp '1,1' --rdg '0,1' --rfg '0,1' --score-min 'L,0,-0.1' --mm --no-unal 2> $LOGS/{}.gtdb.$db.bowtie2.log | samtools view -bS - > $MICROB_OUT/{}.gtdb.$db.bam"
    check_success "Mapping to GTDB chunk $db"
done


# Sort BAM files
log_step "Sorting GTDB BAM files..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "for bam in $MICROB_OUT/{}.gtdb.*.bam; do \
    sorted_bam=\$MICROB_OUT/$(basename \$bam .bam).sorted.bam; \
    samtools sort -@ $THREADS -m 4G -o \$sorted_bam \$bam; \
done"
check_success "Sorting GTDB BAM files"

# Merge BAM files from the 7 chunks
log_step "Merging GTDB BAM files..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "samtools merge -@ $THREADS -f $MICROB_OUT/{}.gtdb.merged.bam $MICROB_OUT/{}.gtdb.*.sorted.bam"
check_success "Merging GTDB BAM files"


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
            --sim_score_low 0.92 \
            --sim_score_high 1.0 \
            --how_many 15 \
            --weight_type 1 \
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
    parallel -j "$THREADSP" "\
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

fi

# ----------------------------------------------
# Step 2: Eukaryotic mapping
# ----------------------------------------------

log_step "Mapping reads to eukaryote database (129 parts) with bowtie2..."
for db in {1..129}; do
    parallel -j "$THREADSP" "\
        bowtie2 --threads $THREADS -k 1000 -t \
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
    bowtie2 --threads $THREADS -k 1000 -t \
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
        bowtie2 --threads $THREADS -k 1000 -t \
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
    bowtie2 --threads $THREADS -k 1000 -t \
        -x $DB_PATH_clean/core_nt.fas.gz \
        -U $EUK_OUT/{}.euk.fastq.gz \
        --no-unal --mm -t 2> $LOGS/{}__core_nt.log \
    | samtools view -bS - \
        > $EUK_OUT/{}.core_nt.bam" \
    :::: "$SAMPLE_LIST"
check_success "Mapping to core NT database"


log_step "Mapping reads to plastid database with bowtie2..."
parallel -j "$THREADSP" "\
    bowtie2 --threads $THREADS -k 1000 -t \
        -x $DB_PATH_clean/refseq_plastid.genomic.fas.gz \
        -U $EUK_OUT/{}.euk.fastq.gz \
        --no-unal --mm -t 2> $LOGS/{}__plastid.log \
    | samtools view -bS - \
        > $EUK_OUT/{}.pla.bam" \
    :::: "$SAMPLE_LIST"
check_success "Mapping to plastid database"


log_step "Mapping finished. Continuing with merging..."


# ------------------------------------
# Step 3: Analysis of bam files
# ------------------------------------
# --- Compress BAM files using metaDMG ---
log_step "Compressing BAM files using metaDMG..."

export EUK_OUT THREADS THREADSP LOGS

parallel -j "$THREADSP" --env EUK_OUT --env LOGS '
  sample={};

  # Find raw BAMs for this sample
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

check_success "Bam files sorted"



# --- Merge ---
log_step "Merging all sorted BAM files..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
  "samtools merge -@ $THREADS -n -f \
     $EUK_OUT/{}.comp.sam.gz $EUK_OUT/{}*.comp.sorted.bam"
check_success "Merging BAM files to sam.gz"


# --- Compress merged SAM to BAM ---
log_step "Compress bam..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
  "/projects/wintherpedersen/apps/metaDMG_14jun24/metaDMG-cpp/misc/compressbam \
     --threads 12 \
     --input $EUK_OUT/{}.comp.sam.gz \
     --output $EUK_OUT/{}.comp.bam \
     > $LOGS/{}__compress_merged.log 2>&1"
check_success "merged sam.gz files with compress bam"


# --- Unicorn or filterBAM ---
if [ "$UNICORN" = true ]; then
    

	log_step "Running unicorn filter..."

	export EUK_OUT LOGS TAX_PATH_NCBI THREADS

	cat "$SAMPLE_LIST" | parallel -j 1 '
	  sample={};
	  outbam="$EUK_OUT/${sample}.comp.filtered.bam";
	  outstat="$EUK_OUT/${sample}.comp.filtered.unicorn.refstats";
	  logfile="$LOGS/${sample}__unicorn_refstats.log";

	  if [[ -s "$outbam" && -s "$outstat" ]]; then
	    echo "[SKIP] $sample: outputs exist and are non-empty" > "$logfile"
	  else
	    /projects/wintherpedersen/apps/unicorn/unicorn refstats \
	      -b "$EUK_OUT/${sample}.comp.bam" \
	      -t "$THREADS" -minreads 3 \
	      --outbam "$outbam" \
	      --outstat "$outstat" \
	      --names "$TAX_PATH_NCBI/taxdump/names.dmp" \
	      --nodes "$TAX_PATH_NCBI/taxdump/nodes.dmp" \
	      > "$logfile" 2>&1
	  fi
	'

	check_success "Unicorn refstats final filtering"
    
	
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
      "/projects/wintherpedersen/apps/unicorn/unicorn bamstats \
        -b $EUK_OUT/{}.comp.unicorn.bam \
        -t $THREADS \
        --outbam $EUK_OUT/{}.comp.filtered.bam \
        --outstat $EUK_OUT/{}.comp.filtered.unicorn.bamstats \
        --printdists $EUK_OUT/{}.comp.filtered.unicorn \
        > $LOGS/{}__unicorn_bamstats.log 2>&1"
    check_success "Unicorn bamstats final filtering"

else
	
    log_step "Sorting merged BAM file for bamfilter..."
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
      "samtools sort -@ $THREADS -m 10G \
         -o $EUK_OUT/{}.sort.comp.bam \
            $EUK_OUT/{}.comp.bam"
    check_success "Sorting BAM file"
    
    log_step "Final filtering with filterBAM..."
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
      "filterBAM filter \
        -m 8G -t 12 -n 3 -A 92 -a 95 -N \
        --bam $EUK_OUT/{}.sort.comp.bam \
        --stats $EUK_OUT/{}.comp.stats.tsv.gz \
        --stats-filtered $EUK_OUT/{}.comp.stats-filtered.tsv.gz \
        --bam-filtered $EUK_OUT/{}.comp.filtered.bam \
        > $LOGS/{}__filterbam.log 2>&1"
    check_success "Final filtering"
fi


# --- Final sorting before metaDMG ---
log_step "Sorting merged BAM file for metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
  "samtools sort -n -@ $THREADS -m 10G \
     -o $EUK_OUT/{}.sort.comp.filtered.bam \
        $EUK_OUT/{}.comp.filtered.bam"
check_success "Sorting BAM file"


# --- metaDMG LCA ---
log_step "Running taxonomic classification with metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
  "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp lca \
     --names $TAX_PATH_NCBI/taxdump/names.dmp \
     --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp \
     --acc2tax <(zcat $TAX_PATH_NCBI/*.acc2taxid.gz /datasets/caeg_dataset/references/phylo_norway/20250127/*.acc2taxid.gz) \
     --sim_score_low 0.95 --sim_score_high 1.0 --how_many 15 --weight_type 1 \
     --fix_ncbi 0 --threads 10 --filtered_acc2tax $EUK_OUT/{}.acc2tax \
     --bam $EUK_OUT/{}.sort.comp.filtered.bam \
     --out_prefix $EUK_OUT/{}.sort.comp.filtered \
     > $LOGS/{}__metadmg_lca.log 2>&1"
check_success "Taxonomic classification"


# --- metaDMG dfit ---
log_step "Running damage estimation with metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
  "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp dfit \
     $EUK_OUT/{}.sort.comp.filtered.bdamage.gz --threads 6 \
     --names $TAX_PATH_NCBI/taxdump/names.dmp \
     --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp \
     --showfits 2 --nopt 10 \
     --nbootstrap 20 --doboot 1 --seed 1234 --lib ds \
     --out_prefix $EUK_OUT/{}.sort.comp.filtered \
     > $LOGS/{}__metadmg_dfit.log 2>&1"
check_success "Damage calculations done"


# --- metaDMG aggregate ---
log_step "Aggregating lca and dfit metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
  "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp aggregate \
     $EUK_OUT/{}.sort.comp.filtered.bdamage.gz \
     --names $TAX_PATH_NCBI/taxdump/names.dmp \
     --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp \
     --lcastat $EUK_OUT/{}.sort.comp.filtered.stat.gz \
     --dfit $EUK_OUT/{}.sort.comp.filtered.dfit.gz \
     --out_prefix $EUK_OUT/{}.sort.comp.filtered.agg \
     > $LOGS/{}__metadmg_agg.log 2>&1"
check_success "Aggregation done."


# --- Unicorn tidstats ---
log_step "Unicorn per taxID statistics..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" \
  "/projects/wintherpedersen/apps/unicorn/unicorn tidstats \
     -b $EUK_OUT/{}.sort.comp.filtered.bam \
     -t $THREADS \
     -o $EUK_OUT/{}.comp.filtered.unicorn.tidstats \
     --names $TAX_PATH_NCBI/taxdump/names.dmp \
     --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp \
     --acc2tax <(zcat $TAX_PATH_NCBI/*.acc2taxid.gz) \
     > $LOGS/{}__unicorn_tidstats.log 2>&1"
check_success "Unicorn refstats final filtering"


echo "Pipeline completed successfully." | tee -a "$LOG_FILE"
