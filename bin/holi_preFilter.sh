#!/bin/bash

# Script: Holi.sh
# Description: Comprehensive automated pipeline for preprocessing, mapping, filtering, and taxonomic classification of FASTQ files.
# Requirements: GNU Parallel, fastp, vsearch, sga, bowtie2, samtools, filterBAM, conda, seqtk, getRtax

module load samtools/1.21
module load seqtk/1.4
module load bowtie2/2.4.2 

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
    DB_PATH_Norwary=$(get_value "DB_PATH_Norwary")
    DB_PATH_bac=$(get_value "DB_PATH_bac")
    TAX_PATH_BAC=$(get_value "TAX_PATH_BAC")
    TAX_PATH_BAC_TSV=$(get_value "TAX_PATH_BAC_TSV")
    THREADS=$(get_value "THREADS")
    OUTPUT_PATH=$(get_value "OUTPUT_PATH")
    RESULT_PATH=$(get_value "RESULT_PATH")
    INPATH=$(get_value "INPATH")
	EUK_PATH=$(get_value "EUK_PATH")
}

load_config "$CONFIG"
echo "[INFO] Loaded config from $CONFIG"

# Redirect all output to logfile
exec > >(tee -i "$LOG_FILE") 2>&1

# Function to check success
check_success() {
    if [ $? -ne 0 ]; then
        echo "[ERROR] $1 failed. Check the logs for details." | tee -a "$LOG_FILE"
        exit 1
    fi
}

# Logging function
log_step() {
    echo "[$(date)] $1" | tee -a "$LOG_FILE"
}

#######################################################################################################################################
## Holi pipeline:

# Check for command-line arguments
SKIP_PREPROCESSING=false
while [[ "$1" != "" ]]; do
    case $1 in
        --skip-preprocessing ) SKIP_PREPROCESSING=true ;;
        * ) SAMPLE_LIST="$1" ;;
    esac
    shift
done

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
      -h '{}.fastp.report.html' -w 1"
    check_success "Trimming and merging reads"

    log_step "Removing duplicates with vsearch..."
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "vsearch \
      --fastx_uniques '{}.ppm.fq' \
      --fastqout '{}.ppm.vs.fq' \
      --minseqlength 30 \
      --strand both"
    check_success "Duplicate removal"

    log_step "Filtering low-complexity reads with SGA..."
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "sga preprocess --dust-threshold=1 -m 30 '{}.ppm.vs.fq' -o '{}.ppm.vs.d4.fq'"
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
# Step 1: Mapping against GTDB
# -------------------------------
log_step "Starting mapping and taxonomic classification..."

#Mapping with Bowtie2 (commented out for now, adjust if needed)
parallel -j 4 "bowtie2 --threads 20 -x '$DB_PATH_bac' -p 8 -k 1000 -D 15 -R 2 -N 1 -L 22 -i S,1,1.15 --np 1 --mp '1,1' --rdg '0,1' --rfg '0,1' --score-min 'L,0,-0.1' --mm --no-unal -U {}.ppm.vs.d4.fq.gz -S '$OUTPUT_PATH'/{}.bam 2> '$OUTPUT_PATH'/{}_bowtie2.log" :::: "$SAMPLE_LIST"

Sorting BAM files
parallel -j 4 "samtools sort -@ 12 -m 8G -o '$OUTPUT_PATH'/{}.sort.bam '$OUTPUT_PATH'/{}.bam" :::: "$SAMPLE_LIST"

# Taxonomic classification
parallel -j 2 "getRTax --bam '$OUTPUT_PATH'/{}.bam -T $TAX_PATH_BAC_TSV -r '{\"domain\":[\"d__Bacteria\", \"d__Archaea\", \"d__Viruses\"]}' --threads 8 --unique --only-read-ids -p '$RESULT_PATH'/{}.bact_reads.txt" :::: "$SAMPLE_LIST"

parallel -j 4 "zcat $RESULT_PATH/{}.bact_reads.txt* > $RESULT_PATH/{}.bact_reads_all.txt" :::: "$SAMPLE_LIST"

# Step 2: Extract sequences for bacterial reads
parallel -j 4 "seqtk subseq {}.ppm.vs.d4.fq.gz $RESULT_PATH/{}.bact_reads_all.txt | gzip > $RESULT_PATH/{}.bact_reads.fq.gz" :::: "$SAMPLE_LIST"

# Step 3: Extract all read IDs
parallel -j 4 "zcat {}.ppm.vs.d4.fq.gz | awk 'NR%4==1 {split(substr(\$0, 2), a, \" \"); print a[1]}' > $RESULT_PATH/{}.all_reads.txt" :::: "$SAMPLE_LIST"

# Step 4: Sort read ID lists
parallel -j 4 "sort $RESULT_PATH/{}.all_reads.txt > $RESULT_PATH/{}.all_reads.sorted.txt" :::: "$SAMPLE_LIST"
parallel -j 4 "sort $RESULT_PATH/{}.bact_reads_all.txt > $RESULT_PATH/{}.bact_reads_all.sorted.txt" :::: "$SAMPLE_LIST"

# Step 5: Identify eukaryotic reads
parallel -j 4 "comm -23 $RESULT_PATH/{}.all_reads.sorted.txt $RESULT_PATH/{}.bact_reads_all.sorted.txt > $RESULT_PATH/{}.euk_reads.txt" :::: "$SAMPLE_LIST"

# Step 6: Extract eukaryotic reads
parallel -j 4 "seqtk subseq {}.ppm.vs.d4.fq.gz $RESULT_PATH/{}.euk_reads.txt | gzip > $RESULT_PATH/{}.euk.fastq.gz" :::: "$SAMPLE_LIST"

######### INTERLUDE FOR PREVIOUSLY MAPPED FILES #############################
#Step X: Remove bacterial reads from BAM
#parallel -j 4 "samtools view -h '$EUK_PATH'/{}.comp.reassign.filtered.bam | grep -v -F -f $RESULT_PATH/{}.bact_reads_all.txt | samtools view -@ 4 -b -o '$RESULT_PATH'/{}.no_bact.bam" :::: "$SAMPLE_LIST"

# ###################### MAPPING READS - PART 2##################################

log_step "Mapping reads to eukaryote database \(250 parts\) with bowtie2..."
for db in {1..250}; do
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "bowtie2 --threads $THREADS -k 1000 -t \
      -x $DB_PATH.$db.fas.gz -U $OUTPUT_PATH/{}.euk.fastq.gz --no-unal --mm -t | \
      samtools view -bS - > $OUTPUT_PATH/{}.euk.$db.bam 2> $OUTPUT_PATH/eukaryota.$db.log.txt"
    check_success "Mapping to eukaryote database part $db"
done

log_step "Mapping reads to mitochondrion database with bowtie2..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "bowtie2 --threads $THREADS -k 1000 -t \
  -x $DB_PATH_clean/refseq_mitochondrion.genomic.fas.gz -U $OUTPUT_PATH/{}.euk.fastq.gz --no-unal --mm -t | \
  samtools view -bS - > $OUTPUT_PATH/{}.mito.bam 2> $OUTPUT_PATH/mitochondrion.log.txt"
check_success "Mapping to mitochondrion database"

log_step "Mapping reads to phylonorwary database \(10 parts\) with bowtie2..."
for db in {1..10}; do
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "bowtie2 --threads $THREADS -k 1000 -t \
      -x $DB_PATH_Norwary.$db-of-10 -U $OUTPUT_PATH/{}.euk.fastq.gz --no-unal --mm -t | \
      samtools view -bS - > $OUTPUT_PATH/{}.phyNor.$db.bam 2> $OUTPUT_PATH/phyloNorwary.$db.log.txt"
    check_success "Mapping to phyloNorwary database part $db"
done

log_step "Mapping reads to core NT database with bowtie2..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "bowtie2 --threads $THREADS -k 1000 -t \
  -x $DB_PATH_clean/core_nt.fas.gz -U $OUTPUT_PATH/{}.euk.fastq.gz --no-unal --mm -t | \
  samtools view -bS - > $OUTPUT_PATH/{}.core_nt.bam 2> $OUTPUT_PATH/core_nt.log.txt"
check_success "Mapping to core NT database"

log_step "Mapping reads to plastid database with bowtie2..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "bowtie2 --threads $THREADS -k 1000 -t \
  -x $DB_PATH_clean/refseq_plastid.genomic.fas.gz -U $OUTPUT_PATH/{}.euk.fastq.gz --no-unal --mm -t | \
  samtools view -bS - > $OUTPUT_PATH/{}.pla.bam 2> $OUTPUT_PATH/plastid.log.txt"
check_success "Mapping to plastid database"


log_step "Mapping finished. Continuing with merging..."

########################## ANALYSIS ######################################

# Now compress the BAM files using metaDMG
log_step "Compressing BAM files using metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "for bam in $OUTPUT_PATH/{}*.bam; do \
  /projects/wintherpedersen/apps/metaDMG_14jun24/metaDMG-cpp/misc/compressbam --threads 12 --input \$bam --output $OUTPUT_PATH/$(basename \$bam .bam).comp.bam; \
done"
check_success "Compressing BAM files"

log_step "Sorting each BAM file before merging..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "for bam in $OUTPUT_PATH/{}*.bam; do \
  sorted_bam=\$OUTPUT_PATH/$(basename \$bam .bam).sorted.bam; \
  samtools sort -n -@ $THREADS -m 4G -o \$sorted_bam \$bam; \
done"
check_success "Bam files sorted"

log_step "Merging all sorted BAM files..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "samtools merge -@ $THREADS -n -f $OUTPUT_PATH/{}.comp.sam.gz $OUTPUT_PATH/{}*.sorted.bam"
check_success "Merging BAM files to sam.gz"

log_step "Compress bam..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/metaDMG_14jun24/metaDMG-cpp/misc/compressbam --threads 12 --input $OUTPUT_PATH/{}.comp.sam.gz --output $OUTPUT_PATH/{}.comp.bam"
check_success "merged sam.gz files with compress bam"

log_step "Sorting merged BAM file for bamfilter..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "samtools sort -@ $THREADS -m 10G -o $OUTPUT_PATH/{}.sort.comp.bam" "$OUTPUT_PATH/{}.comp.bam"
check_success "Sorting BAM file"

log_step "Filtering BAM files with filterBAM..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "filterBAM reassign \
  --bam $OUTPUT_PATH/{}.sort.comp.bam -t 4 -i 0 -A 92 -M 30G -m 5G -n 10 -s 0.0 --squarem-min-improvement 0.001 --squarem-max-step-factor 2.0 \
  -o $OUTPUT_PATH/{}.comp.reassign.bam &> $OUTPUT_PATH/{}.comp.reassign.log.txt"
check_success "filterBAM reassign"

log_step "Final filtering with filterBAM..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "filterBAM filter \
  -e 0.6 -m 8G -t 12 -n 10 -A 92 -a 95 -N \
  --bam $OUTPUT_PATH/{}.comp.reassign.bam \
  --stats $OUTPUT_PATH/{}.comp.reassign.stats.tsv.gz \
  --stats-filtered $OUTPUT_PATH/{}.comp.reassign.stats-filtered.tsv.gz \
  --bam-filtered $OUTPUT_PATH/{}.comp.reassign.filtered.bam"
check_success "Final filtering"

log_step "Sorting merged BAM file for metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "samtools sort -n -@ $THREADS -m 10G -o $OUTPUT_PATH/{}.sort.comp.reassign.filtered.bam" "$OUTPUT_PATH/{}.comp.reassign.filtered.bam"
check_success "Sorting BAM file"

log_step "Running taxonomic classification with metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp lca \
  --names /datasets/caeg_dataset/taxonomy/20250210/names.dmp \
  --nodes /datasets/caeg_dataset/taxonomy/20250210/nodes.dmp \
  --acc2tax /projects/caeg/people/bfj994/hashilan_marsh/newDBall.acc2taxid.gz \
  --sim_score_low 0.95 --sim_score_high 1.0 --how_many 15 --weight_type 0 \
  --fix_ncbi 0 --threads 10 --filtered_acc2tax $OUTPUT_PATH/{}.acc2tax \
  --bam $OUTPUT_PATH/{}.sort.comp.reassign.filtered.bam --out_prefix $OUTPUT_PATH/{}.sort.comp.reassign.filtered"
check_success "Taxonomic classification"

log_step "Running damage estimation with metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp dfit \
	  $OUTPUT_PATH/{}.sort.comp.reassign.filtered.bdamage.gz --threads 6 \
  	  --names /datasets/caeg_dataset/taxonomy/20250210/names.dmp \
  	  --nodes /datasets/caeg_dataset/taxonomy/20250210/nodes.dmp \
      --showfits 2 --nopt 10 \
      --nbootstrap 20 --doboot 1 --seed 1234 --lib ds \
      --out_prefix $OUTPUT_PATH/{}.sort.comp.reassign.filtered"
check_success "Damage calculations done"

log_step "Aggregating lca and dfit metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp aggregate \
	  $OUTPUT_PATH/{}.sort.comp.reassign.filtered.bdamage.gz \
  	  --names /datasets/caeg_dataset/taxonomy/20250210/names.dmp \
  	  --nodes /datasets/caeg_dataset/taxonomy/20250210/nodes.dmp \
      --lcastat $OUTPUT_PATH/{}.sort.comp.reassign.filtered.stat.gz --dfit $OUTPUT_PATH/{}.sort.comp.reassign.filtered.dfit.gz --out_prefix $OUTPUT_PATH/{}.sort.comp.reassign.filtered.agg"
check_success "Aggregation done."

echo "Pipeline completed successfully." | tee -a "$LOG_FILE"
