#!/bin/bash

# Script: Holi.sh
# Description: Comprehensive automated pipeline for preprocessing, mapping, filtering, and taxonomic classification of FASTQ files.
# Requirements: GNU Parallel, fastp, vsearch, sga, bowtie2, samtools, filterBAM, conda

# # load needed tools
# conda activate holi
# module load samtools/1.21

# Log file name
LOG_FILE="holi_newdb_marsh2.log"
# Number of threads for parallel
THREADSP=7
DB_PATH="/datasets/caeg_dataset/references/ncbi/20250205/data/wgs_eukaryota"
DB_PATH_clean="/datasets/caeg_dataset/references/ncbi/20250205/data/"
DB_PATH_Norwary="/datasets/caeg_dataset/references/phylo_norway/20250127/results/shard"
DB_PATH_bac="/datasets/caeg_dataset/references/ncbi/20250205/data/refseq_bacteria.genomic"
THREADS=10
OUTPUT_PATH="/projects/caeg/people/bfj994/hashilan_marsh/out"

# Input path
INPATH="/projects/caeg/scratch/for_nicola/hashilan_marsh"

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

# Check if --skip-preprocessing is passed
SKIP_PREPROCESSING=false
for arg in "$@"; do
  if [ "$arg" == "--skip-preprocessing" ]; then
    SKIP_PREPROCESSING=true
  fi
done

# Check if sample.list is passed as an argument
SAMPLE_LIST=""
for arg in "$@"; do
  if [ -f "$arg" ]; then
    SAMPLE_LIST="$arg"
  fi
done

if [ -z "$SAMPLE_LIST" ]; then
    echo "[ERROR] No sample.list provided. Please provide a sample list file as an argument." | tee -a "$LOG_FILE"
    exit 1
fi

# Check if sample list file exists and is not empty
if [ ! -s "$SAMPLE_LIST" ]; then
    echo "[ERROR] Sample list file ($SAMPLE_LIST) is empty or does not exist." | tee -a "$LOG_FILE"
    exit 1
fi


################### QC #########################
if [ "$SKIP_PREPROCESSING" = false ]; then
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
    log_step "Skipping preprocessing steps."
fi

################### MAPPING ######################
# for db in {1..250}; do
#     cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "
#       bam_file=$OUTPUT_PATH/{}.ppm.vs.d4.euk.$db.bam
#       if [ -s \"\$bam_file\" ]; then
#         echo \"Skipping {} for part $db, BAM file already exists and is not empty.\"
#       else
#         bowtie2 --threads $THREADS -k 1000 -t \
#           -x $DB_PATH.$db.fas.gz -U {}.ppm.vs.d4.fq.gz --no-unal --mm -t | \
#           samtools view -bS - > \"\$bam_file\" 2> $OUTPUT_PATH/eukaryota.$db.log.txt
#       fi"
#     check_success "Mapping to eukaryote database part $db"
# done
#
# log_step "Mapping reads to mitochondrion database with bowtie2..."
# cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "
#   bam_file=$OUTPUT_PATH/{}.ppm.vs.d4.mito.bam
#   if [ -s \"\$bam_file\" ]; then
#     echo \"Skipping {} for mitochondrion, BAM file already exists and is not empty.\"
#   else
#     bowtie2 --threads $THREADS -k 1000 -t \
#       -x $DB_PATH_clean/refseq_mitochondrion.genomic.fas.gz -U {}.ppm.vs.d4.fq.gz --no-unal --mm -t | \
#       samtools view -bS - > \"\$bam_file\" 2> $OUTPUT_PATH/mitochondrion.log.txt
#   fi"
# check_success "Mapping to mitochondrion database"
#
# log_step "Mapping reads to phylonorwary database (10 parts) with bowtie2..."
# for db in {1..10}; do
#     cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "
#       bam_file=$OUTPUT_PATH/{}.ppm.vs.d4.phyNor.$db.bam
#       if [ -s \"\$bam_file\" ]; then
#         echo \"Skipping {} for phylonorwary part $db, BAM file already exists and is not empty.\"
#       else
#         bowtie2 --threads $THREADS -k 1000 -t \
#           -x $DB_PATH_Norwary.$db-of-10 -U {}.ppm.vs.d4.fq.gz --no-unal --mm -t | \
#           samtools view -bS - > \"\$bam_file\" 2> $OUTPUT_PATH/phyloNorwary.$db.log.txt
#       fi"
#     check_success "Mapping to phyloNorwary database part $db"
# done
#
# log_step "Mapping reads to core NT database with bowtie2..."
# cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "
#   bam_file=$OUTPUT_PATH/{}.ppm.vs.d4.core_nt.bam
#   if [ -s \"\$bam_file\" ]; then
#     echo \"Skipping {} for core NT, BAM file already exists and is not empty.\"
#   else
#     bowtie2 --threads $THREADS -k 1000 -t \
#       -x $DB_PATH_clean/core_nt.fas.gz -U {}.ppm.vs.d4.fq.gz --no-unal --mm -t | \
#       samtools view -bS - > \"\$bam_file\" 2> $OUTPUT_PATH/core_nt.log.txt
#   fi"
# check_success "Mapping to core NT database"
#
# log_step "Mapping reads to plastid database with bowtie2..."
# cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "
#   bam_file=$OUTPUT_PATH/{}.ppm.vs.d4.pla.bam
#   if [ -s \"\$bam_file\" ]; then
#     echo \"Skipping {} for plastid, BAM file already exists and is not empty.\"
#   else
#     bowtie2 --threads $THREADS -k 1000 -t \
#       -x \"$DB_PATH_clean/refseq_plastid.genomic.fas.gz\" -U {}.ppm.vs.d4.fq.gz --no-unal --mm -t | \
#       samtools view -bS - > \"\$bam_file\" 2> $OUTPUT_PATH/plastid.log.txt
#   fi"
# check_success "Mapping to plastid database"
#
# log_step "Mapping reads to bacterial database (33 parts) with bowtie2..."
# for db in {1..33}; do
#     cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "
#       bam_file=$OUTPUT_PATH/{}.ppm.vs.d4.bac.$db.bam
#       if [ -s \"\$bam_file\" ]; then
#         echo \"Skipping {} for bacterial part $db, BAM file already exists and is not empty.\"
#       else
#         bowtie2 --threads $THREADS -k 1000 -t \
#           -x $DB_PATH_bac.$db.fas.gz -U {}.ppm.vs.d4.fq.gz --no-unal --mm -t | \
#           samtools view -bS - > \"\$bam_file\" 2> $OUTPUT_PATH/bacterial.$db.log.txt
#       fi"
#     check_success "Mapping to bacterial database part $db"
# done
#
# log_step "Mapping reads to archaea database with bowtie2..."
# cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "
#   bam_file=$OUTPUT_PATH/{}.ppm.vs.d4.arc.bam
#   if [ -s \"\$bam_file\" ]; then
#     echo \"Skipping {} for archaea, BAM file already exists and is not empty.\"
#   else
#     bowtie2 --threads $THREADS -k 1000 -t \
#       -x \"$DB_PATH_clean/refseq_archaea.genomic.fas.gz\" -U {}.ppm.vs.d4.fq.gz --no-unal --mm -t | \
#       samtools view -bS - > \"\$bam_file\" 2> $OUTPUT_PATH/archaea.log.txt
#   fi"
# check_success "Mapping to archaea database"
#
# log_step "Mapping reads to viral database with bowtie2..."
# cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "
#   bam_file=$OUTPUT_PATH/{}.ppm.vs.d4.vir.bam
#   if [ -s \"\$bam_file\" ]; then
#     echo \"Skipping {} for viral, BAM file already exists and is not empty.\"
#   else
#     bowtie2 --threads $THREADS -k 1000 -t \
#       -x \"$DB_PATH_clean/refseq_viral.genomic.fas.gz\" -U {}.ppm.vs.d4.fq.gz --no-unal --mm -t | \
#       samtools view -bS - > \"\$bam_file\" 2> $OUTPUT_PATH/viral.log.txt
#   fi"
# check_success "Mapping to viral database"
#
# log_step "Mapping finished. Continuing with merging..."

# Now compress the BAM files using metaDMG
# log_step "Compressing BAM files using metaDMG..."
# cat "$SAMPLE_LIST" | parallel -j "$THREADSP" '
#   for bam in '"$OUTPUT_PATH"'/{}*.bam; do
#     /projects/wintherpedersen/apps/metaDMG_14jun24/metaDMG-cpp/misc/compressbam --threads 12 --input "$bam" --output "'"$OUTPUT_PATH"'"/$(basename "$bam" .bam).comp.bam;
#   done
# '
# check_success "Compressing BAM files"

# log_step "Sorting each BAM file before merging..."
# cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "for bam in $OUTPUT_PATH/{}*.bam; do \
#   sorted_bam=\$OUTPUT_PATH/$(basename \$bam .bam).sorted.bam; \
#   samtools sort -n -@ $THREADS -m 4G -o \$sorted_bam \$bam; \
# done"
# check_success "Bam files sorted"

log_step "Merging all sorted BAM files..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "samtools merge -@ 5 -n -f $OUTPUT_PATH/{}.comp.sam.gz $OUTPUT_PATH/{}*.comp.bam.sorted.bam"
check_success "Merging BAM files to sam.gz"

log_step "Compress bam..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/metaDMG_14jun24/metaDMG-cpp/misc/compressbam --threads 12 --input $OUTPUT_PATH/{}.comp.sam.gz --output $OUTPUT_PATH/{}.comp.bam"
check_success "merged sam.gz files with compress bam"

log_step "Filtering BAM files with filterBAM..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "filterBAM reassign \
  --bam $OUTPUT_PATH/{}.comp.bam -t 4 -i 0 -A 92 -M 30G -m 5G -n 10 -s 0.75 \
  -o $OUTPUT_PATH/{}.comp.reassign.bam &> $OUTPUT_PATH/{}.comp.reassign.log.txt"
check_success "filterBAM reassign"

log_step "Final filtering with filterBAM..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "filterBAM filter \
  -e 0.6 -m 8G -t 4 -n 10 -A 92 -a 95 -N \
  --bam $OUTPUT_PATH/{}.comp.reassign.bam \
  --stats $OUTPUT_PATH/{}.comp.reassign.stats.tsv.gz \
  --stats-filtered $OUTPUT_PATH/{}.comp.reassign.stats-filtered.tsv.gz \
  --bam-filtered $OUTPUT_PATH/{}.comp.reassign.filtered.bam"
check_success "Final filtering"

log_step "Sorting merged BAM file..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "samtools sort -n -@ $THREADS -m 10G -o $OUTPUT_PATH/{}.sort.comp.reassign.filtered.bam" "$OUTPUT_PATH/{}.comp.reassign.filtered.bam"
check_success "Sorting BAM file"

log_step "Running taxonomic classification with metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp lca \
  --names /datasets/caeg_dataset/taxonomy/20250210/names.dmp \
  --nodes /datasets/caeg_dataset/taxonomy/20250210/nodes.dmp \
  --acc2tax /projects/caeg/people/bfj994/hashilan_marsh/newDBall.acc2taxid.gz \
  --sim_score_low 0.95 --sim_score_high 1.0 --how_many 30 --weight_type 1 \
  --fix_ncbi 0 --threads 10 --filtered_acc2tax $OUTPUT_PATH/{}.acc2tax \
  --bam $OUTPUT_PATH/{}.sort.comp.reassign.filtered.bam --out_prefix $OUTPUT_PATH/{}.sort.comp.reassign.filtered"
check_success "Taxonomic classification"

log_step "Running damage estimation with metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp dfit \
	  $OUTPUT_PATH/{}.sort.comp.reassign.filtered.bdamage.gz --threads 6 \
  	  --names /datasets/caeg_dataset/taxonomy/20250210/names.dmp \
  	  --nodes /datasets/caeg_dataset/taxonomy/20250210/nodes.dmp \
      --showfits 2 --nopt 10 \
      --nbootstrap 20 --doboot 1 --seed 1234 --lib ds
      --out_prefix $OUTPUT_PATH/{}.sort.comp.reassign.filtered" 
check_success "Damage calculations done"

log_step "Aggregating lca and dfit metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp aggregate \
	  $OUTPUT_PATH/{}.sort.comp.reassign.filtered.bdamage.gz
  	  --names /datasets/caeg_dataset/taxonomy/20250210/names.dmp \
  	  --nodes /datasets/caeg_dataset/taxonomy/20250210/nodes.dmp \
      --lcastat $OUTPUT_PATH/{}.sort.comp.reassign.filtered.stat.gz --dfit $OUTPUT_PATH/{}.sort.comp.reassign.filtered.dfit.gz --out_prefix $OUTPUT_PATH/{}.sort.comp.reassign.filtered.agg 
check_success "Aggregation done."

echo "Pipeline completed successfully." | tee -a "$LOG_FILE"
