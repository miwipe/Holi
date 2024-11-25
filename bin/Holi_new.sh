#!/bin/bash

# Script: Holi.sh
# Description: Comprehensive automated pipeline for preprocessing, mapping, filtering, and taxonomic classification of FASTQ files.
# Requirements: GNU Parallel, fastp, vsearch, sga, bowtie2, samtools, filterBAM, metaDMG-cpp, conda

### Add the activation of a conda file - does one exist? 

# Log file name (can be a command line input, or seperated by tool)
LOG_FILE="Holi.log"
# Define number of threads for parallel (can be command line input with default value)
THREADSP=5

# redirect all output to logfile (can be seperated if needed)
exec > >(tee -i "$LOG_FILE") 2>&1

# Funtion to check for output file after every command, to determine if the execusion was successful. 
check_success() {
    if [ $? -ne 0 ]; then
        echo "[ERROR] $1 failed. Check the logs for details." | tee -a "$LOG_FILE"
        exit 1
    fi
}

# Logging function - will report date and time and the log message.
log_step() {
    echo "[$(date)] $1" | tee -a "$LOG_FILE"
}

#######################################################################################################################################
## Holi pipeline:

# Step 1: Generate a sample list
log_step "Generating sample list..."
ll *R1_001.fastq.gz | awk '{print $9}' | cut -f1 -d_ | uniq > sample.list
check_success "Sample list generation"

# Check sample list
if [ ! -s sample.list ]; then
    echo "[ERROR] Sample list is empty or not created." | tee -a "$LOG_FILE"
    exit 1
fi
log_step "Sample list created with $(wc -l < sample.list) samples."

# Step 2: Trimming and merging with fastp
log_step "Trimming and merging reads with fastp..."
cat sample.list | parallel -j $THREADSP 'fastp \
  -i {}*R1*.fastq.gz \
  -I {}*R2*.fastq.gz \
  -m --merged_out {}.ppm.fq \
  -V --detect_adapter_for_pe \
  -D --dup_calc_accuracy 5 \
  -g -x -q 30 -e 25 -l 30 -y -c -p \
  -h {}.fastp.report.html -w 1'
check_success "Trimming and merging reads"

# Step 3: Duplicate removal with vsearch
log_step "Removing duplicates with vsearch..."
cat sample.list | parallel -j $THREADSP 'vsearch \
  --fastx_uniques {}.ppm.fq \
  --fastqout {}.ppm.vs.fq \
  --minseqlength 30 \
  --strand both'
check_success "Duplicate removal"

# Step 4: Low-complexity filtering with SGA
log_step "Filtering low-complexity reads with SGA..."
cat sample.list | parallel -j $THREADSP 'sga --dust {}.ppm.vs.fq > {}.ppm.vs.d4.fq'
check_success "Low-complexity filtering"
cat sample.list | parallel -j $THREADSP 'gzip {}.ppm.vs.d4.fq'
check_success "Compressing filtered files"

# Clean up intermediate files
log_step "Cleaning up intermediate files..."
rm -f *.ppm.fq *.ppm.vs.fq
log_step "Intermediate files removed."

# Part Two: Mapping, filtering, and taxonomic classification
# Mapping with bowtie2
log_step "Mapping reads to vertebrate mammal databases with bowtie2..."
for db in {1..10}; do
    cat sample.list | parallel -j $THREADSP 'bowtie2 --threads 24 -k 1000 -t \
      -x /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.'$db' \
      -U {}.ppm.vs.d4.fq --no-unal --mm -t | samtools view -bS - > {}.vert_mam.'$db'.bam' \
      &> vert_mam.$db.log.txt
    check_success "Mapping to database $db"
done

# Filtering with filterBAM
log_step "Filtering BAM files with filterBAM (step 1)..."
for db in {1..10}; do
    cat sample.list | parallel -j $THREADSP 'filterBAM reassign \
      --bam {}.vert_mam.'$db'.bam -t 4 -i 1 -A 90 -m 20G -n 3 \
      -o {}.vert_mam.'$db'.reassign.bam &> {}.vert_mam.'$db'.reassign.log.txt'
    check_success "filterBAM reassign for database $db"
done

# Mapping with bowtie2
log_step "Mapping reads to invertebrate databases with bowtie2..."
for db in {1..3}; do
    cat sample.list | parallel -j $THREADSP 'bowtie2 --threads 24 -k 1000 -t \
      -x /projects/wintherpedersen/data/refseq_30Aug2022/invertebrate/invert.'$db' \
      -U {}.ppm.vs.fq --no-unal --mm -t | samtools view -bS - > {}.invert.'$db'.bam' \
      &> invert.$db.log.txt
    check_success "Mapping to invertebrate database $db"
done

# Filtering with filterBAM
log_step "Filtering BAM files with filterBAM (step 1) for invertebrates..."
for db in {1..3}; do
    cat sample.list | parallel -j $THREADSP 'filterBAM reassign \
      --bam {}.invert.'$db'.bam -t 4 -i 1 -A 90 -m 20G -n 3 \
      -o {}.invert.'$db'.reassign.bam &> {}.invert.'$db'.reassign.log.txt'
    check_success "filterBAM reassign for invertebrate database $db"
done

log_step "Mapping reads to vertebrate other databases with bowtie2..."
for db in {1..8}; do
    cat sample.list | parallel -j $THREADSP 'bowtie2 --threads 24 -k 1000 -t \
      -x /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_other/vert_other.'$db' \
      -U {}.ppm.vs.fq --no-unal --mm -t | samtools view -bS - > {}.vert_other.'$db'.bam' \
      &> vert_other.$db.log.txt
    check_success "Mapping to vertebrate other database $db"
done

# Filtering with filterBAM for vert_other databases (1 to 8)
log_step "Filtering BAM files with filterBAM (step 1) for vert_other databases..."
for db in {1..8}; do
    cat sample.list | parallel -j $THREADSP 'filterBAM reassign \
      --bam {}.vert_other.'$db'.bam -t 4 -i 1 -A 90 -m 20G -n 3 \
      -o {}.vert_other.'$db'.reassign.bam &> {}.vert_other.'$db'.reassign.log.txt'
    check_success "filterBAM reassign for vertebrate other database $db"
done

# Mapping with bowtie2 for plant databases (1 to 5)
log_step "Mapping reads to plant databases with bowtie2..."
for db in {1..5}; do
    cat sample.list | parallel -j $THREADSP 'bowtie2 --threads 24 -k 1000 -t \
      -x /projects/wintherpedersen/data/refseq_30Aug2022/plant/plant.'$db' \
      -U {}.ppm.vs.fq --no-unal --mm -t | samtools view -bS - > {}.plant.'$db'.bam' \
      &> plant.$db.log.txt
    check_success "Mapping to plant database $db"
done

# Filtering with filterBAM for plant databases (1 to 5)
log_step "Filtering BAM files with filterBAM (step 1) for plant databases..."
for db in {1..5}; do
    cat sample.list | parallel -j $THREADSP 'filterBAM reassign \
      --bam {}.plant.'$db'.bam -t 4 -i 1 -A 90 -m 20G -n 3 \
      -o {}.plant.'$db'.reassign.bam &> {}.plant.'$db'.reassign.log.txt'
    check_success "filterBAM reassign for plant database $db"
done

# Mapping with bowtie2 for archaea_fungi_virus, plastid, mitochondrion, and protozoa databases
log_step "Mapping reads to archaea_fungi_virus, plastid, mitochondrion, and protozoa databases with bowtie2..."
for db in archaea_fungi_virus plastid mitochondrion protozoa; do
    # Mapping for archaea_fungi_virus
    if [ "$db" == "archaea_fungi_virus" ]; then
        cat sample.list | parallel -j $THREADSP 'bowtie2 --threads 24 -k 1000 -t \
          -x /projects/wintherpedersen/data/refseq_30Aug2022/archaea_fungi_virus/archaea_fungi_virus.fa \
          -U {}.ppm.vs.fq --no-unal --mm -t | samtools view -bS - > {}.archaea_fungi_virus.1.bam' \
          &> archaea_fungi_virus.1.log.txt
    fi
    # Mapping for plastid
    if [ "$db" == "plastid" ]; then
        cat sample.list | parallel -j $THREADSP 'bowtie2 --threads 24 -k 1000 -t \
          -x /projects/wintherpedersen/data/refseq_30Aug2022/plastid/plastid.fa \
          -U {}.ppm.vs.fq --no-unal --mm -t | samtools view -bS - > {}.plastid.1.bam' \
          &> plastid.1.log.txt
    fi
    # Mapping for mitochondrion
    if [ "$db" == "mitochondrion" ]; then
        cat sample.list | parallel -j $THREADSP 'bowtie2 --threads 24 -k 1000 -t \
          -x /projects/wintherpedersen/data/refseq_30Aug2022/mitochondrion/mitochondrion.fa \
          -U {}.ppm.vs.fq --no-unal --mm -t | samtools view -bS - > {}.mitochondrion.1.bam' \
          &> mitochondrion.1.log.txt
    fi
    # Mapping for protozoa
    if [ "$db" == "protozoa" ]; then
        cat sample.list | parallel -j $THREADSP 'bowtie2 --threads 24 -k 1000 -t \
          -x /projects/wintherpedersen/data/refseq_30Aug2022/protozoa/protozoa.fa \
          -U {}.ppm.vs.fq --no-unal --mm -t | samtools view -bS - > {}.protozoa.1.bam' \
          &> protozoa.1.log.txt
    fi
    check_success "Mapping to $db database"
done

# Filtering with filterBAM for the new databases
log_step "Filtering BAM files with filterBAM (step 1) for archaea_fungi_virus, plastid, mitochondrion, and protozoa databases..."
for db in archaea_fungi_virus plastid mitochondrion protozoa; do
    # Filtering for archaea_fungi_virus
    if [ "$db" == "archaea_fungi_virus" ]; then
        cat sample.list | parallel -j $THREADSP 'filterBAM reassign \
          --bam {}.archaea_fungi_virus.1.bam -t 4 -i 1 -A 90 -m 20G -n 3 \
          -o {}.archaea_fungi_virus.1.reassign.bam &> {}.archaea_fungi_virus.1.reassign.log.txt'
    fi
    # Filtering for plastid
    if [ "$db" == "plastid" ]; then
        cat sample.list | parallel -j $THREADSP 'filterBAM reassign \
          --bam {}.plastid.1.bam -t 4 -i 1 -A 90 -m 20G -n 3 \
          -o {}.plastid.1.reassign.bam &> {}.plastid.1.reassign.log.txt'
    fi
    # Filtering for mitochondrion
    if [ "$db" == "mitochondrion" ]; then
        cat sample.list | parallel -j $THREADSP 'filterBAM reassign \
          --bam {}.mitochondrion.1.bam -t 4 -i 1 -A 90 -m 20G -n 3 \
          -o {}.mitochondrion.1.reassign.bam &> {}.mitochondrion.1.reassign.log.txt'
    fi
    # Filtering for protozoa
    if [ "$db" == "protozoa" ]; then
        cat sample.list | parallel -j $THREADSP 'filterBAM reassign \
          --bam {}.protozoa.1.bam -t 4 -i 1 -A 90 -m 20G -n 3 \
          -o {}.protozoa.1.reassign.bam &> {}.protozoa.1.reassign.log.txt'
    fi
    check_success "filterBAM reassign for $db database"
done

# Part Seven: Mapping to NCBI nt databases (nt.1 to nt.9)
log_step "Mapping reads to NCBI nt databases with bowtie2..."
for db in {1..9}; do
    cat sample.list | parallel -j $THREADSP 'bowtie2 --threads 24 -k 1000 -t \
      -x /projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.$db \
      -U {}.ppm.vs.fq --no-unal --mm -t | samtools view -bS - > {}.nt.'$db'.bam' \
      &> nt.$db.log.txt
    check_success "Mapping to nt.$db database"
done

# Filtering with bamfilter
log_step "Filtering BAM file for database $db..."
for db in {1..9}; do
    cat sample.list | parallel -j $THREADSP 'bamfilter --input {}.nt.'$db'.bam --minMapQ 30 --maxDepth 1000 \
      --output {}.nt.'$db'reassign.bam' &> {}.nt.'$db'.reassign.log.txt
    check_success "Filtering for nt.$db database"
done

# Part Eight: Mapping to Norway Plant Complete Genomes
log_step "Mapping reads to Norway Plant complete genomes with bowtie2..."
for db in {1..7}; do
    cat sample.list | parallel -j $THREADSP 'bowtie2 --threads 24 -k 1000 -t \
      -x /projects/wintherpedersen/data/norway_plant_ctgenoms/complete_genomes/final_db/norPlantCom.$db \
      -U {}.ppm.vs.fq --no-unal --mm -t | samtools view -bS - > {}.norPlantCom.'$db'.bam' \
      &> norPlantCom.$db.log.txt
    check_success "Mapping to norPlantCom.$db database"
done

# Filtering with bamfilter
log_step "Filtering BAM files Norway Plant complete genomes..."
for db in {1..7}; do
    cat sample.list | parallel -j $THREADSP 'bamfilter --input {}.norPlantCom.'$db'.bam --minMapQ 30 --maxDepth 1000 \
      --output {}.norPlantCom.'$db'.reassign.bam' &> {}.norPlantCom.'$db'.reassign.log.txt
    check_success "Filtering for norPlantCom.$db database"
done

log_step "Merging BAM files per sample..."
cat sample.list | parallel -j $THREADSP 'samtools merge {}.*.reassign.bam {}.comp.reassign.bam -@ 24'
check_success "Merging BAM files"

log_step "Filtering BAM files with filterBAM (step 2)..."
cat sample.list | parallel -j $THREADSP 'filterBAM reassign \
  --bam {}.comp.reassign.bam -t 4 -i 0 -A 92 -m 20G -n 10 \
  -o {}.comp.reassign2.bam &> {}.comp.reassign2.log.txt'
check_success "filterBAM reassign (second pass)"

log_step "Final filtering with filterBAM..."
cat sample.list | parallel -j $THREADSP 'filterBAM filter \
  -e 0.6 -m 8G -t 4 -n 10 -A 92 -a 95 -N \
  --bam {}.comp.reassign2.bam \
  --stats {}.comp.reassign2.stats.tsv.gz \
  --stats-filtered {}.comp.reassign2.stats-filtered.tsv.gz \
  --bam-filtered {}.comp.reassign2.filtered.bam'
check_success "Final filtering"

# Taxonomic classification with metaDMG
log_step "Running taxonomic classification with metaDMG..."
cat sample.list | parallel -j $THREADSP '/projects/lundbeck/people/npl206/programmes/ngsDMG/metaDMG-cpp/metaDMG-cpp lca \
  --names /projects/wintherpedersen/data/ncbi_taxonomy_01Oct2022/names.dmp \
  --nodes /projects/wintherpedersen/data/ncbi_taxonomy_01Oct2022/nodes.dmp \
  --acc2tax /projects/wintherpedersen/data/ncbi_taxonomy_01Oct2022/combined_accession2taxid_20221112.gz \
  --sim_score_low 0.95 --sim_score_high 1.0 --how_many 30 --weight_type 1 \
  --fix_ncbi 0 --threads 4 \
  --bam {}.comp.reassign2.filtered.bam --out_prefix {}.comp.reassign2.filtered'
check_success "Taxonomic classification"

log_step "Extracting DNA damage patterns with metaDMG..."
cat sample.list | parallel -j $THREADSP '/projects/lundbeck/people/npl206/programmes/ngsDMG/metaDMG-cpp/metaDMG-cpp dfit \
  {}.comp.reassign2.filtered.bdamage.gz --threads 6 \
  --names /projects/wintherpedersen/data/ncbi_taxonomy_01Oct2022/names.dmp \
  --nodes /projects/wintherpedersen/data/ncbi_taxonomy_01Oct2022/nodes.dmp --showfits 2 --nopt 10 \
  --nbootstrap 20 --doboot 1 --seed 1234 --lib ds --out {}.comp.reassign2.filtered'
check_success "Extracting DNA damage patterns"

log_step "Merging metaDMG outputs..."
for file in *.bdamage.gz.stat.gz; do
    zcat "$file" | tail -n +2 | awk -v filename="$file" '{print filename "\t" $0}' >> metadmg_data2.tsv
done
zcat *.bdamage.gz.stat.gz | head -1 > metadmg_header.tsv
cat metadmg_header.tsv metadmg_data2.tsv > metadmg_data_final.tsv

log_step "Pipeline completed successfully. Final outputs: metadmg_data_final.tsv, bam-filter_stats.tsv"
