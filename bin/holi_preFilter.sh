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
# Step 1: Mapping against GTDB (7 chunks)
# -------------------------------
log_step "Starting mapping and taxonomic classification against GTDB..."

Map against each of the 7 GTDB chunks
for db in {1..7}; do
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "bowtie2 --threads $THREADS -x $DB_PATH_bac.$db.fas.gz -U {}.ppm.vs.d4.fq.gz -k 1000 -D 15 -R 2 -N 1 -L 22 -i S,1,1.15 --np 1 --mp '1,1' --rdg '0,1' --rfg '0,1' --score-min 'L,0,-0.1' --mm --no-unal \
		 | samtools view -bS - > $MICROB_OUT/{}.gtdb.$db.bam 2> $MICROB_OUT/{}.gtdb.$db.log.txt"
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

	cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp lca \
	  --names '$TAX_PATH_BAC'/names.dmp \
	  --nodes '$TAX_PATH_BAC'/nodes.dmp \
	  --acc2tax '$TAX_PATH_BAC'/../hires-organelles-viruses-smags.acc2tax.gz \
	  --sim_score_low 0.92 --sim_score_high 1.0 --how_many 15 --weight_type 0 \
	  --fix_ncbi 0 --threads 10 \
	  --bam $MICROB_OUT/{}.gtdb.merged.bam --out_prefix $MICROB_OUT/{}"

	cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "zgrep -i -E 'Archaea|virus|bacteria' $MICROB_OUT/{}.lca.gz | cut -f1 > $MICROB_OUT/{}.bact_reads.txt"

	cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "seqtk subseq {}.ppm.vs.d4.fq.gz $MICROB_OUT/{}.bact_reads.txt | gzip > $MICROB_OUT/{}.bact_reads.fq.gz"

	cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "zcat {}.ppm.vs.d4.fq.gz | awk 'NR%4==1 {print substr($0, 2)}' > $MICROB_OUT/{}.all_reads.txt"

	cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "comm -23 <(sort $MICROB_OUT/{}.all_reads.txt) <(sort $MICROB_OUT/{}.bact_reads.txt) > $EUK_OUT/{}.euk_reads.txt"

	cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "seqtk subseq {}.ppm.vs.d4.fq.gz $EUK_OUT/{}.euk_reads.txt | gzip > $EUK_OUT/{}.euk.fastq.gz"

else

	parallel -j "$THREADSP" "getRTax --bam $MICROB_OUT/{}.gtdb.merged.bam -T $TAX_PATH_BAC/taxid.map -r '{\"domain\":[\"d__Bacteria\", \"d__Archaea\", \"d__Viruses\"]}' --threads 8 --unique --only-read-ids -p '$MICROB_OUT'/{}.bact_reads.txt" :::: "$SAMPLE_LIST"

	parallel -j "$THREADSP" "zcat $MICROB_OUT/{}.bact_reads.txt* > $MICROB_OUT/{}.bact_reads_all.txt" :::: "$SAMPLE_LIST"

	parallel -j "$THREADSP" "seqtk subseq {}.ppm.vs.d4.fq.gz $MICROB_OUT/{}.bact_reads_all.txt | gzip > $MICROB_OUT/{}.bact_reads.fq.gz" :::: "$SAMPLE_LIST"

	parallel -j "$THREADSP" "zcat {}.ppm.vs.d4.fq.gz | awk 'NR%4==1 {split(substr(\$0, 2), a, \" \"); print a[1]}' > $MICROB_OUT/{}.all_reads.txt" :::: "$SAMPLE_LIST"

	parallel -j "$THREADSP" "sort $MICROB_OUT/{}.all_reads.txt > $MICROB_OUT/{}.all_reads.sorted.txt" :::: "$SAMPLE_LIST"
	parallel -j "$THREADSP" "sort $MICROB_OUT/{}.bact_reads_all.txt > $MICROB_OUT/{}.bact_reads_all.sorted.txt" :::: "$SAMPLE_LIST"


	parallel -j "$THREADSP" "comm -23 $MICROB_OUT/{}.all_reads.sorted.txt $MICROB_OUT/{}.bact_reads_all.sorted.txt > $EUK_OUT/{}.euk_reads.txt" :::: "$SAMPLE_LIST"


	parallel -j "$THREADSP" "seqtk subseq {}.ppm.vs.d4.fq.gz $EUK_OUT/{}.euk_reads.txt | gzip > $EUK_OUT/{}.euk.fastq.gz" :::: "$SAMPLE_LIST"


fi
######## INTERLUDE FOR PREVIOUSLY MAPPED FILES #############################
Step X: Remove bacterial reads from BAM
parallel -j 4 "samtools view -h '$EUK_OUT'/{}.comp.reassign.filtered.bam | grep -v -F -f $RESULT_PATH/{}.bact_reads_all.txt | samtools view -@ 4 -b -o '$RESULT_PATH'/{}.no_bact.bam" :::: "$SAMPLE_LIST"

###################### MAPPING READS - PART 2##################################

log_step "Mapping reads to eukaryote database \(129 parts\) with bowtie2..."
for db in {1..129}; do
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "bowtie2 --threads $THREADS -k 1000 -t \
      -x $DB_PATH.$db.fas.gz -U $EUK_OUT/{}.euk.fastq.gz --no-unal --mm -t | \
      samtools view -bS - > $EUK_OUT/{}.euk.$db.bam 2> $EUK_OUT/eukaryota.$db.log.txt"
    check_success "Mapping to eukaryote database part $db"
done

log_step "Mapping reads to mitochondrion database with bowtie2..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "bowtie2 --threads $THREADS -k 1000 -t \
  -x $DB_PATH_clean/refseq_mitochondrion.genomic.fas.gz -U $EUK_OUT/{}.euk.fastq.gz --no-unal --mm -t | \
  samtools view -bS - > $EUK_OUT/{}.mito.bam 2> $EUK_OUT/mitochondrion.log.txt"
check_success "Mapping to mitochondrion database"

log_step "Mapping reads to phylonorwary database \(10 parts\) with bowtie2..."
for db in {1..10}; do
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "bowtie2 --threads $THREADS -k 1000 -t \
      -x $DB_PATH_Norwary.$db-of-10 -U $EUK_OUT/{}.euk.fastq.gz --no-unal --mm -t | \
      samtools view -bS - > $EUK_OUT/{}.phyNor.$db.bam 2> $EUK_OUT/phyloNorwary.$db.log.txt"
    check_success "Mapping to phyloNorwary database part $db"
done

log_step "Mapping reads to core NT database with bowtie2..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "bowtie2 --threads $THREADS -k 1000 -t \
  -x $DB_PATH_clean/core_nt.fas.gz -U $EUK_OUT/{}.euk.fastq.gz --no-unal --mm -t | \
  samtools view -bS - > $EUK_OUT/{}.core_nt.bam 2> $EUK_OUT/core_nt.log.txt"
check_success "Mapping to core NT database"

log_step "Mapping reads to plastid database with bowtie2..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "bowtie2 --threads $THREADS -k 1000 -t \
  -x $DB_PATH_clean/refseq_plastid.genomic.fas.gz -U $EUK_OUT/{}.euk.fastq.gz --no-unal --mm -t | \
  samtools view -bS - > $EUK_OUT/{}.pla.bam 2> $EUK_OUT/plastid.log.txt"
check_success "Mapping to plastid database"


log_step "Mapping finished. Continuing with merging..."

########################## ANALYSIS ######################################

# Now compress the BAM files using metaDMG
log_step "Compressing BAM files using metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "
  for bam in $EUK_OUT/{}*.bam; do
    outname=\$(basename \"\$bam\" .bam).comp.bam
    /projects/wintherpedersen/apps/metaDMG_14jun24/metaDMG-cpp/misc/compressbam \
      --threads 12 \
      --input \"\$bam\" \
      --output \"$EUK_OUT/\$outname\"
  done
"
check_success "Compressing BAM files"

log_step "Sorting each BAM file before merging..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "for bam in $EUK_OUT/{}*.bam; do \
  sorted_bam=\$MICROB_OUT/$(basename \$bam .bam).sorted.bam; \
  samtools sort -n -@ $THREADS -m 4G -o \$sorted_bam \$bam; \
done"
check_success "Bam files sorted"

log_step "Merging all sorted BAM files..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "samtools merge -@ $THREADS -n -f $EUK_OUT/{}.comp.sam.gz $EUK_OUT/{}*.comp.bam.sorted.bam"
check_success "Merging BAM files to sam.gz"

log_step "Compress bam..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/metaDMG_14jun24/metaDMG-cpp/misc/compressbam --threads 12 --input $EUK_OUT/{}.comp.sam.gz --output $EUK_OUT/{}.comp.bam"
check_success "merged sam.gz files with compress bam"

log_step "Sorting merged BAM file for bamfilter..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "samtools sort -@ $THREADS -m 10G -o $EUK_OUT/{}.sort.comp.bam" "$EUK_OUT/{}.comp.bam"
check_success "Sorting BAM file"

if [ "$UNICORN" = true ]; then
    
    log_step "Running unicorn reassign..."
    cat "$SAMPLE_LIST" | parallel -j 3 " /projects/wintherpedersen/apps/unicorn/unicorn reassign -b $EUK_OUT/{}.comp.bam -o $EUK_OUT/{}.comp.reassign.bam -t $THREADS &> $EUK_OUT/{}.comp.reassign.unicorn.log.txt"
    check_success "Unicorn reassign"

    log_step "Running unicorn filter..."
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/unicorn/unicorn refstats \
      -b $EUK_OUT/{}.comp.reassign.bam \
      -t $THREADS --outbam $EUK_OUT/{}.comp.reassign.filtered.bam \
	  --outstat $EUK_OUT/{}.comp.reassign.filtered.unicorn.refstats \
	  --names $TAX_PATH_NCBI/taxdump/names.dmp \
	  --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp"
	
    check_success "Unicorn refstats final filtering"
	
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/unicorn/unicorn bamstats \
      -b $EUK_OUT/{}.comp.reassign.unicorn.bam \
      -t $THREADS --outbam $EUK_OUT/{}.comp.reassign.filtered.bam \
	  --outstat $EUK_OUT/{}.comp.reassign.filtered.unicorn.bamstats \
	  --printdists $EUK_OUT/{}.comp.reassign.filtered.unicorn "
	
    check_success "Unicorn bamstats final filtering"

else
    
    log_step "Reassign BAM files with filterBAM..."
    cat "$SAMPLE_LIST" | parallel -j 1 "filterBAM reassign \
      --bam $EUK_OUT/{}.sort.comp.bam -t 4 -i 0 -A 92 -M 30G -m 5G -n 10 -s 0.0 \
      --squarem-min-improvement 0.001 --squarem-max-step-factor 2.0 \
      -o $EUK_OUT/{}.comp.reassign.bam &> $EUK_OUT/{}.comp.reassign.log.txt"
    check_success "filterBAM reassign"

    log_step "Final filtering with filterBAM..."
    cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "filterBAM filter \
      -e 0.6 -m 8G -t 12 -n 10 -A 92 -a 95 -N \
      --bam $EUK_OUT/{}.comp.reassign.bam \
      --stats $EUK_OUT/{}.comp.reassign.stats.tsv.gz \
      --stats-filtered $EUK_OUT/{}.comp.reassign.stats-filtered.tsv.gz \
      --bam-filtered $EUK_OUT/{}.comp.reassign.filtered.bam"
    check_success "Final filtering"
fi


log_step "Sorting merged BAM file for metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "samtools sort -n -@ $THREADS -m 10G -o $EUK_OUT/{}.sort.comp.reassign.filtered.bam" "$EUK_OUT/{}.comp.reassign.filtered.bam"
check_success "Sorting BAM file"

log_step "Running taxonomic classification with metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp lca \
  --names $TAX_PATH_NCBI/taxdump/names.dmp \
  --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp \
  --acc2tax <(zcat $TAX_PATH_NCBI/*.acc2taxid.gz) \
  --sim_score_low 0.95 --sim_score_high 1.0 --how_many 15 --weight_type 0 \
  --fix_ncbi 0 --threads 10 --filtered_acc2tax $EUK_OUT/{}.acc2tax \
  --bam $EUK_OUT/{}.sort.comp.reassign.filtered.bam --out_prefix $EUK_OUT/{}.sort.comp.reassign.filtered"
check_success "Taxonomic classification"

log_step "Running damage estimation with metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp dfit \
	  $EUK_OUT/{}.sort.comp.reassign.filtered.bdamage.gz --threads 6 \
	  --names $TAX_PATH_NCBI/taxdump/names.dmp \
	  --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp \
      --showfits 2 --nopt 10 \
      --nbootstrap 20 --doboot 1 --seed 1234 --lib ds \
      --out_prefix $EUK_OUT/{}.sort.comp.reassign.filtered"
check_success "Damage calculations done"

log_step "Aggregating lca and dfit metaDMG..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/metaDMG_28Nov24/metaDMG-cpp aggregate \
	  $EUK_OUT/{}.sort.comp.reassign.filtered.bdamage.gz \
	 --names $TAX_PATH_NCBI/taxdump/names.dmp \
      --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp \
      --lcastat $EUK_OUT/{}.sort.comp.reassign.filtered.stat.gz --dfit $EUK_OUT/{}.sort.comp.reassign.filtered.dfit.gz --out_prefix $EUK_OUT/{}.sort.comp.reassign.filtered.agg"
check_success "Aggregation done."

log_step "Unicorn per taxID statistics..."
cat "$SAMPLE_LIST" | parallel -j "$THREADSP" "/projects/wintherpedersen/apps/unicorn/unicorn tidstats \
  -b $EUK_OUT/{}.sort.comp.reassign.filtered.bam\
  -t $THREADS \
  --outstat $EUK_OUT/{}.comp.reassign.filtered.unicorn.tidstats \
  --names $TAX_PATH_NCBI/taxdump/names.dmp \
  --nodes $TAX_PATH_NCBI/taxdump/nodes.dmp \
  --acc2tax <(zcat $TAX_PATH_NCBI/*.acc2taxid.gz)"

check_success "Unicorn refstats final filtering" 

echo "Pipeline completed successfully." | tee -a "$LOG_FILE"
