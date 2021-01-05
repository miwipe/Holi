#!/bin/bash
for infile in $(pwd)/*.fq
do
bname=$(basename $infile)
echo $bname
bname2=$(echo $bname | sed 's/.fq*/_holi/')
basepath=$(pwd)/
basefolder=$basepath
echo $basepath
echo $bname2
mkdir $basepath$bname2
cd $basepath$bname2
pwd

## Qaulity check and filtering
echo Step 1. Removing poly A tails
fastq-grep -v "AAAAA$" ../$bname > kmer_$bname
echo Step 2. Removing reverse complemented A tails
fastq-grep -v "^TTTTT" kmer_$bname > kmer2_$bname
echo Step 3. Removing rememnants adapter sequence 1 = AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
fastq-grep -v "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" kmer2_$bname > adap1_kmer2_$bname
echo Step 4. Removing remnants adapter sequence 2 = ATCTCGTATGCCGTCTTCTGCTTG
fastq-grep -v "ATCTCGTATGCCGTCTTCTGCTTG" adap1_kmer2_$bname > adap2_kmer2_$bname

echo Step 6. sga preprocessing
sga preprocess --dust-threshold=4 -m 30 adap2_*.fq -o adap2_kmer2_$bname.pp.fq
echo Step 7. sga indexing
sga index --algorithm=ropebwt --threads=80 adap2_kmer2_$bname.pp.fq
echo Step 8. sga filtering duplicates
sga filter --threads=80  --no-kmer-check adap2_kmer2_$bname.pp.fq -o adap2_kmer2_$bname.pp.rmdup.fq
echo Step 9. Calculating read length distribution and outputting file
cat adap2_kmer2_$bname.pp.rmdup.fq | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > adap2_kmer2_$bname.pp.rmdup.fq.read_length.txt &

ls -lh
rm kmer*
rm adap1*
rm adap2_kmer2_$bname
rm adap2_kmer2_$bname.pp.fq
rm *wt
rm *sai
rm *rmdup.discard.fa

## Mapping remaining reads against downloaded databases
for DB in /willerslev/datasets/ycw/norway_plant_ctgenoms/complete_genomes/final_db/norPlantCom.?
do
echo Mapping adap2_kmer2_$bname.pp.rmdup.fq against $DB
bowtie2 --threads 80 -k 5000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/ncbi_nt_Nov2020/nt.fa.?
do
echo Mapping adap2_kmer2_$bname.pp.rmdup.fq against $DB
bowtie2 --threads 80 -k 5000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/refseq_23dec2020/vertebrate_other/vertebrate_other.?
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 80 -k 5000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/refseq_23dec2020/vertebrate_mammalian/vertebrate_mammalian.?
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 80 -k 5000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/mwp/databases/hs37d5/hs37d5
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 44 -k 5000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/refseq_23dec2020/invertebrate/invertebrate.?
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 80 -k 5000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/refseq_23dec2020/archaea_fungi_virus/archaea_fungi_virus.fa
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 80 -k 5000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/ycw/polar_animals/arctic_animals.fa
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 80 -k 5000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/ycw/polar_animals/arctic_animals_other.fa
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 80 -k 5000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/ycw/polar_animals/arctic_animal_b3.fa
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 80 -k 1000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/refseq_23dec2020/plant/plant.?
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 80 -k 1000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/refseq_23dec2020/protozoa/protozoa.fa
do
echo Mapping $bname.fq against $DB
bowtie2 --threads 80 -k 1000 -x $DB -U adap2_kmer2_$bname.pp.rmdup.fq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

## Merging all alignment files
samtools merge --verbosity 5  $bname.merged.sam.gz *.bam -@ 60

## Sorting the merged sam.gz file
echo Printing header
time samtools view --threads 80  -H $bname.merged.sam.gz | gzip > $bname.merged.Header.sam.gz
echo Printing alignment
time samtools view --threads 80 $bname.merged.sam.gz | gzip > $bname.merged.alignment.sam.gz
echo Sorting alignment file
time /willerslev/software/gz-sort/gz-sort -S 30G -P 10 $bname.merged.alignment.sam.gz $bname.merged.alignment.sort.sam.gz
echo Merging Header and sorted alignment
time zcat $bname.merged.Header.sam.gz $bname.merged.alignment.sort.sam.gz | samtools view -h -o $bname.merged.sort.sam.gz
rm $bname.merged.Header.sam.gz $bname.merged.alignment.sam.gz $bname.merged.alignment.sort.sam.gz

ls -lh *bam

rm *nt.?*
rm *vert_other.?*
rm *vert_mam.?*
rm *vert_mam.??*
rm *invert.?*
rm *norPlantCom*
rm *viral_fungi_archaea*
rm *arctic_animal*
rm $bname.merged.sam.gz

## Running ngsLCA and metaDamage
for file in *.merged.sort.sam.gz
do
nam=/willerslev/edna/ncbi_taxonomy3dec2020/names.dmp
nod=/willerslev/edna/ncbi_taxonomy3dec2020/nodes.dmp
ac2tax=/willerslev/edna/ncbi_taxonomy3dec2020/combined_taxid_accssionNO_20201120.gz
/willerslev/edna/metadamage/metadamage lca -simscorelow 0.95 -simscorehigh 1.0 -names $nam -nodes $nod -acc2tax $ac2tax -bam $file -outnames $file.species -lca_rank species -howmany 15
/willerslev/edna/metadamage/metadamage lca -simscorelow 0.95 -simscorehigh 1.0 -names $nam -nodes $nod -acc2tax $ac2tax -bam $file -outnames $file.genus -lca_rank genus -howmany 15
/willerslev/edna/metadamage/metadamage lca -simscorelow 0.95 -simscorehigh 1.0 -names $nam -nodes $nod -acc2tax $ac2tax -bam $file -outnames $file.family -lca_rank family -howmany 15
done

cd $basepath
done
