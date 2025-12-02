#!/bin/bash -l
#SBATCH -J map

# for i in *.R2.fq.gz; do base=${i%_trimmed.R2.fq.gz*}; sbatch ../1.2.1_map.sh ${base}; done 

conda activate biotools 
#bwa=0.7.17
#samtools=1.7-1
#picard.jar=2.25.7

ref="Path/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="Path/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)
cd $dat

echo bwamem
cd 01_trimmed
bwa mem -t 6 -M ${ref} ${1}_trimmed.R1.fq.gz ${1}_trimmed.R2.fq.gz > ../02_mapped/${1}.sam 

echo samtools
cd ${dat}/02_mapped
samtools view -bhq 6 -@ 12 ${1}.sam > ${1}.bam
samtools sort -m 4G -@ 6 ${1}.bam -o ${1}.sorted.bam

mv ${1}.sorted.bam ./bam_sorted
mv ${1}.sam ./sam
mv ${1}.bam ./bam

rm *.sam
rm *.bam

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) sec to complete this task"
HOUR=3600
echo "It takes $((($ENDTIME - $STARTTIME) / $HOUR)) hr to complete this task"

