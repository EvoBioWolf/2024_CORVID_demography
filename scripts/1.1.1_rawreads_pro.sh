#!/bin/bash -l
#SBATCH -J cutadapt
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=8
#SBATCH --time=5-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/03_fresh/slurms/slurm-%j.out

# sbatch 1.1.1_rawreads_pro.sh
 
conda activate biotools 
#cutadapt=1.18
#fastqc=0.11.9 
#multiqc=0.9.1a0

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat} #path to all folders for this project
mkdir 01_trimmed 02_mapped 03_readgrp 04_markdup 05.1_recal 06_results

#trim adapter using cutadapt
cd rawreads #folder containing all rawreads
for fname in *pair_1*.fq.gz
do
base=${fname%.fq.gz*}
mv ${base}.fq.gz ${base}.R1.fq.gz
done

for fname in *pair_2*.fq.gz
do
base=${fname%.fq.gz*}
mv ${base}.fq.gz ${base}.R2.fq.gz
done

rename "pair_1_Illumina" "Illumina" *
rename "pair_2_Illumina" "Illumina" *

for fname in *.R1.fq.gz
do
base=${fname%.R1.fq.gz*}
cutadapt -j 8 -m 30 -q 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o ../01_trimmed/${base}_trimmed.R1.fq.gz -p ../01_trimmed/${base}_trimmed.R2.fq.gz ${base}.R1.fq.gz ${base}.R2.fq.gz
done

#do after all files have been trimmed
cd ../01_trimmed
mkdir 00_fastqc
for fname in *.fq.gz; do base=${fname%.fq.gz*}; fastqc ./${base}.fq.gz -t 8 --outdir ./00_fastqc; done
multiqc 00_fastqc

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) sec to complete this task"
HOUR=3600
echo "It takes $((($ENDTIME - $STARTTIME) / $HOUR)) hr to complete this task"
