#!/bin/bash -l
#SBATCH -J msmc_ind
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# selected samples only
# sbatch 2.1.2_msmc_ind.sh C.cornix_IRQ641314 cnx6 80.4x

# sbatch 2.1.2_msmc_ind.sh C.cornix_S05 cnx3 29.0x
# sbatch 2.1.2_msmc_ind.sh C.corone_E14 cor1 26.7x
# sbatch 2.1.2_msmc_ind.sh C.corone_D06 cor2 25.4x
# sbatch 2.1.2_msmc_ind.sh C.cornix_A12 cnx4 19.9x
# sbatch 2.1.2_msmc_ind.sh C.orientalis_A18 ori1 19.9x
# sbatch 2.1.2_msmc_ind.sh C.cornix_IRQ641341 cnx6 22.0x
# sbatch 2.1.2_msmc_ind.sh C.corone_FPa01 cor3 17.2x
# sbatch 2.1.2_msmc_ind.sh C.torquatus_X02 pec1 16.5
# sbatch 2.1.2_msmc_ind.sh C.cornix_B04 cnx2 14.7x
# sbatch 2.1.2_msmc_ind.sh C.cornix_P11 cnx3 14.1x
# sbatch 2.1.2_msmc_ind.sh C.cornix_IT12 cnx1 15.7x
# sbatch 2.1.2_msmc_ind.sh C.corone_FPd01 cor3 14.8x

# pip install utils

conda activate biotools 
module load bcftools
#bwa=0.7.17
#samtools=1.7-1
#picard.jar=2.25.7
#bcftools/1.10.2-gcc8

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
scaff="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/scaffold2chr"

echo $(date)
STARTTIME=$(date +%s)

cd $dat
#each chr takes about 2 hrs (28chr+2=30*2=60hrs ~5days for each sample)

#only macrochr are used:
cd $scaff
cat chr1A.* chr1.* chr2.* chr3.* chr4.* chr5.* > macrochr.scaffolds
sed -i "s/scaffold_//g" macrochr.scaffolds
mv macrochr.scaffolds $dat/msmc

cd $dat/msmc
mkdir ${1}
for i in $(cat macrochr.scaffolds)
for i in {0..100}
do 
samtools mpileup -q 20 -Q 20 -C 50 -u -r scaffold_${i} -f ${ref} ${dat}/04_markdup/${1}_markdup_cleaned.ploidy2.bam | \
bcftools call -c -V indels | \
./msmc-tools/bamCaller.py $(echo $(samtools depth -r scaffold_$i ${dat}/04_markdup/${1}_markdup_cleaned.ploidy2.bam | awk '{sum += $3} END {print sum / NR}')) ./${1}/$1_scaffold_${i}_mask.bed.gz | \
gzip -c > ./${1}/$1_scaffold_$i.vcf.gz
done

for i in {101..114}
do
samtools mpileup -q 20 -Q 20 -C 50 -u -r scaffold_${i} -f ${ref} ${dat}/04_markdup/${1}_markdup_cleaned.ploidy2.bam | \
bcftools call -c -V indels | \
./msmc-tools/bamCaller.py $(echo $(samtools depth -r scaffold_$i ${dat}/04_markdup/${1}_markdup_cleaned.ploidy2.bam | awk '{sum += $3} END {print sum / NR}')) ./${1}/$1_scaffold_${i}_mask.bed.gz | \
gzip -c > ./${1}/$1_scaffold_$i.vcf.gz
done

####create input files for msmc####
for i in $(cat macrochr.scaffolds)
for i in {0..100}
do
./msmc-tools/generate_multihetsep.py --mask=./${1}/${1}_scaffold_${i}_mask.bed.gz \
--mask=${dat}/msmc/ref2.5_masked/ref2.5_scaffold_${i}.mask.bed.gz \
./${1}/${1}_scaffold_${i}.vcf.gz > ./${1}/${1}_scaffold_${i}_input.txt
done

for i in {101..114}
do
./msmc-tools/generate_multihetsep.py --mask=./${1}/${1}_scaffold_${i}_mask.bed.gz \
--mask=${dat}/msmc/ref2.5_masked/ref2.5_scaffold_${i}.mask.bed.gz \
./${1}/${1}_scaffold_${i}.vcf.gz > ./${1}/${1}_scaffold_${i}_input.txt
done

#11 scaffs are chrZ (thus filesize=0) and 4 scaffs are chr18
cd ./${1}
mkdir chr18
for i in $(cat ${scaff}/chr18.scaffolds)
do
mv *$i* ./chr18
done

####generate msmc results####
cd ${dat}/msmc
./msmc-tools/getStats $(ls -l ./${1}/${1}_scaffold_*_input.txt | awk '{if ($5 != 0) print $9}') > ./${1}/${1}_getstats.txt

for i in {1..2}
do
./msmc${i} -o ./results/${1}_msmc${i}_output_32segs $(ls -l ./${1}/${1}_scaffold_*_input.txt | awk '{if ($5 != 0) print $9}')
./msmc${i} -p 1*2+15*1+1*2 -o ./results/${1}_msmc${i}_output_19segs $(ls -l ./${1}/${1}_scaffold_*_input.txt | awk '{if ($5 != 0) print $9}')
done

#default -p:1*2+25*1+1*2+1*3 (32segs)
#after multiple runs, 22-20 segments with 2 groups (1*20+1*2) worked the best for MSMC1. but for MSMC2, less segmentation is required

 
ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) sec to complete this task"
HOUR=3600
echo "It takes $((($ENDTIME - $STARTTIME) / $HOUR)) hr to complete this task"

