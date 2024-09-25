#!/bin/bash -l
#SBATCH -J sumstats
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# for i in *.pop; do base=${i%.pop*}; sbatch /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/1.5.5_summarystats.sh 134inds_overlapped_filtered_norepeats_hwe_incinvariants ${base}; done

conda activate biotools
module load bcftools/1.10.2-gcc8
module load vcftools/0.1.14-gcc8
#picard.jar=2.25.7

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/06_results
mkdir sumstats
cd sumstats
mkdir inc_invariants_sites
cd inc_invariants_sites 
#compute persite pi etc including invariant sites
vcftools --gzvcf ${dat}/05.1_recal/overlap/$1.vcf.gz --keep $2.pop --recode --stdout | gzip -c > $2_incinvariants.vcf.gz
vcftools --gzvcf $2_incinvariants.vcf.gz --site-pi --out $2
vcftools --gzvcf $2_incinvariants.vcf.gz --TajimaD 1 --out $2
vcftools --gzvcf $2_incinvariants.vcf.gz --hardy --out $2
vcftools --gzvcf $2_incinvariants.vcf.gz --het --out $2

#window-based over 50kb windows 
vcftools --gzvcf $2.vcf.gz --window-pi 50000 --out $2_50kb
vcftools --gzvcf $2.vcf.gz --TajimaD 50000 --out $2_50kb

vcftools --gzvcf $2.vcf.gz --window-pi 10000 --out $2_10kb
vcftools --gzvcf $2.vcf.gz --TajimaD 10000 --out $2_10kb

for i in *10kb.Tajima.D; do  echo $i; awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' ${i}; done
for i in *10kb.windowed.pi; do  echo $i; awk '{ sum += $5; n++ } END { if (n > 0) print sum / n; }' ${i}; done

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"


