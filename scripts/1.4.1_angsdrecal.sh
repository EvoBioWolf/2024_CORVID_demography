#!/bin/bash -l       
#SBATCH -J angsdREC
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=30
#SBATCH --time=8-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# sbatch 1.4.1_angsdrecal.sh 134inds angsdrecalq30mindepth400genodepth3 poplist_134_angsd 121 400
#final selected filter: angsdrecalq30mindepth400genodepth3 134inds

module load angsd/0.933-gcc8
module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
conda activate biotools

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/05.1_recal

MinDepth=400
MaxDepth=6800
MinDepthInd=3
MinInd=121 #10%
minMapQ=30
minQ=30
genoMinDepth=0

angsd -out ./angsd/$1_$2_min$4 -bam bam_$1.filelist -GL 2 -doMaf 2 -doMajorMinor 1 -nThreads 30 -doCounts 1 -doDepth 1 -SNP_pval 1e-6 -ref ${ref} \
-doPlink 2 -doGeno -4 -doPost 1 -geno_minDepth ${genoMinDepth} \
-setMinDepth ${MinDepth} -setMaxDepth ${MaxDepth} -setMinDepthInd ${MinDepthInd} -minInd ${MinInd} \
-minMapQ ${minMapQ} -minQ ${minQ} 

#-minMaf ${minMAF}

cd angsd 
plink --tfile $1_$2_min$4 --allow-no-sex --allow-extra-chr --make-bed --out $1_$2_min$4
plink --bfile $1_$2_min$4 --allow-extra-chr --recode vcf --out $1_$2_min$4
plink --bfile $1_$2_min$4 --allow-extra-chr --missing
bgzip $1_$2_min$4.vcf
bcftools index $1_$2_min$4.vcf.gz

### Remove repeat regions
vcftools --exclude-bed ${dat}/ref2.5repeats_header.bed --gzvcf $1_$2_min$4.vcf.gz --recode --out $1_$2_min$4_norepeats
vcftools --vcf $1_$2_min$4_norepeats.recode.vcf --singletons --out $1_$2_min$4_norepeats
bcftools view -Oz -o $1_$2_min$4_norepeats.vcf.gz $1_$2_min$4_norepeats.recode.vcf
bcftools index $1_$2_min$4_norepeats.vcf.gz
rm $1_$2_min$4_norepeats.recode.vcf

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
