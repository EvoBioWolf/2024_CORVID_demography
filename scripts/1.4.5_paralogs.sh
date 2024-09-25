#!/bin/bash -l
#SBATCH -J paralogs
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# sbatch 1.4.5_paralogs.sh 134inds overlapped poplist_134

module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
conda activate biotools

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
gatkset="134inds_DP3GQ0Miss10fullinfoQ100"
samtoolsset="samtools_DP3GQ0Miss10Q30"
angsdset="134inds_rescaled_angsdrecalq30mindepth400genodepth3_min121"

cd ${dat}/05.1_recal/overlap

#identify paralogs using cnx3 population
vcftools --gzvcf $1_$2_filtered_norepeats.vcf.gz  --keep ${dat}/06_results/neutral/cnx3.pop --recode --out $1_$2_filtered_norepeats_cnx3
vcftools --vcf $1_$2_filtered_norepeats_cnx3.recode.vcf --hardy --out $1_$2_filtered_norepeats_cnx3
awk 'NR>1 {print $1, $2, $6<0.05}' $1_$2_filtered_norepeats_cnx3.hwe | grep -w 1 | awk '{print $1, $2}' > $1_$2_filtered_norepeats_cnx3_hwedis.pos
grep -v "scaffold_1026\|scaffold_1056\|scaffold_107\|scaffold_1223\|scaffold_246\|scaffold_261\|scaffold_271\|scaffold_305\|scaffold_320\|scaffold_373\|scaffold_458\|scaffold_60\|scaffold_70\|scaffold_750\|scaffold_78\|scaffold_927\|scaffold_971\|scaffold_995" $1_$2_filtered_norepeats_cnx3_hwedis.pos > $1_$2_filtered_norepeats_cnx3_hwedis_nochr18.pos
vcftools --gzvcf $1_$2_filtered_norepeats.vcf.gz --exclude-positions $1_$2_filtered_norepeats_cnx3_hwedis_nochr18.pos --recode --out $1_$2_filtered_norepeats_hwe
bcftools view -Oz -o $1_$2_filtered_norepeats_hwe.vcf.gz $1_$2_filtered_norepeats_hwe.recode.vcf

#LDprune including missing snps
plink --vcf $1_$2_filtered_norepeats_hwe.vcf.gz --allow-no-sex --allow-extra-chr --indep-pairwise 25 10 0.1
plink --vcf $1_$2_filtered_norepeats_hwe.vcf.gz --allow-extra-chr --extract plink.prune.in --make-bed --out $1_$2_filtered_norepeats_hwe_ldpruned
plink --bfile $1_$2_filtered_norepeats_hwe_ldpruned --allow-extra-chr --recode vcf --out $1_$2_filtered_norepeats_hwe_ldpruned
plink --bfile $1_$2_filtered_norepeats_hwe_ldpruned --allow-extra-chr --recode structure --out $1_$2_filtered_norepeats_hwe_ldpruned
bgzip $1_$2_filtered_norepeats_hwe_ldpruned.vcf
bcftools index $1_$2_filtered_norepeats_hwe_ldpruned.vcf.gz
sed -i "s/scaffold_//g" $1_$2_filtered_norepeats_hwe_ldpruned.bim
rm $1_$2_filtered_norepeats_hwe_ldpruned.nosex
rm rm $1_$2_filtered_norepeats_hwe_ldpruned.log
rm plink.nosex
rm plink.prune.in
rm plink.prune.out 

#LDprune and no missing snps
MISS=0
plink --vcf $1_$2_filtered_norepeats_hwe.vcf.gz --allow-extra-chr -allow-no-sex --out $1_$2_filtered_norepeats_hwe_nomissing --recode vcf --geno $MISS
bgzip $1_$2_filtered_norepeats_hwe_nomissing.vcf
bcftools index $1_$2_filtered_norepeats_hwe_nomissing.vcf.gz

#LDprune
plink --vcf $1_$2_filtered_norepeats_hwe_nomissing.vcf.gz --allow-no-sex --allow-extra-chr --indep-pairwise 25 10 0.1
plink --vcf $1_$2_filtered_norepeats_hwe_nomissing.vcf.gz --allow-extra-chr --extract plink.prune.in --make-bed --out $1_$2_filtered_norepeats_hwe_nomissing_ldpruned
plink --bfile $1_$2_filtered_norepeats_hwe_nomissing_ldpruned --allow-extra-chr --recode vcf --out $1_$2_filtered_norepeats_hwe_nomissing_ldpruned
bgzip $1_$2_filtered_norepeats_hwe_nomissing_ldpruned.vcf
bcftools index $1_$2_filtered_norepeats_hwe_nomissing_ldpruned.vcf.gz

conda activate renv
cd ${dat}/05.1_recal/overlap
Rscript ${dat}/pca.R $1_$2_filtered_norepeats_hwe 05.1_recal/overlap 05.1_recal/overlap $3
rm *.gds*
