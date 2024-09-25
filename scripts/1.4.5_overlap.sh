#!/bin/bash -l
#SBATCH -J overlap
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# sbatch 1.4.5_overlap.sh 134inds overlapped poplist_134

module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
conda activate biotools

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
gatkset="134inds_DP3GQ0Miss10fullinfoQ100"
samtoolsset="samtools_DP3GQ0Miss10Q30"
angsdset="134inds_rescaled_angsdrecalq30mindepth400genodepth3_min121"

cd ${dat}/05.1_recal/GATKraw/${gatkset}
echo "GATK_134inds"
bcftools query -f '%CHROM %POS\n' *_norepeats.vcf.gz | awk '{print $1"_"$2}' > ${dat}/05.1_recal/overlap/${gatkset}_pos.txt
	
cd ${dat}/05.1_recal/samtools/134inds/${samtoolsset}
echo "samtools_134inds"
bcftools query -f '%CHROM %POS\n' *_norepeats.vcf.gz | awk '{print $1"_"$2}' > ${dat}/05.1_recal/overlap/${samtoolsset}_pos.txt

cd ${dat}/05.1_recal/angsd
echo "angsds_134inds"
bcftools query -f '%CHROM %POS\n' ${angsdset}_norepeats.vcf.gz | awk '{print $1"_"$2}' > ${dat}/05.1_recal/overlap/${angsdset}_pos.txt

cd ${dat}/05.1_recal/overlap
conda activate renv
Rscript ${dat}/overlap.R

awk 'FS="_" {print $1"_"$2, $3}' overlapping_variants_pos.txt | awk 'NR>1' > overlap_allcallers.pos
vcftools --positions overlap_allcallers.pos --gzvcf ${dat}/05.1_*/samtools/${1}/${samtoolsset}/*_norepeats.vcf.gz --recode --out $1_$2 
vcftools --vcf $1_$2.recode.vcf --singletons --out $1_$2
bcftools view -Oz -o $1_$2.vcf.gz $1_$2.recode.vcf
rm $1_$2.recode.vcf

##### FILTER VCF
MQ=40
QUAL=30
MINDPS=3
MISS=0.1

#Info level
bcftools view --threads 20 -Oz -o ${1}_${2}.V1.vcf.gz --exclude \
        "QUAL < $QUAL || MQ < $MQ" \
        $1_$2.vcf.gz

#Geno level
bcftools filter --threads 20 -Oz -o ${1}_${2}.V2.vcf.gz --set-GTs . --include \
          "FMT/DP >= $MINDPS" \
          ${1}_${2}.V1.vcf.gz

#remove missingness
bcftools view --threads 20 -Oz -o ${1}_${2}_filtered_norepeats.vcf.gz --include \
          "F_MISSING <= $MISS" \
          $1_$2.V2.vcf.gz

#annotate and sort + poplist also sorted alphabetically
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' $1_$2_filtered_norepeats.vcf.gz -Oz -o $1_$2_filtered_norepeats_ID.vcf.gz
java -Xmx20g -Djava.io.tmpdir=$dat -jar $dat/picard.jar SortVcf TMP_DIR=$dat \
      MAX_RECORDS_IN_RAM=10000 \
      I=$1_$2_filtered_norepeats_ID.vcf.gz \
      O=$1_$2_filtered_norepeats_ID_sorted.vcf
rm $1_$2.V*.vcf.gz $1_$2_filtered_norepeats_ID.vcf.gz
mv $1_$2_filtered_norepeats_ID_sorted.vcf $1_$2_filtered_norepeats.vcf
bgzip $1_$2_filtered_norepeats.vcf
bcftools index $1_$2_filtered_norepeats.vcf.gz
vcftools --gzvcf $1_$2_filtered_norepeats.vcf.gz --freq --out  $1_$2_filtered_norepeats
awk 'NR>1 {print $1, $2}' OFS="\t"  $1_$2_filtered_norepeats.frq >  $1_$2_filtered_norepeats.pos

#LDprune including missing snps
plink --vcf $1_$2_filtered_norepeats.vcf.gz --allow-no-sex --allow-extra-chr --indep-pairwise 25 10 0.1
plink --vcf $1_$2_filtered_norepeats.vcf.gz --allow-extra-chr --extract plink.prune.in --make-bed --out $1_$2_filtered_norepeats_ldpruned
plink --bfile $1_$2_filtered_norepeats_ldpruned --allow-extra-chr --recode vcf --out $1_$2_filtered_norepeats_ldpruned
plink --bfile $1_$2_filtered_norepeats_ldpruned --allow-extra-chr --recode structure --out $1_$2_filtered_norepeats_ldpruned
bgzip $1_$2_filtered_norepeats_ldpruned.vcf
bcftools index $1_$2_filtered_norepeats_ldpruned.vcf.gz
sed -i "s/scaffold_//g" $1_$2_filtered_norepeats_ldpruned.bim

#assess missingness across samples
plink --vcf $1_$2_filtered_norepeats.vcf.gz --allow-extra-chr -allow-no-sex --missing
mv plink.imiss $1_$2_filtered_norepeats.imiss
mv plink.lmiss $1_$2_filtered_norepeats.lmiss

#missingness
MISS=0
plink --vcf $1_$2_filtered_norepeats.vcf.gz --allow-extra-chr -allow-no-sex --out $1_$2_filtered_norepeats_nomissing --recode vcf --geno $MISS
bgzip $1_$2_filtered_norepeats_nomissing.vcf
bcftools index $1_$2_filtered_norepeats_nomissing.vcf.gz

#LDprune
plink --vcf $1_$2_filtered_norepeats_nomissing.vcf.gz --allow-no-sex --allow-extra-chr --indep-pairwise 25 10 0.1
plink --vcf $1_$2_filtered_norepeats_nomissing.vcf.gz --allow-extra-chr --extract plink.prune.in --make-bed --out $1_$2_filtered_norepeats_nomissing_ldpruned
plink --bfile $1_$2_filtered_norepeats_nomissing_ldpruned --allow-extra-chr --recode vcf --out $1_$2_filtered_norepeats_nomissing_ldpruned
plink --bfile $1_$2_filtered_norepeats_nomissing_ldpruned --allow-extra-chr --recode structure --out $1_$2_filtered_norepeats_nomissing_ldpruned
bgzip $1_$2_filtered_norepeats_nomissing_ldpruned.vcf
bcftools index $1_$2_filtered_norepeats_nomissing_ldpruned.vcf.gz

# conda activate renv
cd ${dat}/05.1_recal/overlap
Rscript ${dat}/pca.R $1_$2_filtered_norepeats 05.1_recal/overlap 05.1_recal/overlap $3
rm *.gds*
Rscript ${dat}/pca.R $1_$2_filtered_norepeats_ldpruned 05.1_recal/overlap 05.1_recal/overlap $3
rm *.gds*
