#!/bin/bash -l
#SBATCH -J simlinkage
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# sbatch 1.8.1_linkage_jaatha.sh 101 2000 134inds_overlapped_filtered_norepeats_hwe_incinvariants_norepeats cor1_cor2to3_cnx1to3_cnx6 30,30,30,10
# sbatch 1.8.1_linkage_jaatha.sh 101 10000 134inds_overlapped_filtered_norepeats_hwe_incinvariants_norepeats cor1_cor2to3_cnx1to3_cnx6 30,30,30,10
# sbatch 1.8.1_linkage_jaatha.sh 101 1000 134inds_overlapped_filtered_norepeats_hwe_incinvariants_norepeats cor1_cor2to3_cnx1to3_cnx6 30,30,30,10

conda activate py3.8
module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
module load bedtools2

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
scaff="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/scaffold2chr"
gentools="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/06_results/fst/genomics_general"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/crow_demography/linkage_info/rawdata
#removing repeated regions from the intergenic regions
# awk -v OFS="\t" '{print $1, $2-1, $2}' genome_HC_allpaths41687_v2.5_annotated.intergenic_sites > genome_HC_allpaths41687_v2.5_annotated.intergenic_sites.bed
# bedtools intersect -a genome_HC_allpaths41687_v2.5_annotated.intergenic_sites.bed -b ${dat}/ref2.5repeats_header.bed -v > genome_HC_allpaths41687_v2.5_annotated.intergenic_sites.norepeats.bed

# cd $dat/crow_demography/linkage_info/rawdata
# for i in {0..102}
# do
# grep -w scaffold_${i} genome_HC_allpaths41687_v2.5_annotated.intergenic_sites.norepeats.bed > scaffold_${i}.bed
# done

cd $dat/crow_demography/linkage_info
#run python script to generate random position and check that it's intergenic throughout in the 2000bp, i.e. x+1999 = TRUE
#generating 102 blocks and removing scaff 60+78 which are part of chr18 highly differentiated region
python intergenic_block.py $1 $2
mv intergenic_scaff_pos.bed intergenic_scaff_pos_${2}.bed
grep -v "scaffold_60\|scaffold_78" intergenic_scaff_pos_${2}.bed | awk -v OFS="\t" 'BEGIN{print "chrom","chromStart","chromEnd"}'1 > intergenic_scaff_pos_${2}_100blocks.bed

#generate vcf with only the chosen regions
vcftools --bed intergenic_scaff_pos_${2}_100blocks.bed --keep 4pop.txt --gzvcf ${dat}/05.1_*/over*/${3}.vcf.gz --recode --out intergenic_block_${2}
grep -v "##" intergenic_block_${2}.recode.vcf > intergenic_block_${2}.cleaned.vcf
vcftools --min-alleles 2 --vcf intergenic_block_${2}.recode.vcf --recode --out intergenic_block_${2}_variantsitesonly
rm intergenic_block_${2}.recode.vcf

python3 ${dat}/fastsimcoal2/easySFS/easySFS.py -i intergenic_block_${2}_50ind_variantsitesonly.recode.vcf -p ${4}.pop -a -f --unfolded --proj ${5} -o ${4}_${2}
mkdir SFS
cd SFS
mkdir ${2}bp
cd ${2}bp
cp $dat/crow_demography/linkage_info/cor1_cor2to3_cnx1to3_cnx6_${2}/fastsimcoal2/*.obs ./

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
