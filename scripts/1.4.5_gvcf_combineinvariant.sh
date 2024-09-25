#!/bin/bash -l
#SBATCH -J gvcf
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=16
#SBATCH --time=10-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# sbatch 1.4.5_gvcf_combineinvariant.sh 134inds overlapped poplist_134
# combine invariant sites from GATK with the final variant vcf

module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
conda activate biotools

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
gatkset="134inds_DP3GQ0Miss10fullinfoQ100"
samtoolsset="samtools_DP3GQ0Miss10Q30"
angsdset="134inds_rescaled_angsdrecalq30mindepth400genodepth3_min121"


### Merge the filtered vcf and gvcf back together
cd ${dat}/05.1_recal/overlap
bcftools concat --threads 16 --allow-overlaps $1_$2_filtered_norepeats_hwe.vcf.gz ${dat}/05.1_recal/GATKraw/134inds_chrDiploid.V3.gvcf.gz -Oz -o $1_$2_filtered_norepeats_hwe_incinvariants.vcf.gz
bcftools index $1_$2_filtered_norepeats_hwe_incinvariants.vcf.gz
