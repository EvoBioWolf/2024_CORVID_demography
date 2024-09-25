#!/bin/bash -l
#SBATCH -J neutral_invariant
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=1
#SBATCH --time=14-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out
#SBATCH --mem-per-cpu=4763mb

# sbatch 1.6.1_neutral_invariantsites.sh 134inds_chrDiploid_filtered_norepeats #8 days to complete

conda activate py2
module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
#picard.jar=2.25.7

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
neu="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/neutral"
TMPDIR="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/05.1_recal/GATKraw
grep -v "##" ${1}.gvcf.recode.vcf > ${1}_invariant.vcf
grep -v "##FORMAT=<ID=AD" ${1}.gvcf.recode.vcf > ${1}_invariant.vcf
cat header.txt ${1}_invariant.vcf > ${1}_invariant_header.vcf
mv ${1}_invariant_header.vcf ${1}_invariant.vcf
bcftools view --threads 8 -Oz -o ${1}_invariant.vcf.gz ${1}_invariant.vcf
bcftools index ${1}_invariant.vcf.gz

#Generate a vcf summary using vcfSummarizer.py by Robert Williamson
python ${neu}/scripts_rwilliamson/vcfSummarizer.py ${1}_invariant.vcf.gz \
${neu}/genome_HC_allpaths41687_v2.5_annotated.sites -q 0 -d 0 -D 1000 -L 0 -N 0 > ${1}_invariant.summary

#identify synonymous sites on each pop vcf (christen)
awk 'NR>12' $1_invariant.summary | awk '($8==4) {print $1, $2}' > ${1}_invariant_4fold.txt
awk 'NR>12' $1_invariant.summary | awk '($8==0) {print $1, $2}' > ${1}_invariant_intergene.txt
awk 'NR>12' $1_invariant.summary | awk '($8==1) {print $1, $2}' > ${1}_invariant_intron.txt

cat ${1}_invariant_4fold.txt ${1}_invariant_intergene.txt ${1}_invariant_intron.txt | sort | uniq > ${1}_invariant_all_neutral_pos.txt

#remove all scaffolds from chr18 (high linkage disequilibrium)
grep -v "scaffold_1026\|scaffold_1056\|scaffold_107\|scaffold_1223\|scaffold_246\|scaffold_261\|scaffold_271\|scaffold_305\|scaffold_320\|scaffold_373\|scaffold_458\|scaffold_60\|scaffold_70\|scaffold_750\|scaffold_78\|scaffold_927\|scaffold_971\|scaffold_995" ${1}_invariant_all_neutral_pos.txt > ${1}_invariant_all_neutral_nochr18_pos.txt
#1-based

# grep -v "scaffold_78\|scaffold_60" ${1}_invariant_all_neutral_pos.txt > ${1}_invariant_all_neutral_nochr18.pos
#625134484 sites
awk -v OFS="\t" '{print $1, $2}' ${1}_invariant_all_neutral_nochr18.pos > ${1}_invariant_all_neutral_nochr18_pos.txt

#get gatk invariant sites, which are filtered, no repeats and neutral only without chr18 
bcftools view -R ${1}_invariant_all_neutral_nochr18_pos.txt ${1}.gvcf.gz --threads 8 -Oz -o ${1}_neutral.gvcf.gz
bcftools index ${1}_neutral.gvcf.gz
bcftools index -n ${1}_neutral.gvcf.gz

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"


