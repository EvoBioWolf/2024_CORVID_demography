#!/bin/bash -l
#SBATCH -J neutral
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out
#SBATCH --mem-per-cpu=4763mb

# for i in *.pop; do base=${i%.pop*}; sbatch /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/1.6.1_neutral.sh 134inds_overlapped_filtered_norepeats_hwe ${base}; done

conda activate py2
module load vcftools/0.1.14-gcc8
#picard.jar=2.25.7

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
neu="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/neutral"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/06_results
mkdir neutral 
cd neutral 

#keep only varaint sites of each pop
vcftools --gzvcf ${dat}/05.1_recal/overlap/$1.vcf.gz --keep $2.pop --mac 1 --recode --stdout | gzip -c > $2.vcf.gz
zcat $2.vcf.gz | grep -v "INFO=<ID=\|FORMAT=<ID" | gzip -c > $2_edited.vcf.gz
rm $2.vcf.gz
mv $2_edited.vcf.gz $2.vcf.gz

#Generate a vcf summary using vcfSummarizer.py by Robert Williamson
python ${neu}/scripts_rwilliamson/vcfSummarizer.py ./$2.vcf.gz \
${neu}/genome_HC_allpaths41687_v2.5_annotated.sites -q 0 -d 0 -D 1000 -L 0 -N 0 > $2.summary

#identify synonymous sites on each pop vcf (christen)
awk 'NR<=12 || $8==3' $2.summary > ./0fold/${2}_0fold.txt
awk 'NR<=12 || $8==4' $2.summary > ./4fold/${2}_4fold.txt
awk 'NR<=12 || $8==0' $2.summary > ./intergene/${2}_intergene.txt
awk 'NR<=12 || $8==1' $2.summary > ./intron/${2}_intron.txt
awk 'NR<=12 || $8==2' $2.summary > ./exon/${2}_exon.txt

cd $dat/06_results/neutral/intron
awk 'NR>12 {print $1, $2}' ${2}_intron.txt > ${2}.pos
cd $dat/06_results/neutral/intergene
awk 'NR>12 {print $1, $2}' ${2}_intergene.txt > ${2}.pos
cd $dat/06_results/neutral/4fold
awk 'NR>12 {print $1, $2}' ${2}_4fold.txt > ${2}.pos

#once all done, do once
#cd $dat/06_results/neutral/intron
#cat *.pos | sort | uniq > all_pos.txt
#cd $dat/06_results/neutral/intergene
#cat *.pos | sort | uniq > all_pos.txt
#cd $dat/06_results/neutral/4fold   
#cat *.pos | sort | uniq > all_pos.txt
#cd $dat/06_results/neutral
#cat intron/all_pos.txt intergene/all_pos.txt 4fold/all_pos.txt | sort | uniq > neuall_incIRQ_pre.pos

#remove all scaffolds from chr18 (high linkage disequilibrium)
# grep -v "scaffold_1026\|scaffold_1056\|scaffold_107\|scaffold_1223\|scaffold_246\|scaffold_261\|scaffold_271\|scaffold_305\|scaffold_320\|scaffold_373\|scaffold_458\|scaffold_60\|scaffold_70\|scaffold_750\|scaffold_78\|scaffold_927\|scaffold_971\|scaffold_995" neuall_incIRQ_pre.pos > neuall_incIRQ.pos
# cp neuall_incIRQ.pos ${dat}/05.1_recal/overlap/neuall_incIRQ/neuall_incIRQ.pos

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"


