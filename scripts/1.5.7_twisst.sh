#!/bin/bash -l
#SBATCH -J twisst
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=6
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# sbatch 1.5.7_twisst.sh all 134inds_overlapped_filtered_norepeats_hwe_outgroup_biallele 4pop "-g cnx6 -g cnx3 -g cor1 -g cor2 -g O --outgroup O"
# sbatch 1.5.7_twisst.sh all 134inds_overlapped_filtered_norepeats_hwe_outgroup_biallele 3pop "-g cnx3 -g cor1 -g cor2 -g O --outgroup O"
# sbatch 1.5.7_twisst.sh all 134inds_overlapped_filtered_norepeats_hwe_outgroup_biallele 3pop_cnx6 "-g cnx6 -g cor1 -g cor2 -g O --outgroup O"

module load vcftools
module load bcftools
conda activate biotools
#picard.jar=2.25.7
#module load plink2
#using plink1.9 instead

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
scaff="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/scaffold2chr"
gentools="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/06_results/fst/genomics_general"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/06_results/twisst

#subsample for twisst 
vcftools --keep $3.pop --gzvcf ${dat}/05.1_recal/overlap/${2}.vcf.gz --recode --out $1_$3
mv $1_$3.recode.vcf $1_$3.vcf
bgzip $1_$3.vcf
bcftools index $1_$3.vcf.gz
bcftools view --threads 6 -V indels --min-alleles 2 -Oz -o $1_$3_filtered.vcf.gz $1_$3.vcf.gz
bcftools view -r $(echo $(cat chr.txt) | sed 's/ //g') -Oz -o $1_$3_filtered_100scaff.vcf.gz $1_$3.vcf.gz
#Phase with Beagle5.4
java -Xmx12g -jar beagle.jar gt=$1_$3_filtered_100scaff.vcf.gz out=$1_$3_filtered_100scaff_phased impute=true nthreads=6 window=1 overlap=0.1
python ${gentools}/VCF_processing/parseVCF.py -i $1_$3_filtered_100scaff_phased.vcf.gz | gzip > $1_$3_filtered_phased.geno.gz

#NJ tree per 50 SNPs 
python ${gentools}/phyml_sliding_windows.py -T 6 -g $1_$3_filtered_phased.geno.gz --prefix $1_$3.phyml_bionj.w50 -w 50 --windType sites --model GTR --optimise n

awk '{print $1"_A", $2}' $3_groups.tsv > $3_groups_editedA.tsv
awk '{print $1"_B", $2}' $3_groups.tsv > $3_groups_editedB.tsv
cat $3_groups_editedA.tsv $3_groups_editedB.tsv | awk '{print $2, $1}' | sort | awk '{print $2, $1}' OFS="\t" > $3_groups_edited.tsv
rm $3_groups_editedA.tsv $3_groups_editedB.tsv

python twisst.py -t $1_$3.phyml_bionj.w50.trees.gz -w $1_$3.weights.csv.gz --method complete --groupsFile $3_groups_edited.tsv $4
zcat $1_$3.weights.csv.gz | sed "s/outgroup/O/g" >  $1_$3.weights.edited.csv
gzip $1_$3.weights.edited.csv

#to investigates weights of barriers and introgression topology across the genome 
zcat $1_$3.weights.edited.csv.gz | awk 'NR>3' > $1_$3.weights.edited.csv.tmp
paste $1_$3.phyml_bionj.w50.data.tsv $1_$3.weights.edited.csv.tmp > $1_$3_combined.csv.tmp
awk 'NR>1 {print $1, $1, $2, $3, $4, $5, $6, $7, $8, $9, $7/($7+$8+$9), $8/($7+$8+$9), $9/($7+$8+$9)}' $1_$3_combined.csv.tmp > $1_$3_combined.csv
#map scaffolds to chromosome
cd $scaff
for fname in *.scaffolds
do
base=${fname%.scaffolds*}
for i in $(cat $base.scaffolds)
do
echo $i $base
awk -v scaff=$i -v chr=$base '{gsub("^"scaff"$", chr, $1)} 1' ${dat}/06_results/twisst/$1_$3_combined.csv > ${dat}/06_results/twisst/$1_$3_combined.csv.tmp && mv \-v ${dat}/06_results/twisst/$1_$3_combined.csv.tmp ${dat}/06_results/twisst/$1_$3_combined.csv
done
done

cd $dat/06_results/twisst  
nl $1_$3_combined.csv > $1_$3_combined.chr.csv
rm *.tmp

conda activate twisst
Rscript plot_all_$3.R
Rscript twisstntern.R 

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
