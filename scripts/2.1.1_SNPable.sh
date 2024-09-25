#!/bin/bash -l
#SBATCH -J snpable
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=4
#SBATCH --time=4-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

#sbatch 2.0_SNPable.sh ref2.5

conda activate biotools 
#bwa=0.7.17
#samtools=1.7-1
#picard.jar=2.25.7
#bcftools/1.10.2-gcc8

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
echo $(date)
STARTTIME=$(date +%s)

# download SNPable sourcefile, move to msmc folder
# make to compile
# run from the folder

cd $dat
cd $dat/msmc/SNPable

#below is to identify non-uniquely mapping regions using Heng Li's method 
./script/splitfa ${ref} 35 | split -l 20000000
mv x* short_reads

# cd short_reads
for i in x*
do
bwa aln -R 1000000 -O 3 -E 3 ${ref} ${i} > ../bwa_aln/${i}.sai
bwa samse ${ref} ../bwa_aln/${i}.sai ${i} > ../bwa_aln/${i}.sam
gzip ../bwa_aln/${i}.sam
done

cd $dat/msmc/SNPable/bwa_aln
for fname in *.sam.gz
do
base=${fname%.sam.gz*}
gzip -dc ${base}.sam.gz | ../script/gen_raw_mask.pl > ../encoded_mask/rawMask_35_${base}.fa #encoded
../script/gen_mask -l 35 -r 0.5 ../encoded_mask/rawMask_35_${base}.fa > ../encoded_mask/mask_35_50_${base}.fa #encoded
done

cd ${dat}/msmc/SNPable/encoded_mask
cat mask_35_50_*.fa | sed 's/ 35 0.500//g' > ../${1}.mask_35_50_encoded.fa
sed 's/ 35 0.500//g' ../${1}.mask_35_50_encoded.fa > ../${1}.mask_35_50_encoded_edited.fa
rm  ../${1}.mask_35_50_encoded.fa

#modify python script to specify correct path to the mask file and output dir
cd $dat/msmc
python2 ./msmc-tools/makeMappabilityMask.py


ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) sec to complete this task"
HOUR=3600
echo "It takes $((($ENDTIME - $STARTTIME) / $HOUR)) hr to complete this task"
