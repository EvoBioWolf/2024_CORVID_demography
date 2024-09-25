#!/bin/bash -l
#SBATCH -J neusmcpp
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out


# for i in c*.poplist; do base=${i%.poplist*}; sbatch /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/2.3.1_*.sh 134inds_overlapped_filtered_norepeats_hwe_neuall_incIRQ_ldpruned ${base}; done
# for i in o*.poplist; do base=${i%.poplist*}; sbatch /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/2.3.1_*.sh 134inds_overlapped_filtered_norepeats_hwe_neuall_incIRQ_ldpruned ${base}; done
# sbatch /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/2.3.1_*.sh 134inds_overlapped_filtered_norepeats_hwe_neuall_incIRQ_ldpruned pec1

conda activate py2
module load vcftools/0.1.14-gcc8
module load bcftools

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
neu="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/neutral"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/smcpp/vcf

mkdir $2
cd $2
#keep only variant sites of each pop
vcftools --gzvcf ${dat}/05.1_recal/overlap/neuall_incIRQ/$1.vcf.gz --keep ../$2.poplist --mac 1 --recode --out $2
mv $2.recode.vcf $2.vcf
bgzip $2.vcf
bcftools index $2.vcf.gz

#split vcf to individual scaffolds
bcftools index -s $2.vcf.gz | cut -f 1 | \
while read C
do
bcftools view -O z -o ${C}_${2}.vcf.gz ${2}.vcf.gz "${C}"
done

#limit to first 100 scaffolds only due to docker rate limit
mkdir ignore
for i in {101..1400}
do
mv scaffold_${i}_$2.vcf.gz ./ignore
done 
ls scaffold_*_$2.vcf.gz > list

rm chr*_${2}.vcf.gz
rm chr*_${2}.vcf.gz.tbi

#create scaff to chr filelist 
#macrochrom 
cd ../
for fname in *_scaffolds.txt
do
base=${fname%_scaffolds.txt*}
sed "s/134inds_overlapped_filtered_norepeats/${2}/g" ${base}_scaffolds.txt > ./${2}/${base}_${2}_scaffolds.txt
done

#combine scaffolds of the same chr into a single vcf, modified to exclude some scaffs
cd $2
for fname in *_${2}_scaffolds.txt
do
base=${fname%_scaffolds.txt*}
bcftools concat --file-list ${base}_scaffolds.txt --threads 1 -o ${base}.vcf 
bgzip ${base}.vcf
tabix ${base}.vcf.gz
done

conda activate singularity

cd ${dat}/smcpp/vcf
for pop in $(cat run.poplist)
do
cd ./${pop}
for i in {1..5}
do
for fname in chr${i}_${pop}_scaffolds.txt
do
chr=${fname%_${pop}_scaffolds.txt*}
for scaff in $(sed "s/_${pop}.vcf.gz//g" ${fname})
do
singularity run --bind $PWD:$HOME docker://terhorst/smcpp:latest vcf2smc -v ${chr}_${pop}.vcf.gz ./smc_output/${pop}_${chr}_${scaff}.smc.gz ${scaff} ${pop}:$(cat ${dat}/smcpp/vcf/${pop}.poplist | paste -s -d, -)
done
done
done
sleep 21600
done

cd ${dat}/smcpp/vcf
for pop in $(cat run.poplist)
do
cd ./${pop}/smc_output
singularity run --bind $PWD/:$HOME docker://terhorst/smcpp:latest estimate 3.18e-09 $(ls -l ${pop}_chr*.smc.gz | awk '{if ($5 > 1000) print $9}')
cd ${dat}/smcpp/vcf
done

##combine plots, run only after all pops are done
#cp *.json ../results/${3}/${1}.json
#cd ../results/${3}
#singularity run --bind $PWD/:$HOME docker://terhorst/smcpp:latest plot smcpp_all_${3}.pdf *.json -g 5.79 -c

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"


