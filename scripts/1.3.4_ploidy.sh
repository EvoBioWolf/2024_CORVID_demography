#!/bin/bash -l
#SBATCH -J ploidy
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# cd 04_markdup
# for fname in *_markdup_cleaned.bam; do base=${fname%.bam*}; sbatch ../1.3.4_ploidy.sh ${base}; done
# for fname in C.torquatus_X*_markdup_cleaned.bam; do base=${fname%.bam*}; sbatch ../1.3.4_ploidy.sh ${base}; done
# for fname in C.torquatus_X*.rescaled.bam; do base=${fname%.bam*}; sbatch ../1.3.4_ploidy.sh ${base}; done

conda activate biotools 
module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
#picard.jar=2.25.7
#GATK=4.2.1

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
gatk="GATK_4.2.6.1"
kvcf="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/knownsites_filtered.V2.vcf.gz"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/04_markdup

#extract chrM, chrZ and chrW from bam
#chrZ of the 100 males should be ploidy2 but for now lumped with ploidy1 
samtools view -h $1.bam -@ 8 $(cat ${dat}/chrHaploid.scaffolds.list) -o $1.ploidy1.bam
samtools view -h $1.bam -@ 8 $(cat ${dat}/chrDiploid.scaffolds.list) -o $1.ploidy2.bam
samtools index -@ 8 $1.ploidy1.bam
samtools index -@ 8 $1.ploidy2.bam

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
