#!/bin/bash -l
#SBATCH -J mapdam
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

#cd 04*
#for i in C.torquatus_X*_cleaned.bam; do base=${i%.bam*}; sbatch ../1.3.2_mapdamage.sh ${base}; done

conda activate mapdamage

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/04_markdup

mapDamage -i $1.bam -r ${ref}
mapDamage -i $1.bam -r ${ref} -d results_$1 --plot-only

mv ./results_$1/$1.rescaled.bam ./ 
java -Xmx4g -jar $dat/picard.jar BuildBamIndex INPUT=$1.rescaled.bam

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) sec to complete this task"
HOUR=3600
echo "It takes $((($ENDTIME - $STARTTIME) / $HOUR)) hr to complete this task"
