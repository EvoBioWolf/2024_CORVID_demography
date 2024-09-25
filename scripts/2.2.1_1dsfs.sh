#!/bin/bash -l
#SBATCH -J 1dSFS
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/slurms/slurm-%j-%x.out

# for i in *.poplist
# do
# base=${i%.poplist*}
# sbatch /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/2.2.1_1dsfs.sh 05.1_recal/overlap/neuall_incIRQ/134inds_overlapped_filtered_norepeats_hwe_neuall_incIRQ_nomissing_ldpruned ${base} $(wc -l $i | sed "s/$i//g")
# done

module load vcftools/0.1.14-gcc8
module load bcftools/1.10.2-gcc8
conda activate easySFS

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/stairway

nseq=$(echo $(( 2*$3 )))

#folded
python3 ${dat}/fastsimcoal2/easySFS/easySFS.py -i ${dat}/${1}.vcf.gz -p ${2}.poplist -a -f -o $2_folded --prefix $2_folded --proj ${nseq}

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
