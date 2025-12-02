#!/bin/bash -l
#SBATCH -J stairwayrescaled
# for i in *.poplist
# do
# base=${i%.poplist*}
# sbatch Path/04_fresh2/2.2.2_stairway_rescaled.sh ${base}
# done

conda activate biotools 

ref="Path/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="Path/04_fresh2"

echo $(date)
STARTTIME=$(date +%s)


cd $dat/stairway

#edit blueprint
cd stairway_scaled

java -cp stairway_plot_es Stairbuilder ${1}_folded.blueprint
bash ${1}_folded.blueprint.sh
bash ${1}_folded.blueprint.plot.sh

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
