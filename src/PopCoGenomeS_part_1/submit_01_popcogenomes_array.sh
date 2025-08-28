#!/bin/bash
#
#SBATCH --job-name=submit_popcogenomes_01
#SBATCH --cpus-per-task=4
#SBATCH --mem=10GB
#SBATCH --nice=0
#SBATCH --partition=basic
#SBATCH --output=log/%x-%A_%a.out
#SBATCH --error=log/%x-%A_%a.err
#SBATCH --time=10:00:00

#usage
#no need for config file since the folders are numbered
#sbatch -a 1-number of chunks submit_01_popcogenomes_array.sh 

n=${SLURM_ARRAY_TASK_ID}
export n

window_size=500
export window_size 

#define directory where the script is called
HOMEFOLDER=`pwd`

#define and construct workspace on /tmp/$USER
THISWORKFOLDER=$TMPDIR

#move genomes and popcogenomes code
cp -R ${HOMEFOLDER}/PopCoGenomeS_part_1 ${THISWORKFOLDER}
cp -R ${HOMEFOLDER}/example_genomes ${THISWORKFOLDER}

cd ${THISWORKFOLDER}/PopCoGenomeS_part_1

#run Popcogenomes
sh 01_PopCOGenomeS_array.sh


rsync -a --no-p ${THISWORKFOLDER}/PopCoGenomeS_part_1/./output_${n} ${HOMEFOLDER}





