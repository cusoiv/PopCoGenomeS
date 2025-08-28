#!/bin/bash
#
#SBATCH --job-name=submit_02_popcogenomes
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --nice=0
#SBATCH --partition=basic
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.err
#SBATCH --time=10:00:00

#usage
#no need for config file since the folders are numbered
#sbatch -a 1-number of chunks ../Scripts/00_popcogent.sh 

window_size=500
export window_size 

#define directory where the script is called
HOMEFOLDER=`pwd`

#define and construct workspace on /tmp/$USER
THISWORKFOLDER=$TMPDIR

#move previous output and popcogenomes code
cp -R ${HOMEFOLDER}/PopCoGenomeS_part_1 ${THISWORKFOLDER}
cp -R ${HOMEFOLDER}/output_* ${THISWORKFOLDER}/PopCoGenomeS_part_1

cd ${THISWORKFOLDER}/PopCoGenomeS_part_1

#run clustering step
sh 02_PopCOGenomeS_cluster.sh

#rsync back to local folder
rsync -a --no-p ${THISWORKFOLDER}/PopCoGenomeS_part_1/./output ${HOMEFOLDER}





