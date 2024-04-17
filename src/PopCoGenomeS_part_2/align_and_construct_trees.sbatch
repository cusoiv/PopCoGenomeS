#!/bin/bash
#
#SBATCH --job-name=align_and_construct_trees
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB
#SBATCH --nice=0
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

#usage
#sbatch -a 1-number of chunks align_and_construct_trees.sbatch <CONFIGFILE>

PHY_CONFIG=./phybreak_config.sh
source ${PHY_CONFIG}

# The configfile is the list of clonal frames
CONFIGFILE=$1
OUTFOLDER_2=`awk '{if (NR=='${SLURM_ARRAY_TASK_ID}') print $1}' $CONFIGFILE`

#define directory where the script is called
HOMEFOLDER=`pwd`

#define and construct workspace
THISWORKFOLDER=${project_dir_1}/${OUTFOLDER_2}
mkdir -p $THISWORKFOLDER
mkdir -p $THISWORKFOLDER/genomes
cp -R ${path_to_phybreak}/align_and_construct_trees $THISWORKFOLDER
cp ${PHY_CONFIG} $THISWORKFOLDER/align_and_construct_trees

project_dir=$THISWORKFOLDER
export project_dir
export OUTFOLDER_2

# make the output directory
mkdir -p ${output_dir}
echo ${output_dir}

cd ${project_dir}/align_and_construct_trees

#Run prep code
conda activate PopCoGenomeS
sh 00_prepare_phybreak.sh
echo Done!

#Generate multiple sequence alignment 
python phybreak1.generate_maf.py
python phybreak2.maf_to_fasta.py

#Run the tree building step
conda activate PopCoGenomeS_R
Rscript fasta_to_phy.R $basename
cd ${project_dir}/align
#This tree uses a GTR+G+I model
phyml -i ${basename}.core.fasta.phy -m 'GTR' -t 'e' -a 'e' -f 'm' -v 'e' > ${basename}.core.phy_phyml_stat.txt

#copy tree back to output folder
cp ${basename}.core.fasta.phy_phyml_tree.txt $output_dir/${OUTFOLDER_2}.core.fasta.phy_phyml_tree.txt


