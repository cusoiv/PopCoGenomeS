#!/bin/bash

THISDIR=${BASH_SOURCE[0]}
THISDIR=${THISDIR%/PopCoGenomeS.sh}
echo $THISDIR

configfile=./config.sh
source ${configfile}

module load conda
conda activate ${path_to_PopCoGenomeS}

export window_size

#Get pairwise length bias and recombined fraction data
python get_alignment_and_length_bias_select.py --genome_dir ${genome_dir} --genome_ext ${genome_ext} --alignment_dir ${alignment_dir} --base_name ${basename} --window_size ${window_size} --final_output_dir ${final_output_dir} --num_threads ${num_threads} ${keep_alignments}


#Cluster into clusters with >clonal_cutoff clonal fraction

conda activate ${path_to_PopCoGenomeS_R}
Rscript cluster.R ${basename} ${clonal_cutoff} ${window_size} ${final_output_dir}


