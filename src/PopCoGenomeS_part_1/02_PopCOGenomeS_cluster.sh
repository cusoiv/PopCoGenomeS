#!/bin/bash

THISDIR=${BASH_SOURCE[0]}
THISDIR=${THISDIR%/PopCoGenomeS.sh}
echo $THISDIR

configfile=./config.sh
source ${configfile}

export window_size
awk 'FNR==1 && NR!=1 { next } { print }' ${output_dirs}/${basename}.length_bias_${window_size}.txt > ${final_output_dir}/${basename}.length_bias_${window_size}.txt

#Cluster into clusters with >clonal_cutoff clonal fraction
module load conda
conda activate ${path_to_PopCoGenomeS_R}
Rscript cluster.R ${basename} ${clonal_cutoff} ${window_size} ${final_output_dir}


