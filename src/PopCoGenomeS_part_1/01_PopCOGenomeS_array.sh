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
python get_alignment_and_length_bias_select.py --genome_dir ${genome_dir} --genome_ext ${genome_ext} --alignment_dir ${alignment_dir} --base_name ${basename} --window_size ${window_size} --final_output_dir ${output_dir} --num_threads ${num_threads} --chunks ${num_chunks} --chunk_id ${chunk_id} ${keep_alignments}

