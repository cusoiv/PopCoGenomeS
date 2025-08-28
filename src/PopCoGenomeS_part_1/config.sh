# Base name for final output files ust a prefix to identify your outputs.
basename='example'

# Output directory for the final output files.

# This will create the directory if it does not already exist.
# you only need to edit the parameter output_dir. The final_output_dir will be directly used as the output folder if you are running on single machine format, while outputs will be in output_dir for the 01 step if you are running in slurm array format. The 02 step will concatenate all outputs in the output_dir to generate the main clusters.

final_output_dir=./output # make sure it doesn't end in a / if you are running array mode

output_dir="${final_output_dir}${n:+_$n}"
output_dirs=${final_output_dir}_*
mkdir -p ${output_dir}
mkdir -p ${final_output_dir}

# Path to genome files. (make sure it ends in a "/")
genome_dir=../example_genomes/

# Genome file filename extension.
genome_ext=.fna

# Are you running on a single machine? Please specify the number of threads to run.
# This can, at maximum, be the number of logical cores your machine has.
num_threads=4

# Are you running on slurm and want to split the job between different nodes? 
#Please specifify how many different array jobs you want to split your job into.
num_chunks=4
#Please specific how many cores each array job uses.
num_threads=4
#This is the n th chunk that is run and should be same as the array task number
chunk_id=${n}

# Whether to keep alignments after length bias is calculated. 
# Alignment files can be 10MB each and thus a run on 100 genomes can take up on the order of 50 GB of space if alignment files are not discarded. 
# If you want to keep alignments, set to --keep_alignments. Otherwise leave as ''.
#keep_alignments=--keep_alignments
keep_alignments=''

# Directory for output alignments. Must provide absolute path.
alignment_dir=$(pwd)/proc/
mkdir -p ${alignment_dir}

# window size to use
window_size=500

# cutoff for a cluster of strains to be considered as vertically inherited, ranging from 0-1
clonal_cutoff=0.5

# Path to the PopCoGenomeS conda environment
path_to_PopCoGenomeS='/lisc/home/user/yu/envs/popcogenomes'

# Path to the PopCoGenomeS_R conda environment
path_to_PopCoGenomeS_R='/lisc/home/user/yu/envs/popcogenomes_R'
