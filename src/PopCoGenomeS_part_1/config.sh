# Base name for final output files ust a prefix to identify your outputs.
base_name='example'

# Output directory for the final output files.

# This will create the directory if it does not already exist.
final_output_dir=./output/
mkdir -p ${final_output_dir}

# Path to genome files. (make sure it ends in a "/")
genome_dir=./genomes/

# Genome file filename extension.
genome_ext=.fna

# Are you running on a single machine? Please specify the number of threads to run.
# This can, at maximum, be the number of logical cores your machine has.
num_threads=1

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

# hmm cutoff for performing hmm
hmm_cutoff=0.1

# cutoff for cluster of strains to be considered to share a clonal frame, ranging from 0-1
clonal_cutoff=0.5
