# PopCoGenomeS part 2
This part performs tree construction for each vertically-inherited cluster and eventually outputs a list of sweeps according to the 5X rule.
# Usage

## Step 1 
### Performs genome alignment for each vertically-inherited genome cluster, and constructs a GTR+G+I tree for each vertically-inherited genome cluster.

1. editing the phybreak_config.sh file according to the instructions in the file. Also make sure that phybreak_config.sh is in the same folder as 
align_and_construct_trees.sbatch or align_and_construct_trees.sh.

2. if you have access to a slurm computing cluster, you can run
sbatch -a 1-${number of vertically-inherited genome clusters} align_and_construct_trees.sbatch ${basename}_cf_size_3.list

The ${basename}_cf_size_3.list is the output from the previous part 1 that consists of the names of all vertically-inherited genome clusters.

if you don't have access to a slurm computing cluster, you can run this as :

while read line; do   
	sh align_and_construct_trees.sh ${line}  
done<${basename}_cf_size_3.list

## Step 2 
### Applies the 5X rule to all vertically-inherited genome clusters and its subclusters

This requires the output trees from Step 1, as well as the ${basename}_cf_size_3.list and ${basename}_length_bias_${window_size}.txt files to all be present in the same folder

Rscript find_sweeps.R ${basename} ${path_to_folder_with_files}

This should return the output ${basename}_sweeps.txt
