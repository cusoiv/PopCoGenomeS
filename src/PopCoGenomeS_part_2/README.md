# PopCoGenomeS part 2
This part performs tree construction for each vertically-inherited cluster and eventually outputs a list of sweeps according to the 5X rule.
# Usage
The first part performs genome alignment for each vertically-inherited genome cluster, and constructs a GTR+G+I tree for each vertically-inherited genome cluster.

1. editing the phybreak_config.sh file according to the instructions in the file. Also make sure that phybreak_config.sh is in the same folder as 
align_and_construct_trees.sbatch or align_and_construct_trees.sh.

2. if you have access to a slurm computing cluster, you can run
sbatch -a 1-${number of vertically-inherited genome clusters} align_and_construct_trees.sbatch ${basename}_cf_size_3.list

The ${basename}_cf_size_3.list is the output from the previous part 1 that consists of the names of all vertically-inherited genome clusters.

if you don't have access to a slurm computing cluster, you can run this as :

while read line; do
        sh align_and_construct_trees.sh $line
<${basename}_cf_size_3.list

3. 
