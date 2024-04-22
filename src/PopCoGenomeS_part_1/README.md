# Part 1
This part calculates the vertically-inherited fraction for each pair of genomes and then generates clusters of genomes that are vertically inherited at a user's defined cutoff. 

# Usage
Edit config.sh according to the instructions in the file, and then run PopCoGenomeS.sh

# Input
A directory of genomes in fasta format.

# Expected outputs
This should generate four files:
1. ${basename}_length_bias_${window_size}.txt  This is the data file containing information for all pairwise genomes
   What each column means in this file:
   
   
3. ${basename}_length_bias_filtered.txt  This is the filtered datafile for all pairwise genomes, which makes three major changes
   a. Pairwise genomes with the mean divergence of the recombined fraction being no more than 2.5X of the vertically-inherited fraction are omitted 
   b. Pairwise genomes that diverge by <1500 SNPs are counted as 100% vertically inherited with divergence being 10^-5.
   c. We cross-checked whether the MLE-estimated recombination fraction exceeded twice that determined by the HMM. If such a discrepancy occurred, we  substituted the MLE-estimated parameters with those derived from the HMM, and marked the change as 'C' in the type column.
4. ${basename}_${clonal_cutoff}.cluster.tab.txt This is a table that assigns all the isolate genomes to a corresponding vertically-inherited cluster 
5. ${basename}_cf_size_3.list The list of all vertically-inherited genome clusters with >=3 isolate genomes to be inputed to Part 2
   


