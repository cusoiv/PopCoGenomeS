# Part 1
This part calculates the vertically-inherited fraction for each pair of genomes and then generates clusters of genomes that are vertically inherited at a user's defined cutoff. 

# Usage
Edit config.sh according to the instructions in the file, and then run PopCoGenomeS.sh

# Input
A directory of genomes in fasta format.

# Expected outputs
This should generate five groups of files:
## 1. ${basename}.length_bias_${window_size}.txt  
   This is the data file containing information for all pairwise genomes.  
   What each column means in this file:  
   Strain 1 : genome name 1 in focal pair of genomes  
   Strain 2 : genome name 2 in focal pair of genomes  
   Initial divergence raw: Divergence between genome 1 and genome 2, directly calculated from mugsy alignment  
   Initial divergence: Divergence between genome 1 and genome 2, after preliminary removal of regions with more than expected SNPs  
   Alignment size raw: Size of genome alignment, directly calculated from musgy alignment  
   Alignment size: Size of genome alignment, after preliminary removal of regions with more than expected SNPs  
   Genome 1 size: Size of genome 1  
   Genome 2 size: Size of genome 2  
   Observed SSD: Length bias calculated as in previous PopCOGenT paper  
   SSD 95 CI low: 95% confidence interval of length bias, lower end  
   SSD 95 CI high: 95% confidence interval of length bias, higher end  
   mu_div: Divergence of the vertically inherited fraction, estimated by maximum likelihood estimation (MLE)  
   mnb_div: Divergence of the recombined fraction, estimated by MLE  
   sim_fr: Percent genome recombined, estimated by MLE  
   nb_size: alpha or size factor of negative binomial distribution, estimated by maximum likelihood estimation  
   hmm_fr: Percent genome recombined, estimated by hidden markov models  
   hmm_mu: Divergence of the vertically inherited fraction, estimated by HMM  
   hmm_mnb: Divergence of the recombined fraction, estimated by HMM  
   
## 2. ${basename}.length_bias.filtered.txt  
   This is the filtered datafile for all pairwise genomes, which makes three major changes compared to the complete length bias file:  
   
   a. Pairwise genomes with the mean divergence of the recombined fraction being no more than 2.5X of the vertically-inherited fraction are omitted.     
   b. Pairwise genomes that diverge by <1500 SNPs are counted as 100% vertically inherited with divergence being 10^-5.   
   c. We cross-checked whether the MLE-estimated recombination fraction exceeded twice that determined by the HMM. If such a discrepancy occurred, we  substituted the MLE-estimated parameters with those derived from the HMM, and marked the change as 'C' in the type column.  
   
   Compared to the previous length bias file, there are three extra columns:  
   
   totalR: Percent genome recombined, after correction  
   div: Divergence of the vertically inherited fraction, after correction  
   type: Whether there was a discrepency in the HMM and MLE based estimations  

## 3. ${basename}_${clonal_cutoff}.cluster.tab.txt  
   This is a table that assigns all the isolate genomes to a corresponding vertically-inherited cluster.  
## 4. ${basename}_cf_size_3.list      
   The list of all vertically-inherited genome clusters with >=3 isolate genomes to be inputed to Part 2.  
## 5. ${basename}_XXX.txt
   These are the strains in each of the vertically-inherited genome clusters in the previous list.  
   


