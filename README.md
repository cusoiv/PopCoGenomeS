# PopCoGenomeS
This is the github repository for the tool PopCoGenomeS
# Purpose
Identifying microbial Populations as Clusters Of Genome Sweeps (PopCoGenomeS). This is built on the basis of the first part of the tool PopCOGenT, but focuses on:
1. Estimation of the recombined fractions and vertically inherited fractions of pairwise genomes.
2. Clustering genomes into (mostly) vertically inherited genome clusters.
3. Identifying microbial populations as genome-wide selective sweeps.

PopCoGenomeS contrasts with the tool PopCOGenT in that PopCOGenT identifies populations as units of gene flow, thus assuming high amounts of recombination within populations, and very low vertical inheritance within populations; meanwhile, PopCoGenomeS identifies genome clusters where the majority of the genome is vertically inherited.

For more detailed information on PopCOGenT, see:
https://github.com/philarevalo/PopCOGenT/tree/master

# Suitability
We suggest applying PopCoGenomeS to a collection of genomes from the same species. Since the identification of genome-wide selective sweeps could be confounded by repeated sampling of the same person/host, we suggest that the group of genomes to apply this tool be dereplicated according to host ids. Also, as this tool relies on pairwise genome alignments, we caution against applying the tool to tens of thousands of genomes simultaneously. It is suggested to not apply the tool to >300 genomes at once. If you want to quickly check if there are possible genome-wide sweeps in a very large dataset, we suggest randomly subsampling the collection of genomes down to ~200 as a start. 

# Dependencies
## System requirement: 
A linux-based system. The typical installation time for the conda installation should take between 30mins - 1h. 

## Setting up work environment:

The required python and R packages can be installed by creating two conda environments with the included PopCOGenomeS.yml & PopCOGenomeS_R.yml files as follows:
    
    conda env create -f PopCoGenomeS.yml
    conda env create -f PopCoGenomeS_R.yml

The first conda environment consists of all python modules involved, as well as the tools mugsy and seqkit.
The second conda enviroment consists of all R libraries involved, as well as the tool phyml.

The PopCoGenomeS.yml and PopCoGenomeS_R.yml files reflect the versions we used for our publication and at this time we cannot guarantee forward compatability with new versions of the dependencies.

We have come to notice that sometimes there is a problem with the mugsy that it is not using the perl version of the conda environment but the system perl, which can sometime lead to issues. If you run into this issue, edit the nucmer and mummerplot files (scripts invoked by mugsy) in the conda environment (under bin/MUMmer3.20) by changing '#!/usr/bin/perl' to '#!/usr/bin/env perl' in those two scripts. 


# Usage

This pipeline is suitable with isolate genomes with N50>50K. We provide a workflow consisting of two modules, which allows the user to start from a folder of genomes (in fasta format) to a list of identified genome-wide sweeps. We provide an example folder and an example output folder for you to test and run the pipeline. 

Part 1: Estimation of the recombined fractions and vertically inherited fractions of pairwise genomes & Clustering genomes into (mostly) vertically inherited genome clusters.

This part is in /src/PopCOGenomeS_part_1. The file to run is PopCOGenomeS.sh. Please follow the instructions to edit the config.sh file as requested. Depending on the number of cores assigned, the run time of the genomes in the example folder will vary. But it should take a few hours to run on a 4 core machine.

For larger number of genomes (i.e. >100), we also provide a slurm array based version that will allow you to run the code on multiple nodes on a slurm controlled cluster. Please see the instructions in the README of /src/PopCOGenomeS_part_1.

Part 2: Identifying microbial populations as genome-wide selective sweeps.

This part is in /src/PopCOGenomeS_part_2. Please follow the instructions to edit the phybreak_config.sh file as requested. This should take less than 30 mins for the example dataset.

