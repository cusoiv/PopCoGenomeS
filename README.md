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

# Dependencies
## System requirement: 
A linux-based system 

## Setting up work environment:

The required python and R packages can be installed by creating two conda environments with the included PopCOGenomeS.yml & PopCOGenomeS_R.yml files as follows:
    
    conda env create -f PopCoGenomeS.yml
    conda env create -f PopCoGenomeS_R.yml

The first conda environment consists of all python modules involved, as well as the tools mugsy and seqkit.
The second conda enviroment consists of all R libraries involved, as well as the tool phyml.

The PopCoGenomeS.yml and PopCoGenomeS_R.yml files reflect the versions we used for our publication and at this time we cannot guarantee forward compatability with new versions of the dependencies.

We have come to notice that sometimes there is a problem with the mugsy that it is not using the perl version of the conda environment but the system perl, which can sometime lead to issues. If you run into this issue, edit the nucmer and mummerplot files (scripts invoked by mugsy) in the conda environment (under bin/MUMmer3.20) by changing '#!/usr/bin/perl' to '#!/usr/bin/env perl' in those two scripts. 


# Usage

This pipeline is suitable with isolate genomes with N50>50K. We provide a workflow consisting of two modules, which allows the user to start from a folder of genomes (in fasta format) to a list of identified genome-wide sweeps.

Part 1: Estimation of the recombined fractions and vertically inherited fractions of pairwise genomes & Clustering genomes into (mostly) vertically inherited genome clusters.

This part is in /src/PopCOGenomeS_part_1. The file to run is PopCOGenomeS.sh. Please follow the instructions to edit the config.sh file as requested.

Part 2: Identifying microbial populations as genome-wide selective sweeps.

This part is in /src/PopCOGenomeS_part_2. Please follow the instructions to edit the phybreak_config.sh file as requested.

