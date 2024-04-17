library(phylotools, lib.loc="~/Rlibs_4")
args = commandArgs(trailingOnly=TRUE)
fasta=read.fasta(file = paste0("../align/",args[1],".core.fasta"))
dat2phylip(fasta, outfile = paste0("../align/",args[1],".core.fasta.phy"))







