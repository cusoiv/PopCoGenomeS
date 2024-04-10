library(dplyr)
library(micropan)

#setwd("C:/Users/Xiaoqian/Desktop/pop_gen/PopcoGenomeS/")
args = commandArgs(trailingOnly=TRUE)
basename=args[1]
clonal_cutoff=args[2]
window_size=args[3]
final_output_dir=args[4]

#Read in length bias file
lbfile=read.table(paste0(final_output_dir,basename,'.length_bias_',window_size,'.txt'),sep='\t',header = T,stringsAsFactors = F)

#Two types of filtering

lbfile_filter=lbfile %>% 
  filter(mnb_div/mu_div>2.5 | Initial.divergence.iter1 < 1500/Alignment.size.iter1) %>%
  mutate(totalR=sim_fr,div=mu_div) %>%
  mutate(totalR=case_when(sim_fr/hmm_fr>2 ~ hmm_fr, TRUE ~ sim_fr),
         div=case_when(sim_fr/hmm_fr>2 ~ hmm_mu, TRUE ~ mu_div),
         type=case_when(sim_fr/hmm_fr>2 ~ 'C', TRUE ~ 'NC')) %>%
  mutate(totalR=case_when(Initial.divergence.iter1 < 1500/Alignment.size.iter1  ~ 0, TRUE ~ totalR),
         div=case_when(Initial.divergence.iter1 < 1500/Alignment.size.iter1 ~ 10^-5, TRUE ~ div))

lbfile_filter_1=lbfile %>% 
  mutate(totalR=case_when(sim_fr/hmm_fr>=2 & (mnb_div/mu_div>2.5 | Initial.divergence.iter1 < 1500/Alignment.size.iter1) ~ hmm_fr, 
                          TRUE ~ sim_fr),
         div=case_when(sim_fr/hmm_fr>2 & (mnb_div/mu_div>2.5 | Initial.divergence.iter1 < 1500/Alignment.size.iter1) ~ hmm_mu, 
                       sim_fr/hmm_fr<2 & (mnb_div/mu_div>2.5 | Initial.divergence.iter1 < 1500/Alignment.size.iter1) ~ mu_div,
                       TRUE ~ Initial.divergence.iter1),
         type=case_when(sim_fr/hmm_fr>2 ~ 'C', TRUE ~ 'NC')) %>%
  mutate(totalR=case_when(Initial.divergence.iter1 < 1500/Alignment.size.iter1 ~ 0, TRUE ~ totalR),
         div=case_when(Initial.divergence.iter1 < 1500/Alignment.size.iter1 ~ 10^-5, TRUE ~ div)) 

#Select and prepare for clustering
lbfile_filter_select=lbfile_filter[,c("Strain.1","Strain.2","totalR")]
names(lbfile_filter_select)=c('Dbase','Query','Distance')
lbfile_filter_select$Distance=lbfile_filter_select$Distance/100

#cluster with average clustering
if (min(lbfile_filter_select$Distance)<=0.5){
  z=as.data.frame(bClust(lbfile_filter_select,linkage = "average",threshold = clonal_cutoff))
  names(z)='cluster'
  z$strain=rownames(z)
  z2=split(z,z$cluster)
  z3=Filter(function(x) nrow(x) > 1, z2)
  lbfile_filter_list_cl=z3
  lbfile_filter_list_all_cl=z2
}

#this step is to clean the clusters so that some of the points that come after the decline is not included
#basically, points within the cluster that have pairwise divergence larger than the mean divergence of points outside of the cluster are omitted


if (!is.null(lbfile_filter_list_cl)){
    lb=lbfile_filter_1
    total_strains=unique(c(lb$Strain.1,lb$Strain.2))
    lb_in_list=lapply(lbfile_filter_list_cl,function(x) {
      lb_in=lb %>% filter(Strain.1 %in% x$strain & Strain.2 %in% x$strain) %>% select(div,totalR,Strain.1,Strain.2); lb_in})
	print (lb_in_list)
    lb_not_in_list=lapply(lbfile_filter_list_cl,function(x) {
      lb_not_in=lb %>% filter(Strain.1 %in% x$strain & ! (Strain.2 %in% x$strain) |
                                !(Strain.1 %in% x$strain) &  Strain.2 %in% x$strain) %>% select(div,totalR,Strain.1,Strain.2);lb_not_in})
	print (lb_not_in_list)
    lb_in_list_clean=mapply(function(x,y) {if (nrow(y)>0) {x=x %>% filter(x$div < mean(y$div))};x},
                            lb_in_list,lb_not_in_list, SIMPLIFY = F)
    lbfile_filter_list_cl_all_clean=mapply(function(x,z){x %>% filter(x$strain %in% unique(c(z$Strain.1,z$Strain.2)))},
                   lbfile_filter_list_cl,lb_in_list_clean, SIMPLIFY = F)
    lbfile_filter_list_cl_all_clean_bind=bind_rows(lbfile_filter_list_cl_all_clean)
    omitted_strains=total_strains[!(total_strains %in% lbfile_filter_list_cl_all_clean_bind$strain)]
    #print (omitted_strains)
    #print (max(lbfile_filter_list_cl_all_clean_bind$cluster))
    if (length(omitted_strains)>0) {
    omitted_strains_df=data.frame(cluster=(max(lbfile_filter_list_cl_all_clean_bind$cluster)+1):(max(lbfile_filter_list_cl_all_clean_bind$cluster)+length(omitted_strains)),
                                  strain=omitted_strains)
    lbfile_filter_list_cl_all_clean_bind=bind_rows(lbfile_filter_list_cl_all_clean_bind,omitted_strains_df)
    cluster_mapping = data.frame(cluster = unique(lbfile_filter_list_cl_all_clean_bind$cluster)) %>%
      mutate(new_cluster = row_number())
    lbfile_filter_list_cl_all_clean_bind=left_join(lbfile_filter_list_cl_all_clean_bind, cluster_mapping, by = "cluster")
    lbfile_filter_list_cl_all_clean_bind = lbfile_filter_list_cl_all_clean_bind %>%
      select(-cluster) %>%
      rename(cluster = new_cluster)}else{
	lbfile_filter_list_cl_all_clean_bind=bind_rows(lbfile_filter_list_all_cl) %>% select(strain,cluster)}
    lbfile_filter_list_cl_all_clean_list=split(lbfile_filter_list_cl_all_clean_bind,f=lbfile_filter_list_cl_all_clean_bind$cluster)
    lbfile_filter_list_cl_clean_list=Filter(function(x) nrow(x) >= 3,  lbfile_filter_list_cl_all_clean_list)
    
}

# write out files

#1. The filtered length bias file 
colnames(lbfile_filter) <- gsub("\\.", " ", colnames(lbfile_filter))
write.table(lbfile_filter,paste0(final_output_dir, basename,'.length_bias.filtered.txt'), sep="\t", quote = FALSE, row.names = FALSE)
#2. The cluster file
write.table(lbfile_filter_list_cl_all_clean_bind,paste0(final_output_dir,basename,'_',clonal_cutoff,'.cluster.tab.txt'), sep="\t", quote = FALSE, row.names = FALSE)
#3. Write out clonal frames with >= 3 strains
sweep_names=c()
for (i in 1:length(lbfile_filter_list_cl_clean_list)){
  x=as.numeric(names(lbfile_filter_list_cl_clean_list)[i])
  write.table(lbfile_filter_list_cl_clean_list[[i]]$strain,paste0(final_output_dir,basename,'_cf_',sprintf("%03d", x),'.txt'),
              quote = FALSE, row.names = FALSE,col.names = FALSE)
  sweep_names=c(sweep_names,paste0(basename,'_cf_',sprintf("%03d", x)))
}
#4. Write out list of sweeps
write.table(sweep_names,paste0(final_output_dir, basename,'_cf_size_3.list'),quote = FALSE, row.names = FALSE,col.names = FALSE)



