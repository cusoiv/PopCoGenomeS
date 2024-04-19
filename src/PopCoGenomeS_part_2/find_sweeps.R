library(phangorn)
library(castor)
library(dplyr)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
basename=args[1]
working_dir=args[2]
setwd(working_dir)
folder_list_all = read.table(paste0(basename,'_cf_size_3.list'))
tree_list=lapply(folder_list_all$V1, function(x) paste0(x,'.core.fasta.phy_phyml_tree.txt'))

#Read in length bias file
lbfile=read.table(paste0(basename,'.length_bias_500.txt'),sep='\t',header = T)

#filter
lbfile_1=lbfile%>% 
  mutate(totalR=case_when(sim_fr/hmm_fr>=2 & (mnb_div/mu_div>2.5 | Initial.divergence.iter1 < 1500/Alignment.size.iter1) ~ hmm_fr, 
                          TRUE ~ sim_fr),
         div=case_when(sim_fr/hmm_fr>2 & (mnb_div/mu_div>2.5 | Initial.divergence.iter1 < 1500/Alignment.size.iter1) ~ hmm_mu, 
                       sim_fr/hmm_fr<2 & (mnb_div/mu_div>2.5 | Initial.divergence.iter1 < 1500/Alignment.size.iter1) ~ mu_div,
                       TRUE ~ Initial.divergence.iter1),
         type=case_when(sim_fr/hmm_fr>2 ~ 'C', TRUE ~ 'NC')) %>%
  mutate(totalR=case_when(Initial.divergence.iter1 < 1500/Alignment.size.iter1 ~ 0, TRUE ~ totalR),
         div=case_when(Initial.divergence.iter1 < 1500/Alignment.size.iter1 ~ 10^-5, TRUE ~ div)) 

#loop through trees for all vertically-inherited clusters 
cc_list={}
cc_sister_list={}
for (i in 1:nrow(folder_list_all)){
    tree_file=tree_list[[i]]
      if (file.size(tree_file) != 0L){

        t=read_tree(file=tree_file)
        t=midpoint(t)
        t$node.label=1:t$Nnode+Ntip(t)
        clade_list=get_clade_list(t)$clades
        ntips=length(t$tip.label)
        ratio_df=data.frame(node=ntips+1,ratio=NA, cfratio=NA)
        dist_list=list()
        for (n in 2:t$Nnode){ #because 1 is always root, and we are only looking at nodes with sister nodes
          
          # find last ancestral node
          an=get_ancestral_nodes(t,ntips+n, Nsplits=1)
          
          # find sister node
          sisters=clade_list[ntips+an,c(2,3)]
          sister_n=sisters[sisters!=(ntips+n)]
          
          # if sister node is a node:
          if (sister_n>ntips){
            
            # find all tip to node distance for sister node
            sister_st=get_subtree_at_node(t, sister_n-ntips)$subtree
            all_distances_sister_st = get_all_distances_to_root(sister_st)
            tip_distances_sister_st=all_distances_sister_st[1:Ntip(sister_st)]
            tips_sister_st=sister_st$tip.label
            
            # find all tip to node distance for node
            st=get_subtree_at_node(t, n)$subtree
            all_distances_st = get_all_distances_to_root(st)
            tip_distances_st=all_distances_st[1:Ntip(st)]
            
            # find all pairwise distances between tips for node
            all_tip_pairs_st_sq=get_all_pairwise_distances(st,only_clades = st$tip.label)
            all_tip_pairs_st=all_tip_pairs_st_sq[lower.tri(all_tip_pairs_st_sq)]
            tips_st=st$tip.label
            all_tip_pairs_st_sq=as.data.frame(all_tip_pairs_st_sq)
            names(all_tip_pairs_st_sq)=tips_st
            all_tip_pairs_st_sq$id=tips_st
            all_tip_pairs_st_melt=melt(all_tip_pairs_st_sq) 
            names(all_tip_pairs_st_melt)=c('Strain.1','Strain.2','dist')
            all_tip_pairs_st_melt=all_tip_pairs_st_melt %>% filter(Strain.1!=Strain.2)
            
            # all combinations of tip distances 
            all_diff_clade_dist=expand.grid(tip_distances_sister_st,tip_distances_st)
            all_diff_clade_tips=expand.grid(tips_sister_st,tips_st)
            names(all_diff_clade_tips)=c('Strain.1','Strain.2')
            ad=get_pairwise_distances(t,sister_n,ntips+n)  #add in distance between nodes
            all_diff_clade_dist$dist=all_diff_clade_dist$Var1+all_diff_clade_dist$Var2+ad
            
            all_diff_clade_dist_2=cbind(all_diff_clade_tips,all_diff_clade_dist$dist)
            names(all_diff_clade_dist_2)=c('Strain.1','Strain.2','dist')
            all_diff_clade_dist_1=all_diff_clade_dist_2
            all_diff_clade_dist_1$Strain.1=all_diff_clade_dist_2$Strain.2
            all_diff_clade_dist_1$Strain.2=all_diff_clade_dist_2$Strain.1
            all_diff_clade_dist_full=rbind(all_diff_clade_dist_2,all_diff_clade_dist_1)
            all_diff_clade_dist_full$cc="diff_cc"
            all_tip_pairs_st_melt$cc="same_cc"
            all_clade_dist_full=rbind(all_diff_clade_dist_full,all_tip_pairs_st_melt)
            all_clade_dist_full=inner_join(all_clade_dist_full, lbfile_1, 
                                           by=c("Strain.1","Strain.2"))
            all_clade_dist_full_cfsummary=all_clade_dist_full %>% group_by(cc) %>% summarise(meandiv=mean(div))
          }else{
            # find sister tip
            tips_sister_st=t$tip.label[sister_n]
            
            # find all tip to node distance for node
            st=get_subtree_at_node(t, n)$subtree
            all_distances_st = get_all_distances_to_root(st)
            tip_distances_st=all_distances_st[1:Ntip(st)]
            
            # find all pairwise distances between tips for node
            all_tip_pairs_st_sq=get_all_pairwise_distances(st,only_clades = st$tip.label)
            all_tip_pairs_st=all_tip_pairs_st_sq[lower.tri(all_tip_pairs_st_sq)]
            tips_st=st$tip.label
            all_tip_pairs_st_sq=as.data.frame(all_tip_pairs_st_sq)
            names(all_tip_pairs_st_sq)=tips_st
            all_tip_pairs_st_sq$id=tips_st
            all_tip_pairs_st_melt=melt(all_tip_pairs_st_sq) 
            names(all_tip_pairs_st_melt)=c('Strain.1','Strain.2','dist')
            all_tip_pairs_st_melt=all_tip_pairs_st_melt %>% filter(Strain.1!=Strain.2)
            
            # all combinations of tip distances 
            ad=get_pairwise_distances(t,sister_n,ntips+n) #sister tip to node distance
            all_diff_clade_dist=expand.grid(ad,tip_distances_st)
            all_diff_clade_tips=expand.grid(tips_sister_st,tips_st)
            all_diff_clade_dist$dist=all_diff_clade_dist$Var1+all_diff_clade_dist$Var2
            
            all_diff_clade_dist_2=cbind(all_diff_clade_tips,all_diff_clade_dist$dist)
            names(all_diff_clade_dist_2)=c('Strain.1','Strain.2','dist')
            all_diff_clade_dist_1=all_diff_clade_dist_2
            all_diff_clade_dist_1$Strain.1=all_diff_clade_dist_2$Strain.2
            all_diff_clade_dist_1$Strain.2=all_diff_clade_dist_2$Strain.1
            all_diff_clade_dist_full=rbind(all_diff_clade_dist_2,all_diff_clade_dist_1)
            all_diff_clade_dist_full$cc="diff_cc"
            all_tip_pairs_st_melt$cc="same_cc"
            all_clade_dist_full=rbind(all_diff_clade_dist_full,all_tip_pairs_st_melt)
            all_clade_dist_full=inner_join(all_clade_dist_full, lbfile_1, 
                                           by=c("Strain.1","Strain.2"))
            all_clade_dist_full_cfsummary=all_clade_dist_full %>% group_by(cc) %>% summarise(meandiv=mean(div))
          }
          
          
          # find mean distance ratio
          ratio_df[n,]=c(ntips+n,mean(all_diff_clade_dist$dist)/mean(all_tip_pairs_st),
                         all_clade_dist_full_cfsummary$meandiv[1]/all_clade_dist_full_cfsummary$meandiv[2])
          
          dist_list[[n]]=all_clade_dist_full
        }
        
        t$node.label=ratio_df$cfratio>=5
        t$node.label[t$node.label==FALSE]=NA
        t$node.label[t$node.label==TRUE]='S'
        ratio_df$sweep=t$node.label
        sweep_node=ratio_df$node[ratio_df$cfratio>=5]
        sweep_node=sweep_node[!is.na(sweep_node)]
        
        
        an=get_ancestral_nodes(t,ntips+n, Nsplits=1)
        # find sister node
        sisters=clade_list[ntips+an,c(2,3)]
        sister_n=sisters[sisters!=(ntips+n)]
        sweep_tips=lapply(sweep_node, function(x){data.frame(sweep_genome=
          get_subtree_at_node(t, x-Ntip(t))$subtree$tip.label)
        })
        sister_sweep_tips=lapply(sweep_node, function(x) {
          an=get_ancestral_nodes(t,x, Nsplits=1) # find last ancestral node
          sisters=clade_list[ntips+an,c(2,3)]
          sister_n=sisters[sisters!=x] # find sister node
          if (sister_n>ntips){
            data.frame(sister_of_sweep=get_subtree_at_node(t, sister_n-ntips)$subtree$tip.label)
          }else{
            data.frame(sister_of_sweep=t$tip.label[sister_n])
          }
        })
      #Is this clonal frame directly a sweep?
      strains=t$tip.label
      lb_in=lbfile_1 %>% filter(Strain.1 %in% strains & Strain.2 %in% strains) %>% select(div,totalR,Strain.1,Strain.2)

      lb_not_in=lbfile_1 %>% filter(Strain.1 %in% strains & ! (Strain.2 %in% strains) |
                                  !(Strain.1 %in% strains) &  Strain.2 %in% strains) %>% select(div,totalR,Strain.1,Strain.2)
      ratio=min(lb_not_in$div)/mean(lb_in$div)
      
      }
        if (length(sweep_tips)>0){
          sweep_tips_2=Filter(function(x) nrow(x) >= 3,sweep_tips)
          cc_list[[i]]=bind_rows(sweep_tips_2,.id='sweep_id')
          names(cc_list)[i]=folder_list_all$V1[i]
        }
        if (ratio >=5){
          sweep_df=data.frame(sweep_id='cf',sweep_genome=strains)
          cc_list[[i]]=bind_rows(cc_list[[i]],sweep_df)
          names(cc_list)[i]=folder_list_all$V1[i]
        }
  }

cc_list_df=bind_rows(cc_list,.id='vertical_cluster_name')

write.table(cc_list_df,paste0(basename,'_sweeps.txt'),quote = F,row.names = F,sep='\t')


