library(matrixStats)
library(WGCNA)
library(ggplot2)
library(ggridges)
library(SummarizedExperiment)
source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")
source("/users/ywang/10_25_2019/functions/functions_2d_quantile_diagnal_only.R")

source("/users/ywang/10_25_2019/functions/functions_evaluation_Aug_5.R")
#library(Rfast)#Sort/Rank
library(WGCNA)
library(matrixStats)
library("sva")
library(SummarizedExperiment)

source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")
source("/users/ywang/10_25_2019/functions/functions_evaluation_Aug_5.R")

dir_data=paste0("/users/ywang/10_25_2019/scRNA/data/GSE45719_RAW/")
dir_inpu = paste0("/users/ywang/10_25_2019/scRNA/data/GSE45719_RAW/")
dir_output="/users/ywang/10_25_2019/scRNA/data/"
dir_fig="/users/ywang/10_25_2019/scRNA/fig/"

#########################################################################################################
#########################################################################################################

######load data and filtration
# load counts data
load(paste0(dir_output,"rpkm_midblast.RData"))#rpkm_midblast
# dim(rpkm_midblast )
# [1] 22958    60
logrpkm_all=log2(rpkm_midblast+0.5)#not_sure
keep=which(rowMedians(data.matrix(logrpkm_all))>0)
log2rpkm_kp=logrpkm_all[keep,]
# nrow(log2rpkm_kp)
# [1] 6915
# set.seed(1)

# # 
n_PC=16
if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
ave_logrpkm_kp=rowMeans(log2rpkm_kp)[order(rowMeans(log2rpkm_kp))]


##### step 1, asssign outer bins
group_loc=get_grps(ave_logrpkm_kp,ngrp=60,size_grp=400)
##### step 2, asssign inner bins
group_loc_adj=get_grps_inner(group_loc)
##### step 3, get rank for each running bin
group_loc_ref_label=cut(1:length(ave_logrpkm_kp),10)
group_loc_ref=split(1:length(ave_logrpkm_kp),group_loc_ref_label)
cor_ref=corr_ori[group_loc_ref[[9]],group_loc_ref[[9]]]
#cor_ref=corr_ori[group_loc[[18]],group_loc[[18]]]
rank_bin=get_bin_rank2(corr_ori,group_loc,group_loc_adj,cor_ref) #6s for 20*20 bins, 4000 genes; 2min for all genes,ngrp=20,size_grp=1000
remove(corr_ori)
##### step 4, transform rank to cor_est
cor_est=est_cor2(rank_bin,cor_ref)#6s for 20*20 bins, 4000 genes; 3.5min for all genes, 20*20 bins,siez_grp=1000
remove(rank_bin)
save(cor_est,file=paste0(dir_output,"cor_est_",n_PC,"_PCs.RData"))




  
