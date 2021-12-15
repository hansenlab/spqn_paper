
source("/functions/functions_2d_quantile_no_surf_July_18.R")

source("/functions/functions_evaluation_Aug_5.R")
library(WGCNA)
library(matrixStats)
library("sva")
library(SummarizedExperiment)

# choose tissue name and output folder name
i=1


dir_figure = paste0("/corplot_adj_4_svapcs/fig2/")
load(paste0(dir_figure,"nsv_",i,"_th_tissue.RData"))
n_PC=nsv

list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
               "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")
name_tissue=list_tissues[[i]]
dir_outpu = paste0("/users/ywang/10_25_2019/corplot_ori_svapcs/data/")
dir_data="/users/ywang/10_25_2019/data/"


# load and input raw counts, normalize and filter low expression genes(average(log2cpm)>0)
name_var= paste0("gtex.",name_tissue)
load(paste0(dir_data, name_var,".rda" ))
data_raw=eval(as.name(name_var))
counts_raw=assay(data_raw)
counts_raw=counts_raw[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA")),]



list_filt=rpkm_filt(counts_raw,rowData(data_raw)$gene_length[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))])
log2rpkm_kp=list_filt$counts_log2rpkm_keep

# mod=matrix(1,nrow=dim(log2rpkm_kp)[2],ncol=1)
# colnames(mod)="Intercept"
# nsv=num.sv(log2rpkm_kp,mod, method = "be") ## num.sv requires data matrix with features(genes) in the rows and samples in the column
# print(paste("Number of PCs estimated to be removed:", nsv))
# save(nsv,file=paste0(dir_figure,"nsv_",i,"_th_tissue.RData"))

# gtex.rse.sub <- variable.selection.average(gtex.rse, num.variables)

# pc.estimate = function(rse.object, less.pc = F, alt.pc = F, ...){
#   dat <- t(SummarizedExperiment::assay(rse.object, 1))
#   
#   ## determine number of principal components to adjust for
#   mod = matrix(1, nrow = dim(dat)[1], ncol = 1)
#   colnames(mod) = "Intercept"
#   n.pc <- num.sv(t(dat), mod,method="be")
#   n.pc
# }


# 
if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
set.seed(1)



corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
ave_logrpkm_kp=rowMeans(log2rpkm_kp)





########################################################################
#quantile normzalization

##### step 1, asssign outer bins
group_loc=get_grps(ave_logrpkm_kp,ngrp=60,size_grp=400)
group_loc=get_grps(ave_logrpkm_kp,ngrp=12,size_grp=30)
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
save(cor_est,file=paste0(dir_outpu,"cor_est_",n_PC,"_PCs_",i,"th_tissue.RData"))



















