library(matrixStats)
library(WGCNA)
library(ggplot2)
library(ggridges)
library(SummarizedExperiment)

source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")
source("/users/ywang/10_25_2019/functions/functions_evaluation_Aug_5.R")


dir_data="/users/ywang/10_25_2019/data/"
# dir_output="/users/ywang/10_25_2019/boxplot_9tissues/data/"
# dir_input="/users/ywang/10_25_2019/corplot_ori_svapcs/data/"

# dir_input="/users/ywang/July_28_2020/regulatome/output/
dir_output="/users/ywang/July_28_2020/VST/output_4pc/"
setwd(dir_output)
dir_data="/dcl01/hansen/data/meanCoexp/"
n_PC=4

#########################################################################################################
#########################################################################################################

for(i in 1:9){
  print(i)
#   ######load data and filtration
#   # load counts data
#   
  list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
                 "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")
  #tissue_list=c( "adiposesub", "adrenalgl","arteryti","braince","brainco","breastma",
  #                 "colontr","esophagusmu","heartleft")
  name_tissue=list_tissues[[i]]

  #load and input raw counts, normalize and filter low expression genes(average(log2cpm)>0)
  name_var= paste0("gtex.",name_tissue)
  load(paste0(dir_data, name_var,".rda" ))
  data_raw=eval(as.name(name_var))
  counts_raw=assay(data_raw)
  counts_raw=counts_raw[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA")),]

  list_filt=rpkm_filt(counts_raw,rowData(data_raw)$gene_length[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))])
  log2rpkm_kp=list_filt$counts_log2rpkm_keep

  
  ####### get cor - without vst 
  if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
  corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
  ave_logrpkm_kp=rowMeans(log2rpkm_kp)[order(rowMeans(log2rpkm_kp))]
  
  list_IQR_ori = cal_m_sd_IQR(corr_ori,ave_logrpkm_kp)
  save(list_IQR_ori,file=paste0(dir_output,"list_IQR_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))
  
  
#   #########################################################################################################
#   #########################################################################################################
  ### svaPCs

  # load(paste0("/users/ywang/10_25_2019/corplot_adj_4_svapcs/fig2/","nsv_",i,"_th_tissue.RData"))
  # n_PC=nsv
  # 
  # # calculate cor_ori
  # if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
  # corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
  # ave_logrpm_kp=rowMeans(log2rpkm_kp)[order(rowMeans(log2rpkm_kp))]
  # #########################################################################################################
  # ##calculate IQR - svaPCs
  # list_IQR_ori = cal_m_sd_IQR(corr_ori,ave_logrpm_kp)
  # remove(corr_ori)
  # save(list_IQR_ori,file=paste0(dir_output,"list_IQR_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))
  # #
  # load(paste0(dir_input,"cor_est_",n_PC,"_PCs_",i,"th_tissue.RData"))
  # list_IQR_est = cal_m_sd_IQR(cor_est,ave_logrpm_kp)
  # remove(cor_est)
  # save(list_IQR_est,file=paste0(dir_output,"list_IQR_est_",n_PC,"_PCs_",i,"th_tissue.RData"))

    
#   #########################################################################################################
#   #########################################################################################################
#   ### 4PCs
#   
  # calculate cor_ori - 4PCs

  
  ##########################################################################################################
  ##########################################################################################################
  
  ####### get cor -vst 
  # load vst-normalized counts
  load(paste0("counts_raw_keep_vst_",name_tissue,".RData"))
  # save(counts_raw_keep_vst,file=paste0("counts_raw_keep_vst_",name_tissue,".RData"))

  list_filt_vst=rpkm_filt(counts_raw_keep_vst,rowData(data_raw)$gene_length[rownames(counts_raw_keep_vst)])
  log2rpkm_kp_vst=list_filt_vst$counts_log2rpkm_keep
  log2rpkm_kp_rmPCt_vst = removePCs(log2rpkm_kp_vst,n_PC)#remove 4PC from the expressed genes in this tissue
  cor_ori_vst=cor(t(log2rpkm_kp_rmPCt_vst[order(rowMeans(log2rpkm_kp)),]))
  list_IQR_ori_vst = cal_m_sd_IQR(cor_ori_vst,ave_logrpkm_kp)
  save(list_IQR_ori_vst,file=paste0(dir_output,"list_IQR_ori_",n_PC,"_PCs_vst_",i,"th_tissue.RData"))
  # 
  ##########################################################################################################
  ##########################################################################################################
  
  # load(paste0(dir_input,"cor_est_",n_PC,"_PCs_",i,"th_tissue.RData"))
  # list_IQR_est = cal_m_sd_IQR(cor_est,ave_logrpm_kp)
  # remove(cor_est)
  # save(list_IQR_est,file=paste0(dir_output,"list_IQR_est_",n_PC,"_PCs_",i,"th_tissue.RData"))
# 
#   #########################################################################################################
#   #########################################################################################################
#   ### 0PCs
# 
#   # calculate cor_ori - 0PC
#   n_PC=0
# 
#   # calculate cor_ori
#   if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
#   corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
#   ave_logrpkm_kp=rowMeans(log2rpkm_kp)[order(rowMeans(log2rpkm_kp))]
#   #########################################################################################################
#   ##calculate IQR - svaPCs
#   list_IQR_ori = cal_m_sd_IQR(corr_ori,ave_logrpkm_kp)
#   remove(corr_ori)
#   save(list_IQR_ori,file=paste0(dir_output,"list_IQR_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))
# 
    # print(i)
# 
}


#########################################################################################################
#########################################################################################################
#########################################################################################################
#######format IQR for ggplot input
#########################################################################################################
####### 100 IQRs - sva_PCs
#IQR_compare_sva_PC_trend.pdf

# f_g_m2=function(grp_median){
#   grp_median_m=array(dim=c(10,10))
#   for(i in 1:10){
#     for(j in 1:10){
#       grp_median_m[i,j]=grp_median[min(i,j)]
#     }
#   }
#   grp_median_m
# }
# 
# IQR_ori_all=IQR_adj_all=wid_box=grp_mean_all=nPC_all=c()
# for(i in 1:9){
#   #n_PC=list_nsv[i]
# 
#   load(paste0("/users/ywang/10_25_2019/corplot_adj_4_svapcs/fig2/","nsv_",i,"_th_tissue.RData"))
#   n_PC=nsv
#   load(paste0(dir_output,"list_IQR_est_",n_PC,"_PCs_",i,"th_tissue.RData"))
#   load(paste0(dir_output,"list_IQR_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))
# 
#   IQR_grp_ori=list_IQR_ori$IQR_cor_mat_m
#   IQR_grp_adj=list_IQR_est$IQR_cor_mat_m
# 
#   IQR_ori_all=c(IQR_ori_all,as.vector(IQR_grp_ori))
#   IQR_adj_all=c(IQR_adj_all,as.vector(IQR_grp_adj))
# 
#   grp_mean=f_g_m2(list_IQR_ori$grp_mean)
#   grp_mean_all=c(grp_mean_all,as.vector(grp_mean))
#   #wid_box=c(wid_box,IQR(IQR_grp_ori))
# }
# 
# IQR_all=c(IQR_ori_all,IQR_adj_all)
# group=c(rep("original",900),rep("adjusted",900))
# tissue=rep(rep(list_tissues,each=100),2)
# d_quantile=data.frame(IQR_all,group,tissue)
# d_quantile$tissue=as.character(d_quantile$tissue)
# d_quantile$group <- factor(d_quantile$group, levels = c("original","adjusted"))
# d_quantile$tissue=as.character(d_quantile$tissue)
# d_quantile$IQR_all=IQR_all
# d_quantile$expression_level <- grp_mean_all
# col_=rep(c("#00BFC4","#F8766D"),length(grp_mean_all)/2)
# d_quantile$col=col_
# #fill_=rep(c("#00BFC4","#F8766D"),length(grp_mean_all)/2)
# ylim_=range(IQR_all)
# save(d_quantile,file=paste0(dir_output,"d_quantile_svaPCs.RData"))




#########################################################################################################
#######plot 100 IQRs - 4_PCs
# f_g_m2=function(grp_median){
#   grp_median_m=array(dim=c(10,10))
#   for(i in 1:10){
#     for(j in 1:10){
#       grp_median_m[i,j]=grp_median[min(i,j)]
#     }
#   }
#   grp_median_m
# }
# 
# list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
#                "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")
# tissue_list=c( "adiposesub", "adrenalgl","arteryti","braince","brainco","breastma",
#                  "colontr","esophagusmu","heartleft")
# IQR_ori_all=IQR_adj_all=wid_box=grp_mean_all=nPC_all=c()
# for(i in 1:9){
#   n_PC=4
#   load(paste0(dir_output,"list_IQR_est_",n_PC,"_PCs_",i,"th_tissue.RData"))
#   load(paste0(dir_output,"list_IQR_ori_",n_PC,"_PCs_",i,"th_tissue.RData"))
# 
#   IQR_grp_ori=list_IQR_ori$IQR_cor_mat_m
#   IQR_grp_adj=list_IQR_est$IQR_cor_mat_m
# 
#   IQR_ori_all=c(IQR_ori_all,as.vector(IQR_grp_ori))
#   IQR_adj_all=c(IQR_adj_all,as.vector(IQR_grp_adj))
# 
#   grp_mean=f_g_m2(list_IQR_ori$grp_mean)
#   grp_mean_all=c(grp_mean_all,as.vector(grp_mean))
#   #wid_box=c(wid_box,IQR(IQR_grp_ori))
# }
# 
# IQR_all=c(IQR_ori_all,IQR_adj_all)
# group=c(rep("original",900),rep("adjusted",900))
# tissue=rep(rep(list_tissues,each=100),2)
# d_quantile=data.frame(IQR_all,group,tissue)
# d_quantile$tissue=as.character(d_quantile$tissue)
# d_quantile$group <- factor(d_quantile$group, levels = c("original","adjusted"))
# d_quantile$tissue=as.character(d_quantile$tissue)
# d_quantile$IQR_all=IQR_all
# d_quantile$expression_level <- grp_mean_all
# col_=rep(c("#00BFC4","#F8766D"),length(grp_mean_all)/2)
# d_quantile$col=col_
# #fill_=rep(c("#00BFC4","#F8766D"),length(grp_mean_all)/2)
# ylim_=range(IQR_all)
# save(d_quantile,file=paste0(dir_output,"d_quantile_4PCs.RData"))
# 
# 
# 
# 
