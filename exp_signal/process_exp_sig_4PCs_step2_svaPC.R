library(matrixStats)
library(WGCNA)
library(SummarizedExperiment)
source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")
source("/users/ywang/10_25_2019/functions/functions_evaluation_Aug_5.R")
dir_data="/users/ywang/10_25_2019/data/"
dir_input="/users/ywang/10_25_2019/corplot_ori_svapcs/data/"
dir_output="/users/ywang/10_25_2019/exp_signal_4PC/data/"
  
  
list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
               "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")

#########################################################################################################
### 4PCs
# n_PC=4
# for(i in 1:9){
#   ######load data and filtration
#   # load counts data
# 
#   name_tissue=list_tissues[[i]]
#   
#   # load and input raw counts, normalize and filter low expression genes(average(log2cpm)>0)
#   name_var= paste0("gtex.",name_tissue)
#   load(paste0(dir_data, name_var,".rda" ))
#   data_raw=eval(as.name(name_var))
#   counts_raw=assay(data_raw)
#   counts_raw=counts_raw[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA")),]
#   
#   list_filt=rpkm_filt(counts_raw,rowData(data_raw)$gene_length[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))])
#   log2rpkm_kp=list_filt$counts_log2rpkm_keep
# 
#   if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
#   set.seed(1)
#   
#   corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
#   ave_logrpkm_kp=rowMeans(log2rpkm_kp)
#   
#   
#   ind_up=which(upper.tri(corr_ori),arr.ind = T)
#   n_up=nrow(ind_up)  
#   exp_min=ave_logrpkm_kp[ind_up[,1]]
#   exp_max=ave_logrpkm_kp[ind_up[,2]]
#   
#   order_cor_ori=order(abs(corr_ori)[ind_up],decreasing =TRUE)
#   remove(corr_ori)
#   
#   load(paste0(dir_input,"cor_est_",n_PC,"_PCs_",i,"th_tissue.RData"))
#   
#   order_cor_est=order(abs(cor_est)[ind_up],decreasing=TRUE)
#   remove(cor_est)
#   
#   list_nsignal=c(round(n_up/500000),round(n_up/100000),round(n_up/10000),round(n_up/2000),round(n_up/1000),round(n_up/500),round(n_up/100),round(n_up/20))
#   percent_list=c(1/500000*100,1/100000*100,1/10000*100,1/2000*100,1/1000*100,1/500*100,1/100*100,1/20*100)
#   percent_list2=c("0_0002","0_001","0_01","0_05", "0_1","0_2", "1", "5")
# 
#   k = length(list_nsignal)-1 # 1%
#   nsignal=list_nsignal[k]
#   percent_=percent_list2[k]
#   loc_signal_ori=order_cor_ori[1:nsignal]
#   exp_signal_m1m2_ori=ave_logrpkm_kp[ind_up[loc_signal_ori,]]
#   exp_signal_min_ori=exp_min[loc_signal_ori]
#   exp_signal_max_ori=exp_max[loc_signal_ori]
#   
#   loc_signal_est=order_cor_est[1:nsignal]
#   exp_signal_m1m2_est=ave_logrpkm_kp[ind_up[loc_signal_est,]]
#   exp_signal_min_est=exp_min[loc_signal_est]
#   exp_signal_max_est=exp_max[loc_signal_est]
#   
#   save(exp_signal_max_ori,file=paste0(dir_output,"exp_signal_max_ori_",n_PC,"PCs_",i,".RData"))  
#   save(exp_signal_min_ori,file=paste0(dir_output,"exp_signal_min_ori_",n_PC,"PCs_",i,".RData"))  
#   save(exp_signal_m1m2_ori,file=paste0(dir_output,"exp_signal_m1m2_ori_",n_PC,"PCs_",i,".RData"))  
#   
#   save(exp_signal_min_est,file=paste0(dir_output,"exp_signal_min_est_",n_PC,"PCs_",i,".RData"))  
#   save(exp_signal_max_est,file=paste0(dir_output,"exp_signal_max_est_",n_PC,"PCs_",i,".RData"))  
#   save(exp_signal_m1m2_est,file=paste0(dir_output,"exp_signal_m1m2_est_",n_PC,"PCs_",i,".RData"))  
# }
# 
# 
# for(i in 1:9){
#   name_tissue=list_tissues[[i]]
#   
#   # load and input raw counts, normalize and filter low expression genes(average(log2cpm)>0)
#   name_var= paste0("gtex.",name_tissue)
#   load(paste0(dir_data, name_var,".rda" ))
#   data_raw=eval(as.name(name_var))
#   counts_raw=assay(data_raw)
#   counts_raw=counts_raw[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA")),]
#   
#   list_filt=rpkm_filt(counts_raw,rowData(data_raw)$gene_length[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))])
#   log2rpkm_kp=list_filt$counts_log2rpkm_keep
#   
#   if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
#   set.seed(1)
#   
#   corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
#   ave_logrpkm_kp=rowMeans(log2rpkm_kp)
#   
#   ind_up=which(upper.tri(corr_ori),arr.ind = T)
#   
#     
#   n_up=nrow(ind_up)  
#   exp_min=ave_logrpkm_kp[ind_up[,1]]
#   exp_max=ave_logrpkm_kp[ind_up[,2]]
# 
#   load(paste0(dir_output,"exp_signal_max_ori_",n_PC,"PCs_",i,".RData"))  
#   load(paste0(dir_output,"exp_signal_min_ori_",n_PC,"PCs_",i,".RData"))  
#   load(paste0(dir_output,"exp_signal_m1m2_ori_",n_PC,"PCs_",i,".RData"))  
# 
#   load(paste0(dir_output,"exp_signal_min_est_",n_PC,"PCs_",i,".RData"))  
#   load(paste0(dir_output,"exp_signal_max_est_",n_PC,"PCs_",i,".RData"))  
#   load(paste0(dir_output,"exp_signal_m1m2_est_",n_PC,"PCs_",i,".RData"))  
#   
#   nsignal=round(length(exp_signal_max_ori)/10) #0.1%
# 
#   exp_signal_max_ori=exp_signal_max_ori[1:nsignal]
#   exp_signal_min_ori=exp_signal_min_ori[1:nsignal]
#   exp_signal_max_est=exp_signal_max_est[1:nsignal]
#   exp_signal_min_est=exp_signal_min_est[1:nsignal]
# 
#   d1=density((exp_signal_max_ori+exp_signal_min_ori)/2)
#   d2=density((exp_signal_max_est+exp_signal_min_est)/2)
#   d0=density((exp_max+exp_min)/2)
#   
#   percent_="0_1"
# 
#   save(d1,file=paste0(dir_output,percent_,"_percent_",n_PC,"PCs_",i,"_d1.RData"))
#   save(d0,file=paste0(dir_output,percent_,"_percent_",n_PC,"PCs_",i,"_d0.RData"))
#   save(d2,file=paste0(dir_output,percent_,"_percent_",n_PC,"PCs_",i,"_d2.RData"))
# }

#########################################################################################################
### sva_PCs
for(i in 1:9){
  ######load data and filtration
  # load counts data

  name_tissue=list_tissues[[i]]
  dir_figure = paste0("/users/ywang/10_25_2019/corplot_adj_4_svapcs/fig2/")
  load(paste0(dir_figure,"nsv_",i,"_th_tissue.RData"))
  n_PC=nsv

  # load and input raw counts, normalize and filter low expression genes(average(log2cpm)>0)
  name_var= paste0("gtex.",name_tissue)
  load(paste0(dir_data, name_var,".rda" ))
  data_raw=eval(as.name(name_var))
  counts_raw=assay(data_raw)
  counts_raw=counts_raw[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA")),]

  list_filt=rpkm_filt(counts_raw,rowData(data_raw)$gene_length[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))])
  log2rpkm_kp=list_filt$counts_log2rpkm_keep

  if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
  set.seed(1)

  corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
  ave_logrpkm_kp=rowMeans(log2rpkm_kp)[order(rowMeans(log2rpkm_kp))]

  gene_type_filt1=rowData(data_raw)$gene_type[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))]
  gene_type_filt12=gene_type_filt1[list_filt$keep][order(rowMeans(log2rpkm_kp))]
  
  corr_ori=corr_ori[which(gene_type_filt12=="protein_coding"),which(gene_type_filt12=="protein_coding")]
  ave_logrpkm_kp=ave_logrpkm_kp[which(gene_type_filt12=="protein_coding")]
  
  ind_up=which(upper.tri(corr_ori),arr.ind = T)
  n_up=nrow(ind_up)
  exp_min=ave_logrpkm_kp[ind_up[,1]]
  exp_max=ave_logrpkm_kp[ind_up[,2]]

  order_cor_ori=order(abs(corr_ori)[ind_up],decreasing=TRUE)
  remove(corr_ori)

  load(paste0(dir_input,"cor_est_",n_PC,"_PCs_",i,"th_tissue.RData"))

  cor_est=cor_est[which(gene_type_filt12=="protein_coding"),which(gene_type_filt12=="protein_coding")]

  order_cor_est=order(abs(cor_est)[ind_up],decreasing=TRUE)
  remove(cor_est)

  list_nsignal=c(round(n_up/500000),round(n_up/100000),round(n_up/10000),round(n_up/2000),round(n_up/1000),round(n_up/500),round(n_up/100),round(n_up/20))
  percent_list=c(1/500000*100,1/100000*100,1/10000*100,1/2000*100,1/1000*100,1/500*100,1/100*100,1/20*100)
  percent_list2=c("0_0002","0_001","0_01","0_05", "0_1","0_2", "1", "5")

  k = length(list_nsignal)-3 # 0.1%
  nsignal=list_nsignal[k]
  percent_=percent_list2[k]
  loc_signal_ori=order_cor_ori[1:nsignal]
  exp_signal_m1m2_ori=ave_logrpkm_kp[ind_up[loc_signal_ori,]]
  exp_signal_min_ori=exp_min[loc_signal_ori]
  exp_signal_max_ori=exp_max[loc_signal_ori]

  loc_signal_est=order_cor_est[1:nsignal]
  exp_signal_m1m2_est=ave_logrpkm_kp[ind_up[loc_signal_est,]]
  exp_signal_min_est=exp_min[loc_signal_est]
  exp_signal_max_est=exp_max[loc_signal_est]

  save(exp_signal_max_ori,file=paste0(dir_output,"exp_signal_max_ori_",n_PC,"PCs_",i,".RData"))
  save(exp_signal_min_ori,file=paste0(dir_output,"exp_signal_min_ori_",n_PC,"PCs_",i,".RData"))
  save(exp_signal_m1m2_ori,file=paste0(dir_output,"exp_signal_m1m2_ori_",n_PC,"PCs_",i,".RData"))

  save(exp_signal_min_est,file=paste0(dir_output,"exp_signal_min_est_",n_PC,"PCs_",i,".RData"))
  save(exp_signal_max_est,file=paste0(dir_output,"exp_signal_max_est_",n_PC,"PCs_",i,".RData"))
  save(exp_signal_m1m2_est,file=paste0(dir_output,"exp_signal_m1m2_est_",n_PC,"PCs_",i,".RData"))
}



for(i in 1:9){
  name_tissue=list_tissues[[i]]
  dir_figure = paste0("/users/ywang/10_25_2019/corplot_adj_4_svapcs/fig2/")
  load(paste0(dir_figure,"nsv_",i,"_th_tissue.RData"))
  n_PC=nsv
  
  # load and input raw counts, normalize and filter low expression genes(average(log2cpm)>0)
  name_var= paste0("gtex.",name_tissue)
  load(paste0(dir_data, name_var,".rda" ))
  data_raw=eval(as.name(name_var))
  counts_raw=assay(data_raw)
  counts_raw=counts_raw[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA")),]
  
  list_filt=rpkm_filt(counts_raw,rowData(data_raw)$gene_length[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))])
  log2rpkm_kp=list_filt$counts_log2rpkm_keep
  
  if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
  set.seed(1)
  
  corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
  ave_logrpkm_kp=rowMeans(log2rpkm_kp)[order(rowMeans(log2rpkm_kp))]
  gene_type_filt1=rowData(data_raw)$gene_type[which(rowData(data_raw)$gene_type %in% c("protein_coding", "lincRNA"))]
  gene_type_filt12=gene_type_filt1[list_filt$keep][order(rowMeans(log2rpkm_kp))]
  
  corr_ori=corr_ori[which(gene_type_filt12=="protein_coding"),which(gene_type_filt12=="protein_coding")]
  ave_logrpkm_kp=ave_logrpkm_kp[which(gene_type_filt12=="protein_coding")]
  
  ind_up=which(upper.tri(corr_ori),arr.ind = T)
  
  n_up=nrow(ind_up)  
  
  
  list_nsignal=c(round(n_up/500000),round(n_up/100000),round(n_up/10000),round(n_up/2000),round(n_up/1000),round(n_up/500),round(n_up/100),round(n_up/20))
  percent_list=c(1/500000*100,1/100000*100,1/10000*100,1/2000*100,1/1000*100,1/500*100,1/100*100,1/20*100)
  percent_list2=c("0_0002","0_001","0_01","0_05", "0_1","0_2", "1", "5")
  k = length(list_nsignal)-3 # 0.1%
  nsignal=list_nsignal[k]
  
  
  
  exp_min=ave_logrpkm_kp[ind_up[,1]]
  exp_max=ave_logrpkm_kp[ind_up[,2]]
  
  load(paste0(dir_output,"exp_signal_max_ori_",n_PC,"PCs_",i,".RData"))  
  load(paste0(dir_output,"exp_signal_min_ori_",n_PC,"PCs_",i,".RData"))  
  load(paste0(dir_output,"exp_signal_m1m2_ori_",n_PC,"PCs_",i,".RData"))  
  
  load(paste0(dir_output,"exp_signal_min_est_",n_PC,"PCs_",i,".RData"))  
  load(paste0(dir_output,"exp_signal_max_est_",n_PC,"PCs_",i,".RData"))  
  load(paste0(dir_output,"exp_signal_m1m2_est_",n_PC,"PCs_",i,".RData"))  
  
  nsignal=round(length(exp_signal_max_ori)) #0.1%
  
  exp_signal_max_ori=exp_signal_max_ori[1:nsignal]
  exp_signal_min_ori=exp_signal_min_ori[1:nsignal]
  exp_signal_max_est=exp_signal_max_est[1:nsignal]
  exp_signal_min_est=exp_signal_min_est[1:nsignal]
  
  d1=density((exp_signal_max_ori+exp_signal_min_ori)/2)
  d2=density((exp_signal_max_est+exp_signal_min_est)/2)
  d0=density((exp_max+exp_min)/2)
  
  percent_=percent_list2[k]
  
  save(d0,file=paste0(dir_output,percent_,"_percent_",n_PC,"PCs_",i,"_d0.RData"))
  save(d1,file=paste0(dir_output,percent_,"_percent_",n_PC,"PCs_",i,"_d1.RData"))
  save(d2,file=paste0(dir_output,percent_,"_percent_",n_PC,"PCs_",i,"_d2.RData"))
}













