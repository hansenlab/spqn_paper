library(matrixStats)
library(WGCNA)
library(ggplot2)
library(ggridges)
library(SummarizedExperiment)

source("/users/ywang/10_25_2019/functions/functions_2d_quantile_no_surf_July_18.R")
source("/users/ywang/10_25_2019/functions/functions_evaluation_Aug_5.R")

dir_data=paste0("/users/ywang/10_25_2019/scRNA/data/GSE45719_RAW/")
dir_inpu = paste0("/users/ywang/10_25_2019/scRNA/data/GSE45719_RAW/")
dir_output="/users/ywang/10_25_2019/scRNA/data/"
dir_fig="/users/ywang/10_25_2019/scRNA/fig/"
dir.create(dir_fig)
  
list_tissues=c("Adipose_Subcutaneous","Adrenal_Gland","Artery_Tibial","Brain_Cerebellum",
               "Brain_Cortex", "Breast_Mammary","Colon_Transverse","Esophagus_Mucosa","Heart_Left_Ventricle")

#########################################################################################################
### 4PCs
# n_PC=4
#   ######load data and filtration
#   # load counts data
# i=1
# load(paste0(dir_output,"rpkm_midblast.RData"))#rpkm_midblast
# # dim(rpkm_midblast )
# # [1] 22958    60
# logrpkm_all=log2(rpkm_midblast+0.5)#not_sure
# keep=which(rowMedians(data.matrix(logrpkm_all))>0)
# log2rpkm_kp=logrpkm_all[keep,]
# 
# 
#   if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
# 
#   corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
#   ave_logrpkm_kp=rowMeans(log2rpkm_kp)[order(rowMeans(log2rpkm_kp))]
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
#   # load(paste0(dir_input,"cor_est_",n_PC,"_PCs_",i,"th_tissue.RData"))
#   #
#   # order_cor_est=order(abs(cor_est)[ind_up],decreasing=TRUE)
#   # remove(cor_est)
# 
#   list_nsignal=c(round(n_up/500000),round(n_up/100000),round(n_up/10000),round(n_up/2000),round(n_up/1000),round(n_up/500),round(n_up/100),round(n_up/20))
#   percent_list=c(1/500000*100,1/100000*100,1/10000*100,1/2000*100,1/1000*100,1/500*100,1/100*100,1/20*100)
#   percent_list2=c("0_0002","0_001","0_01","0_05", "0_1","0_2", "1", "5")
# 
#   k = length(list_nsignal)-3 # 0.1%
#   nsignal=list_nsignal[k]
#   percent_=percent_list2[k]
#   loc_signal_ori=order_cor_ori[1:nsignal]
#   exp_signal_m1m2_ori=ave_logrpkm_kp[ind_up[loc_signal_ori,]]
#   exp_signal_min_ori=exp_min[loc_signal_ori]
#   exp_signal_max_ori=exp_max[loc_signal_ori]
# 
#   # loc_signal_est=order_cor_est[1:nsignal]
#   # exp_signal_m1m2_est=ave_logrpkm_kp[ind_up[loc_signal_est,]]
#   # exp_signal_min_est=exp_min[loc_signal_est]
#   # exp_signal_max_est=exp_max[loc_signal_est]
# 
#   save(exp_signal_max_ori,file=paste0(dir_output,"exp_signal_max_ori_",n_PC,"PCs_",i,".RData"))
#   save(exp_signal_min_ori,file=paste0(dir_output,"exp_signal_min_ori_",n_PC,"PCs_",i,".RData"))
#   save(exp_signal_m1m2_ori,file=paste0(dir_output,"exp_signal_m1m2_ori_",n_PC,"PCs_",i,".RData"))
# #
# #   save(exp_signal_min_est,file=paste0(dir_output,"exp_signal_min_est_",n_PC,"PCs_",i,".RData"))
# #   save(exp_signal_max_est,file=paste0(dir_output,"exp_signal_max_est_",n_PC,"PCs_",i,".RData"))
# #   save(exp_signal_m1m2_est,file=paste0(dir_output,"exp_signal_m1m2_est_",n_PC,"PCs_",i,".RData"))
# # 
# 
# 
# for(i in 1){
#   n_PC=4
#   ######load data and filtration
#   # load counts data
#   i=1
#   load(paste0(dir_output,"rpkm_midblast.RData"))#rpkm_midblast
#   # dim(rpkm_midblast )
#   # [1] 22958    60
#   logrpkm_all=log2(rpkm_midblast+0.5)#not_sure
#   keep=which(rowMedians(data.matrix(logrpkm_all))>0)
#   log2rpkm_kp=logrpkm_all[keep,]
#   
#   
#   if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
#   
#   corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
#   ave_logrpkm_kp=rowMeans(log2rpkm_kp)[order(rowMeans(log2rpkm_kp))]
#   
#   
#   ind_up=which(upper.tri(corr_ori),arr.ind = T)
#   
#     
#   n_up=nrow(ind_up)  
#   
#   list_nsignal=c(round(n_up/500000),round(n_up/100000),round(n_up/10000),round(n_up/2000),round(n_up/1000),round(n_up/500),round(n_up/100),round(n_up/20))
#   percent_list=c(1/500000*100,1/100000*100,1/10000*100,1/2000*100,1/1000*100,1/500*100,1/100*100,1/20*100)
#   percent_list2=c("0_0002","0_001","0_01","0_05", "0_1","0_2", "1", "5")
#   k = length(list_nsignal)-3 # 0.1%
#   nsignal=list_nsignal[k]
#   
#   exp_min=ave_logrpkm_kp[ind_up[,1]]
#   exp_max=ave_logrpkm_kp[ind_up[,2]]
# 
#   load(paste0(dir_output,"exp_signal_max_ori_",n_PC,"PCs_",i,".RData"))  
#   load(paste0(dir_output,"exp_signal_min_ori_",n_PC,"PCs_",i,".RData"))  
#   load(paste0(dir_output,"exp_signal_m1m2_ori_",n_PC,"PCs_",i,".RData"))  
# # 
# #   load(paste0(dir_output,"exp_signal_min_est_",n_PC,"PCs_",i,".RData"))  
# #   load(paste0(dir_output,"exp_signal_max_est_",n_PC,"PCs_",i,".RData"))  
# #   load(paste0(dir_output,"exp_signal_m1m2_est_",n_PC,"PCs_",i,".RData"))  
#   
#   nsignal=round(length(exp_signal_max_ori)) #0.1%
# 
#   exp_signal_max_ori=exp_signal_max_ori[1:nsignal]
#   exp_signal_min_ori=exp_signal_min_ori[1:nsignal]
#   # exp_signal_max_est=exp_signal_max_est[1:nsignal]
#   # exp_signal_min_est=exp_signal_min_est[1:nsignal]
# 
#   d1=density((exp_signal_max_ori+exp_signal_min_ori)/2)
#   # d2=density((exp_signal_max_est+exp_signal_min_est)/2)
#   d0=density((exp_max+exp_min)/2)
#   
#   percent_=percent_list2[k]
#   
#   save(d1,file=paste0(dir_output,percent_,"_percent_",n_PC,"PCs_",i,"_d1.RData"))
#   save(d0,file=paste0(dir_output,percent_,"_percent_",n_PC,"PCs_",i,"_d0.RData"))
#   # save(d2,file=paste0(dir_output,percent_,"_percent_",n_PC,"PCs_",i,"_d2.RData"))
# }

#########################################################################################################
### sva_PCs
  n_PC=16
  ######load data and filtration
  # load counts data
  i=1
  load(paste0(dir_output,"rpkm_midblast.RData"))#rpkm_midblast
  # dim(rpkm_midblast )
  # [1] 22958    60
  logrpkm_all=log2(rpkm_midblast+0.5)#not_sure
  keep=which(rowMedians(data.matrix(logrpkm_all))>0)
  log2rpkm_kp=logrpkm_all[keep,]
  
  
  if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
  
  corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
  ave_logrpkm_kp=rowMeans(log2rpkm_kp)[order(rowMeans(log2rpkm_kp))]
  
  
  ind_up=which(upper.tri(corr_ori),arr.ind = T)
  n_up=nrow(ind_up)
  exp_min=ave_logrpkm_kp[ind_up[,1]]
  exp_max=ave_logrpkm_kp[ind_up[,2]]
  
  order_cor_ori=order(abs(corr_ori)[ind_up],decreasing =TRUE)
  remove(corr_ori)
  
  load(paste0(dir_input,"cor_est_",n_PC,"_PCs.RData"))

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
  #
    save(exp_signal_min_est,file=paste0(dir_output,"exp_signal_min_est_",n_PC,"PCs_",i,".RData"))
    save(exp_signal_max_est,file=paste0(dir_output,"exp_signal_max_est_",n_PC,"PCs_",i,".RData"))
    save(exp_signal_m1m2_est,file=paste0(dir_output,"exp_signal_m1m2_est_",n_PC,"PCs_",i,".RData"))
  # 
  
  
  for(i in 1){
    n_PC=16
    ######load data and filtration
    # load counts data
    i=1
    load(paste0(dir_output,"rpkm_midblast.RData"))#rpkm_midblast
    # dim(rpkm_midblast )
    # [1] 22958    60
    logrpkm_all=log2(rpkm_midblast+0.5)#not_sure
    keep=which(rowMedians(data.matrix(logrpkm_all))>0)
    log2rpkm_kp=logrpkm_all[keep,]
    
    
    if(n_PC==0){log2rpkm_kp_rmPC=log2rpkm_kp}else{log2rpkm_kp_rmPC = removePCs(log2rpkm_kp,n_PC)}
    
    corr_ori=cor(t(log2rpkm_kp_rmPC[order(rowMeans(log2rpkm_kp)),]))
    ave_logrpkm_kp=rowMeans(log2rpkm_kp)[order(rowMeans(log2rpkm_kp))]
    
    
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
    # 
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
    
    save(d1,file=paste0(dir_output,percent_,"_percent_",n_PC,"PCs_",i,"_d1.RData"))
    save(d0,file=paste0(dir_output,percent_,"_percent_",n_PC,"PCs_",i,"_d0.RData"))
    save(d2,file=paste0(dir_output,percent_,"_percent_",n_PC,"PCs_",i,"_d2.RData"))
  }